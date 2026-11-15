// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "gpu_tile.h"
#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#include "optimotu.h"

#include <algorithm>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <utility>

#ifdef OPTIMOTU_CUDA
#include <cuda_runtime.h>
#endif

#ifdef __linux__
#include <pthread.h>
#include <sched.h>
#endif

std::size_t GpuTileSpec::n_pairs() const {
  if (kind == Kind::Triangle) {
    if (n < 2) {
      return 0;
    }
    return static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1) / 2;
  }
  const int ni = i1 - i0;
  const int nj = j1 - j0;
  if (ni <= 0 || nj <= 0) {
    return 0;
  }
  return static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
}

namespace {

std::size_t gpu_pair_index(const GpuTileSpec &spec, int i0, int j0) {
  if (spec.kind == GpuTileSpec::Kind::Triangle) {
    int i = i0 - spec.offset;
    int j = j0 - spec.offset;
    if (i < j) {
      std::swap(i, j);
    }
    return static_cast<std::size_t>(i) * static_cast<std::size_t>(i - 1) / 2 +
      static_cast<std::size_t>(j);
  }
  const int nj = spec.j1 - spec.j0;
  const int row = i0 - spec.i0;
  const int col = j0 - spec.j0;
  return static_cast<std::size_t>(row) * static_cast<std::size_t>(nj) +
    static_cast<std::size_t>(col);
}

}  // namespace

double gpu_result_distance(
  const GpuTileSpec &spec,
  const GpuTileResult &result,
  int i0,
  int j0
) {
  return result.distances[gpu_pair_index(spec, i0, j0)];
}

int gpu_result_score(
  const GpuTileSpec &spec,
  const GpuTileResult &result,
  int i0,
  int j0
) {
  return result.scores[gpu_pair_index(spec, i0, j0)];
}

GpuTileSpec gpu_tile_from_generator(const PairGenerator &pg) {
  GpuTileSpec spec;
  if (auto *diag = dynamic_cast<const AllPairGenerator *>(&pg)) {
    spec.kind = GpuTileSpec::Kind::Triangle;
    spec.offset = static_cast<int>(diag->seq_offset());
    spec.n = static_cast<int>(diag->seq_count());
    return spec;
  }
  if (auto *bip = dynamic_cast<const BipartitePairGenerator *>(&pg)) {
    spec.kind = GpuTileSpec::Kind::Rectangle;
    spec.i0 = static_cast<int>(bip->seq_begin_i());
    spec.i1 = static_cast<int>(bip->seq_end_i());
    spec.j0 = static_cast<int>(bip->seq_begin_j());
    spec.j1 = static_cast<int>(bip->seq_end_j());
    return spec;
  }
  OPTIMOTU_STOP("GPU tiles require an all-pairs or bipartite generator");
  return spec;
}

#ifndef OPTIMOTU_CUDA

bool gpu_runtime_available() {
  return false;
}

int gpu_device_count() {
  return 0;
}

std::vector<int> gpu_selected_devices() {
  return {};
}

std::string gpu_device_name(int) {
  return "";
}

void gpu_run_specs(
  const std::vector<GpuTileSpec> &,
  std::size_t,
  const std::function<std::unique_ptr<GpuDistanceKernel>()> &,
  const std::function<
    void(std::size_t, const GpuTileSpec &, const GpuTileResult &)> &
) {
  OPTIMOTU_STOP("this build of optimotu was compiled without CUDA");
}

void gpu_run_tiles(
  const std::vector<std::unique_ptr<PairGenerator>> &,
  std::size_t,
  const std::function<std::unique_ptr<GpuDistanceKernel>()> &,
  const std::function<
    void(std::size_t, PairGenerator &, const GpuTileResult &)> &
) {
  OPTIMOTU_STOP("this build of optimotu was compiled without CUDA");
}

#else

namespace {

void check_cuda(cudaError_t st, const char *what) {
  if (st != cudaSuccess) {
    OPTIMOTU_STOP(
      std::string(what) + ": " + cudaGetErrorString(st)
    );
  }
}

std::vector<int> parse_cpulist(const std::string &text) {
  std::vector<int> out;
  std::stringstream ss(text);
  std::string part;
  while (std::getline(ss, part, ',')) {
    if (part.empty()) {
      continue;
    }
    const auto dash = part.find('-');
    if (dash == std::string::npos) {
      out.push_back(std::atoi(part.c_str()));
      continue;
    }
    const int a = std::atoi(part.substr(0, dash).c_str());
    const int b = std::atoi(part.substr(dash + 1).c_str());
    for (int c = a; c <= b; ++c) {
      out.push_back(c);
    }
  }
  return out;
}

int numa_from_sysfs(int device) {
  cudaDeviceProp prop{};
  if (cudaGetDeviceProperties(&prop, device) != cudaSuccess) {
    return -1;
  }
  char path[128];
  std::snprintf(
    path,
    sizeof(path),
    "/sys/bus/pci/devices/%04x:%02x:%02x.0/numa_node",
    prop.pciDomainID,
    prop.pciBusID,
    prop.pciDeviceID
  );
  std::ifstream in(path);
  int n = -1;
  if (in) {
    in >> n;
  }
  return n;
}

int host_numa_id(int device) {
  int numa = -1;
#if defined(cudaDevAttrHostNumaId)
  if (cudaDeviceGetAttribute(&numa, cudaDevAttrHostNumaId, device) !=
      cudaSuccess) {
    numa = -1;
  }
#endif
  if (numa >= 0) {
    return numa;
  }
  return numa_from_sysfs(device);
}

std::vector<int> cpus_for_device(int device) {
  const int numa = host_numa_id(device);
  if (numa < 0) {
    return {};
  }
  std::ifstream in(
    "/sys/devices/system/node/node" + std::to_string(numa) + "/cpulist"
  );
  if (!in) {
    return {};
  }
  std::string text;
  in >> text;
  return parse_cpulist(text);
}

void pin_current_thread(const std::vector<int> &cpus) {
#ifdef __linux__
  if (cpus.empty()) {
    return;
  }
  cpu_set_t set;
  CPU_ZERO(&set);
  for (int c : cpus) {
    if (c >= 0 && c < CPU_SETSIZE) {
      CPU_SET(c, &set);
    }
  }
  pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
#else
  (void)cpus;
#endif
}

std::vector<int> parse_device_env(int n_dev) {
  const char *env = std::getenv("OPTIMOTU_CUDA_DEVICES");
  if (env == nullptr || env[0] == '\0') {
    std::vector<int> all(static_cast<std::size_t>(n_dev));
    for (int i = 0; i < n_dev; ++i) {
      all[static_cast<std::size_t>(i)] = i;
    }
    return all;
  }
  std::vector<int> out;
  std::stringstream ss(env);
  std::string part;
  while (std::getline(ss, part, ',')) {
    if (part.empty()) {
      continue;
    }
    const int d = std::atoi(part.c_str());
    if (d < 0 || d >= n_dev) {
      OPTIMOTU_STOP(
        "OPTIMOTU_CUDA_DEVICES lists a device that is not present"
      );
    }
    out.push_back(d);
  }
  if (out.empty()) {
    OPTIMOTU_STOP("OPTIMOTU_CUDA_DEVICES is empty");
  }
  return out;
}

struct TileJob {
  std::size_t index = 0;
  GpuTileSpec spec;
  GpuTileResult result;
  std::mutex mu;
  std::condition_variable cv;
  bool done = false;
  std::string error;
};

void pump_loop(
  int device,
  GpuDistanceKernel *kernel,
  const std::vector<int> &cpus,
  std::queue<TileJob *> *queue,
  std::mutex *qmu,
  std::condition_variable *qcv,
  int *stop
) {
  pin_current_thread(cpus);
  check_cuda(cudaSetDevice(device), "cudaSetDevice");
  kernel->bind_device(device);
  while (true) {
    TileJob *job = nullptr;
    {
      std::unique_lock<std::mutex> lock(*qmu);
      qcv->wait(lock, [&] { return *stop != 0 || !queue->empty(); });
      if (*stop != 0 && queue->empty()) {
        return;
      }
      job = queue->front();
      queue->pop();
    }
    try {
      const GpuTileSpec spec = job->spec;
      const std::size_t np = spec.n_pairs();
      job->result.kind = kernel->kind();
      if (kernel->kind() == GpuKernelKind::Distance) {
        job->result.distances.assign(np, 1.0);
        if (np > 0) {
          kernel->fill_distances(spec, job->result.distances.data());
        }
      } else {
        job->result.scores.assign(np, 0);
        if (np > 0) {
          kernel->fill_scores(spec, job->result.scores.data());
        }
      }
    } catch (const std::exception &e) {
      job->error = e.what();
    }
    {
      std::lock_guard<std::mutex> lock(job->mu);
      job->done = true;
    }
    job->cv.notify_one();
  }
}

}  // namespace

bool gpu_runtime_available() {
  int n = 0;
  const cudaError_t st = cudaGetDeviceCount(&n);
  return st == cudaSuccess && n > 0;
}

int gpu_device_count() {
  int n = 0;
  if (cudaGetDeviceCount(&n) != cudaSuccess) {
    return 0;
  }
  return n;
}

std::vector<int> gpu_selected_devices() {
  const int n = gpu_device_count();
  if (n <= 0) {
    return {};
  }
  return parse_device_env(n);
}

std::string gpu_device_name(int device) {
  cudaDeviceProp prop{};
  if (cudaGetDeviceProperties(&prop, device) != cudaSuccess) {
    return "";
  }
  return prop.name;
}

void gpu_run_specs(
  const std::vector<GpuTileSpec> &specs,
  std::size_t threads,
  const std::function<std::unique_ptr<GpuDistanceKernel>()> &make_kernel,
  const std::function<
    void(std::size_t, const GpuTileSpec &, const GpuTileResult &)> &consume
) {
  const auto devices = gpu_selected_devices();
  if (devices.empty()) {
    OPTIMOTU_STOP("no CUDA device is available");
  }
  const std::size_t n_tiles = specs.size();
  if (n_tiles == 0) {
    return;
  }
  const std::size_t n_dev = devices.size();
  std::vector<std::unique_ptr<GpuDistanceKernel>> kernels(n_dev);
  std::vector<std::queue<TileJob *>> queues(n_dev);
  std::vector<std::mutex> qmu(n_dev);
  std::vector<std::condition_variable> qcv(n_dev);
  std::vector<int> stop(n_dev, 0);
  std::vector<std::thread> pumps;
  pumps.reserve(n_dev);
  for (std::size_t d = 0; d < n_dev; ++d) {
    kernels[d] = make_kernel();
    const std::vector<int> cpus = cpus_for_device(devices[d]);
    pumps.emplace_back(
      pump_loop,
      devices[d],
      kernels[d].get(),
      cpus,
      &queues[d],
      &qmu[d],
      &qcv[d],
      &stop[d]
    );
  }

  std::vector<std::vector<std::size_t>> by_dev(n_dev);
  for (std::size_t t = 0; t < n_tiles; ++t) {
    by_dev[t % n_dev].push_back(t);
  }
  std::size_t workers_per = std::max<std::size_t>(1, threads / n_dev);
  if (threads < n_dev) {
    workers_per = 1;
  }

  std::vector<std::thread> workers;
  std::string worker_error;
  std::mutex err_mu;
  for (std::size_t d = 0; d < n_dev; ++d) {
    const std::size_t n_w = std::min(
      workers_per,
      std::max<std::size_t>(1, by_dev[d].size())
    );
    const std::vector<int> cpus = cpus_for_device(devices[d]);
    for (std::size_t w = 0; w < n_w; ++w) {
      workers.emplace_back([&, d, w, n_w, cpus]() {
        pin_current_thread(cpus);
        const auto &mine = by_dev[d];
        for (std::size_t k = w; k < mine.size(); k += n_w) {
          const std::size_t idx = mine[k];
          TileJob job;
          job.index = idx;
          job.spec = specs[idx];
          {
            std::lock_guard<std::mutex> lock(qmu[d]);
            queues[d].push(&job);
          }
          qcv[d].notify_one();
          {
            std::unique_lock<std::mutex> lock(job.mu);
            job.cv.wait(lock, [&] { return job.done; });
          }
          if (!job.error.empty()) {
            std::lock_guard<std::mutex> lock(err_mu);
            if (worker_error.empty()) {
              worker_error = job.error;
            }
            return;
          }
          try {
            consume(idx, specs[idx], job.result);
          } catch (const std::exception &e) {
            std::lock_guard<std::mutex> lock(err_mu);
            if (worker_error.empty()) {
              worker_error = e.what();
            }
            return;
          }
        }
      });
    }
  }
  for (auto &th : workers) {
    th.join();
  }
  for (std::size_t d = 0; d < n_dev; ++d) {
    {
      std::lock_guard<std::mutex> lock(qmu[d]);
      stop[d] = 1;
    }
    qcv[d].notify_all();
  }
  for (auto &th : pumps) {
    th.join();
  }
  if (!worker_error.empty()) {
    OPTIMOTU_STOP(worker_error);
  }
}

void gpu_run_tiles(
  const std::vector<std::unique_ptr<PairGenerator>> &generators,
  std::size_t threads,
  const std::function<std::unique_ptr<GpuDistanceKernel>()> &make_kernel,
  const std::function<
    void(std::size_t, PairGenerator &, const GpuTileResult &)> &consume
) {
  std::vector<GpuTileSpec> specs;
  specs.reserve(generators.size());
  for (const auto &g : generators) {
    specs.push_back(gpu_tile_from_generator(*g));
  }
  gpu_run_specs(
    specs,
    threads,
    make_kernel,
    [&](std::size_t i, const GpuTileSpec &, const GpuTileResult &res) {
      consume(i, *generators[i], res);
    }
  );
}

#endif

#ifdef OPTIMOTU_R
#include <Rcpp.h>

//' Whether a CUDA device is usable
//'
//' @name optimotu_cuda
//' @return `TRUE` when this build linked CUDA and at least one device
//'   is visible (after `OPTIMOTU_CUDA_DEVICES` if that is set).
//' @export
// [[Rcpp::export]]
bool optimotu_cuda() {
  return gpu_runtime_available() && !gpu_selected_devices().empty();
}

//' Whether this build of optimotu was compiled with CUDA
//'
//' @name optimotu_cuda_build
//' @return `TRUE` when `.cu` kernels were compiled (`OPTIMOTU_CUDA`).
//' @keywords internal
//' @export
// [[Rcpp::export]]
bool optimotu_cuda_build() {
#ifdef OPTIMOTU_CUDA
  return true;
#else
  return false;
#endif
}

//' Number of selected CUDA devices
//'
//' @name optimotu_cuda_device_count
//' @return integer count of devices `optimotu_cuda()` would use
//' @keywords internal
//' @export
// [[Rcpp::export]]
int optimotu_cuda_device_count() {
  return static_cast<int>(gpu_selected_devices().size());
}
#endif
