#include "AlignClusterWorker.h"

size_t AlignClusterWorker::prealigned() {
  return _prealigned;
}

size_t AlignClusterWorker::aligned() {
  return _aligned;
}

std::uint8_t AlignClusterWorker::nthreads() {
  return threads;
}
