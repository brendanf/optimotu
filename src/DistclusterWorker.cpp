// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "DistClusterWorker.h"

size_t DistClusterWorker::prealigned() {
  return _prealigned;
}

size_t DistClusterWorker::aligned() {
  return _aligned;
}

std::size_t DistClusterWorker::n_threads() {
  return threads;
}

ClusterAlgorithm *DistClusterWorker::tile_algo(PairGenerator *pg)
{
  if (threads == 1)
  {
    return &clust_algo;
  }
  return clust_algo.make_child(pg);
}

void DistClusterWorker::finish_tile(ClusterAlgorithm *my_algo)
{
  if (!my_algo)
  {
    return;
  }
  if (threads == 1)
  {
    return;
  }
  my_algo->finalize();
  my_algo->merge_into_parent();
  clust_algo.release_child(my_algo);
}
