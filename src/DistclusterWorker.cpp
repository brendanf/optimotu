// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "DistClusterWorker.h"

size_t DistClusterWorker::prealigned() {
  return _prealigned;
}

size_t DistClusterWorker::aligned() {
  return _aligned;
}

std::uint8_t DistClusterWorker::n_threads() {
  return threads;
}
