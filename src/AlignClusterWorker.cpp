// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "AlignClusterWorker.h"

size_t AlignClusterWorker::prealigned() {
  return _prealigned;
}

size_t AlignClusterWorker::aligned() {
  return _aligned;
}

std::uint8_t AlignClusterWorker::n_threads() {
  return threads;
}
