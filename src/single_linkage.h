#ifndef OPTIMOTU_SINGLE_LINKAGE_H_INCLUDED
#define OPTIMOTU_SINGLE_LINKAGE_H_INCLUDED

#include "optimotu.h"
#include <cstdint>
#include <limits>
#include <cstdint>

// type of "j", which is the index of a cluster
typedef std::uint32_t j_t;
// type of "d", which is the index corresponding to a distance
typedef std::int32_t d_t;

static const j_t NO_CLUST = std::numeric_limits<std::uint32_t>::max();
static const d_t NO_DIST = std::numeric_limits<std::int32_t>::max();
#endif //OPTIMOTU_SINGLELINKAGE_H_INCLUDED
