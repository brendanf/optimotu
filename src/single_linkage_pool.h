#ifndef _SINGLE_LINKAGE_POOL_
#define _SINGLE_LINKAGE_POOL_
#include <cstdint>
#include <limits>
// type of "j", which is the index of a cluster
typedef std::uint32_t j_t;
// type of "d", which is the index corresponding to a distance
typedef std::int32_t d_t;

static const j_t NO_CLUST = std::numeric_limits<std::uint32_t>::max();
static const d_t NO_DIST = std::numeric_limits<std::int32_t>::max();

#endif
