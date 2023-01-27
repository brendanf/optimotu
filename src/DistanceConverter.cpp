#include "DistanceConverter.hpp"
#include <cmath>

d_t UniformDistanceConverter::convert(double dist) const {
   d_t i = std::max((int) ceilf((dist - dmin) / dstep), 0);
   return i;
}

UniformDistanceConverter::UniformDistanceConverter(
   const double dmin,
   const double dstep
): dmin(dmin), dstep(dstep) {}

d_t ArrayDistanceConverter::convert(double dist) const {
   if (dist > max_threshold) return thresholds.size();
   auto thresh_i = std::lower_bound(thresholds.begin(), thresholds.end(), dist);
   d_t i = std::distance(thresholds.begin(), thresh_i);
   return i;
}

ArrayDistanceConverter::ArrayDistanceConverter(std::vector<double> thresholds):
      thresholds(thresholds),
      max_threshold(thresholds[thresholds.size() - 1]) {}
