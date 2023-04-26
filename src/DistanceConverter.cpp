#include "DistanceConverter.h"
#include <cmath>

d_t UniformDistanceConverter::convert(double dist) const {
   d_t i = std::max((int) ceilf((dist - dmin) / dstep), 0);
   return i;
}

double UniformDistanceConverter::inverse(d_t d) const {
   double dist = dmin + d * dstep;
   return dist;
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

double ArrayDistanceConverter::inverse(d_t d) const {
   if (d >= (d_t)thresholds.size()) return max_threshold;
   if (d < 0) return thresholds[0];
   return thresholds[d];
}

ArrayDistanceConverter::ArrayDistanceConverter(std::vector<double> thresholds):
      thresholds(thresholds),
      max_threshold(thresholds[thresholds.size() - 1]) {}

d_t CachedDistanceConverter::convert(double dist) const {
   if (dist > max_threshold) return thresholds.size();
   size_t i = 0;
   if (dist > thresholds[0]) {
      i = (size_t)round((dist - thresholds[0])/precision);
   }
   return cache[i];
}

std::vector<d_t> init_cache(
      std::vector<double> thresholds,
      double max_threshold,
      double precision
) {
   size_t cache_size = (size_t)ceil((max_threshold - thresholds[0]) / precision);
   std::vector<d_t> cache;
   cache.reserve(cache_size);
   d_t i = 0;
   for (size_t j = 0; j <= cache_size; j++) {
      double t = thresholds[0] + precision * j;
      while(t - thresholds[i] > sqrt(std::numeric_limits<double>::epsilon())) i++;
      cache.push_back(i);
   }
   return cache;
}

CachedDistanceConverter::CachedDistanceConverter(
   std::vector<double> thresholds,
   double precision
):
   ArrayDistanceConverter(thresholds),
   precision(precision),
   cache(init_cache(thresholds, max_threshold, precision)) {}

