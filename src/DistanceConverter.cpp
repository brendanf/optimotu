#include "DistanceConverter.h"
#include <cmath>

d_t UniformDistanceConverter::convert(double dist) const {
   d_t i = std::max((int) ceil((dist - dmin) / dstep), 0);
   return i;
}

double UniformDistanceConverter::inverse(d_t d) const {
   double dist = dmin + d * dstep;
   return dist;
}

UniformDistanceConverter::UniformDistanceConverter(
   const double dmin,
   const double dmax,
   const double dstep
): DistanceConverter(std::ceil((dmax - dmin)/dstep) + 1, dmax), dmin(dmin),
dstep(dstep) {}

d_t ArrayDistanceConverter::convert(double dist) const {
   if (dist > dmax) return thresholds.size();
   auto thresh_i = std::lower_bound(thresholds.begin(), thresholds.end(), dist);
   d_t i = std::distance(thresholds.begin(), thresh_i);
   return i;
}

double ArrayDistanceConverter::inverse(d_t d) const {
   if (d >= (d_t)thresholds.size()) return dmax;
   if (d < 0) return thresholds[0];
   return thresholds[d];
}

ArrayDistanceConverter::ArrayDistanceConverter(std::vector<double> thresholds):
  DistanceConverter(thresholds.size() - 1, thresholds[thresholds.size() - 1]),
      thresholds(thresholds) {}

d_t CachedDistanceConverter::convert(double dist) const {
   if (dist > dmax) return thresholds.size();
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
   cache(init_cache(thresholds, dmax, precision)) {}

