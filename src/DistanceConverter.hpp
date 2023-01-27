#ifndef _DISTANCE_CONVERTER_
#define _DISTANCE_CONVERTER_

#include <vector>
#include "single_linkage_pool.hpp"

class DistanceConverter
{
public:
   virtual d_t convert(double dist) const = 0;
};

class UniformDistanceConverter : public DistanceConverter
{
   const double dmin, dstep;
public:
   d_t convert(double dist) const override;
   UniformDistanceConverter(const double dmin, const double dstep);
};

class ArrayDistanceConverter : public DistanceConverter
{
   const std::vector<double> thresholds;
   const double max_threshold;
public:
   d_t convert(double dist) const override;
   ArrayDistanceConverter(std::vector<double> thresholds);
};

#endif
