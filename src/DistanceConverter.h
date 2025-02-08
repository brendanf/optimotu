// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_DISTANCECONVERTER_H_INCLUDED
#define OPTIMOTU_DISTANCECONVERTER_H_INCLUDED

#include <vector>
#include "single_linkage.h"

class DistanceConverter
{
public:
  const d_t m;
  const double dmax;
  DistanceConverter(const d_t m, const double dmax) : m(m), dmax(dmax) {};

  virtual d_t convert(double dist) const = 0;
  virtual double inverse(d_t d) const = 0;
  virtual ~DistanceConverter() = default;
};

class UniformDistanceConverter : public DistanceConverter
{
protected:
   const double dmin, dstep;
public:
   d_t convert(double dist) const override;
   double inverse(d_t d) const override;
   UniformDistanceConverter(const double dmin, const double dmax, const double dstep);
};

class ArrayDistanceConverter : public DistanceConverter
{
protected:
   const std::vector<double> thresholds;
public:
   d_t convert(double dist) const override;
   double inverse(d_t d) const override;
   ArrayDistanceConverter(std::vector<double> thresholds);
};

class CachedDistanceConverter : public ArrayDistanceConverter
{
   const double precision;
   const std::vector<d_t> cache;
public:
   d_t convert(double dist) const override;
   CachedDistanceConverter(std::vector<double> thresholds, double precision);
};

#endif //OPTIMOTU_DISTANCECONVERTER_H_INCLUDED
