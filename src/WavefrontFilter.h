#ifndef OPTIMOTU_WAVEFRONTFILTER_H_INCLUDED
#define OPTIMOTU_WAVEFRONTFILTER_H_INCLUDED

#include <memory>
#include <bindings/cpp/WFAligner.hpp>
#include "PrealignFilter.h"
#include "pad_strings.h"

class WavefrontFilter : public PrealignFilter {
private:
  wfa::WFAlignerEdit aligner;
  WavefrontFilter(const WavefrontFilter &obj);
public:
  WavefrontFilter(const std::vector<std::string> & seq);
  bool operator()(
      const std::size_t i,
      const std::size_t j,
      const int max_ed,
      const int min_k,
      const int max_k
  ) override;
  std::unique_ptr<PrealignFilter> copy() const override;
};

#endif // OPTIMOTU_WAVEFRONTFILTER_H_INCLUDED
