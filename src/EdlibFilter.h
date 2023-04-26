#ifndef OPTIMOTU_EDLIBFILTER_H_INCLUDED
#define OPTIMOTU_EDLIBFILTER_H_INCLUDED

#include <memory>
#include <edlib.h>
#include "PrealignFilter.h"
#include "pad_strings.h"

class EdlibFilter : public PrealignFilter {
private:
  EdlibAlignConfig aligner = edlibDefaultAlignConfig();
public:
  EdlibFilter(const std::vector<std::string> & seq);
  bool operator()(
      const std::size_t i,
      const std::size_t j,
      const int max_ed,
      const int min_k,
      const int max_k
  ) override;
  std::unique_ptr<PrealignFilter> copy() const override;
};

#endif // OPTIMOTU_EDLIBFILTER_H_INCLUDED
