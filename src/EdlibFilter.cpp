#include "EdlibFilter.h"

EdlibFilter::EdlibFilter(const std::vector<std::string> & seq):
  PrealignFilter(seq){}
bool EdlibFilter::operator()(
    const std::size_t i,
    const std::size_t j,
    const int max_ed,
    const int min_k,
    const int max_k
) {
  aligner.k = max_ed;
  auto result = edlibAlign(seq[i].c_str(), seq[i].size(),
                           seq[j].c_str(), seq[j].size(),
                           aligner);
  int status = result.status;
  int ed = result.editDistance;
  edlibFreeAlignResult(result);
  return (status == EDLIB_STATUS_OK && ed != -1);
}
std::unique_ptr<PrealignFilter> EdlibFilter::copy() const {
  return std::make_unique<EdlibFilter>(seq);
}
