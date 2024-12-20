#include "WavefrontFilter.h"

WavefrontFilter::WavefrontFilter(const std::vector<std::string> & seq):
  PrealignFilter(seq), aligner(wfa::WFAligner::Score) {}
bool WavefrontFilter::operator()(
    const std::size_t i,
    const std::size_t j,
    const int max_ed,
    const int min_k,
    const int max_k
) {
  aligner.setMaxAlignmentSteps(max_ed);
  aligner.setHeuristicBandedStatic(min_k, max_k);
  aligner.alignEnd2End(seq[i], seq[j]);
  auto status = aligner.getAlignmentStatus();
  return status == wfa::WFAligner::AlignmentStatus::StatusAlgCompleted;
}
std::unique_ptr<PrealignFilter> WavefrontFilter::copy() const {
  return std::make_unique<WavefrontFilter>(seq);
}
