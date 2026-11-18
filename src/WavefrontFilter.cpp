#include "WavefrontFilter.h"
#include "wfa_fatal.h"

WavefrontFilter::WavefrontFilter(const std::vector<std::string> & seq):
  PrealignFilter(seq), aligner(wfa::WFAligner::Score) {}
bool WavefrontFilter::operator()(
    const std::size_t i,
    const std::size_t j,
    const int max_ed,
    const int min_k,
    const int max_k
) {
  if (optimotu_wfa_fatal_pending() || aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return false;
  }
  aligner.setMaxAlignmentSteps(max_ed);
  aligner.setHeuristicBandedStatic(min_k, max_k);
  auto status = aligner.alignEnd2End(seq[i], seq[j]);
  if (aligner.fatal() || status == wfa::WFAligner::StatusOOM)
  {
    optimotu_wfa_note_fatal();
    return false;
  }
  return status == wfa::WFAligner::AlignmentStatus::StatusAlgCompleted;
}
std::unique_ptr<PrealignFilter> WavefrontFilter::copy() const {
  return std::make_unique<WavefrontFilter>(seq);
}
