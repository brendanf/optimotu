#include <cstdio>

/*
 * Opened from C++ so the C force-include does not rewrite this
 * translation unit's own stderr declaration.
 */
extern "C" {
FILE *optimotu_wfa_null = nullptr;
}

namespace {

struct OptimotuWfaNullInit {
  OptimotuWfaNullInit() {
#ifdef _WIN32
    optimotu_wfa_null = std::fopen("nul", "w");
#else
    optimotu_wfa_null = std::fopen("/dev/null", "w");
#endif
  }
  ~OptimotuWfaNullInit() {
    if (optimotu_wfa_null != nullptr) {
      std::fclose(optimotu_wfa_null);
      optimotu_wfa_null = nullptr;
    }
  }
};

OptimotuWfaNullInit optimotu_wfa_null_init;

}  // namespace
