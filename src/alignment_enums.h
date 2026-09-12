#ifndef OPTIMOTU_ALIGNMENT_ENUMS_H_INCLUDED
#define OPTIMOTU_ALIGNMENT_ENUMS_H_INCLUDED

enum class AlignmentSpan {
  GLOBAL = 0,
  EXTEND = 1
};

// Dependent false for discarded if-constexpr else branches.
// `static_assert(span != span)` triggers -Wtautological-compare on
// GCC/Windows even when the branch is not instantiated.
template <auto>
inline constexpr bool always_false_v = false;

#endif //OPTIMOTU_ALIGNMENT_ENUMS_H_INCLUDED
