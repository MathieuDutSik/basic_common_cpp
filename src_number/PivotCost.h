// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_PIVOTCOST_H_
#define SRC_NUMBER_PIVOTCOST_H_

// clang-format off
#include <cstddef>
// clang-format on

// Pivot selection is a classic, type-dependent problem in Gaussian elimination.
// The elimination kernels select, among the non-zero candidates, the pivot with
// the most preferable "cost", where the cost and the preference direction are
// specialized per numeric type:
//   * floating point : the pivot of largest magnitude (numerical stability).
//   * exact rational : the one with the smallest numerator + denominator size
//                      (keeps intermediate fractions small).
//   * jet<T, N>      : the one of smallest order, i.e. non-zero constant term
//                      first (correctness: only such a pivot is invertible;
//                      the inverse of a jet with a leading t^k, k >= 1, would be
//                      a Laurent series). See jet_number.h.
//
// A numeric type participates by providing:
//   f_cost_pivot(x)                 -> its pivot cost (the "delegated type")
//   is_preferable_pivot(cx, cy)     -> true iff cost cx is a better pivot than cy
//
// Types without a specialization fall back to the trivial cost below, for which
// no candidate is ever preferred, so the elimination keeps the first non-zero
// pivot it encounters -- reproducing the historical behaviour exactly.

struct PivotCostTrivial {};

template <typename T> PivotCostTrivial f_cost_pivot(T const &) {
  return PivotCostTrivial{};
}
inline bool is_preferable_pivot(PivotCostTrivial const &,
                                PivotCostTrivial const &) {
  return false;
}

// Floating point: cost is the magnitude, larger is preferred.
inline double f_cost_pivot(double const &x) { return x < 0 ? -x : x; }
inline float f_cost_pivot(float const &x) { return x < 0 ? -x : x; }
inline bool is_preferable_pivot(double const &x, double const &y) {
  return x > y;
}
inline bool is_preferable_pivot(float const &x, float const &y) {
  return x > y;
}

// Integer / rational size cost (number of machine words): smaller is preferred.
inline bool is_preferable_pivot(size_t const &x, size_t const &y) {
  return x < y;
}

// clang-format off
#endif  // SRC_NUMBER_PIVOTCOST_H_
// clang-format on
