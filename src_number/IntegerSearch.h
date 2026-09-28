// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_INTEGERSEARCH_H_
#define SRC_NUMBER_INTEGERSEARCH_H_

// clang-format off
#include <cmath>
// clang-format on

/*
  The exact search behind the floor, the ceiling and the nearest integer of
  the real algebraic types (QuadField, RealField), whose elements no double
  represents exactly. It is not meant to be called outside of src_number:
  the callers use UniversalFloorScalarInteger, UniversalCeilScalarInteger and
  UniversalNearestScalarInteger.

  It returns the largest n in Z with pred(n), for a predicate that holds on
  every integer below some threshold and on none above it. The three
  roundings of x are each one such search, with no arithmetic on x:

      floor(x)   = largest n with n <= x
      ceil(x)    = 1 + largest n with n < x
      nearest(x) = largest n with n - 1/2 < x   (a tie y + 1/2 goes to y)

  The answer is fixed by pred alone, so it does not depend on how the search
  starts. start_d is a floating point estimate of the answer: when it is
  finite and small enough for the integers to be exact in a double, the
  search starts there and a couple of evaluations of pred correct it.
  Otherwise the threshold is bracketed by powers of two and the bracket is
  halved, both logarithmic in the value, which is what makes the search safe
  on an element no double can locate.
 */
template <typename Tint, typename Fpred>
Tint helper_largest_integer_satisfying(double const &start_d,
                                       Fpred const &pred) {
  if (std::isfinite(start_d) && std::abs(start_d) < 9e15) {
    Tint n(static_cast<long>(start_d));
    while (!pred(n)) {
      n -= 1;
    }
    while (true) {
      Tint n_next = n + 1;
      if (!pred(n_next)) {
        return n;
      }
      n = n_next;
    }
  }
  // The bracket: pred(lo) holds and pred(hi) does not.
  Tint lo(-1);
  Tint hi(1);
  while (!pred(lo)) {
    lo *= 2;
  }
  while (pred(hi)) {
    hi *= 2;
  }
  while (hi - lo > 1) {
    // Strictly between lo and hi since hi - lo >= 2, whichever way the
    // division rounds.
    Tint mid = (lo + hi) / 2;
    if (pred(mid)) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

// clang-format off
#endif  // SRC_NUMBER_INTEGERSEARCH_H_
// clang-format on
