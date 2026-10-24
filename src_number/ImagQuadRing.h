// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_IMAGQUADRING_H_
#define SRC_NUMBER_IMAGQUADRING_H_
// The number types (UniversalFloorScalarInteger) are those of the
// NumberTheory headers included by the caller, as for fractions.h.
// clang-format off
#include "Temp_common.h"
#include <cmath>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  The ring of integers O_K = Z[omega] of an imaginary quadratic field
  K = Q(omega), omega^2 = s omega - p, with s in {0, 1}:
  - D = 4m, m = -p, s = 0 for O_K = Z[sqrt(-p)] (discriminant -4p),
  - s = 1, p = (1 - D) / 4 for the discriminant D = 1 mod 4.
  The parameters s and p are given at runtime (unlike QuadField, which is
  real and of fixed d), the same as RealImagSum and RealImagProd of the
  Hermitian T-spaces.

  An element a + b omega is the pair (a, b) of the type T, which must be a
  field (mpq_class, fmpq_class, ...): the elements of O_K are the pairs of
  integers, the others are the elements of K. The norm is
    N(a + b omega) = a^2 + s a b + p b^2.

  The ring is norm-Euclidean iff the covering radius of O_K, the largest
  distance from a point of C to O_K, is < 1. Its square is computed
  exactly: (1 + p) / 4 for s = 0 and (|D| + 1)^2 / (16 |D|) for s = 1, which
  gives the norm-Euclidean rings for D = -3, -4, -7, -8, -11 only. In the
  other rings the Euclidean division fails for some pairs: the functions
  using it return an empty optional then.
 */
template <typename T> struct ImagQuadRing {
  using Elt = std::pair<T, T>;
  T s;
  T p;
  ImagQuadRing(T const &_s, T const &_p) : s(_s), p(_p) {
    if (s != 0 && s != 1) {
      std::cerr << "ImagQuadRing: s should be 0 or 1, s=" << s << "\n";
      throw TerminalException{1};
    }
    if (discriminant() >= 0) {
      std::cerr << "ImagQuadRing: the discriminant should be negative, s=" << s
                << " p=" << p << "\n";
      throw TerminalException{1};
    }
  }
  // The ring of discriminant D (D < 0, D = 0 or 1 mod 4).
  static ImagQuadRing from_discriminant(int const &D) {
    int res = ((D % 4) + 4) % 4;
    if (D >= 0 || (res != 0 && res != 1)) {
      std::cerr << "ImagQuadRing: D=" << D << " is not a negative discriminant\n";
      throw TerminalException{1};
    }
    if (res == 0) {
      return ImagQuadRing(T(0), T(-D / 4));
    }
    return ImagQuadRing(T(1), T((1 - D) / 4));
  }
  T discriminant() const { return s * s - 4 * p; }
  // Arithmetic
  static Elt zero() { return {T(0), T(0)}; }
  static Elt one() { return {T(1), T(0)}; }
  static bool is_zero(Elt const &x) { return x.first == 0 && x.second == 0; }
  static Elt add(Elt const &x, Elt const &y) {
    return {x.first + y.first, x.second + y.second};
  }
  static Elt sub(Elt const &x, Elt const &y) {
    return {x.first - y.first, x.second - y.second};
  }
  static Elt neg(Elt const &x) { return {-x.first, -x.second}; }
  Elt mul(Elt const &x, Elt const &y) const {
    T bd = x.second * y.second;
    return {x.first * y.first - p * bd,
            x.first * y.second + x.second * y.first + s * bd};
  }
  Elt conj(Elt const &x) const { return {x.first + s * x.second, -x.second}; }
  T norm(Elt const &x) const {
    return x.first * x.first + s * x.first * x.second +
           p * x.second * x.second;
  }
  T trace(Elt const &x) const { return 2 * x.first + s * x.second; }
  Elt inverse(Elt const &x) const {
    T nrm = norm(x);
    if (nrm == 0) {
      std::cerr << "ImagQuadRing: inverse of zero\n";
      throw TerminalException{1};
    }
    Elt c = conj(x);
    return {c.first / nrm, c.second / nrm};
  }
  Elt div(Elt const &x, Elt const &y) const { return mul(x, inverse(y)); }
  static bool is_integral(Elt const &x) {
    return IsInteger(x.first) && IsInteger(x.second);
  }
  // Whether x divides y in O_K (x, y integral).
  bool divides(Elt const &x, Elt const &y) const {
    if (is_zero(x)) {
      return is_zero(y);
    }
    return is_integral(div(y, x));
  }
  // The units: the elements of norm 1, 2 of them except for D = -3, -4.
  std::vector<Elt> units() const {
    std::vector<Elt> l_unit;
    for (int a = -1; a <= 1; a++) {
      for (int b = -1; b <= 1; b++) {
        Elt u{T(a), T(b)};
        if (norm(u) == 1) {
          l_unit.push_back(u);
        }
      }
    }
    return l_unit;
  }
  // The square of the covering radius of O_K in C, for the norm.
  T covering_radius_squared() const {
    if (s == 0) {
      return (1 + p) / 4;
    }
    T absD = -discriminant();
    return (absD + 1) * (absD + 1) / (16 * absD);
  }
  // A point of K at which the distance to O_K is the covering radius: the
  // circumcenter of the triangle (0, 1, omega).
  Elt deep_hole() const {
    if (s == 0) {
      return {T(1) / T(2), T(1) / T(2)};
    }
    T absD = -discriminant();
    return {(absD + 1) / (4 * absD), (absD - 1) / (2 * absD)};
  }
  bool is_norm_euclidean() const { return covering_radius_squared() < 1; }
  /*
    The elements q of O_K nearest to x (for the norm): all of them when
    there are ties. For a fixed coordinate b of q, N(x - q) is a quadratic
    function of a minimized at t = x_a + s (x_b - b) / 2, so a is floor(t)
    or floor(t) + 1. The coordinate b differs from x_b by at most
    2 R / sqrt(|D|) < 1 (R the covering radius), so the range below covers
    it with a margin.
   */
  std::vector<Elt> nearest_integers(Elt const &x) const {
    T fb = UniversalFloorScalarInteger<T, T>(x.second);
    std::optional<T> best;
    std::vector<Elt> l_best;
    for (int k = -1; k <= 2; k++) {
      T b = fb + T(k);
      T t = x.first + s * (x.second - b) / 2;
      T fa = UniversalFloorScalarInteger<T, T>(t);
      for (int l = 0; l <= 1; l++) {
        Elt q{fa + T(l), b};
        T val = norm(sub(x, q));
        if (!best || val < *best) {
          best = val;
          l_best.clear();
        }
        if (val == *best) {
          l_best.push_back(q);
        }
      }
    }
    return l_best;
  }
  // The smallest value of N(x - q) for q in O_K.
  T distance_squared(Elt const &x) const {
    return norm(sub(x, nearest_integers(x)[0]));
  }
  // x = q y + r with N(r) < N(y), q the nearest integer to x / y. Empty if
  // there is none (the ring is then not norm-Euclidean).
  std::optional<std::pair<Elt, Elt>> euclidean_division(Elt const &x,
                                                        Elt const &y) const {
    Elt q = nearest_integers(div(x, y))[0];
    Elt r = sub(x, mul(q, y));
    if (norm(r) >= norm(y)) {
      return {};
    }
    return std::pair<Elt, Elt>{q, r};
  }
  // The gcd g with u x + v y = g, as (g, u, v), by the extended Euclidean
  // algorithm. Empty if a division fails.
  std::optional<std::vector<Elt>> extended_gcd(Elt const &x,
                                               Elt const &y) const {
    Elt r0 = x, r1 = y;
    Elt u0 = one(), u1 = zero();
    Elt v0 = zero(), v1 = one();
    while (!is_zero(r1)) {
      std::optional<std::pair<Elt, Elt>> opt = euclidean_division(r0, r1);
      if (!opt) {
        return {};
      }
      Elt const &q = opt->first;
      Elt r2 = opt->second;
      Elt u2 = sub(u0, mul(q, u1));
      Elt v2 = sub(v0, mul(q, v1));
      r0 = r1;
      r1 = r2;
      u0 = u1;
      u1 = u2;
      v0 = v1;
      v1 = v2;
    }
    return std::vector<Elt>{r0, u0, v0};
  }
  std::optional<Elt> gcd(Elt const &x, Elt const &y) const {
    std::optional<std::vector<Elt>> opt = extended_gcd(x, y);
    if (!opt) {
      return {};
    }
    return (*opt)[0];
  }
  /*
    The ideal generated by the elements l_x of O_K as a lattice of Z^2:
    its Hermite basis (e1, f1), (0, e2) with e1, e2 > 0 and 0 <= f1 < e2,
    returned as {e1, f1, e2}; the ideal is nonzero. Its norm is e1 e2.
   */
  std::vector<T> ideal_lattice(std::vector<Elt> const &l_x) const {
    Elt v1 = zero();
    T e2(0);
    auto f_gcd = [](T a, T b) -> T {
      a = T_abs(a);
      b = T_abs(b);
      while (b != 0) {
        T r = a - UniversalFloorScalarInteger<T, T>(a / b) * b;
        a = b;
        b = r;
      }
      return a;
    };
    for (auto &x : l_x) {
      for (Elt w : {x, mul(x, Elt{T(0), T(1)})}) {
        while (w.first != 0) {
          if (v1.first == 0) {
            std::swap(v1, w);
            continue;
          }
          T q = UniversalFloorScalarInteger<T, T>(w.first / v1.first);
          w = sub(w, Elt{q * v1.first, q * v1.second});
          if (w.first != 0) {
            std::swap(v1, w);
          }
        }
        e2 = f_gcd(e2, w.second);
      }
    }
    if (v1.first < 0) {
      v1 = neg(v1);
    }
    if (v1.first == 0 || e2 == 0) {
      std::cerr << "ImagQuadRing: ideal_lattice, the ideal should be nonzero\n";
      throw TerminalException{1};
    }
    T f1 = v1.second - UniversalFloorScalarInteger<T, T>(v1.second / e2) * e2;
    return {v1.first, f1, e2};
  }
  // Whether x belongs to the lattice {e1, f1, e2} of ideal_lattice.
  static bool in_lattice(std::vector<T> const &L, Elt const &x) {
    if (!IsInteger(x.first / L[0])) {
      return false;
    }
    T k = x.first / L[0];
    return IsInteger((x.second - k * L[1]) / L[2]);
  }
  /*
    A generator of the ideal generated by l_x, if it is principal (always
    for a class number 1): an element of the ideal whose norm is the norm
    of the ideal, found among the elements of that norm. With
    N(a + b omega) = (a + s b / 2)^2 + (p - s^2 / 4) b^2 they have
    |b| <= sqrt(N / (p - s^2 / 4)).
   */
  std::optional<Elt> principal_generator(std::vector<Elt> const &l_x) const {
    std::vector<T> L = ideal_lattice(l_x);
    T nrm = L[0] * L[2];
    T c = p - s * s / 4;
    double b_max = std::sqrt(UniversalScalarConversion<double, T>(nrm / c)) + 1;
    int64_t b_bound = static_cast<int64_t>(b_max);
    for (int64_t b_i = -b_bound; b_i <= b_bound; b_i++) {
      T b = UniversalScalarConversion<T, int64_t>(b_i);
      T rest = nrm - c * b * b;
      if (rest < 0) {
        continue;
      }
      double sq = std::sqrt(UniversalScalarConversion<double, T>(rest));
      T center = -s * b / 2;
      int64_t a0 = static_cast<int64_t>(std::floor(UniversalScalarConversion<double, T>(center) - sq)) - 1;
      int64_t a1 = static_cast<int64_t>(std::ceil(UniversalScalarConversion<double, T>(center) + sq)) + 1;
      for (int64_t a_i = a0; a_i <= a1; a_i++) {
        Elt x{UniversalScalarConversion<T, int64_t>(a_i), b};
        if (norm(x) == nrm && in_lattice(L, x)) {
          return x;
        }
      }
    }
    return {};
  }
  // The gcd of the elements l_x: by the Euclidean algorithm when it
  // applies, otherwise a generator of the ideal they generate (empty if
  // it is not principal).
  std::optional<Elt> gcd(std::vector<Elt> const &l_x) const {
    Elt g = zero();
    bool euclid_ok = true;
    for (auto &x : l_x) {
      std::optional<Elt> opt = gcd(g, x);
      if (!opt) {
        euclid_ok = false;
        break;
      }
      g = *opt;
    }
    if (euclid_ok) {
      return g;
    }
    return principal_generator(l_x);
  }
  // The largest of the unit multiples of x in the lexicographic order of
  // the pairs.
  Elt canonical_associate(Elt const &x) const {
    Elt best = x;
    for (auto &u : units()) {
      Elt xu = mul(u, x);
      if (best < xu) {
        best = xu;
      }
    }
    return best;
  }
};

#endif // SRC_NUMBER_IMAGQUADRING_H_
