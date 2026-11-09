// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_POLYNOMIALREALROOTS_H_
#define SRC_NUMBER_POLYNOMIALREALROOTS_H_
// clang-format off
#include "Basic_functions.h"
#include <algorithm>
#include <utility>
#include <vector>
// clang-format on

/*
  Approximations of real numbers with exact rational arithmetic, generic
  over the rational types (mpq_class, fmpq_class, boost::multiprecision
  rationals): no floating point, so the error bounds are rigorous.

  - real_roots_approximation(poly, eps): the real roots, in increasing
    order, of the squarefree polynomial poly = sum_i poly[i] x^i with
    rational coefficients, each given by a rational within eps of it. The
    roots are isolated by the Sturm sequence f, f', -rem(f, f'), ... and
    bisection of (-B, B], B = 1 + max |c_i / c_d| (Cauchy), then each is
    refined by bisection on the sign of f.
  - square_root_approximation(x, eps): for x >= 0, a rational within eps of
    the square root of x, by bisection.
 */

template <typename Tq>
Tq polynomial_evaluation(std::vector<Tq> const &poly, Tq const &x) {
  Tq res(0);
  for (size_t i = poly.size(); i-- > 0;) {
    res *= x;
    res += poly[i];
  }
  return res;
}

template <typename Tq> void polynomial_trim(std::vector<Tq> &poly) {
  while (!poly.empty() && poly.back() == 0) {
    poly.pop_back();
  }
}

// the remainder of a by b, b nonzero
template <typename Tq>
std::vector<Tq> polynomial_remainder(std::vector<Tq> a,
                                     std::vector<Tq> const &b) {
  polynomial_trim(a);
  size_t db = b.size() - 1;
  while (a.size() >= b.size()) {
    Tq c = a.back() / b.back();
    size_t shift = a.size() - 1 - db;
    for (size_t j = 0; j <= db; j++) {
      Tq t = c * b[j];
      a[shift + j] -= t;
    }
    a.pop_back();
    polynomial_trim(a);
  }
  return a;
}

template <typename Tq>
std::vector<std::vector<Tq>> sturm_sequence(std::vector<Tq> const &poly) {
  std::vector<std::vector<Tq>> seq{poly};
  std::vector<Tq> der;
  for (size_t i = 1; i < poly.size(); i++) {
    Tq coef(0);
    for (size_t k = 0; k < i; k++) {
      coef += poly[i];
    }
    der.push_back(coef);
  }
  polynomial_trim(der);
  if (der.empty()) {
    return seq;
  }
  seq.push_back(der);
  while (true) {
    std::vector<Tq> r =
        polynomial_remainder(seq[seq.size() - 2], seq[seq.size() - 1]);
    if (r.empty()) {
      break;
    }
    for (auto &x : r) {
      x = -x;
    }
    seq.push_back(r);
  }
  return seq;
}

template <typename Tq>
size_t sturm_sign_changes(std::vector<std::vector<Tq>> const &seq,
                          Tq const &x) {
  size_t n_change = 0;
  int prev = 0;
  for (auto const &p : seq) {
    Tq v = polynomial_evaluation(p, x);
    int s = (v > 0) ? 1 : ((v < 0) ? -1 : 0);
    if (s == 0) {
      continue;
    }
    if (prev != 0 && s != prev) {
      n_change++;
    }
    prev = s;
  }
  return n_change;
}

template <typename Tq>
std::vector<Tq> real_roots_approximation(std::vector<Tq> poly,
                                         Tq const &eps) {
  polynomial_trim(poly);
  std::vector<Tq> roots;
  if (poly.size() <= 1) {
    return roots;
  }
  std::vector<std::vector<Tq>> seq = sturm_sequence(poly);
  Tq B(1);
  {
    Tq lead = poly.back();
    Tq mx(0);
    for (size_t i = 0; i + 1 < poly.size(); i++) {
      Tq v = poly[i] / lead;
      if (v < 0) {
        v = -v;
      }
      if (v > mx) {
        mx = v;
      }
    }
    B += mx;
  }
  // the number of roots in (a, b] is V(a) - V(b)
  std::vector<std::pair<Tq, Tq>> todo{{-B, B}};
  std::vector<std::pair<Tq, Tq>> isolated;
  while (!todo.empty()) {
    auto [a, b] = todo.back();
    todo.pop_back();
    size_t va = sturm_sign_changes(seq, a);
    size_t vb = sturm_sign_changes(seq, b);
    size_t n = va - vb;
    if (n == 0) {
      continue;
    }
    if (n == 1) {
      isolated.push_back({a, b});
      continue;
    }
    Tq mid = a + b;
    mid /= 2;
    todo.push_back({a, mid});
    todo.push_back({mid, b});
  }
  for (auto [a, b] : isolated) {
    // the root is in (a, b]; f(b) = 0 gives it exactly
    if (polynomial_evaluation(poly, b) == 0) {
      roots.push_back(b);
      continue;
    }
    while (true) {
      Tq w = b - a;
      if (!(w > eps)) {
        break;
      }
      Tq mid = a + b;
      mid /= 2;
      Tq fm = polynomial_evaluation(poly, mid);
      if (fm == 0) {
        a = mid;
        b = mid;
        break;
      }
      size_t vm = sturm_sign_changes(seq, mid);
      size_t vb = sturm_sign_changes(seq, b);
      if (vm - vb == 1) {
        a = mid;
      } else {
        b = mid;
      }
    }
    Tq c = a + b;
    c /= 2;
    roots.push_back(c);
  }
  std::sort(roots.begin(), roots.end());
  return roots;
}

template <typename Tq>
Tq square_root_approximation(Tq const &x, Tq const &eps) {
  Tq a(0);
  Tq b = (x > 1) ? x : Tq(1);
  while (true) {
    Tq w = b - a;
    if (!(w > eps)) {
      break;
    }
    Tq mid = a + b;
    mid /= 2;
    Tq sq = mid * mid;
    if (sq > x) {
      b = mid;
    } else {
      a = mid;
    }
  }
  Tq c = a + b;
  c /= 2;
  return c;
}

// clang-format off
#endif  // SRC_NUMBER_POLYNOMIALREALROOTS_H_
// clang-format on
