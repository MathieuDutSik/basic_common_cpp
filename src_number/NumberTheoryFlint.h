// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYFLINT_H_
#define SRC_NUMBER_NUMBERTHEORYFLINT_H_

// The flint arithmetic: fmpz_class / fmpq_class are RAII wrappers around the
// fmpz_t / fmpq_t types of flint (https://flintlib.org), exposing the same
// interface surface as mpz_class / mpq_class of gmpxx so that the generic
// matrix code accepts them unchanged. The point of flint over gmp is the
// small-integer optimization: an fmpz is a single word holding either the
// value itself (when it fits in 62 bits) or a pointer to an mpz, so the
// typical small entries of a matrix computation never touch the allocator.
//
// This header is self-contained and is only included when building with
// ENABLE_FLINT_SUPPORT.

// clang-format off
#include "BasicNumberTypes.h"
#include "NumberTheoryTryInt.h"
#include "ResidueQuotient.h"
#include "Temp_common.h"
#include "TypeConversion.h"
#include "TemplateTraits.h"
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_free.hpp>
#include <limits>
#include <string>
#include <utility>
// clang-format on

class fmpq_class;

// Deferred products (light expression templates): a * b between two
// fmpz_class (resp. fmpq_class) returns a proxy holding the two operands,
// and the consuming operation picks the fused flint call:
//   acc += a * b        -> fmpz_addmul
//   acc -= a * b        -> fmpz_submul
//   x    = a * b        -> fmpz_mul, no temporary
//   x = a*b + c, a*b - c*d, ...  -> mul/addmul/submul chains
// Anything else converts the proxy to the concrete class (one temporary,
// the behavior before the proxies existed), so the generic code and Eigen
// keep compiling unchanged. Same caveat as the gmpxx expression templates:
// the proxy references its operands, so `auto p = a * b;` kept beyond the
// full expression dangles -- write the target type instead of auto.
class fmpz_product;
class fmpq_product;

class fmpz_class {
private:
  fmpz_t a;

public:
  ~fmpz_class() { fmpz_clear(a); }
  fmpz_class() { fmpz_init(a); }
  fmpz_class(fmpz_class const &x) {
    fmpz_init(a);
    fmpz_set(a, x.a);
  }
  // fmpz_init does not allocate, so a moved-from value is a valid zero and
  // the swap hands the possibly heap-backed limbs over without a copy.
  fmpz_class(fmpz_class &&x) noexcept {
    fmpz_init(a);
    fmpz_swap(a, x.a);
  }
  fmpz_class(int const &u) {
    fmpz_init(a);
    fmpz_set_si(a, u);
  }
  fmpz_class(long const &u) {
    fmpz_init(a);
    fmpz_set_si(a, u);
  }
  fmpz_class(std::string const &str) {
    fmpz_init(a);
    if (fmpz_set_str(a, str.c_str(), 10) != 0) {
      std::cerr << "fmpz_class: failed to parse the string " << str << "\n";
      throw TerminalException{1};
    }
  }
  fmpz_class &operator=(fmpz_class const &x) {
    fmpz_set(a, x.a);
    return *this;
  }
  fmpz_class &operator=(fmpz_class &&x) noexcept {
    fmpz_swap(a, x.a);
    return *this;
  }
  fmpz_class &operator=(int const &u) {
    fmpz_set_si(a, u);
    return *this;
  }
  fmpz_class &operator=(long const &u) {
    fmpz_set_si(a, u);
    return *this;
  }
  // Access to the underlying fmpz for the free functions of this header.
  fmpz *get_fmpz_t() { return a; }
  const fmpz *get_fmpz_t() const { return a; }
  double get_d() const { return fmpz_get_d(a); }
  long get_si() const { return fmpz_get_si(a); }
  bool fits_slong_p() const { return fmpz_fits_si(a) != 0; }
  std::string get_str() const {
    char *str = fmpz_get_str(NULL, 10, a);
    std::string eStr = str;
    flint_free(str);
    return eStr;
  }
  //
  // Arithmetic operators
  //
  fmpz_class &operator+=(fmpz_class const &x) {
    fmpz_add(a, a, x.a);
    return *this;
  }
  fmpz_class &operator-=(fmpz_class const &x) {
    fmpz_sub(a, a, x.a);
    return *this;
  }
  fmpz_class &operator*=(fmpz_class const &x) {
    fmpz_mul(a, a, x.a);
    return *this;
  }
  // Division truncates towards zero, matching mpz_class.
  fmpz_class &operator/=(fmpz_class const &x) {
    fmpz_tdiv_q(a, a, x.a);
    return *this;
  }
  fmpz_class &operator++() {
    fmpz_add_si(a, a, 1);
    return *this;
  }
  fmpz_class &operator--() {
    fmpz_sub_si(a, a, 1);
    return *this;
  }
  fmpz_class operator++(int) {
    fmpz_class ret(*this);
    fmpz_add_si(a, a, 1);
    return ret;
  }
  fmpz_class operator--(int) {
    fmpz_class ret(*this);
    fmpz_sub_si(a, a, 1);
    return ret;
  }
  friend fmpz_class operator+(fmpz_class const &x, fmpz_class const &y) {
    fmpz_class z;
    fmpz_add(z.a, x.a, y.a);
    return z;
  }
  friend fmpz_class operator+(fmpz_class const &x, int const &y) {
    fmpz_class z;
    fmpz_add_si(z.a, x.a, y);
    return z;
  }
  friend fmpz_class operator+(int const &x, fmpz_class const &y) {
    fmpz_class z;
    fmpz_add_si(z.a, y.a, x);
    return z;
  }
  friend fmpz_class operator-(fmpz_class const &x, fmpz_class const &y) {
    fmpz_class z;
    fmpz_sub(z.a, x.a, y.a);
    return z;
  }
  friend fmpz_class operator-(fmpz_class const &x, int const &y) {
    fmpz_class z;
    fmpz_sub_si(z.a, x.a, y);
    return z;
  }
  friend fmpz_class operator-(int const &x, fmpz_class const &y) {
    fmpz_class z;
    fmpz_sub_si(z.a, y.a, x);
    fmpz_neg(z.a, z.a);
    return z;
  }
  friend fmpz_class operator-(fmpz_class const &x) {
    fmpz_class z;
    fmpz_neg(z.a, x.a);
    return z;
  }
  friend fmpz_product operator*(fmpz_class const &x, fmpz_class const &y);
  fmpz_class(fmpz_product const &p);
  fmpz_class &operator=(fmpz_product const &p);
  fmpz_class &operator+=(fmpz_product const &p);
  fmpz_class &operator-=(fmpz_product const &p);
  friend fmpz_class operator*(fmpz_class const &x, int const &y) {
    fmpz_class z;
    fmpz_mul_si(z.a, x.a, y);
    return z;
  }
  friend fmpz_class operator*(int const &x, fmpz_class const &y) {
    fmpz_class z;
    fmpz_mul_si(z.a, y.a, x);
    return z;
  }
  friend fmpz_class operator/(fmpz_class const &x, fmpz_class const &y) {
    fmpz_class z;
    fmpz_tdiv_q(z.a, x.a, y.a);
    return z;
  }
  friend fmpz_class operator/(fmpz_class const &x, int const &y) {
    fmpz_class z, b;
    b = y;
    fmpz_tdiv_q(z.a, x.a, b.a);
    return z;
  }
  friend fmpz_class operator%(fmpz_class const &x, fmpz_class const &y) {
    fmpz_class z, q;
    fmpz_tdiv_qr(q.a, z.a, x.a, y.a);
    return z;
  }
  //
  // Comparison operators
  //
  friend bool operator==(fmpz_class const &x, fmpz_class const &y) {
    return fmpz_equal(x.a, y.a) != 0;
  }
  friend bool operator==(fmpz_class const &x, int const &y) {
    return fmpz_cmp_si(x.a, y) == 0;
  }
  friend bool operator==(int const &x, fmpz_class const &y) {
    return fmpz_cmp_si(y.a, x) == 0;
  }
  friend bool operator!=(fmpz_class const &x, fmpz_class const &y) {
    return fmpz_equal(x.a, y.a) == 0;
  }
  friend bool operator!=(fmpz_class const &x, int const &y) {
    return fmpz_cmp_si(x.a, y) != 0;
  }
  friend bool operator!=(int const &x, fmpz_class const &y) {
    return fmpz_cmp_si(y.a, x) != 0;
  }
  friend bool operator<(fmpz_class const &x, fmpz_class const &y) {
    return fmpz_cmp(x.a, y.a) < 0;
  }
  friend bool operator<(fmpz_class const &x, int const &y) {
    return fmpz_cmp_si(x.a, y) < 0;
  }
  friend bool operator<(int const &x, fmpz_class const &y) {
    return fmpz_cmp_si(y.a, x) > 0;
  }
  friend bool operator>(fmpz_class const &x, fmpz_class const &y) {
    return fmpz_cmp(x.a, y.a) > 0;
  }
  friend bool operator>(fmpz_class const &x, int const &y) {
    return fmpz_cmp_si(x.a, y) > 0;
  }
  friend bool operator>(int const &x, fmpz_class const &y) {
    return fmpz_cmp_si(y.a, x) < 0;
  }
  friend bool operator<=(fmpz_class const &x, fmpz_class const &y) {
    return fmpz_cmp(x.a, y.a) <= 0;
  }
  friend bool operator<=(fmpz_class const &x, int const &y) {
    return fmpz_cmp_si(x.a, y) <= 0;
  }
  friend bool operator<=(int const &x, fmpz_class const &y) {
    return fmpz_cmp_si(y.a, x) >= 0;
  }
  friend bool operator>=(fmpz_class const &x, fmpz_class const &y) {
    return fmpz_cmp(x.a, y.a) >= 0;
  }
  friend bool operator>=(fmpz_class const &x, int const &y) {
    return fmpz_cmp_si(x.a, y) >= 0;
  }
  friend bool operator>=(int const &x, fmpz_class const &y) {
    return fmpz_cmp_si(y.a, x) <= 0;
  }
  //
  // Input / output
  //
  friend std::ostream &operator<<(std::ostream &os, fmpz_class const &v) {
    return os << v.get_str();
  }
  // Character-wise parse of [+-]?digits, stopping at the first non-matching
  // character, so that separators such as commas are left in the stream the
  // way the gmpxx extractor does.
  friend std::istream &operator>>(std::istream &is, fmpz_class &v) {
    is >> std::ws;
    std::string s;
    int c = is.peek();
    if (c == '-' || c == '+') {
      s += static_cast<char>(is.get());
      c = is.peek();
    }
    while (c != std::char_traits<char>::eof() && isdigit(c)) {
      s += static_cast<char>(is.get());
      c = is.peek();
    }
    if (s.empty() || fmpz_set_str(v.a, s.c_str(), 10) != 0)
      is.setstate(std::ios::failbit);
    return is;
  }
};

class fmpq_class {
private:
  fmpq_t a;

public:
  ~fmpq_class() { fmpq_clear(a); }
  fmpq_class() { fmpq_init(a); }
  fmpq_class(fmpq_class const &x) {
    fmpq_init(a);
    fmpq_set(a, x.a);
  }
  fmpq_class(fmpq_class &&x) noexcept {
    fmpq_init(a);
    fmpq_swap(a, x.a);
  }
  fmpq_class(int const &u) {
    fmpq_init(a);
    fmpq_set_si(a, u, 1);
  }
  fmpq_class(long const &u) {
    fmpq_init(a);
    fmpq_set_si(a, u, 1);
  }
  fmpq_class(fmpz_class const &u) {
    fmpq_init(a);
    fmpz_set(fmpq_numref(a), u.get_fmpz_t());
  }
  fmpq_class(fmpz_class const &num, fmpz_class const &den) {
    fmpq_init(a);
    fmpz_set(fmpq_numref(a), num.get_fmpz_t());
    fmpz_set(fmpq_denref(a), den.get_fmpz_t());
    fmpq_canonicalise(a);
  }
  fmpq_class &operator=(fmpq_class const &x) {
    fmpq_set(a, x.a);
    return *this;
  }
  fmpq_class &operator=(fmpq_class &&x) noexcept {
    fmpq_swap(a, x.a);
    return *this;
  }
  fmpq_class &operator=(int const &u) {
    fmpq_set_si(a, u, 1);
    return *this;
  }
  fmpq_class &operator=(long const &u) {
    fmpq_set_si(a, u, 1);
    return *this;
  }
  fmpq_class &operator=(fmpz_class const &u) {
    fmpz_set(fmpq_numref(a), u.get_fmpz_t());
    fmpz_one(fmpq_denref(a));
    return *this;
  }
  fmpq *get_fmpq_t() { return a; }
  const fmpq *get_fmpq_t() const { return a; }
  double get_d() const { return fmpq_get_d(a); }
  fmpz_class get_num() const {
    fmpz_class z;
    fmpz_set(z.get_fmpz_t(), fmpq_numref(a));
    return z;
  }
  fmpz_class get_den() const {
    fmpz_class z;
    fmpz_set(z.get_fmpz_t(), fmpq_denref(a));
    return z;
  }
  std::string get_str() const {
    char *str = fmpq_get_str(NULL, 10, a);
    std::string eStr = str;
    flint_free(str);
    return eStr;
  }
  //
  // Arithmetic operators
  //
  fmpq_class &operator+=(fmpq_class const &x) {
    fmpq_add(a, a, x.a);
    return *this;
  }
  fmpq_class &operator-=(fmpq_class const &x) {
    fmpq_sub(a, a, x.a);
    return *this;
  }
  fmpq_class &operator*=(fmpq_class const &x) {
    fmpq_mul(a, a, x.a);
    return *this;
  }
  fmpq_class &operator/=(fmpq_class const &x) {
    fmpq_div(a, a, x.a);
    return *this;
  }
  friend fmpq_class operator+(fmpq_class const &x, fmpq_class const &y) {
    fmpq_class z;
    fmpq_add(z.a, x.a, y.a);
    return z;
  }
  friend fmpq_class operator+(fmpq_class const &x, int const &y) {
    fmpq_class z;
    fmpq_add_si(z.a, x.a, y);
    return z;
  }
  friend fmpq_class operator+(int const &x, fmpq_class const &y) {
    fmpq_class z;
    fmpq_add_si(z.a, y.a, x);
    return z;
  }
  friend fmpq_class operator-(fmpq_class const &x, fmpq_class const &y) {
    fmpq_class z;
    fmpq_sub(z.a, x.a, y.a);
    return z;
  }
  friend fmpq_class operator-(fmpq_class const &x, int const &y) {
    fmpq_class z;
    fmpq_sub_si(z.a, x.a, y);
    return z;
  }
  friend fmpq_class operator-(int const &x, fmpq_class const &y) {
    fmpq_class z;
    fmpq_sub_si(z.a, y.a, x);
    fmpq_neg(z.a, z.a);
    return z;
  }
  friend fmpq_class operator-(fmpq_class const &x) {
    fmpq_class z;
    fmpq_neg(z.a, x.a);
    return z;
  }
  friend fmpq_product operator*(fmpq_class const &x, fmpq_class const &y);
  fmpq_class(fmpq_product const &p);
  fmpq_class &operator=(fmpq_product const &p);
  fmpq_class &operator+=(fmpq_product const &p);
  fmpq_class &operator-=(fmpq_product const &p);
  friend fmpq_class operator*(fmpq_class const &x, int const &y) {
    fmpq_class z, b;
    b = y;
    fmpq_mul(z.a, x.a, b.a);
    return z;
  }
  friend fmpq_class operator*(int const &x, fmpq_class const &y) {
    fmpq_class z, b;
    b = x;
    fmpq_mul(z.a, b.a, y.a);
    return z;
  }
  friend fmpq_class operator/(fmpq_class const &x, fmpq_class const &y) {
    fmpq_class z;
    fmpq_div(z.a, x.a, y.a);
    return z;
  }
  friend fmpq_class operator/(fmpq_class const &x, int const &y) {
    fmpq_class z, b;
    b = y;
    fmpq_div(z.a, x.a, b.a);
    return z;
  }
  friend fmpq_class operator/(int const &x, fmpq_class const &y) {
    fmpq_class z, b;
    b = x;
    fmpq_div(z.a, b.a, y.a);
    return z;
  }
  //
  // Comparison operators
  //
  friend bool operator==(fmpq_class const &x, fmpq_class const &y) {
    return fmpq_equal(x.a, y.a) != 0;
  }
  friend bool operator==(fmpq_class const &x, int const &y) {
    return fmpq_cmp_si(x.a, y) == 0;
  }
  friend bool operator==(int const &x, fmpq_class const &y) {
    return fmpq_cmp_si(y.a, x) == 0;
  }
  friend bool operator!=(fmpq_class const &x, fmpq_class const &y) {
    return fmpq_equal(x.a, y.a) == 0;
  }
  friend bool operator!=(fmpq_class const &x, int const &y) {
    return fmpq_cmp_si(x.a, y) != 0;
  }
  friend bool operator!=(int const &x, fmpq_class const &y) {
    return fmpq_cmp_si(y.a, x) != 0;
  }
  friend bool operator<(fmpq_class const &x, fmpq_class const &y) {
    return fmpq_cmp(x.a, y.a) < 0;
  }
  friend bool operator<(fmpq_class const &x, int const &y) {
    return fmpq_cmp_si(x.a, y) < 0;
  }
  friend bool operator<(int const &x, fmpq_class const &y) {
    return fmpq_cmp_si(y.a, x) > 0;
  }
  friend bool operator>(fmpq_class const &x, fmpq_class const &y) {
    return fmpq_cmp(x.a, y.a) > 0;
  }
  friend bool operator>(fmpq_class const &x, int const &y) {
    return fmpq_cmp_si(x.a, y) > 0;
  }
  friend bool operator>(int const &x, fmpq_class const &y) {
    return fmpq_cmp_si(y.a, x) < 0;
  }
  friend bool operator<=(fmpq_class const &x, fmpq_class const &y) {
    return fmpq_cmp(x.a, y.a) <= 0;
  }
  friend bool operator<=(fmpq_class const &x, int const &y) {
    return fmpq_cmp_si(x.a, y) <= 0;
  }
  friend bool operator<=(int const &x, fmpq_class const &y) {
    return fmpq_cmp_si(y.a, x) >= 0;
  }
  friend bool operator>=(fmpq_class const &x, fmpq_class const &y) {
    return fmpq_cmp(x.a, y.a) >= 0;
  }
  friend bool operator>=(fmpq_class const &x, int const &y) {
    return fmpq_cmp_si(x.a, y) >= 0;
  }
  friend bool operator>=(int const &x, fmpq_class const &y) {
    return fmpq_cmp_si(y.a, x) <= 0;
  }
  //
  // Input / output
  //
  friend std::ostream &operator<<(std::ostream &os, fmpq_class const &v) {
    return os << v.get_str();
  }
  // Parses [+-]?digits(/digits)? the way the gmpxx extractor does.
  friend std::istream &operator>>(std::istream &is, fmpq_class &v) {
    is >> std::ws;
    std::string s;
    auto read_integer = [&]() {
      int c = is.peek();
      if (c == '-' || c == '+') {
        s += static_cast<char>(is.get());
        c = is.peek();
      }
      while (c != std::char_traits<char>::eof() && isdigit(c)) {
        s += static_cast<char>(is.get());
        c = is.peek();
      }
    };
    read_integer();
    if (is.peek() == '/') {
      s += static_cast<char>(is.get());
      read_integer();
    }
    if (s.empty() || fmpq_set_str(v.a, s.c_str(), 10) != 0) {
      is.setstate(std::ios::failbit);
      return is;
    }
    fmpq_canonicalise(v.a);
    return is;
  }
};

// The deferred product proxies. They hold plain pointers to the operands:
// valid within the full expression that created them, which is the only
// place they are supposed to live.

class fmpz_product {
  friend class fmpz_class;
  friend fmpz_product operator*(fmpz_class const &x, fmpz_class const &y);
  const fmpz *x;
  const fmpz *y;
  fmpz_product(const fmpz *x, const fmpz *y) : x(x), y(y) {}
  // The conversion to fmpz_class goes through the converting constructor
  // of fmpz_class; a conversion operator here as well would make the
  // conversion ambiguous.
};

inline fmpz_product operator*(fmpz_class const &x, fmpz_class const &y) {
  return {x.get_fmpz_t(), y.get_fmpz_t()};
}

inline fmpz_class::fmpz_class(fmpz_product const &p) {
  fmpz_init(a);
  fmpz_mul(a, p.x, p.y);
}

inline fmpz_class &fmpz_class::operator=(fmpz_product const &p) {
  fmpz_mul(a, p.x, p.y);
  return *this;
}

inline fmpz_class &fmpz_class::operator+=(fmpz_product const &p) {
  // acc += acc * y needs the product materialized first.
  if (p.x == a || p.y == a) {
    fmpz_class tmp(p);
    fmpz_add(a, a, tmp.get_fmpz_t());
  } else {
    fmpz_addmul(a, p.x, p.y);
  }
  return *this;
}

inline fmpz_class &fmpz_class::operator-=(fmpz_product const &p) {
  if (p.x == a || p.y == a) {
    fmpz_class tmp(p);
    fmpz_sub(a, a, tmp.get_fmpz_t());
  } else {
    fmpz_submul(a, p.x, p.y);
  }
  return *this;
}

// The sums and differences involving a product go through the fused calls.
// The overloads taking int exist to keep expressions like a*b + 1
// unambiguous (both operator+(T, int) and operator+(product, T) would
// otherwise be equally good).

inline fmpz_class operator+(fmpz_class const &c, fmpz_product const &p) {
  fmpz_class z(c);
  z += p;
  return z;
}

inline fmpz_class operator+(fmpz_product const &p, fmpz_class const &c) {
  fmpz_class z(c);
  z += p;
  return z;
}

inline fmpz_class operator+(fmpz_product const &p, fmpz_product const &q) {
  fmpz_class z(p);
  z += q;
  return z;
}

inline fmpz_class operator+(fmpz_product const &p, int const &c) {
  fmpz_class z(p);
  fmpz_add_si(z.get_fmpz_t(), z.get_fmpz_t(), c);
  return z;
}

inline fmpz_class operator+(int const &c, fmpz_product const &p) {
  return p + c;
}

inline fmpz_class operator-(fmpz_class const &c, fmpz_product const &p) {
  fmpz_class z(c);
  z -= p;
  return z;
}

inline fmpz_class operator-(fmpz_product const &p, fmpz_class const &c) {
  fmpz_class z(p);
  z -= c;
  return z;
}

inline fmpz_class operator-(fmpz_product const &p, fmpz_product const &q) {
  fmpz_class z(p);
  z -= q;
  return z;
}

inline fmpz_class operator-(fmpz_product const &p, int const &c) {
  fmpz_class z(p);
  fmpz_sub_si(z.get_fmpz_t(), z.get_fmpz_t(), c);
  return z;
}

inline fmpz_class operator-(int const &c, fmpz_product const &p) {
  fmpz_class z(p);
  fmpz_sub_si(z.get_fmpz_t(), z.get_fmpz_t(), c);
  fmpz_neg(z.get_fmpz_t(), z.get_fmpz_t());
  return z;
}

inline fmpz_class operator-(fmpz_product const &p) {
  fmpz_class z(p);
  fmpz_neg(z.get_fmpz_t(), z.get_fmpz_t());
  return z;
}

// Comparisons between two products (the cross-multiplication pattern
// a*d < c*b): the hidden friends of fmpz_class are not found by ADL when
// both arguments are proxies, so these exist explicitly.

inline bool operator==(fmpz_product const &p, fmpz_product const &q) {
  return fmpz_class(p) == fmpz_class(q);
}
inline bool operator!=(fmpz_product const &p, fmpz_product const &q) {
  return fmpz_class(p) != fmpz_class(q);
}
inline bool operator<(fmpz_product const &p, fmpz_product const &q) {
  return fmpz_class(p) < fmpz_class(q);
}
inline bool operator<=(fmpz_product const &p, fmpz_product const &q) {
  return fmpz_class(p) <= fmpz_class(q);
}
inline bool operator>(fmpz_product const &p, fmpz_product const &q) {
  return fmpz_class(p) > fmpz_class(q);
}
inline bool operator>=(fmpz_product const &p, fmpz_product const &q) {
  return fmpz_class(p) >= fmpz_class(q);
}

class fmpq_product {
  friend class fmpq_class;
  friend fmpq_product operator*(fmpq_class const &x, fmpq_class const &y);
  const fmpq *x;
  const fmpq *y;
  fmpq_product(const fmpq *x, const fmpq *y) : x(x), y(y) {}
};

inline fmpq_product operator*(fmpq_class const &x, fmpq_class const &y) {
  return {x.get_fmpq_t(), y.get_fmpq_t()};
}

inline fmpq_class::fmpq_class(fmpq_product const &p) {
  fmpq_init(a);
  fmpq_mul(a, p.x, p.y);
}

inline fmpq_class &fmpq_class::operator=(fmpq_product const &p) {
  fmpq_mul(a, p.x, p.y);
  return *this;
}

inline fmpq_class &fmpq_class::operator+=(fmpq_product const &p) {
  if (p.x == a || p.y == a) {
    fmpq_class tmp(p);
    fmpq_add(a, a, tmp.get_fmpq_t());
  } else {
    fmpq_addmul(a, p.x, p.y);
  }
  return *this;
}

inline fmpq_class &fmpq_class::operator-=(fmpq_product const &p) {
  if (p.x == a || p.y == a) {
    fmpq_class tmp(p);
    fmpq_sub(a, a, tmp.get_fmpq_t());
  } else {
    fmpq_submul(a, p.x, p.y);
  }
  return *this;
}

inline fmpq_class operator+(fmpq_class const &c, fmpq_product const &p) {
  fmpq_class z(c);
  z += p;
  return z;
}

inline fmpq_class operator+(fmpq_product const &p, fmpq_class const &c) {
  fmpq_class z(c);
  z += p;
  return z;
}

inline fmpq_class operator+(fmpq_product const &p, fmpq_product const &q) {
  fmpq_class z(p);
  z += q;
  return z;
}

inline fmpq_class operator+(fmpq_product const &p, int const &c) {
  fmpq_class z(p);
  fmpq_add_si(z.get_fmpq_t(), z.get_fmpq_t(), c);
  return z;
}

inline fmpq_class operator+(int const &c, fmpq_product const &p) {
  return p + c;
}

inline fmpq_class operator-(fmpq_class const &c, fmpq_product const &p) {
  fmpq_class z(c);
  z -= p;
  return z;
}

inline fmpq_class operator-(fmpq_product const &p, fmpq_class const &c) {
  fmpq_class z(p);
  z -= c;
  return z;
}

inline fmpq_class operator-(fmpq_product const &p, fmpq_product const &q) {
  fmpq_class z(p);
  z -= q;
  return z;
}

inline fmpq_class operator-(fmpq_product const &p, int const &c) {
  fmpq_class z(p);
  fmpq_sub_si(z.get_fmpq_t(), z.get_fmpq_t(), c);
  return z;
}

inline fmpq_class operator-(int const &c, fmpq_product const &p) {
  fmpq_class z(p);
  fmpq_sub_si(z.get_fmpq_t(), z.get_fmpq_t(), c);
  fmpq_neg(z.get_fmpq_t(), z.get_fmpq_t());
  return z;
}

inline fmpq_class operator-(fmpq_product const &p) {
  fmpq_class z(p);
  fmpq_neg(z.get_fmpq_t(), z.get_fmpq_t());
  return z;
}

inline bool operator==(fmpq_product const &p, fmpq_product const &q) {
  return fmpq_class(p) == fmpq_class(q);
}
inline bool operator!=(fmpq_product const &p, fmpq_product const &q) {
  return fmpq_class(p) != fmpq_class(q);
}
inline bool operator<(fmpq_product const &p, fmpq_product const &q) {
  return fmpq_class(p) < fmpq_class(q);
}
inline bool operator<=(fmpq_product const &p, fmpq_product const &q) {
  return fmpq_class(p) <= fmpq_class(q);
}
inline bool operator>(fmpq_product const &p, fmpq_product const &q) {
  return fmpq_class(p) > fmpq_class(q);
}
inline bool operator>=(fmpq_product const &p, fmpq_product const &q) {
  return fmpq_class(p) >= fmpq_class(q);
}

// The sgn / abs free functions that gmpxx provides for its types.

inline int sgn(fmpz_class const &x) { return fmpz_sgn(x.get_fmpz_t()); }

inline int sgn(fmpq_class const &x) { return fmpq_sgn(x.get_fmpq_t()); }

inline fmpz_class abs(fmpz_class const &x) {
  fmpz_class z;
  fmpz_abs(z.get_fmpz_t(), x.get_fmpz_t());
  return z;
}

inline fmpq_class abs(fmpq_class const &x) {
  fmpq_class z;
  fmpq_abs(z.get_fmpq_t(), x.get_fmpq_t());
  return z;
}

// The fast conversion into a try-type, see NumberTheoryGmp.h for the design.

template <typename Ttry = TryInt64>
inline Ttry ConvertToTryInt64(fmpz_class const &val,
                              [[maybe_unused]] fmpz_class &scratch) {
  if (!fmpz_fits_si(val.get_fmpz_t()))
    throw TryIntException{1};
  return Ttry(static_cast<int64_t>(fmpz_get_si(val.get_fmpz_t())));
}

// get_bit

inline size_t get_bit(fmpz_class const &v) {
  return fmpz_sizeinbase(v.get_fmpz_t(), 2);
}

// Traits: same choices as for mpz_class / mpq_class.

template <> struct is_fmpz_class<fmpz_class> {
  static const bool value = true;
};

template <> struct is_fmpq_class<fmpq_class> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<fmpz_class> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<fmpq_class> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<fmpz_class> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<fmpq_class> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<fmpz_class> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<fmpq_class> {
  static const bool value = true;
};

template <> struct is_ring_field<fmpz_class> {
  static const bool value = false;
};

template <> struct is_ring_field<fmpq_class> {
  static const bool value = true;
};

template <> struct use_bareiss_for_determinants<fmpq_class> {
  static const bool value = true;
};

template <> struct use_hnf_mod_D<fmpz_class> {
  static const bool value = true;
};

template <> struct is_totally_ordered<fmpz_class> {
  static const bool value = true;
};

template <> struct is_totally_ordered<fmpq_class> {
  static const bool value = true;
};

template <> struct is_exact_arithmetic<fmpz_class> {
  static const bool value = true;
};

template <> struct is_exact_arithmetic<fmpq_class> {
  static const bool value = true;
};

// The compound acc += a * b materializes a temporary for the product; flint
// has native fused calls, used through the AddMul / SubMul specializations.
template <> struct is_fma_prefered<fmpz_class> {
  static const bool value = false;
};

template <> struct is_fma_prefered<fmpq_class> {
  static const bool value = false;
};

template <>
inline void AddMul(fmpz_class &acc, fmpz_class const &a, fmpz_class const &b) {
  fmpz_addmul(acc.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t());
}

template <>
inline void SubMul(fmpz_class &acc, fmpz_class const &a, fmpz_class const &b) {
  fmpz_submul(acc.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t());
}

template <>
inline void AddMul(fmpq_class &acc, fmpq_class const &a, fmpq_class const &b) {
  fmpq_addmul(acc.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t());
}

template <>
inline void SubMul(fmpq_class &acc, fmpq_class const &a, fmpq_class const &b) {
  fmpq_submul(acc.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t());
}

template <> struct underlying_ring<fmpz_class> {
  typedef fmpz_class ring_type;
};

template <> struct underlying_ring<fmpq_class> {
  typedef fmpz_class ring_type;
};

template <> struct overlying_field<fmpz_class> {
  typedef fmpq_class field_type;
};

template <> struct overlying_field<fmpq_class> {
  typedef fmpq_class field_type;
};

template <> struct underlying_totally_ordered_ring<fmpz_class> {
  typedef fmpz_class real_type;
};

template <> struct underlying_totally_ordered_ring<fmpq_class> {
  typedef fmpq_class real_type;
};

// hash functionality

namespace std {
template <> struct hash<fmpz_class> {
  std::size_t operator()(const fmpz_class &val) const {
    const fmpz *a = val.get_fmpz_t();
    if (fmpz_fits_si(a))
      return static_cast<std::size_t>(fmpz_get_si(a));
    // Largest prime below 2^64 keeps the big values spread out.
    return fmpz_fdiv_ui(a, UWORD(18446744073709551557));
  }
};
template <> struct hash<fmpq_class> {
  std::size_t operator()(const fmpq_class &val) const {
    size_t hash1 = std::hash<fmpz_class>()(val.get_den());
    size_t hash2 = std::hash<fmpz_class>()(val.get_num());
    return hash1 + (hash2 << 6) + (hash2 >> 2);
  }
};
// clang-format off
}  // namespace std
// clang-format on

// std::format support, through the operator<<

template <>
struct std::formatter<fmpz_class> : ostream_formatter<fmpz_class> {};

template <>
struct std::formatter<fmpq_class> : ostream_formatter<fmpq_class> {};

// Quotient / remainder, same conventions as the mpz_class / mpq_class code in
// NumberTheoryGmp.h (transliterated from it).

inline void ResInt_Kernel(fmpz_class const &a, fmpz_class const &b,
                          fmpz_class &res) {
  fmpz_class q;
  fmpz_cdiv_qr(q.get_fmpz_t(), res.get_fmpz_t(), a.get_fmpz_t(),
               b.get_fmpz_t());
  if (b > 0 && res != 0) {
    if (b < 0)
      res -= b;
    else
      res += b;
  }
}

inline void ResInt_Kernel(fmpq_class const &a, fmpq_class const &b,
                          fmpq_class &res) {
  fmpz_class a_den = a.get_den();
  fmpz_class b_den = b.get_den();
  fmpz_class eGcd;
  fmpz_gcd(eGcd.get_fmpz_t(), a_den.get_fmpz_t(), b_den.get_fmpz_t());
  fmpz_class eLCM_z = a_den * b_den / eGcd;
  fmpq_class eLCM(eLCM_z);
  fmpq_class aProd = a * eLCM;
  fmpq_class bProd = b * eLCM;
  fmpz_class a_num = aProd.get_num();
  fmpz_class b_num = bProd.get_num();
  fmpz_class b_num_pos = abs(b_num);
  fmpz_class res_z;
  fmpz_mod(res_z.get_fmpz_t(), a_num.get_fmpz_t(), b_num_pos.get_fmpz_t());
  res = fmpq_class(res_z) / eLCM;
}

inline void QUO_INT(stc<fmpz_class> const &a, stc<fmpz_class> const &b,
                    fmpz_class &q) {
  fmpz_cdiv_q(q.get_fmpz_t(), a.val.get_fmpz_t(), b.val.get_fmpz_t());
  if (b.val > 0 && b.val * q != a.val) {
    if (b.val > 0)
      --q;
    else
      ++q;
  }
}

inline void QUO_INT(stc<fmpq_class> const &a, stc<fmpq_class> const &b,
                    fmpq_class &q) {
  fmpq_class res = ResInt(a.val, b.val);
  q = (a.val - res) / b.val;
}

#include "QuoIntFcts.h"

inline fmpz_class CanonicalizationUnit(fmpz_class const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline fmpq_class CanonicalizationUnit(fmpq_class const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline fmpz_class T_NormGen(fmpz_class const &x) { return abs(x); }

inline fmpq_class T_NormGen(fmpq_class const &x) { return abs(x); }

// Pivot cost: the number of machine words, as for gmp (see PivotCost.h).

inline size_t f_cost_pivot(fmpz_class const &x) {
  return fmpz_size(x.get_fmpz_t());
}

inline size_t f_cost_pivot(fmpq_class const &x) {
  return fmpz_size(fmpq_numref(x.get_fmpq_t())) +
         fmpz_size(fmpq_denref(x.get_fmpq_t()));
}

inline bool IsInteger(fmpq_class const &x) {
  return fmpz_is_one(fmpq_denref(x.get_fmpq_t()));
}

//

inline fmpq_class GetDenominator(fmpq_class const &x) {
  return fmpq_class(x.get_den());
}

inline fmpz_class GetDenominator([[maybe_unused]] fmpz_class const &x) {
  return 1;
}

inline fmpq_class GetNumerator(fmpq_class const &x) {
  return fmpq_class(x.get_num());
}

inline fmpz_class GetNumerator(fmpz_class const &x) { return x; }

inline fmpz_class GetDenominator_z(fmpq_class const &x) { return x.get_den(); }

inline fmpz_class GetDenominator_z([[maybe_unused]] fmpz_class const &x) {
  return 1;
}

inline fmpz_class GetNumerator_z(fmpq_class const &x) { return x.get_num(); }

inline fmpz_class GetNumerator_z(fmpz_class const &x) { return x; }

//

inline void ScalingInteger_Kernel(stc<fmpq_class> const &x, fmpz_class &x_ret) {
  x_ret = x.val.get_den();
}

inline void ScalingInteger_Kernel([[maybe_unused]] stc<fmpz_class> const &x,
                                  fmpz_class &x_ret) {
  x_ret = 1;
}

// gcd / lcm / extended gcd

inline fmpz_class KernelGcdPair(fmpz_class const &a, fmpz_class const &b) {
  fmpz_class eGCD;
  fmpz_gcd(eGCD.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t());
  return eGCD;
}

inline PairGCD_dot<fmpz_class> ComputePairGcdDot(fmpz_class const &m,
                                                 fmpz_class const &n) {
  fmpz_class eGCD;
  if (n == 0 && m == 0) {
    eGCD = 0;
    return {0, 0, eGCD};
  }
  fmpz_class s, t;
  // The canonical Bezout coefficients satisfy the same size bounds as the
  // ones of mpz_gcdext, which the HNF code relies on to control growth.
  fmpz_xgcd_canonical_bezout(eGCD.get_fmpz_t(), s.get_fmpz_t(), t.get_fmpz_t(),
                             m.get_fmpz_t(), n.get_fmpz_t());
  return {s, t, eGCD};
}

inline fmpz_class KernelLCMpair(fmpz_class const &a, fmpz_class const &b) {
  fmpz_class eLCM;
  fmpz_lcm(eLCM.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t());
  return eLCM;
}

// TYPE_CONVERSION: fmpq_class as input

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, fmpq_class &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, double &a2) {
  a2 = a1.val.get_d();
}

inline void Termination_fmpq_not_integer(stc<fmpq_class> const &a1) {
  if (!IsInteger(a1.val)) {
    std::string str =
        "a1=" + a1.val.get_str() + " is not an integer";
    throw ConversionException{str};
  }
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, fmpz_class &a2) {
  Termination_fmpq_not_integer(a1);
  a2 = a1.val.get_num();
}

// fmpz_class as source for the bounded integer types, with the same overflow
// checking contract as mpz_class_to_small_integer. slong is 64-bit on the
// supported platforms, so the fits_si test covers every narrower target.
template <typename Tout>
inline void fmpz_class_to_small_integer(fmpz_class const &val, Tout &out) {
  static_assert(std::is_integral_v<Tout>,
                "fmpz_class_to_small_integer is for the integral types");
  auto fail = [&]() {
    std::string str = "value=" + val.get_str() +
                      " does not fit in the destination integer type";
    throw ConversionException{str};
  };
  if constexpr (std::is_signed_v<Tout>) {
    if (!fmpz_fits_si(val.get_fmpz_t()))
      fail();
    slong e_val = fmpz_get_si(val.get_fmpz_t());
    if (e_val < static_cast<slong>(std::numeric_limits<Tout>::min()) ||
        e_val > static_cast<slong>(std::numeric_limits<Tout>::max()))
      fail();
    out = static_cast<Tout>(e_val);
  } else {
    if (fmpz_sgn(val.get_fmpz_t()) < 0)
      fail();
    if (fmpz_cmp_ui(val.get_fmpz_t(),
                    static_cast<ulong>(std::numeric_limits<Tout>::max())) > 0)
      fail();
    out = static_cast<Tout>(fmpz_get_ui(val.get_fmpz_t()));
  }
}

template <typename Tout>
inline void fmpq_class_to_small_integer(stc<fmpq_class> const &a1, Tout &out) {
  Termination_fmpq_not_integer(a1);
  fmpz_class a1_z = a1.val.get_num();
  fmpz_class_to_small_integer(a1_z, out);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, int8_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, uint8_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, int16_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, uint16_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, int32_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, uint32_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, int64_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, uint64_t &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

template <typename T>
  requires (std::is_same_v<T, long>
            && !std::is_same_v<long, int64_t>
            && !std::is_same_v<long, int32_t>)
inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, T &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, T &a2) {
  fmpq_class_to_small_integer(a1, a2);
}

// fmpz_class as input

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, fmpz_class &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, fmpq_class &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, double &a2) {
  a2 = a1.val.get_d();
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, int8_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, uint8_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, int16_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, uint16_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, int32_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, uint32_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, int64_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, uint64_t &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, T &a2) {
  fmpz_class_to_small_integer(a1.val, a2);
}

// The small integer types as input

inline void TYPE_CONVERSION(stc<int8_t> const &a1, fmpz_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, fmpq_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, fmpz_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, fmpq_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, fmpz_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, fmpq_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, fmpz_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, fmpq_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, fmpz_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, fmpq_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, fmpz_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, fmpq_class &a2) {
  a2 = static_cast<long>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, fmpz_class &a2) {
  fmpz_set_si(a2.get_fmpz_t(), static_cast<slong>(a1.val));
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, fmpq_class &a2) {
  fmpq_set_si(a2.get_fmpq_t(), static_cast<slong>(a1.val), 1);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, fmpz_class &a2) {
  fmpz_set_ui(a2.get_fmpz_t(), static_cast<ulong>(a1.val));
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, fmpq_class &a2) {
  fmpz_set_ui(fmpq_numref(a2.get_fmpq_t()), static_cast<ulong>(a1.val));
  fmpz_one(fmpq_denref(a2.get_fmpq_t()));
}

template <typename T>
  requires (std::is_same_v<T, long>
            && !std::is_same_v<long, int64_t>
            && !std::is_same_v<long, int32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, fmpz_class &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, long>
            && !std::is_same_v<long, int64_t>
            && !std::is_same_v<long, int32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, fmpq_class &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, fmpz_class &a2) {
  fmpz_set_ui(a2.get_fmpz_t(), static_cast<ulong>(a1.val));
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, fmpq_class &a2) {
  fmpz_set_ui(fmpq_numref(a2.get_fmpq_t()), static_cast<ulong>(a1.val));
  fmpz_one(fmpq_denref(a2.get_fmpq_t()));
}

// double as input, same limited contract as for mpz_class

inline void TYPE_CONVERSION(stc<double> const &a1, fmpz_class &a2) {
  int64_t a1_b = static_cast<int64_t>(a1.val);
  fmpz_set_si(a2.get_fmpz_t(), static_cast<slong>(a1_b));
}

#ifdef SRC_NUMBER_NUMBERTHEORYGMP_H_
// Conversions between the flint and gmp types, available when both headers
// are included (NumberTheory.h always comes first in the .cpp files).

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, mpz_class &a2) {
  fmpz_get_mpz(a2.get_mpz_t(), a1.val.get_fmpz_t());
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, fmpz_class &a2) {
  fmpz_set_mpz(a2.get_fmpz_t(), a1.val.get_mpz_t());
}

inline void TYPE_CONVERSION(stc<fmpq_class> const &a1, mpq_class &a2) {
  fmpq_get_mpq(a2.get_mpq_t(), a1.val.get_fmpq_t());
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, fmpq_class &a2) {
  fmpq_set_mpq(a2.get_fmpq_t(), a1.val.get_mpq_t());
}

inline void TYPE_CONVERSION(stc<fmpz_class> const &a1, mpq_class &a2) {
  mpz_class a2_z;
  fmpz_get_mpz(a2_z.get_mpz_t(), a1.val.get_fmpz_t());
  a2 = a2_z;
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, fmpq_class &a2) {
  fmpz_class a2_z;
  fmpz_set_mpz(a2_z.get_fmpz_t(), a1.val.get_mpz_t());
  a2 = a2_z;
}
#endif

// square root

inline bool universal_square_root(fmpz_class &ret, fmpz_class const &val) {
  if (sgn(val) < 0)
    return false;
  fmpz_sqrt(ret.get_fmpz_t(), val.get_fmpz_t());
  return ret * ret == val;
}

inline bool universal_square_root(fmpq_class &ret, fmpq_class const &val) {
  fmpz_class ret_num, ret_den;
  if (!universal_square_root(ret_num, val.get_num()))
    return false;
  if (!universal_square_root(ret_den, val.get_den()))
    return false;
  ret = fmpq_class(ret_num, ret_den);
  return true;
}

inline void set_to_infinity(fmpz_class &x) {
  fmpz_set_ui(x.get_fmpz_t(),
              static_cast<ulong>(std::numeric_limits<uint64_t>::max()));
}

inline void set_to_infinity(fmpq_class &x) {
  fmpz_set_ui(fmpq_numref(x.get_fmpq_t()),
              static_cast<ulong>(std::numeric_limits<uint64_t>::max()));
  fmpz_one(fmpq_denref(x.get_fmpq_t()));
}

//
// Nearest integer and similar stuff.
//

inline fmpq_class FractionalPart(fmpq_class const &x) {
  fmpz_class eNum = x.get_num();
  fmpz_class eDen = x.get_den();
  fmpz_class res;
  fmpz_mod(res.get_fmpz_t(), eNum.get_fmpz_t(), eDen.get_fmpz_t());
  return fmpq_class(res, eDen);
}

inline fmpq_class Floor_fmpq(fmpq_class const &x) {
  fmpz_class q;
  fmpz_fdiv_q(q.get_fmpz_t(), fmpq_numref(x.get_fmpq_t()),
              fmpq_denref(x.get_fmpq_t()));
  return fmpq_class(q);
}

inline fmpq_class Ceil_fmpq(fmpq_class const &x) {
  fmpz_class q;
  fmpz_cdiv_q(q.get_fmpz_t(), fmpq_numref(x.get_fmpq_t()),
              fmpq_denref(x.get_fmpq_t()));
  return fmpq_class(q);
}

inline void FloorInteger(fmpq_class const &xI, fmpq_class &xO) {
  xO = Floor_fmpq(xI);
}

inline void FloorInteger(fmpq_class const &xI, fmpz_class &xO) {
  fmpz_fdiv_q(xO.get_fmpz_t(), fmpq_numref(xI.get_fmpq_t()),
              fmpq_denref(xI.get_fmpq_t()));
}

inline void FloorInteger(fmpq_class const &xI, int &xO) {
  fmpz_class xO_z;
  FloorInteger(xI, xO_z);
  xO = static_cast<int>(xO_z.get_si());
}

inline void FloorInteger(fmpq_class const &xI, long &xO) {
  fmpz_class xO_z;
  FloorInteger(xI, xO_z);
  xO = xO_z.get_si();
}

inline void CeilInteger(fmpq_class const &xI, fmpq_class &xO) {
  xO = Ceil_fmpq(xI);
}

inline void CeilInteger(fmpq_class const &xI, fmpz_class &xO) {
  fmpz_cdiv_q(xO.get_fmpz_t(), fmpq_numref(xI.get_fmpq_t()),
              fmpq_denref(xI.get_fmpq_t()));
}

inline void CeilInteger(fmpq_class const &xI, int &xO) {
  fmpz_class xO_z;
  CeilInteger(xI, xO_z);
  xO = static_cast<int>(xO_z.get_si());
}

inline void CeilInteger(fmpq_class const &xI, long &xO) {
  fmpz_class xO_z;
  CeilInteger(xI, xO_z);
  xO = xO_z.get_si();
}

// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y.
inline fmpq_class NearestInteger_rni(fmpq_class const &x) {
  fmpq_class eFrac = FractionalPart(x);
  fmpq_class eDiff1 = eFrac;
  fmpq_class eDiff2 = 1 - eFrac;
  fmpq_class RetVal = x - eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return RetVal + 1;
  }
}

inline void NearestInteger(fmpq_class const &xI, fmpq_class &xO) {
  xO = NearestInteger_rni(xI);
}

inline void NearestInteger(fmpq_class const &xI, fmpz_class &xO) {
  fmpq_class xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_num();
}

inline void NearestInteger(int const &xI, fmpq_class &xO) { xO = xI; }

inline void NearestInteger(long const &xI, fmpq_class &xO) { xO = xI; }

inline void NearestInteger(fmpq_class const &xI, int &xO) {
  fmpq_class xO_q = NearestInteger_rni(xI);
  xO = static_cast<int>(xO_q.get_num().get_si());
}

inline void NearestInteger(fmpq_class const &xI, long &xO) {
  fmpq_class xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_num().get_si();
}

namespace boost::serialization {

// fmpq_class

template <class Archive>
inline void load(Archive &ar, fmpq_class &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("fmpq", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, fmpq_class const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str = val.get_str();
  ar &make_nvp("fmpq", str);
}

template <class Archive>
inline void serialize(Archive &ar, fmpq_class &val,
                      [[maybe_unused]] const unsigned int version) {
  split_free(ar, val, version);
}

// fmpz_class

template <class Archive>
inline void load(Archive &ar, fmpz_class &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("fmpz", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, fmpz_class const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str = val.get_str();
  ar &make_nvp("fmpz", str);
}

template <class Archive>
inline void serialize(Archive &ar, fmpz_class &val,
                      const unsigned int version) {
  split_free(ar, val, version);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYFLINT_H_
// clang-format on
