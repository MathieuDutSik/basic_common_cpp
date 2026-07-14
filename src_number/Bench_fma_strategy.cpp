// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Microbenchmark backing the is_fma_prefered<T> trait. For each numerical type
// it times the three ways of forming a multiply-accumulate  acc (+/-)= a * b
// inside a hot loop:
//
//   fused   : acc += a * b;                     (the natural / direct form)
//   scratch : prod = a * b; acc += prod;        (one reused scratch temporary)
//   inplace : prod = a; prod *= b; acc += prod; (compound ops, no binary *)
//
// and prints is_fma_prefered<T>::value alongside, so the trait can be checked
// against the measured winner: is_fma_prefered should be true exactly when the
// `fused` column is (approximately) the fastest.
//
// Summary of the measurements this reproduces (small operands, -O2):
//   mpz_class      scratch   (gmpxx does not fuse += product)
//   boost mpz_int  fused     (its expression templates DO fuse)
//   boost cpp_int  scratch
//   mpq_class      scratch
//   mpq_rational   scratch
//   cpp_rational   scratch (marginal)
//   Rational<_>    ~neutral  -> scratch
//   QuadField      scratch   (operator=(prod) uses fewer temporaries)
//   RealField      fused     (operator+=(prod) accumulates in place)
//   jet            fused     (operator+=(prod) convolves in place, one pass)
//   native ints    neutral   -> fused
//
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "rational.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryRealField.h"
#include "jet_number.h"
#include "approximation.h"
// clang-format on
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

static auto now() { return std::chrono::steady_clock::now(); }
static double ms(std::chrono::steady_clock::duration d) {
  return std::chrono::duration<double, std::milli>(d).count();
}

// The three FMA code paths, timed back to back on the same data.
template <typename T>
void bench(std::string const &name, std::vector<T> const &a,
           std::vector<T> const &b, long outer) {
  int m = a.size();
  auto t0 = now();
  T acc1(0);
  for (long it = 0; it < outer; it++) {
    T x(0);
    for (int i = 0; i < m; i++)
      x += a[i] * b[i];
    acc1 += x;
  }
  auto t1 = now();
  T acc2(0), prod;
  for (long it = 0; it < outer; it++) {
    T x(0);
    for (int i = 0; i < m; i++) {
      prod = a[i] * b[i];
      x += prod;
    }
    acc2 += x;
  }
  auto t2 = now();
  T acc3(0), p2;
  for (long it = 0; it < outer; it++) {
    T x(0);
    for (int i = 0; i < m; i++) {
      p2 = a[i];
      p2 *= b[i];
      x += p2;
    }
    acc3 += x;
  }
  auto t3 = now();
  double f = ms(t1 - t0), s = ms(t2 - t1), ip = ms(t3 - t2);
  std::string best = (f <= s && f <= ip) ? "fused" : (s <= ip ? "scratch" : "inplace");
  std::cout << name << ": fused=" << f << " scratch=" << s << " inplace=" << ip
            << " ms  best=" << best
            << "  is_fma_prefered=" << (is_fma_prefered<T>::value ? "true" : "false")
            << "  (eq=" << ((acc1 == acc2 && acc2 == acc3) ? 1 : 0) << ")\n";
}

template <typename T> std::vector<T> mk_int(int m) {
  std::vector<T> v(m);
  for (int i = 0; i < m; i++)
    v[i] = T(i * 7 + 3);
  return v;
}

int main() {
  try {
    using Tq = mpq_class;
    int m = 64;
    long outer = 60000;

    std::cout << "--- integer types ---\n";
    bench<mpz_class>("mpz_class    ", mk_int<mpz_class>(m), mk_int<mpz_class>(m), outer);
    bench<boost::multiprecision::mpz_int>("boost mpz_int", mk_int<boost::multiprecision::mpz_int>(m), mk_int<boost::multiprecision::mpz_int>(m), outer);
    bench<boost::multiprecision::cpp_int>("boost cpp_int", mk_int<boost::multiprecision::cpp_int>(m), mk_int<boost::multiprecision::cpp_int>(m), outer);

    std::cout << "--- rational types ---\n";
    std::vector<mpq_class> qa(m), qb(m);
    std::vector<boost::multiprecision::mpq_rational> ba(m), bb(m);
    std::vector<boost::multiprecision::cpp_rational> ca(m), cb(m);
    std::vector<Rational<mpz_class>> rza(m), rzb(m);
    std::vector<Rational<long>> rla(m), rlb(m);
    for (int i = 0; i < m; i++) {
      qa[i] = mpq_class(i * 7 + 3, i + 2); qa[i].canonicalize();
      qb[i] = mpq_class(i * 5 + 11, i + 3); qb[i].canonicalize();
      ba[i] = boost::multiprecision::mpq_rational(i * 7 + 3) / (i + 2);
      bb[i] = boost::multiprecision::mpq_rational(i * 5 + 11) / (i + 3);
      ca[i] = boost::multiprecision::cpp_rational(i * 7 + 3) / (i + 2);
      cb[i] = boost::multiprecision::cpp_rational(i * 5 + 11) / (i + 3);
      rza[i] = Rational<mpz_class>(mpz_class(i * 7 + 3), mpz_class(i + 2));
      rzb[i] = Rational<mpz_class>(mpz_class(i * 5 + 11), mpz_class(i + 3));
      rla[i] = Rational<long>(i * 7 + 3, i + 2);
      rlb[i] = Rational<long>(i * 5 + 11, i + 3);
    }
    bench<mpq_class>("mpq_class    ", qa, qb, outer);
    bench<boost::multiprecision::mpq_rational>("boost mpq_rat", ba, bb, outer);
    bench<boost::multiprecision::cpp_rational>("boost cpp_rat", ca, cb, outer);
    bench<Rational<mpz_class>>("Rational<mpz>", rza, rzb, outer);
    bench<Rational<long>>("Rational<long>", rla, rlb, outer);

    std::cout << "--- algebraic / composite types ---\n";
    std::vector<QuadField<mpq_class, 3>> qfa(m), qfb(m);
    for (int i = 0; i < m; i++) {
      qfa[i] = QuadField<mpq_class, 3>(mpq_class(i * 7 + 3), mpq_class(i + 1));
      qfb[i] = QuadField<mpq_class, 3>(mpq_class(i * 5 + 11), mpq_class(i + 2));
    }
    bench<QuadField<mpq_class, 3>>("QuadField<mpq,3>", qfa, qfb, outer);

    // RealField needs a registered field (2cos(2pi/7), X^3 + X^2 - 2X - 1).
    std::string eFile = "Examples/RealAlgebraicField/CubicFieldDisc_49";
    for (int lev = 0; lev <= 10 && !IsExistingFile(eFile); lev++)
      eFile = "../" + eFile;
    if (IsExistingFile(eFile)) {
      HelperClassRealField<Tq> hcrf(eFile);
      int const idx = 1;
      insert_helper_real_algebraic_field(idx, hcrf);
      using Trf = RealField<idx>;
      std::vector<Trf> rfa(m), rfb(m);
      for (int i = 0; i < m; i++) {
        rfa[i] = Trf(std::vector<Tq>{Tq(i * 7 + 3, i + 2), Tq(i + 1, 3), Tq(2, i + 5)});
        rfb[i] = Trf(std::vector<Tq>{Tq(i * 5 + 11, i + 3), Tq(i + 2, 5), Tq(3, i + 7)});
      }
      bench<Trf>("RealField<cubic>", rfa, rfb, outer);
    } else {
      std::cout << "RealField<cubic>: skipped (field data file not found)\n";
    }

    using Tj = jet<mpq_class, 2>;
    std::vector<Tj> ja(m), jb(m);
    for (int i = 0; i < m; i++) {
      Tj u, v;
      u[0] = Tq(i * 7 + 3, i + 2); u[1] = Tq(i + 1, 3); u[2] = Tq(2, i + 5);
      v[0] = Tq(i * 5 + 11, i + 3); v[1] = Tq(i + 2, 5); v[2] = Tq(3, i + 7);
      ja[i] = u; jb[i] = v;
    }
    bench<Tj>("jet<mpq,2>   ", ja, jb, outer);

    std::cout << "Normal termination of Bench_fma_strategy\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Bench_fma_strategy\n";
    exit(e.eVal);
  }
}
