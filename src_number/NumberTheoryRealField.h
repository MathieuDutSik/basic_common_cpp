// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYREALFIELD_H_
#define SRC_NUMBER_NUMBERTHEORYREALFIELD_H_

// clang-format off
#include "MAT_MatrixSolutionMat.h"
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "Temp_common.h"
#include "InputOutput.h"
#include <boost/container/small_vector.hpp>
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

// For general real field.
// We need to use more sophisticated and slower algorithms than quadratic
// fields. (A) linear algebra is needed. (B) analysis is needed for deciding
// signs.
//
// An element is stored as an integer polynomial together with a single
// denominator:
//     elt = (num[0] + num[1] y + ... + num[deg-1] y^{deg-1}) / den
// where y = scal * x is a rescaling of the field generator x chosen so that
// the minimal polynomial of y is monic with integer coefficients. Products
// and reductions then stay in integer arithmetic and the gcd normalization
// is done once per operation instead of once per rational coefficient
// operation, which removes most of the gcd cost that dominates the
// vector-of-rationals representation.
// The external interface (input files, string I/O, serialization, the
// vector<T> constructor) remains expressed over the powers of x; the change
// of basis is diagonal since y^i = scal^i x^i.

#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
double threshold_real_alg_check = 0.0001;
#endif

#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
using Trat_real_field = boost::multiprecision::mpq_rational;
using Tint_real_field = boost::multiprecision::mpz_int;
#else
using Trat_real_field = mpq_class;
using Tint_real_field = mpz_class;
#endif

// Inline storage capacity of the numerator coefficients: elements of fields
// of degree up to REALFIELD_INLINE_DEG are stored without heap allocation
// while larger degrees transparently fall back to heap storage, exactly as
// with std::vector. This is a performance knob, not a limit on the degree;
// computations on fields of low degree can tighten it, for example with
// -DREALFIELD_INLINE_DEG=4.
#ifndef REALFIELD_INLINE_DEG
#define REALFIELD_INLINE_DEG 8
#endif
using Tvec_real_field =
    boost::container::small_vector<Tint_real_field, REALFIELD_INLINE_DEG>;
using Tconv_real_field =
    boost::container::small_vector<Tint_real_field,
                                   2 * REALFIELD_INLINE_DEG - 1>;

template <typename Tfield> struct HelperClassRealField {
private:
  using T = Tfield;
  using Tz = Tint_real_field;
  using Tvec = Tvec_real_field;
  using Tconv = Tconv_real_field;
  // One level of the approximant ladder: the powers of the lower / upper
  // bound of y as integers over the single denominator den shared by both
  // lists. A shared denominator is required because a bound evaluation mixes
  // lower and upper powers depending on the coefficient signs.
  struct ApproximantLevel {
    std::vector<Tz> pow_low;
    std::vector<Tz> pow_upp;
    Tz den;
  };
  void Initialize(std::vector<T> const &Pminimal, double const &_val_double,
                  std::vector<std::pair<T, T>> const &l_approx) {
    val_double = _val_double;
    if (val_double < 0) {
      std::cerr << "NTRF: We require that the value val is positive. val_double="
                << val_double << "\n";
      std::cerr << "NTRF: This is arbitrary of us, but this is our choice and easy "
                   "for you to correct\n";
      throw TerminalException{1};
    }
#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
    double sum = 0;
    double expo = 1;
    for (size_t i = 0; i < Pminimal.size(); i++) {
      sum += UniversalScalarConversion<double, T>(Pminimal[i]) * expo;
      expo *= val_double;
    }
    if (T_abs(sum) > threshold_real_alg_check) {
      std::cerr << "NTRF: Error in Initialize\n";
      std::cerr << "NTRF: sum=" << sum << " val_double=" << val_double << "\n";
      std::cerr << "NTRF: Pminimal =";
      for (auto &eVal : Pminimal)
        std::cerr << " " << eVal;
      std::cerr << "\n";
      throw TerminalException{1};
    }
#endif
    // Finding the expression of X^deg
    deg = Pminimal.size() - 1;
    for (int u = 0; u < deg; u++) {
      T val = -Pminimal[u] / Pminimal[deg];
      ExprXdeg.push_back(val);
    }
    // The rescaling y = scal * x. With x^deg = sum_i ExprXdeg[i] x^i we get
    // y^deg = sum_i ExprXdeg[i] scal^{deg-i} y^i and taking scal to be the
    // lcm of the denominators of the ExprXdeg[i] makes every coefficient
    // ExprXdeg[i] scal^{deg-i} an integer since deg - i >= 1.
    scal = 1;
    for (int u = 0; u < deg; u++) {
      Tz eDen = GetDenominator_z(ExprXdeg[u]);
      scal = KernelLCMpair(scal, eDen);
    }
    pow_scal.push_back(Tz(1));
    for (int u = 1; u < deg; u++)
      pow_scal.push_back(pow_scal[u - 1] * scal);
    std::vector<Tz> ExprYdeg(deg);
    Tz pow = scal;
    for (int u = deg - 1; u >= 0; u--) {
      T val = ExprXdeg[u] * T(pow);
#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
      if (GetDenominator_z(val) != 1) {
        std::cerr << "NTRF: The rescaled coefficient should be an integer\n";
        throw TerminalException{1};
      }
#endif
      ExprYdeg[u] = GetNumerator_z(val);
      pow *= scal;
    }
    // The expressions of y^{deg+k} for k=0..deg-2 over 1, y, ..., y^{deg-1},
    // obtained by iterated multiplication by y.
    ExprYpow.push_back(ExprYdeg);
    for (int k = 1; k < deg - 1; k++) {
      std::vector<Tz> const &prev = ExprYpow[k - 1];
      std::vector<Tz> next(deg);
      Tz const &carry = prev[deg - 1];
      next[0] = carry * ExprYdeg[0];
      for (int j = 1; j < deg; j++)
        next[j] = prev[j - 1] + carry * ExprYdeg[j];
      ExprYpow.push_back(next);
    }
    val_double_y = val_double * UniversalScalarConversion<double, Tz>(scal);
    T scal_T(scal);
    for (auto &e_approx : l_approx) {
      T val_low = e_approx.first;
      T val_upp = e_approx.second;
      if (val_low < 0) {
        std::cerr << "NTRF: We require that the value val_low is positive. val_low="
                  << val_low << "\n";
        std::cerr << "NTRF: This is arbitrary of us, but this is our choice and easy "
                     "for you to correct\n";
        throw TerminalException{1};
      }
      if (val_upp <= val_low) {
        std::cerr << "NTRF: The upper bound should be higher than the lower bound.\n";
        std::cerr << "NTRF: val_low=" << val_low << " val_upp=" << val_upp << "\n";
        throw TerminalException{1};
      }
      // The bounds are on x but the sign determination works on the numerator
      // polynomial in y, so the stored powers are those of y = scal * x.
      // For the bounds pl/ql and pu/qu the powers (pl/ql)^i and (pu/qu)^i
      // are stored as integers over the shared denominator
      // den = ql^{deg-1} qu^{deg-1}, so that the bound evaluations and their
      // sign tests stay in integer arithmetic.
      T y_low = val_low * scal_T;
      T y_upp = val_upp * scal_T;
      ApproximantLevel level;
      Tz pl = GetNumerator_z(y_low);
      Tz ql = GetDenominator_z(y_low);
      Tz pu = GetNumerator_z(y_upp);
      Tz qu = GetDenominator_z(y_upp);
      std::vector<Tz> pow_pl(deg), pow_ql(deg), pow_pu(deg), pow_qu(deg);
      pow_pl[0] = 1;
      pow_ql[0] = 1;
      pow_pu[0] = 1;
      pow_qu[0] = 1;
      for (int i = 1; i < deg; i++) {
        pow_pl[i] = pow_pl[i - 1] * pl;
        pow_ql[i] = pow_ql[i - 1] * ql;
        pow_pu[i] = pow_pu[i - 1] * pu;
        pow_qu[i] = pow_qu[i - 1] * qu;
      }
      level.den = pow_ql[deg - 1] * pow_qu[deg - 1];
      level.pow_low.resize(deg - 1);
      level.pow_upp.resize(deg - 1);
      for (int i = 1; i < deg; i++) {
        level.pow_low[i - 1] =
            pow_pl[i] * pow_ql[deg - 1 - i] * pow_qu[deg - 1];
        level.pow_upp[i - 1] =
            pow_pu[i] * pow_qu[deg - 1 - i] * pow_ql[deg - 1];
      }
      SequenceApproximant.push_back(level);
    }
  }

public:
  HelperClassRealField(std::vector<T> const &Pminimal,
                       double const &_val_double,
                       std::vector<std::pair<T, T>> const &l_approx) {
    Initialize(Pminimal, _val_double, l_approx);
  }
  HelperClassRealField(std::string const &eFile) {
    if (!IsExistingFile(eFile)) {
      std::cerr << "NTRF: HelperClassRealField constructor error. eFile=" << eFile << "\n";
      std::cerr << "NTRF: does not exist\n";
      throw TerminalException{1};
    }
    std::ifstream is(eFile);
    is >> deg;
    std::vector<T> Pminimal;
    for (int u = 0; u <= deg; u++) {
      T val;
      is >> val;
      Pminimal.push_back(val);
    }
    //
    double _val_double;
    is >> _val_double;
    //
    std::vector<std::pair<T, T>> l_approx;
    size_t n_approx;
    is >> n_approx;
    for (size_t u = 0; u < n_approx; u++) {
      T val_low, val_upp;
      is >> val_low;
      is >> val_upp;
      l_approx.push_back({val_low, val_upp});
    }
    Initialize(Pminimal, _val_double, l_approx);
  }
  void normalize(Tvec &num, Tz &den) const {
    Tz g = den;
    for (int u = 0; u < deg; u++) {
      if (g == 1)
        return;
      g = KernelGcdPair(g, num[u]);
    }
    if (g == 1)
      return;
    for (int u = 0; u < deg; u++)
      num[u] /= g;
    den /= g;
  }
  // Product computed into a caller provided buffer. The buffer keeps its
  // 2 deg - 1 size across calls so that the limb storage of its entries is
  // reused (no allocation in steady state); only the first deg entries are
  // meaningful on return.
  void ComputeProductInto(Tconv &conv, Tvec const &a, Tvec const &b) const {
    // Schoolbook convolution into degrees 0..2 deg - 2 followed by the
    // reduction of the upper part with the precomputed integer rows.
    size_t conv_len = 2 * deg - 1;
    if (conv.size() != conv_len)
      conv.resize(conv_len);
    for (size_t k = 0; k < conv_len; k++)
      conv[k] = 0;
    for (int i = 0; i < deg; i++)
      for (int j = 0; j < deg; j++)
        AddMul(conv[i + j], a[i], b[j]);
    for (int k = deg - 2; k >= 0; k--) {
      Tz const &val = conv[deg + k];
      if (val != 0) {
        std::vector<Tz> const &row = ExprYpow[k];
        for (int j = 0; j < deg; j++)
          AddMul(conv[j], val, row[j]);
      }
    }
#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
    double result_d = evaluate_as_double(conv, Tz(1));
    double a_d = evaluate_as_double(a, Tz(1));
    double b_d = evaluate_as_double(b, Tz(1));
    if (T_abs(result_d - a_d * b_d) >
        threshold_real_alg_check * (1 + T_abs(a_d * b_d))) {
      std::cerr << "Error in ComputeProduct\n";
      throw TerminalException{1};
    }
#endif
  }
  Tvec ComputeProduct(Tvec const &a, Tvec const &b) const {
    Tconv conv;
    ComputeProductInto(conv, a, b);
    Tvec res(deg);
    for (int u = 0; u < deg; u++)
      res[u] = std::move(conv[u]);
    return res;
  }
  // The quotient of the two numerator polynomials: a / b = qnum / qden with
  // qnum integral and qden > 0. Solved via the linear system M(b) sol = a
  // where M(b) is the matrix of the multiplication by b.
  std::pair<Tvec, Tz> FindQuotient(Tvec const &a, Tvec const &b) const {
    MyMatrix<T> M(deg, deg);
    Tvec row = b;
    for (int i_row = 0; i_row < deg; i_row++) {
      for (int i_col = 0; i_col < deg; i_col++)
        M(i_row, i_col) = T(row[i_col]);
      if (i_row < deg - 1) {
        // Multiplication of the row by y.
        Tz carry = row[deg - 1];
        for (int j = deg - 1; j > 0; j--)
          row[j] = row[j - 1];
        row[0] = 0;
        if (carry != 0) {
          std::vector<Tz> const &red = ExprYpow[0];
          for (int j = 0; j < deg; j++)
            AddMul(row[j], carry, red[j]);
        }
      }
    }
    MyVector<T> w(deg);
    for (int i = 0; i < deg; i++)
      w(i) = T(a[i]);
    std::optional<MyVector<T>> opt = SolutionMat(M, w);
    if (!opt) {
      std::cerr << "NTRF: Failed to solve the linear system\n";
      throw TerminalException{1};
    }
    MyVector<T> const &eSol = *opt;
    Tz qden(1);
    for (int u = 0; u < deg; u++)
      qden = KernelLCMpair(qden, GetDenominator_z(eSol(u)));
    Tvec qnum(deg);
    for (int u = 0; u < deg; u++)
      qnum[u] = GetNumerator_z(eSol(u)) * (qden / GetDenominator_z(eSol(u)));
#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
    double sol_d = evaluate_as_double(qnum, qden);
    double a_d = evaluate_as_double(a, Tz(1));
    double b_d = evaluate_as_double(b, Tz(1));
    if (T_abs(sol_d * b_d - a_d) >
        threshold_real_alg_check * (1 + T_abs(a_d))) {
      std::cerr << "Error in FindQuotient\n";
      std::cerr << "sol_d=" << sol_d << " a_d=" << a_d << " b_d=" << b_d
                << "\n";
      throw TerminalException{1};
    }
#endif
    return {std::move(qnum), std::move(qden)};
  }
  bool IsStrictlyPositive(Tvec const &x) const {
    // x is the numerator polynomial, assumed to be non-zero; the denominator
    // is positive and does not affect the sign. The bound evaluations are
    // integers over the positive common denominators of the level, so the
    // sign tests stay in integer arithmetic.
    for (auto &level : SequenceApproximant) {
      Tz val_low = x[0] * level.den;
      Tz val_upp = val_low;
      for (int i = 1; i < deg; i++) {
        if (x[i] > 0) {
          AddMul(val_low, x[i], level.pow_low[i - 1]);
          AddMul(val_upp, x[i], level.pow_upp[i - 1]);
        }
        if (x[i] < 0) {
          AddMul(val_low, x[i], level.pow_upp[i - 1]);
          AddMul(val_upp, x[i], level.pow_low[i - 1]);
        }
      }
#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
      if (val_low > val_upp) {
        std::cerr << "The ordering of values is not respected\n";
        throw TerminalException{1};
      }
#endif
      if (val_upp <= 0) {
#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
        double val_upp_d =
            UniversalScalarConversion<double, T>(T(val_upp) / T(level.den));
        if (val_upp_d > threshold_real_alg_check) {
          std::cerr << "Error in IsStrictlyPositive (it is negative)\n";
          throw TerminalException{1};
        }
#endif
        return false;
      }
      if (val_low >= 0) {
#ifdef SANITY_CHECK_REAL_ALG_NUMERIC
        double val_low_d =
            UniversalScalarConversion<double, T>(T(val_low) / T(level.den));
        if (val_low_d < -threshold_real_alg_check) {
          std::cerr << "Error in IsStrictlyPositive (it is positive)\n";
          throw TerminalException{1};
        }
#endif
        return true;
      }
    }
    std::cerr << "x =";
    for (auto &eVal : x)
      std::cerr << " " << eVal;
    std::cerr << "\n";
    std::cerr << "Failed to find an approximant that allows to conclude, "
                 "please produce better approximants\n";
    throw TerminalException{1};
  }
  template <typename Tvect>
  double evaluate_as_double(Tvect const &num, Tz const &den) const {
    double ret_val = 0;
    double pow_double = 1.0;
    for (int i = 0; i < deg; i++) {
      double coeff = UniversalScalarConversion<double, Tz>(num[i]);
      ret_val += coeff * pow_double;
      pow_double *= val_double_y;
    }
    return ret_val / UniversalScalarConversion<double, Tz>(den);
  }
  // Conversion to the coefficients over the powers of x: the change of basis
  // is diagonal since y^i = scal^i x^i.
  std::vector<T> get_x_basis(Tvec const &num, Tz const &den) const {
    std::vector<T> V(deg);
    for (int u = 0; u < deg; u++)
      V[u] = T(num[u] * pow_scal[u]) / T(den);
    return V;
  }
  void set_from_x_basis(std::vector<T> const &V, Tvec &num, Tz &den) const {
    std::vector<T> b(deg);
    for (int u = 0; u < deg; u++)
      b[u] = V[u] / T(pow_scal[u]);
    den = 1;
    for (int u = 0; u < deg; u++)
      den = KernelLCMpair(den, GetDenominator_z(b[u]));
    num.resize(deg);
    for (int u = 0; u < deg; u++)
      num[u] = GetNumerator_z(b[u]) * (den / GetDenominator_z(b[u]));
  }
  int deg;
  std::vector<T> ExprXdeg;

private:
  Tz scal;
  std::vector<Tz> pow_scal;
  std::vector<std::vector<Tz>> ExprYpow;
  double val_double;
  double val_double_y;
  std::vector<ApproximantLevel> SequenceApproximant;
};

std::map<int, HelperClassRealField<Trat_real_field>> list_helper;

void insert_helper_real_algebraic_field(
    int i_field, HelperClassRealField<Trat_real_field> const &hcrf) {
  list_helper.emplace(i_field, hcrf);
}

void print_all_helpers(int val) {
  for (auto &kv : list_helper) {
    std::cerr << "val=" << val << " key=" << kv.first
              << " kv.second.deg=" << kv.second.deg
              << " |kv.second.ExprXdeg|=" << kv.second.ExprXdeg.size() << "\n";
  }
}

// Zero test for std::vector like containers (std::vector, small_vector); the
// MyVector case is handled by IsZeroVector in MAT_MatrixFund.h.
template <typename Tvect> bool IsZeroStdVector(Tvect const &V) {
  for (auto &val : V)
    if (val != 0)
      return false;
  return true;
}

template <int i_field> class RealField;

// Lazy product of two RealField elements -- a minimal expression template (see
// the analogous RatProd / QuadProd). `a * b` returns this proxy; the fast sinks
// evaluate it directly into their own storage:
//   prod  = a * b;   -> RealField::operator=(RealProd)   (move, reuse storage)
//   acc  += a * b;   -> RealField::operator+=(RealProd)  (no wrapper copy)
//   acc  -= a * b;   -> RealField::operator-=(RealProd)
//   RealField r=a*b; -> RealField(RealProd)              (fresh, moved)
// Every other use materializes it into a RealField through the operators after
// the class, so results are identical to the eager version. Holds references:
// consume within the same full-expression, do not bind with `auto`.
template <int i_field> struct RealProd {
  RealField<i_field> const &x;
  RealField<i_field> const &y;
};

template <int i_field> class RealField {
public:
  using T = Trat_real_field;
  using Tz = Tint_real_field;
  using Tvec = Tvec_real_field;
  using Tresidual = T;

private:
  // The element is (num[0] + num[1] y + ... + num[deg-1] y^{deg-1}) / den
  // with den > 0 and gcd(den, num[0], ..., num[deg-1]) = 1.
  Tvec num;
  Tz den;

  static HelperClassRealField<T> const &get_hcrf() {
    static HelperClassRealField<T> const &hcrf = list_helper.at(i_field);
    return hcrf;
  }
  void normalize() { get_hcrf().normalize(num, den); }
  // this += (or -=) onum / oden, followed by the normalization. Templated on
  // the container so that both elements (Tvec) and the convolution scratch
  // buffer (Tconv_real_field, of which only the first deg entries are read)
  // can be merged in.
  template <typename Tvect>
  void axpy_merge(Tvect const &onum, Tz const &oden, bool negate) {
    size_t len = num.size();
    if (den == oden) {
      if (negate) {
        for (size_t u = 0; u < len; u++)
          num[u] -= onum[u];
      } else {
        for (size_t u = 0; u < len; u++)
          num[u] += onum[u];
      }
    } else {
      Tz g = KernelGcdPair(den, oden);
      Tz m_this = oden / g;
      Tz m_o = den / g;
      // Split as in-place multiply then addmul/submul so that no temporary
      // integer is created per coefficient.
      if (negate) {
        for (size_t u = 0; u < len; u++) {
          num[u] *= m_this;
          SubMul(num[u], onum[u], m_o);
        }
      } else {
        for (size_t u = 0; u < len; u++) {
          num[u] *= m_this;
          AddMul(num[u], onum[u], m_o);
        }
      }
      den *= m_this;
    }
    normalize();
  }
  // The numerator polynomial of x - y, whose sign is the sign of x - y since
  // the denominators are positive. No normalization is needed for sign or
  // zero tests.
  static Tvec diff_numerator(RealField<i_field> const &x,
                             RealField<i_field> const &y) {
    size_t len = x.num.size();
    Tvec V(len);
    if (x.den == y.den) {
      for (size_t u = 0; u < len; u++)
        V[u] = x.num[u] - y.num[u];
    } else {
      for (size_t u = 0; u < len; u++)
        V[u] = x.num[u] * y.den - y.num[u] * x.den;
    }
    return V;
  }
  // The numerator polynomial of x - y for y integer.
  static Tvec diff_numerator_int(RealField<i_field> const &x, int const &y) {
    Tvec V = x.num;
    Tz y_z(y);
    SubMul(V[0], y_z, x.den);
    return V;
  }
  RealField(Tvec &&_num, Tz &&_den)
      : num(std::move(_num)), den(std::move(_den)) {}

public:
  Tvec const &get_num() const { return num; }
  Tz const &get_den() const { return den; }
  std::vector<T> get_x_basis() const {
    return get_hcrf().get_x_basis(num, den);
  }
  static size_t get_deg() { return get_hcrf().deg; }

  // Note: We are putting "int" as argument here because we want to do the
  // comparison with the stuff like x > 0 or x = 1. For the type "rational<T>"
  // we had to forbid that because this lead to erroneous conversion of say
  // int64_t to int with catastrophic loss of precision. But for the
  // QuadField<T> the loss of precision does not occur because T is typically
  // mpq_class. or some other type that does not convert to integers easily.
  // And at the same time the natural conversion of int to int64_t allows the
  // comparison x > 0 and equality set x = 1 to work despite the lack of a
  // operator=(int const& u)

  // Constructor
  RealField() : num(get_deg(), Tz(0)), den(1) {}
  RealField(int const &u) : num(get_deg(), Tz(0)), den(1) { num[0] = u; }
  RealField(T const &u) : num(get_deg(), Tz(0)), den(GetDenominator_z(u)) {
    num[0] = GetNumerator_z(u);
  }
  // Constructor from the coefficients over the powers of x.
  RealField(std::vector<T> const &V) {
    get_hcrf().set_from_x_basis(V, num, den);
  }
  // Construct from a lazy product a*b.
  RealField(RealProd<i_field> const &e) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    num = hcrf.ComputeProduct(e.x.num, e.y.num);
    den = e.x.den * e.y.den;
    normalize();
  }
  // assignment operator from int
  RealField<i_field> &operator=(int const &val) {
    size_t len = num.size();
    num[0] = val;
    for (size_t u = 1; u < len; u++)
      num[u] = 0;
    den = 1;
    return *this;
  }
  // Assign from a lazy product a*b. Aliasing-safe: the product is fully
  // computed into the scratch buffer from the operands before this->num is
  // overwritten, and the existing limb storage of this->num is reused.
  RealField<i_field> &operator=(RealProd<i_field> const &e) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    static thread_local Tconv_real_field conv;
    hcrf.ComputeProductInto(conv, e.x.num, e.y.num);
    Tz pd = e.x.den * e.y.den;
    size_t len = num.size();
    for (size_t u = 0; u < len; u++)
      num[u] = conv[u];
    den = std::move(pd);
    normalize();
    return *this;
  }
  //
  // Arithmetic operators below:
  void operator+=(RealField<i_field> const &x) {
    axpy_merge(x.num, x.den, false);
  }
  // Fused accumulate of a lazy product: this += a*b. The product is left
  // unnormalized and a single normalization runs after the merge. The
  // convolution goes into a thread-local scratch buffer whose limb storage
  // survives across the accumulations of a dot product.
  void operator+=(RealProd<i_field> const &e) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    static thread_local Tconv_real_field conv;
    hcrf.ComputeProductInto(conv, e.x.num, e.y.num);
    axpy_merge(conv, e.x.den * e.y.den, false);
  }
  void operator-=(RealField<i_field> const &x) {
    axpy_merge(x.num, x.den, true);
  }
  // Fused subtract of a lazy product: this -= a*b.
  void operator-=(RealProd<i_field> const &e) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    static thread_local Tconv_real_field conv;
    hcrf.ComputeProductInto(conv, e.x.num, e.y.num);
    axpy_merge(conv, e.x.den * e.y.den, true);
  }
  void operator/=(RealField<i_field> const &x) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    std::pair<Tvec, Tz> quot = hcrf.FindQuotient(num, x.num);
    Tz db = x.den;
    size_t len = num.size();
    for (size_t u = 0; u < len; u++)
      num[u] = quot.first[u] * db;
    den *= quot.second;
    normalize();
  }
  friend RealField<i_field> operator+(RealField<i_field> const &x,
                                      RealField<i_field> const &y) {
    RealField<i_field> res = x;
    res.axpy_merge(y.num, y.den, false);
    return res;
  }
  friend RealField<i_field> operator-(RealField<i_field> const &x,
                                      RealField<i_field> const &y) {
    RealField<i_field> res = x;
    res.axpy_merge(y.num, y.den, true);
    return res;
  }
  friend RealField<i_field> operator-(RealField<i_field> const &x,
                                      int const &y) {
    // gcd(num[0] - y den, den) = gcd(num[0], den) so the canonical form is
    // preserved without normalization.
    RealField<i_field> res = x;
    Tz y_z(y);
    SubMul(res.num[0], y_z, res.den);
    return res;
  }
  friend RealField<i_field> operator-(RealField<i_field> const &x) {
    RealField<i_field> res = x;
    size_t len = res.num.size();
    for (size_t u = 0; u < len; u++)
      res.num[u] = -res.num[u];
    return res;
  }
  friend RealField<i_field> operator/(int const &x,
                                      RealField<i_field> const &y) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    Tvec xnum(y.num.size(), Tz(0));
    xnum[0] = x;
    std::pair<Tvec, Tz> quot = hcrf.FindQuotient(xnum, y.num);
    size_t len = y.num.size();
    for (size_t u = 0; u < len; u++)
      quot.first[u] *= y.den;
    RealField<i_field> res(std::move(quot.first), std::move(quot.second));
    res.normalize();
    return res;
  }
  friend RealField<i_field> operator/(RealField<i_field> const &x,
                                      RealField<i_field> const &y) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    std::pair<Tvec, Tz> quot = hcrf.FindQuotient(x.num, y.num);
    size_t len = x.num.size();
    for (size_t u = 0; u < len; u++)
      quot.first[u] *= y.den;
    RealField<i_field> res(std::move(quot.first), Tz(quot.second * x.den));
    res.normalize();
    return res;
  }
  double get_d() const { return get_hcrf().evaluate_as_double(num, den); }
  void operator*=(RealField<i_field> const &x) {
    HelperClassRealField<T> const &hcrf = get_hcrf();
    static thread_local Tconv_real_field conv;
    hcrf.ComputeProductInto(conv, num, x.num);
    size_t len = num.size();
    for (size_t u = 0; u < len; u++)
      num[u] = conv[u];
    den *= x.den;
    normalize();
  }
  // Lazy: returns a RealProd proxy (see above), evaluated in place by the
  // consumer. Mixed int*RealField stays eager below.
  friend RealProd<i_field> operator*(RealField<i_field> const &x,
                                     RealField<i_field> const &y) {
    return RealProd<i_field>{x, y};
  }
  friend RealField<i_field> operator*(int const &x,
                                      RealField<i_field> const &y) {
    RealField<i_field> res = y;
    size_t len = res.num.size();
    for (size_t u = 0; u < len; u++)
      res.num[u] *= x;
    res.normalize();
    return res;
  }
  friend std::ostream &operator<<(std::ostream &os,
                                  RealField<i_field> const &v) {
    std::vector<T> V = v.get_x_basis();
    WriteVectorFromRealAlgebraicString(os, V);
    return os;
  }
  friend std::istream &operator>>(std::istream &is, RealField<i_field> &v) {
    size_t deg = get_deg();
    std::vector<T> V = ReadVectorFromRealAlgebraicString<T>(is, deg);
    v = RealField<i_field>(V);
    return is;
  }
  friend bool operator==(RealField<i_field> const &x,
                         RealField<i_field> const &y) {
    // Both sides are in canonical form.
    if (x.den != y.den)
      return false;
    size_t deg = x.num.size();
    for (size_t u = 0; u < deg; u++) {
      if (x.num[u] != y.num[u]) {
        return false;
      }
    }
    return true;
  }
  friend bool operator!=(RealField<i_field> const &x,
                         RealField<i_field> const &y) {
    return !(x == y);
  }
  friend bool operator!=(RealField<i_field> const &x, int const &y) {
    size_t deg = x.num.size();
    for (size_t u = 1; u < deg; u++) {
      if (x.num[u] != 0) {
        return true;
      }
    }
    return x.num[0] != y * x.den;
  }
  friend bool IsNonNegative(RealField<i_field> const &x) {
    if (IsZeroStdVector(x.num))
      return true;
    return get_hcrf().IsStrictlyPositive(x.num);
  }
  friend bool operator>=(RealField<i_field> const &x,
                         RealField<i_field> const &y) {
    Tvec V = diff_numerator(x, y);
    if (IsZeroStdVector(V))
      return true;
    return get_hcrf().IsStrictlyPositive(V);
  }
  friend bool operator>=(RealField<i_field> const &x, int const &y) {
    Tvec V = diff_numerator_int(x, y);
    if (IsZeroStdVector(V))
      return true;
    return get_hcrf().IsStrictlyPositive(V);
  }
  friend bool operator<=(RealField<i_field> const &x,
                         RealField<i_field> const &y) {
    return y >= x;
  }
  friend bool operator<=(RealField<i_field> const &x, int const &y) {
    Tvec V = diff_numerator_int(x, y);
    if (IsZeroStdVector(V))
      return true;
    for (auto &val : V)
      val = -val;
    return get_hcrf().IsStrictlyPositive(V);
  }
  friend bool operator>(RealField<i_field> const &x,
                        RealField<i_field> const &y) {
    Tvec V = diff_numerator(x, y);
    if (IsZeroStdVector(V)) {
      return false;
    }
    return get_hcrf().IsStrictlyPositive(V);
  }
  friend bool operator>(RealField<i_field> const &x, int const &y) {
    Tvec V = diff_numerator_int(x, y);
    if (IsZeroStdVector(V)) {
      return false;
    }
    return get_hcrf().IsStrictlyPositive(V);
  }
  friend bool operator<(RealField<i_field> const &x,
                        RealField<i_field> const &y) {
    Tvec V = diff_numerator(y, x);
    if (IsZeroStdVector(V)) {
      return false;
    }
    return get_hcrf().IsStrictlyPositive(V);
  }
  friend bool operator<(RealField<i_field> const &x, int const &y) {
    Tvec V = diff_numerator_int(x, y);
    if (IsZeroStdVector(V)) {
      return false;
    }
    for (auto &val : V)
      val = -val;
    return get_hcrf().IsStrictlyPositive(V);
  }
};

// ---------------------------------------------------------------------------
// RealProd (the lazy a*b proxy) as a first-class value. Every use other than the
// in-place sinks above materializes the proxy into a RealField and delegates to
// the ordinary RealField operators, so results are identical to the eager
// implementation. Arithmetic operators return RealField explicitly so that a
// RealProd produced on the right-hand side is materialized before the operand
// temporaries die.
// ---------------------------------------------------------------------------
template <int i_field>
inline RealField<i_field> const &real_eval(RealField<i_field> const &x) {
  return x;
}
template <int i_field>
inline RealField<i_field> real_eval(RealProd<i_field> const &e) {
  return RealField<i_field>(e);
}

#define REALFIELD_REALPROD_ARITH(OP)                                           \
  template <int i_field>                                                       \
  inline RealField<i_field> operator OP(RealProd<i_field> const &a,            \
                                        RealProd<i_field> const &b) {          \
    return real_eval(a) OP real_eval(b);                                       \
  }                                                                            \
  template <int i_field>                                                       \
  inline RealField<i_field> operator OP(RealProd<i_field> const &a,            \
                                        RealField<i_field> const &b) {         \
    return real_eval(a) OP b;                                                  \
  }                                                                            \
  template <int i_field>                                                       \
  inline RealField<i_field> operator OP(RealField<i_field> const &a,           \
                                        RealProd<i_field> const &b) {          \
    return a OP real_eval(b);                                                  \
  }
REALFIELD_REALPROD_ARITH(+)
REALFIELD_REALPROD_ARITH(-)
REALFIELD_REALPROD_ARITH(*)
REALFIELD_REALPROD_ARITH(/)
#undef REALFIELD_REALPROD_ARITH

#define REALFIELD_REALPROD_CMP(OP)                                             \
  template <int i_field>                                                       \
  inline bool operator OP(RealProd<i_field> const &a,                          \
                          RealProd<i_field> const &b) {                        \
    return real_eval(a) OP real_eval(b);                                       \
  }                                                                            \
  template <int i_field>                                                       \
  inline bool operator OP(RealProd<i_field> const &a,                          \
                          RealField<i_field> const &b) {                       \
    return real_eval(a) OP b;                                                  \
  }                                                                            \
  template <int i_field>                                                       \
  inline bool operator OP(RealField<i_field> const &a,                         \
                          RealProd<i_field> const &b) {                        \
    return a OP real_eval(b);                                                  \
  }                                                                            \
  template <int i_field>                                                       \
  inline bool operator OP(RealProd<i_field> const &a, int const &b) {          \
    return real_eval(a) OP b;                                                  \
  }
REALFIELD_REALPROD_CMP(==)
REALFIELD_REALPROD_CMP(!=)
REALFIELD_REALPROD_CMP(<)
REALFIELD_REALPROD_CMP(>)
REALFIELD_REALPROD_CMP(<=)
REALFIELD_REALPROD_CMP(>=)
#undef REALFIELD_REALPROD_CMP

template <int i_field>
inline RealField<i_field> operator-(RealProd<i_field> const &e) {
  return -real_eval(e);
}
template <int i_field>
inline bool IsNonNegative(RealProd<i_field> const &e) {
  return IsNonNegative(RealField<i_field>(e));
}
template <int i_field>
inline std::ostream &operator<<(std::ostream &os, RealProd<i_field> const &e) {
  return os << RealField<i_field>(e);
}

// For this construction we cannot hope to handle rings and fields nicely

template <int i_field> struct overlying_field<RealField<i_field>> {
  typedef RealField<i_field> field_type;
};

// Note that the underlying ring is not unique, there are many possibiliies
// actually but we can represent only one in our scheme.
template <int i_field> struct underlying_ring<RealField<i_field>> {
  typedef RealField<i_field> ring_type;
};

template <int i_field>
inline void TYPE_CONVERSION(stc<RealField<i_field>> const &eQ, double &eD) {
  eD = eQ.val.get_d();
}

template <int i_field>
inline void TYPE_CONVERSION(stc<RealField<i_field>> const &eQ,
                            RealField<i_field> &eD) {
  eD = eQ.val;
}

template <int i_field> struct is_totally_ordered<RealField<i_field>> {
  static const bool value = true;
};

template <int i_field> struct is_ring_field<RealField<i_field>> {
  static const bool value = true;
};

// Exact real-algebraic field with heavy arithmetic (polynomial reduction per
// operation): the same reasoning as QuadField makes Bareiss preferable to
// classical Gaussian elimination (see use_bareiss_for_determinants).
template <int i_field>
struct use_bareiss_for_determinants<RealField<i_field>> {
  static const bool value = true;
};

// Fraction-free LU inverse measured 23% slower on the G553 dual description
// benchmark (dimension 7 systems over a degree 4 field): with the integer
// polynomial representation the classical elimination normalizes cheaply and
// fraction-free growth does not pay off at these sizes.
template <int i_field> struct use_fraction_free_lu<RealField<i_field>> {
  static const bool value = false;
};

// FMA form (see is_fma_prefered). The direct/fused form is fastest for RealField
// (measured): operator+=(RealProd) accumulates the product in place, while the
// scratch form move-assigns the product vector into a temporary first.
template <int i_field> struct is_fma_prefered<RealField<i_field>> {
  static const bool value = true;
};

template <int i_field> struct is_exact_arithmetic<RealField<i_field>> {
  static const bool value = true;
};

// Hashing function

template <int i_field> struct is_implementation_of_Z<RealField<i_field>> {
  static const bool value = false;
};

template <int i_field> struct is_implementation_of_Q<RealField<i_field>> {
  static const bool value = false;
};

// Hashing function

namespace std {
template <int i_field> struct hash<RealField<i_field>> {
  std::size_t operator()(const RealField<i_field> &x) const {
    auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    size_t seed = 0x9e2479b9;
    for (auto &val : x.get_num()) {
      size_t e_hash = std::hash<Tint_real_field>()(val);
      combine_hash(seed, e_hash);
    }
    combine_hash(seed, std::hash<Tint_real_field>()(x.get_den()));
    return seed;
  }
};
// clang-format off
}  // namespace std
// clang-format on

// Local typing info

template <typename T> struct is_real_algebraic_field {};

template <int i_field> struct is_real_algebraic_field<RealField<i_field>> {
  static const bool value = true;
};

// Some functionality

template <int i_field> bool IsInteger(RealField<i_field> const &x) {
  Tvec_real_field const &num = x.get_num();
  size_t len = num.size();
  for (size_t u = 1; u < len; u++)
    if (num[u] != 0)
      return false;
  // The element is num[0] / den in canonical form.
  return x.get_den() == 1;
}

// The conversion tools (int)

template <typename T2, int i_field>
requires (!is_real_algebraic_field<T2>::value)
inline void TYPE_CONVERSION(stc<RealField<i_field>> const &x1, T2 &x2) {
  Tvec_real_field const &num = x1.val.get_num();
  size_t len = num.size();
  for (size_t u = 1; u < len; u++) {
    if (num[u] != 0) {
      std::string str = "Conversion error for real algebraic field";
      throw ConversionException{str};
    }
  }
  Trat_real_field val =
      Trat_real_field(num[0]) / Trat_real_field(x1.val.get_den());
  stc<Trat_real_field> a1{val};
  TYPE_CONVERSION(a1, x2);
}

// Serialization stuff. The archive contains the coefficients over the powers
// of x as rationals, so the format is independent of the internal
// representation.

namespace boost::serialization {

template <class Archive, int i_field>
inline void serialize(Archive &ar, RealField<i_field> &val,
                      [[maybe_unused]] const unsigned int version) {
  size_t deg = RealField<i_field>::get_deg();
  if constexpr (Archive::is_saving::value) {
    std::vector<Trat_real_field> V = val.get_x_basis();
    for (auto &e_val : V)
      ar &make_nvp("realfield_seq", e_val);
  } else {
    std::vector<Trat_real_field> V(deg);
    for (auto &e_val : V)
      ar &make_nvp("realfield_seq", e_val);
    val = RealField<i_field>(V);
  }
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// Turning into something rational

template <typename Tring, int i_field>
void ScalingInteger_Kernel(stc<RealField<i_field>> const &x, Tring &x_res) {
  std::vector<Trat_real_field> V = x.val.get_x_basis();
  std::vector<Tring> Vd;
  for (auto &val : V)
    Vd.push_back(GetDenominator_z(val));
  x_res = LCMlist(Vd);
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYREALFIELD_H_
// clang-format on
