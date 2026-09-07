// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Test of RealRing, the order Z[x] spanned by the powers of the generator of a
// real algebraic field. The field used is the cubic field of discriminant 49:
// the generator is 2*cos(2*pi/7), of minimal polynomial X^3 + X^2 - 2X - 1,
// which is monic, so Z[x] is a ring.
//
// The points checked are:
//   --- the ring arithmetic agrees with the field arithmetic,
//   --- the exact division succeeds when the quotient is in the ring,
//   --- it throws when the quotient is not in the ring,
//   --- a non-monic description of the same field is rejected,
//   --- the conversions ring <-> field and the string round trip,
//   --- the boost serialization, on which the MPI dual description relies.

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "MAT_Matrix.h"
// clang-format on
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <sstream>
#include <string>
#include <vector>

int const idx_monic = 1;
int const idx_non_monic = 2;
using Tfield = RealField<idx_monic>;
using Tring = underlying_ring<Tfield>::ring_type;
using Tz = Tint_real_field;

static int n_error = 0;

static void check(bool test, std::string const &name) {
  if (test) {
    std::cerr << "PASS: " << name << "\n";
  } else {
    std::cerr << "FAIL: " << name << "\n";
    n_error++;
  }
}

static Tring MakeRing(int a, int b, int c) {
  std::vector<Tz> V{Tz(a), Tz(b), Tz(c)};
  return Tring(V);
}

static Tfield MakeField(int a, int b, int c) {
  std::vector<mpq_class> V{mpq_class(a), mpq_class(b), mpq_class(c)};
  return Tfield(V);
}

int main() {
  try {
    using T_rat = mpq_class;
    std::string eFile = "CI_tests/RealAlgebraicField/CubicFieldDisc_49";
    bool found = false;
    for (int level = 0; level <= 10; level++) {
      if (FILE_IsExistingFile(eFile)) {
        found = true;
        break;
      }
      eFile = "../" + eFile;
    }
    if (!found) {
      std::cerr << "Failed to find RealAlgebraicField test data after checking "
                   "paths from CI_tests/ up to 10 parent levels\n";
      throw TerminalException{1};
    }
    HelperClassRealField<T_rat> hcrf(eFile);
    insert_helper_real_algebraic_field(idx_monic, hcrf);
    check(hcrf.is_monic(), "X^3 + X^2 - 2X - 1 is recognized as monic");
    std::cerr << "STEP 1: field registered, degree " << hcrf.deg << "\n";
    //
    // The arithmetic agrees with the field arithmetic.
    Tring a = MakeRing(1, 2, 4);
    Tring b = MakeRing(1, 3, 1);
    Tfield af = MakeField(1, 2, 4);
    Tfield bf = MakeField(1, 3, 1);
    auto same = [](Tring const &x, Tfield const &y) -> bool {
      return UniversalScalarConversion<Tfield, Tring>(x) == y;
    };
    check(same(a + b, af + bf), "addition agrees with the field");
    check(same(a - b, af - bf), "subtraction agrees with the field");
    check(same(a * b, af * bf), "multiplication agrees with the field");
    check(same(-a, -af), "negation agrees with the field");
    check(same(3 * a, 3 * af), "scalar multiplication agrees with the field");
    std::cerr << "STEP 2: ring arithmetic checked against the field\n";
    //
    // The ordering agrees with the field ordering.
    check((a > b) == (af > bf), "the comparison a > b agrees with the field");
    check((a < b) == (af < bf), "the comparison a < b agrees with the field");
    check(IsNonNegative(a * a), "a square is non-negative");
    check(MakeRing(0, 0, 0) == Tring(0), "the zero element compares equal");
    check(MakeRing(5, 0, 0) == Tring(Tz(5)), "an integer element compares equal");
    check(IsInteger(MakeRing(5, 0, 0)), "IsInteger on a rational integer");
    check(!IsInteger(a), "IsInteger is false on a genuine ring element");
    std::cerr << "STEP 3: ordering and predicates checked\n";
    //
    // The exact division. a * b is divisible by b, with quotient a.
    Tring prod = a * b;
    Tring quot = prod / b;
    check(quot == a, "the division (a*b)/b returns a");
    check(same(quot, af), "the division agrees with the field");
    check(prod / a == b, "the division (a*b)/a returns b");
    check(MakeRing(6, 12, 24) / MakeRing(3, 0, 0) == MakeRing(2, 4, 8),
          "division by a rational integer");
    std::cerr << "STEP 4: exact divisions checked\n";
    //
    // b happens to be a unit of Z[x]: its norm is +-1, so 1/b lies in the ring
    // and the division goes through.
    Tring b_inv = Tring(1) / b;
    check(b_inv * b == Tring(1), "the inverse of a unit is found in the ring");
    std::cerr << "INFO: 1/b=" << b_inv << "\n";
    // A division whose result is not in the ring has to throw. 1/2 is not in
    // Z[x], and neither is a/2 for the a above.
    auto division_throws = [](Tring const &x, Tring const &y) -> bool {
      try {
        Tring bad = x / y;
        std::cerr << "INFO: the division did not throw, it gave " << bad
                  << "\n";
        return false;
      } catch (TerminalException const &e) {
        return true;
      }
    };
    check(division_throws(Tring(1), Tring(2)),
          "the quotient 1/2, outside the ring, throws");
    check(division_throws(a, Tring(2)),
          "the quotient a/2, outside the ring, throws");
    check(division_throws(a, MakeRing(0, 2, 0)),
          "the quotient a/(2x), outside the ring, throws");
    std::cerr << "STEP 5: the non-representable divisions throw\n";
    //
    // The conversions and the string round trip.
    Tfield a_field = UniversalScalarConversion<Tfield, Tring>(a);
    Tring a_back = UniversalScalarConversion<Tring, Tfield>(a_field);
    check(a_back == a, "the round trip ring -> field -> ring");
    // UniversalScalarConversion turns the ConversionException raised by the
    // conversion into a TerminalException.
    bool conv_thrown = false;
    try {
      Tfield half = MakeField(1, 0, 0) / MakeField(2, 0, 0);
      Tring half_ring = UniversalScalarConversion<Tring, Tfield>(half);
      std::cerr << "INFO: the conversion did not throw, it gave " << half_ring
                << "\n";
    } catch (TerminalException const &e) {
      conv_thrown = true;
    }
    check(conv_thrown, "converting 1/2 to the ring throws");
    //
    std::ostringstream os;
    os << a;
    std::string str = os.str();
    Tring a_read;
    std::istringstream(str) >> a_read;
    check(a_read == a, "the round trip ring -> string -> ring");
    std::cerr << "INFO: a=" << str << "\n";
    std::cerr << "STEP 6: conversions and string round trip checked\n";
    //
    // The boost serialization. POLY_MPI_DualDesc moves matrices between ranks,
    // so both the scalar and the matrix round trip are exercised.
    {
      std::ostringstream oss;
      {
        boost::archive::text_oarchive oa(oss);
        oa << a;
      }
      Tring a_ser;
      {
        std::istringstream iss(oss.str());
        boost::archive::text_iarchive ia(iss);
        ia >> a_ser;
      }
      check(a_ser == a, "the round trip ring -> archive -> ring");
      //
      MyMatrix<Tring> Mser(2, 2);
      Mser(0, 0) = a;
      Mser(0, 1) = MakeRing(5, 0, 0);
      Mser(1, 0) = -a;
      Mser(1, 1) = a * b;
      std::ostringstream oss2;
      {
        boost::archive::text_oarchive oa(oss2);
        oa << Mser;
      }
      MyMatrix<Tring> Nser;
      {
        std::istringstream iss(oss2.str());
        boost::archive::text_iarchive ia(iss);
        ia >> Nser;
      }
      bool mat_ok = (Nser.rows() == 2 && Nser.cols() == 2);
      if (mat_ok) {
        for (int i = 0; i < 2; i++)
          for (int j = 0; j < 2; j++)
            if (Mser(i, j) != Nser(i, j))
              mat_ok = false;
      }
      check(mat_ok, "the round trip MyMatrix<ring> -> archive -> MyMatrix<ring>");
    }
    std::cerr << "STEP 7: serialization checked\n";
    //
    // The matrix operations run over the ring, with the determinant matching
    // the one computed over the field.
    int n = 4;
    MyMatrix<Tring> Mr(n, n);
    MyMatrix<Tfield> Mf(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        int c0 = 1 + ((i * 7 + j * 3) % 5);
        int c1 = ((i * 5 + j * 11) % 7) - 3;
        int c2 = ((i * 13 + j * 2) % 5) - 2;
        Mr(i, j) = MakeRing(c0, c1, c2);
        Mf(i, j) = MakeField(c0, c1, c2);
      }
    Tring det_r = DeterminantMat(Mr);
    Tfield det_f = DeterminantMat(Mf);
    check(same(det_r, det_f), "the determinant over the ring matches the field");
    std::cerr << "INFO: det=" << det_r << "\n";
    MyMatrix<Tring> Pr = Mr * Mr;
    MyMatrix<Tfield> Pf = Mf * Mf;
    bool prod_ok = true;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (!same(Pr(i, j), Pf(i, j)))
          prod_ok = false;
    check(prod_ok, "the matrix product over the ring matches the field");
    std::cerr << "STEP 8: matrix operations checked\n";
    //
    // A non-monic description of the very same field must be rejected. The
    // polynomial 2X^3 + 2X^2 - 4X - 2 has the same root as X^3 + X^2 - 2X - 1.
    std::vector<T_rat> Pmin_non_monic{-2, -4, 2, 2};
    std::vector<std::pair<T_rat, T_rat>> l_approx;
    // Reuse the bracketing of the monic description: it is the same number.
    {
      std::ifstream is(eFile);
      int deg_read;
      is >> deg_read;
      for (int u = 0; u <= deg_read; u++) {
        T_rat val;
        is >> val;
      }
      double val_double_read;
      is >> val_double_read;
      size_t n_approx;
      is >> n_approx;
      for (size_t u = 0; u < n_approx; u++) {
        T_rat val_low, val_upp;
        is >> val_low >> val_upp;
        l_approx.push_back({val_low, val_upp});
      }
      HelperClassRealField<T_rat> hcrf_nm(Pmin_non_monic, val_double_read,
                                          l_approx);
      check(!hcrf_nm.is_monic(),
            "2X^3 + 2X^2 - 4X - 2 is recognized as non-monic");
      insert_helper_real_algebraic_field(idx_non_monic, hcrf_nm);
    }
    // The field still works over a non-monic description.
    using Tfield_nm = RealField<idx_non_monic>;
    Tfield_nm x_nm(std::vector<T_rat>{1, 2, 4});
    std::cerr << "INFO: over the non-monic description, x_nm=" << x_nm << "\n";
    check(x_nm != 0, "RealField works over a non-monic description");
    // The ring does not.
    using Tring_nm = underlying_ring<Tfield_nm>::ring_type;
    bool monic_thrown = false;
    try {
      Tring_nm y_nm = Tring_nm(1);
      std::cerr << "INFO: the ring was built, it gave " << y_nm << "\n";
    } catch (TerminalException const &e) {
      monic_thrown = true;
    }
    check(monic_thrown, "RealRing over a non-monic description throws");
    std::cerr << "STEP 9: the monic requirement is enforced\n";
    //
    std::cerr << "n_error=" << n_error << "\n";
    if (n_error > 0) {
      std::cerr << "Erroneous termination of Test_RealRing\n";
      return 1;
    }
    std::cerr << "Normal termination of Test_RealRing\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
}
