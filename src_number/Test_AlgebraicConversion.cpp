// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Test of the conversions between the algebraic types (QuadField, RealField,
// RealRing) and the types that are not algebraic, and of the rounding out of
// those fields.
//
// What is covered:
//   --- out of the algebraic types: an element that is a rational number
//       converts, one that is not is refused. This direction was written but
//       unreachable: it is guarded by is_quad_field / is_real_algebraic_field,
//       whose primary templates carried no value member, so the guard was a
//       substitution failure on every type that does not specialize them and
//       the conversion was dropped from the overload set instead of chosen.
//   --- into the algebraic types from the types that are not: the scalar
//       becomes the constant coefficient. A fraction enters a field and is
//       refused by a ring.
//   --- NearestInteger out of the fields into a plain integer type, which is
//       what an LLL reduction over such a field needs: the lattice is Z^n
//       whatever the field the form takes its values in, so the basis
//       transformation stays in Z.
//   --- the rounding of a value that no double can locate, where the start
//       comes from an exact search rather than from the estimate.
//   --- FloorInteger and CeilInteger out of RealField, and the commutation of
//       the floor, the ceiling and the nearest integer with the translations
//       by Z, which also pins the tie of the nearest integer, y + 1/2 going to
//       y, on both sides of the range a double locates.

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryRealField.h"
// clang-format on
#include <string>
#include <vector>

int const idx_field = 1;
using T_rat = mpq_class;
using T_int = mpz_class;

static int n_error = 0;

static void check(bool test, std::string const &name) {
  if (test) {
    std::cerr << "PASS: " << name << "\n";
  } else {
    std::cerr << "FAIL: " << name << "\n";
    n_error++;
  }
}

// The rounding is correct when it lands on an integer no further than 1/2
// from the value. That is the property the size reduction relies on and it
// does not depend on how the tie is broken.
template <typename Tfield, typename Tint>
static bool IsNearestInteger(Tfield const &x, Tint const &n) {
  Tfield n_f = UniversalScalarConversion<Tfield, Tint>(n);
  Tfield two(2);
  return T_abs(x - n_f) * two <= Tfield(1);
}

template <int d> static void process_quad(std::string const &name) {
  using T = QuadField<T_rat, d>;
  using Tring = QuadField<T_int, d>;
  std::string pre = name + ": ";
  // Out of the field.
  T seven(T_rat(7), T_rat(0));
  check(UniversalScalarConversion<T_int, T>(seven) == T_int(7),
        pre + "7 leaves the field as an integer");
  check(UniversalScalarConversion<T_rat, T>(seven) == T_rat(7),
        pre + "7 leaves the field as a rational");
  T irr(T_rat(1), T_rat(1));
  check(!UniversalScalarConversionCheck<T_int, T>(irr).has_value(),
        pre + "1 + sqrt(d) is refused on the way out");
  // Into the field.
  check(UniversalScalarConversion<T, T_int>(T_int(-13)) ==
            T(T_rat(-13), T_rat(0)),
        pre + "an integer enters the field");
  check(UniversalScalarConversion<T, T_rat>(T_rat(3, 4)) ==
            T(T_rat(3, 4), T_rat(0)),
        pre + "a fraction enters the field");
  // Rounding out of the field.
  T sq(T_rat(0), T_rat(1));
  T sq3 = sq * T(T_rat(3));
  check(IsNearestInteger(sq, UniversalNearestScalarInteger<T_int, T>(sq)),
        pre + "sqrt(d) rounds to within 1/2");
  check(IsNearestInteger(sq3, UniversalNearestScalarInteger<T_int, T>(sq3)),
        pre + "3 sqrt(d) rounds to within 1/2");
  check(IsNearestInteger(-sq3, UniversalNearestScalarInteger<T_int, T>(-sq3)),
        pre + "-3 sqrt(d) rounds to within 1/2");
  check(UniversalNearestScalarInteger<T_int, T>(seven) == T_int(7),
        pre + "an integer rounds to itself");
  T half(T_rat(5, 2), T_rat(0));
  check(IsNearestInteger(half, UniversalNearestScalarInteger<T_int, T>(half)),
        pre + "an exact half rounds to within 1/2");
  // The ring typed rounding and the integer typed one agree.
  Tring r = UniversalNearestScalarInteger<Tring, T>(sq3);
  check(r == Tring(UniversalNearestScalarInteger<T_int, T>(sq3), T_int(0)),
        pre + "the ring typed and integer typed roundings agree");
  // Past the range where a double locates the integer, so the start comes
  // from the exact truncation. A walk from zero would not come back here.
  T_rat big(1);
  for (int i = 0; i < 40; i++) {
    big *= 1024;
  }
  T big_f(big, T_rat(0));
  check(UniversalNearestScalarInteger<T_int, T>(big_f) ==
            UniversalScalarConversion<T_int, T_rat>(big),
        pre + "2^400 rounds to itself");
  T big_irr = big_f + sq;
  check(IsNearestInteger(big_irr,
                         UniversalNearestScalarInteger<T_int, T>(big_irr)),
        pre + "2^400 + sqrt(d) rounds to within 1/2");
  // Floor and ceiling. On these small values the double is exact enough to
  // give the expected answer.
  auto floor_z = [](T const &v) -> T_int {
    return UniversalFloorScalarInteger<T_int, T>(v);
  };
  auto ceil_z = [](T const &v) -> T_int {
    return UniversalCeilScalarInteger<T_int, T>(v);
  };
  auto nearest_z = [](T const &v) -> T_int {
    return UniversalNearestScalarInteger<T_int, T>(v);
  };
  auto to_z = [](double v) -> T_int { return T_int(static_cast<long>(v)); };
  double sq_d = std::sqrt(static_cast<double>(d));
  check(floor_z(sq) == to_z(std::floor(sq_d)) &&
            ceil_z(sq) == to_z(std::ceil(sq_d)),
        pre + "floor and ceil of sqrt(d)");
  check(floor_z(-sq3) == to_z(std::floor(-3 * sq_d)) &&
            ceil_z(-sq3) == to_z(std::ceil(-3 * sq_d)),
        pre + "floor and ceil of -3 sqrt(d)");
  check(floor_z(seven) == 7 && ceil_z(seven) == 7,
        pre + "an integer is its own floor and ceil");
  // The three output types agree.
  check(UniversalFloorScalarInteger<Tring, T>(sq3) ==
            Tring(floor_z(sq3), T_int(0)),
        pre + "the ring typed and integer typed floors agree");
  check(UniversalFloorScalarInteger<T, T>(sq3) ==
            T(T_rat(floor_z(sq3)), T_rat(0)),
        pre + "the field typed and integer typed floors agree");
  check(UniversalCeilScalarInteger<Tring, T>(sq3) ==
            Tring(ceil_z(sq3), T_int(0)),
        pre + "the ring typed and integer typed ceils agree");
  check(UniversalNearestScalarInteger<T, T>(sq3) ==
            T(T_rat(nearest_z(sq3)), T_rat(0)),
        pre + "the field typed and integer typed nearest agree");
  // Floor(v + n) = Floor(v) + n, and the same for the ceiling and the nearest
  // integer, on both sides of the range where a double locates the integer.
  // On the halves it pins the tie as well: y + 1/2 goes to y, whatever the
  // magnitude.
  T one_half(T_rat(1, 2), T_rat(0));
  std::vector<T_rat> l_shift{T_rat(0), T_rat(1), T_rat(-1), T_rat(1000003),
                             T_rat(-1000003), big, -big};
  std::vector<T> l_val{sq, -sq, sq3, one_half, -one_half, half,
                       T(T_rat(-5, 2), T_rat(0)), T(T_rat(7, 3), T_rat(1))};
  bool is_covariant = true;
  for (auto &shift : l_shift) {
    T_int shift_z = UniversalScalarConversion<T_int, T_rat>(shift);
    T shift_f(shift, T_rat(0));
    for (auto &v : l_val) {
      T w = v + shift_f;
      if (floor_z(w) != floor_z(v) + shift_z ||
          ceil_z(w) != ceil_z(v) + shift_z ||
          nearest_z(w) != nearest_z(v) + shift_z) {
        std::cerr << "shift=" << shift << " v=" << v << "\n";
        is_covariant = false;
      }
    }
  }
  check(is_covariant, pre + "floor, ceil and nearest commute with Z");
  check(nearest_z(one_half) == 0 && nearest_z(-one_half) == -1 &&
            nearest_z(half) == 2 && nearest_z(-half) == -3,
        pre + "a tie y + 1/2 goes to y");
  check(nearest_z(big_f + one_half) ==
            UniversalScalarConversion<T_int, T_rat>(big),
        pre + "2^400 + 1/2 goes to 2^400, as the small ties do");
  check(floor_z(big_irr) ==
            UniversalScalarConversion<T_int, T_rat>(big) + to_z(sq_d),
        pre + "the floor of 2^400 + sqrt(d) is 2^400 + floor(sqrt(d))");
}

static void process_real(std::string const &eFile) {
  // The rationals and integers the field is built on.
  using T_rat = Trat_real_field;
  using T_int = Tint_real_field;
  using T = RealField<idx_field>;
  using Tring = RealRing<idx_field>;
  HelperClassRealField<T_rat> hcrf(eFile);
  insert_helper_real_algebraic_field(idx_field, hcrf);
  std::string pre = "RealField: ";
  // Out of the field.
  T seven(T_rat(7));
  check(UniversalScalarConversion<T_int, T>(seven) == T_int(7),
        pre + "7 leaves the field as an integer");
  check(UniversalScalarConversion<T_rat, T>(seven) == T_rat(7),
        pre + "7 leaves the field as a rational");
  // The generator itself, which is irrational.
  size_t deg = static_cast<size_t>(hcrf.deg);
  std::vector<T_rat> Vx(deg, T_rat(0));
  Vx[1] = T_rat(1);
  T x(Vx);
  check(!UniversalScalarConversionCheck<T_int, T>(x).has_value(),
        pre + "the generator is refused on the way out");
  check(!UniversalScalarConversionCheck<T_int, Tring>(
             UniversalScalarConversion<Tring, T>(x))
             .has_value(),
        pre + "the generator is refused on the way out of the ring");
  // Into the field and into the ring.
  check(UniversalScalarConversion<T, T_int>(T_int(-13)) == T(T_rat(-13)),
        pre + "an integer enters the field");
  check(UniversalScalarConversion<T, T_rat>(T_rat(3, 4)) == T(T_rat(3, 4)),
        pre + "a fraction enters the field");
  check(UniversalScalarConversion<Tring, T_int>(T_int(-13)) ==
            Tring(T_int(-13)),
        pre + "an integer enters the ring");
  check(!UniversalScalarConversionCheck<Tring, T_rat>(T_rat(3, 4)).has_value(),
        pre + "a fraction is refused by the ring");
  // Rounding out of the field.
  T x3 = x * T(T_rat(3));
  check(IsNearestInteger(x, UniversalNearestScalarInteger<T_int, T>(x)),
        pre + "the generator rounds to within 1/2");
  check(IsNearestInteger(x3, UniversalNearestScalarInteger<T_int, T>(x3)),
        pre + "3 x rounds to within 1/2");
  check(IsNearestInteger(-x3, UniversalNearestScalarInteger<T_int, T>(-x3)),
        pre + "-3 x rounds to within 1/2");
  check(UniversalNearestScalarInteger<T_int, T>(seven) == T_int(7),
        pre + "an integer rounds to itself");
  Tring r = UniversalNearestScalarInteger<Tring, T>(x3);
  check(r == Tring(UniversalNearestScalarInteger<T_int, T>(x3)),
        pre + "the ring typed and integer typed roundings agree");
  // Past the range where a double locates the integer.
  T_rat big(1);
  for (int i = 0; i < 40; i++) {
    big *= 1024;
  }
  T big_f(big);
  check(UniversalNearestScalarInteger<T_int, T>(big_f) ==
            UniversalScalarConversion<T_int, T_rat>(big),
        pre + "2^400 rounds to itself");
  T big_irr = big_f + x;
  check(IsNearestInteger(big_irr,
                         UniversalNearestScalarInteger<T_int, T>(big_irr)),
        pre + "2^400 + x rounds to within 1/2");
  // Floor and ceiling, on exact values: x = 2 cos(2 pi / 7) = 1.2469...
  auto floor_z = [](T const &v) -> T_int {
    return UniversalFloorScalarInteger<T_int, T>(v);
  };
  auto ceil_z = [](T const &v) -> T_int {
    return UniversalCeilScalarInteger<T_int, T>(v);
  };
  auto nearest_z = [](T const &v) -> T_int {
    return UniversalNearestScalarInteger<T_int, T>(v);
  };
  check(floor_z(x) == 1 && ceil_z(x) == 2, pre + "floor and ceil of x");
  check(floor_z(-x) == -2 && ceil_z(-x) == -1, pre + "floor and ceil of -x");
  check(floor_z(x3) == 3 && ceil_z(x3) == 4, pre + "floor and ceil of 3 x");
  check(floor_z(seven) == 7 && ceil_z(seven) == 7,
        pre + "an integer is its own floor and ceil");
  T mseven(T_rat(-7));
  check(floor_z(mseven) == -7 && ceil_z(mseven) == -7,
        pre + "a negative integer is its own floor and ceil");
  T half(T_rat(1, 2));
  check(floor_z(half) == 0 && ceil_z(half) == 1, pre + "floor and ceil of 1/2");
  check(floor_z(-half) == -1 && ceil_z(-half) == 0,
        pre + "floor and ceil of -1/2");
  // The three output types agree.
  check(UniversalFloorScalarInteger<Tring, T>(x3) == Tring(floor_z(x3)),
        pre + "the ring typed and integer typed floors agree");
  check(UniversalFloorScalarInteger<T, T>(x3) == T(T_rat(floor_z(x3))),
        pre + "the field typed and integer typed floors agree");
  check(UniversalCeilScalarInteger<Tring, T>(x3) == Tring(ceil_z(x3)),
        pre + "the ring typed and integer typed ceils agree");
  // Floor(v + n) = Floor(v) + n, and the same for the ceiling and the nearest
  // integer, on both sides of the range where a double locates the integer.
  // On the halves it pins the tie as well: y + 1/2 goes to y, whatever the
  // magnitude.
  std::vector<T_rat> l_shift{T_rat(0), T_rat(1), T_rat(-1), T_rat(1000003),
                             T_rat(-1000003), big, -big};
  std::vector<T> l_val{x, -x, x3, half, -half, T(T_rat(3, 2)),
                       T(T_rat(-5, 2)), T(T_rat(7, 3))};
  bool is_covariant = true;
  for (auto &shift : l_shift) {
    T_int shift_z = UniversalScalarConversion<T_int, T_rat>(shift);
    T shift_f(shift);
    for (auto &v : l_val) {
      T w = v + shift_f;
      if (floor_z(w) != floor_z(v) + shift_z ||
          ceil_z(w) != ceil_z(v) + shift_z ||
          nearest_z(w) != nearest_z(v) + shift_z) {
        std::cerr << "shift=" << shift << " v=" << v << "\n";
        is_covariant = false;
      }
    }
  }
  check(is_covariant, pre + "floor, ceil and nearest commute with Z");
  check(nearest_z(half) == 0 && nearest_z(-half) == -1 &&
            nearest_z(T(T_rat(5, 2))) == 2 &&
            nearest_z(T(T_rat(-5, 2))) == -3,
        pre + "a tie y + 1/2 goes to y");
  check(nearest_z(big_f + half) == UniversalScalarConversion<T_int, T_rat>(big),
        pre + "2^400 + 1/2 goes to 2^400, as the small ties do");
  check(floor_z(big_irr) == UniversalScalarConversion<T_int, T_rat>(big) + 1,
        pre + "the floor of 2^400 + x is 2^400 + 1");
}

int main() {
  try {
    process_quad<2>("Qsqrt2");
    process_quad<5>("Qsqrt5");
    //
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
    process_real(eFile);
    //
    if (n_error > 0) {
      std::cerr << "Test_AlgebraicConversion: " << n_error << " error(s)\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of Test_AlgebraicConversion\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_AlgebraicConversion\n";
    exit(e.eVal);
  }
}
