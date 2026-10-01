// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "factorizations.h"
// clang-format on

// Every computed value is compared with its exact expected value: a test
// that only prints is passed by a wrong arithmetic.
template <typename Trat> void process(std::string const &name) {
  using T = QuadField<Trat, 5>;
  std::cerr << "name=" << name << "\n";
  auto check = [&](std::string const &context, T const &result,
                   T const &expected) -> void {
    if (result != expected) {
      std::cerr << "name=" << name << " " << context << " result=" << result
                << " expected=" << expected << "\n";
      throw TerminalException{1};
    }
  };
  auto q = [](int num, int den) -> Trat { return Trat(num) / Trat(den); };
  T x;
  T y = UniversalScalarConversion<T, T>(x);
  check("UniversalScalarConversion", y, x);
  T near = UniversalNearestScalarInteger<T, T>(x);
  check("UniversalNearestScalarInteger", near, x);
  //
  // The powers of the golden ratio phi = (1 + sqrt(5)) / 2: phi^n is
  // (L_n + F_n sqrt(5)) / 2 with the Lucas and Fibonacci numbers.
  T phi(q(1, 2), q(1, 2));
  T pow(1);
  for (int i = 0; i < 10; i++)
    pow *= phi;
  check("phi^10", pow, T(q(123, 2), q(55, 2)));
  check("phi^2 = phi + 1", T(phi * phi), phi + T(1));
  //
  // The product assigned to one of its operands, z = z * w, z = w * z and
  // z = z * z, the form of the in place rescalings M(i,j) = M(i,j) * c.
  T u(q(-3, 4), q(-1, 4));
  T w(Trat(2), q(1, 3));
  T z = u;
  z = z * T(1);
  check("z = z * 1", z, u);
  T expected_uw(q(-23, 12), q(-3, 4));
  z = u;
  z = z * w;
  check("z = z * w", z, expected_uw);
  z = u;
  z = w * z;
  check("z = w * z", z, expected_uw);
  z = u;
  z = z * z;
  check("z = z * z", z, T(q(7, 8), q(3, 8)));
  //
  // The icosahedron: its 12 vertices are the cyclic permutations of
  // (0, +-1, +-phi). Its symmetry group is transitive on the vertices, so
  // every row of the Gram matrix has the same values: the norm phi + 2 on the
  // diagonal, phi for the 5 neighbors, -phi for the 5 vertices at distance 2
  // and -(phi + 2) for the antipode.
  MyMatrix<T> M(12, 3);
  int pos = 0;
  for (int s1 : {1, -1})
    for (int s2 : {1, -1}) {
      T a = T(s1);
      T b = T(s2) * phi;
      for (int shift = 0; shift < 3; shift++) {
        M(pos, shift) = T(0);
        M(pos, (shift + 1) % 3) = a;
        M(pos, (shift + 2) % 3) = b;
        pos++;
      }
    }
  auto check_gram = [&](std::string const &context, MyMatrix<T> const &G,
                        T const &scale) -> void {
    std::vector<std::pair<T, int>> pattern{{-phi - T(2), 1},
                                           {-phi, 5},
                                           {phi, 5},
                                           {phi + T(2), 1}};
    for (int i = 0; i < 12; i++) {
      check(context + " diagonal", G(i, i), (phi + T(2)) * scale);
      for (auto const &[val, mult] : pattern) {
        T val_scaled = val * scale;
        int count = 0;
        for (int j = 0; j < 12; j++)
          if (G(i, j) == val_scaled)
            count++;
        if (count != mult) {
          std::cerr << "name=" << name << " " << context << " row " << i
                    << " has " << count << " entries " << val_scaled
                    << " expected " << mult << "\n";
          throw TerminalException{1};
        }
      }
    }
  };
  // The Gram matrix by the matrix product and by the accumulation of the
  // lazy products, the two forms the matrix code uses.
  MyMatrix<T> G = M * M.transpose();
  check_gram("M * M^T", G, T(1));
  MyMatrix<T> Gacc(12, 12);
  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 12; j++) {
      T sum(0);
      for (int k = 0; k < 3; k++)
        sum += M(i, k) * M(j, k);
      Gacc(i, j) = sum;
    }
  check_gram("accumulated", Gacc, T(1));
  // Rescaling the vertices in place, by 1 and by phi, scales the Gram matrix
  // by the square.
  for (auto const &[c_name, c] : {std::pair<std::string, T>{"1", T(1)},
                                  std::pair<std::string, T>{"phi", phi}}) {
    MyMatrix<T> Mc = M;
    for (int i = 0; i < 12; i++)
      for (int k = 0; k < 3; k++)
        Mc(i, k) = Mc(i, k) * c;
    MyMatrix<T> Gc = Mc * Mc.transpose();
    check_gram("rescaled by " + c_name, Gc, c * c);
  }
}

int main() {
  try {
    process<mpq_class>("mpq_class");
#ifdef ENABLE_FLINT_SUPPORT
    process<fmpq_class>("fmpq_class");
#endif
    process<Rational<SafeInt64>>("Rational<SafeInt64>");
    std::cerr << "All the QuadField checks passed\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
