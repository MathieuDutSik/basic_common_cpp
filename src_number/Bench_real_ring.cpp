// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Benchmark of the real algebraic arithmetic run over the field RealField
// against the same computation run over its underlying ring RealRing, on the
// cubic field of discriminant 49 (the generator is 2*cos(2*pi/7), of minimal
// polynomial X^3 + X^2 - 2X - 1, which is monic so the powers of the generator
// span a ring).
//
// The point of the ring is that its elements carry no denominator, hence no
// gcd normalization after every addition and multiplication, which is what
// dominates the field arithmetic. Everything timed here has integral input
// and integral output, so both types compute the very same values and the
// comparison is exact.

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "MAT_Matrix.h"
// clang-format on
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

static auto now() { return std::chrono::steady_clock::now(); }
static double ms(std::chrono::steady_clock::duration d) {
  return std::chrono::duration<double, std::milli>(d).count();
}

int const idx_discriminant_49 = 1;
using Tfield = RealField<idx_discriminant_49>;
using Tring = underlying_ring<Tfield>::ring_type;

// A deterministic pseudo-random sequence, so that the field run and the ring
// run see exactly the same matrix.
struct SmallRandom {
  uint64_t state;
  explicit SmallRandom(uint64_t seed) : state(seed) {}
  int next(int modulo) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<int>((state >> 33) % modulo);
  }
};

// The same matrix over the two types, with entries in Z[alpha].
template <typename T>
MyMatrix<T> BuildMatrix(int n, int deg, int spread, uint64_t seed) {
  SmallRandom rnd(seed);
  MyMatrix<T> M(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      std::vector<typename T::Tresidual> V(deg);
      for (int u = 0; u < deg; u++)
        V[u] = rnd.next(2 * spread + 1) - spread;
      M(i, j) = T(V);
    }
  return M;
}

template <typename T>
void run(std::string const &name, int n, int deg, int spread, long n_prod,
         long n_det, std::vector<double> &out_det) {
  MyMatrix<T> A = BuildMatrix<T>(n, deg, spread, 987654321ULL);
  MyMatrix<T> B = BuildMatrix<T>(n, deg, spread, 123456789ULL);
  //
  auto t0 = now();
  MyMatrix<T> C = A;
  for (long it = 0; it < n_prod; it++)
    C = A * B;
  double t_prod = ms(now() - t0);
  //
  t0 = now();
  T det(0);
  for (long it = 0; it < n_det; it++)
    det = DeterminantMat(A);
  double t_det = ms(now() - t0);
  //
  // A dot product loop, the innermost kernel of every matrix routine.
  t0 = now();
  T acc(0);
  for (long it = 0; it < n_prod; it++)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        acc += A(i, j) * B(j, i);
  double t_dot = ms(now() - t0);
  //
  // The comparison operator, which drives every sign test.
  t0 = now();
  long n_pos = 0;
  for (long it = 0; it < n_prod; it++)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (A(i, j) > B(i, j))
          n_pos++;
  double t_cmp = ms(now() - t0);
  //
  std::cout << "  " << name << "  product=" << t_prod << " ms  determinant="
            << t_det << " ms  dot=" << t_dot << " ms  compare=" << t_cmp
            << " ms\n";
  out_det.push_back(UniversalScalarConversion<double, T>(det));
  std::cerr << "INFO: " << name << " det=" << det << " trace_acc=" << acc
            << " n_pos=" << n_pos << "\n";
}

int main(int argc, char *argv[]) {
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
    insert_helper_real_algebraic_field(idx_discriminant_49, hcrf);
    int deg = hcrf.deg;
    std::cout << "Cubic field of discriminant 49, degree " << deg
              << ", monic=" << hcrf.is_monic() << "\n";
    //
    int n = 8;
    long n_prod = 40;
    long n_det = 40;
    if (argc == 4) {
      n = ParseScalar<int>(std::string(argv[1]));
      n_prod = ParseScalar<long>(std::string(argv[2]));
      n_det = ParseScalar<long>(std::string(argv[3]));
    }
    std::cout << "Matrices of size " << n << " x " << n << ", " << n_prod
              << " products, " << n_det << " determinants\n";
    std::vector<double> l_det;
    run<Tfield>("RealField", n, deg, 3, n_prod, n_det, l_det);
    run<Tring>("RealRing ", n, deg, 3, n_prod, n_det, l_det);
    //
    // Both runs must produce the same determinant.
    MyMatrix<Tfield> Af = BuildMatrix<Tfield>(n, deg, 3, 987654321ULL);
    MyMatrix<Tring> Ar = BuildMatrix<Tring>(n, deg, 3, 987654321ULL);
    Tfield det_f = DeterminantMat(Af);
    Tring det_r = DeterminantMat(Ar);
    Tfield det_r_f = UniversalScalarConversion<Tfield, Tring>(det_r);
    if (det_f != det_r_f) {
      std::cerr << "The field and the ring disagree on the determinant\n";
      std::cerr << "det_f=" << det_f << " det_r=" << det_r << "\n";
      throw TerminalException{1};
    }
    std::cout << "The two runs agree on the determinant: " << det_f << "\n";
    std::cout << "Normal termination of Bench_real_ring\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Bench_real_ring\n";
    exit(e.eVal);
  }
}
