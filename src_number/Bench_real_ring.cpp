// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Benchmark of the matrix operations over a real algebraic field, run over
// the field RealField and over its underlying ring RealRing. The field is
// given by a description file, by default the cubic field of discriminant 49
// of CI_tests/RealAlgebraicField (the generator is 2*cos(2*pi/7), of minimal
// polynomial X^3 + X^2 - 2X - 1). Its minimal polynomial has to be monic, so
// that the powers of the generator span a ring.
//
// The point of the ring is that its elements carry no denominator, hence no
// gcd normalization after every addition and multiplication, which is what
// dominates the field arithmetic. The operations shared by both types have
// integral input and integral output, so both compute the very same values;
// the program checks that they agree on the determinant.
//
// The integers and rationals the elements are built on are GMP by default,
// flint when compiled with -DENABLE_FLINT_SUPPORT -DREALFIELD_USE_FLINT
// (make ENABLE_FLINT_SUPPORT=1 REALFIELD_USE_FLINT=1). Every result is folded
// into a checksum, printed on stderr, so that two builds can be compared
// value for value.
//
// The operations:
//   product      C = A * B
//   determinant  Bareiss, coefficient growth
//   dot          the inner kernel of every matrix routine
//   compare      the sign test, which drives every comparison
//   inverse      field only
//   nullspace    field only, of a matrix of rank n/2 with n + n/2 columns
//   rank         field only, of the same matrix
// Each operation is run several times and the minimum is reported, which is
// the standard way of removing scheduling noise from a benchmark.

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

int const idx_bench_field = 1;
using Tfield = RealField<idx_bench_field>;
using Tring = underlying_ring<Tfield>::ring_type;

// A deterministic pseudo-random sequence, so that the field run and the ring
// run, and the GMP and flint builds, see exactly the same matrices.
struct SmallRandom {
  uint64_t state;
  explicit SmallRandom(uint64_t seed) : state(seed) {}
  int next(int modulo) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<int>((state >> 33) % modulo);
  }
};

// A matrix with entries in Z[alpha], of coefficients in [-spread, spread].
template <typename T>
MyMatrix<T> BuildMatrix(int n_row, int n_col, int deg, int spread,
                        uint64_t seed) {
  SmallRandom rnd(seed);
  MyMatrix<T> M(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++) {
      std::vector<typename T::Tresidual> V(deg);
      for (int u = 0; u < deg; u++)
        V[u] = rnd.next(2 * spread + 1) - spread;
      M(i, j) = T(V);
    }
  return M;
}

template <typename F> double time_best(int reps, F f) {
  double best = -1;
  for (int i = 0; i < reps; i++) {
    auto t0 = now();
    f();
    double dur = ms(now() - t0);
    if (best < 0 || dur < best)
      best = dur;
  }
  return best;
}

// Folds a value into the checksum through its printed form, which is the
// same whatever the integers it is built on.
template <typename T> void fold(size_t &seed, T const &x) {
  std::ostringstream os;
  os << x;
  size_t new_hash = std::hash<std::string>()(os.str());
  seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
template <typename T> void fold(size_t &seed, MyMatrix<T> const &M) {
  for (int i = 0; i < M.rows(); i++)
    for (int j = 0; j < M.cols(); j++)
      fold(seed, M(i, j));
}

template <typename T>
T run(std::string const &name, int n, int deg, int reps) {
  int const spread = 3;
  MyMatrix<T> A = BuildMatrix<T>(n, n, deg, spread, 987654321ULL);
  MyMatrix<T> B = BuildMatrix<T>(n, n, deg, spread, 123456789ULL);
  size_t checksum = 0;
  std::cout << "  " << name;
  auto report = [&](std::string const &op, double dur) {
    std::cout << "  " << op << "=" << dur << " ms";
  };
  //
  MyMatrix<T> C;
  report("product", time_best(reps, [&]() { C = A * B; }));
  fold(checksum, C);
  //
  T det(0);
  report("determinant", time_best(reps, [&]() { det = DeterminantMat(A); }));
  fold(checksum, det);
  //
  T acc(0);
  report("dot", time_best(reps, [&]() {
           acc = 0;
           for (int i = 0; i < n; i++)
             for (int j = 0; j < n; j++)
               acc += A(i, j) * B(j, i);
         }));
  fold(checksum, acc);
  //
  long n_pos = 0;
  report("compare", time_best(reps, [&]() {
           n_pos = 0;
           for (int i = 0; i < n; i++)
             for (int j = 0; j < n; j++)
               if (A(i, j) > B(i, j))
                 n_pos++;
         }));
  fold(checksum, n_pos);
  //
  if constexpr (is_ring_field<T>::value) {
    MyMatrix<T> Inv;
    report("inverse", time_best(reps, [&]() { Inv = Inverse(A); }));
    fold(checksum, Inv);
    //
    int r = n / 2;
    MyMatrix<T> L = BuildMatrix<T>(n, r, deg, spread, 192837465ULL);
    MyMatrix<T> R = BuildMatrix<T>(r, n + r, deg, spread, 564738291ULL);
    MyMatrix<T> D = L * R;
    MyMatrix<T> Ker;
    report("nullspace", time_best(reps, [&]() { Ker = NullspaceMat(D); }));
    fold(checksum, Ker);
    int rnk = 0;
    report("rank", time_best(reps, [&]() { rnk = RankMat(D); }));
    fold(checksum, rnk);
  }
  std::cout << "\n";
  std::cerr << "CHECKSUM " << name << " " << checksum << "\n";
  return det;
}

std::string find_default_description() {
  std::string eFile = "CI_tests/RealAlgebraicField/CubicFieldDisc_49";
  for (int level = 0; level <= 10; level++) {
    if (FILE_IsExistingFile(eFile))
      return eFile;
    eFile = "../" + eFile;
  }
  std::cerr << "Failed to find RealAlgebraicField test data after checking "
               "paths from CI_tests/ up to 10 parent levels\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 1 && argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "Bench_real_ring\n";
      std::cerr << "    or\n";
      std::cerr << "Bench_real_ring [FileDesc] [n] [reps]\n";
      std::cerr << "with FileDesc the description of a real algebraic field "
                   "of monic minimal polynomial, n the size of the matrices "
                   "and reps the number of runs of each operation\n";
      return -1;
    }
    std::string eFile;
    int n = 8;
    int reps = 5;
    if (argc == 4) {
      eFile = argv[1];
      n = ParseScalar<int>(std::string(argv[2]));
      reps = ParseScalar<int>(std::string(argv[3]));
    } else {
      eFile = find_default_description();
    }
    HelperClassRealField<Trat_real_field> hcrf(eFile);
    if (!hcrf.is_monic()) {
      std::cerr << "The minimal polynomial of " << eFile
                << " is not monic, so there is no ring to compare with\n";
      throw TerminalException{1};
    }
    insert_helper_real_algebraic_field(idx_bench_field, hcrf);
    int deg = hcrf.deg;
#ifdef REALFIELD_USE_FLINT
    std::string arith = "flint";
#else
    std::string arith = "gmp";
#endif
    std::cout << "Field " << eFile << " of degree " << deg << ", matrices "
              << n << " x " << n << ", best of " << reps << " runs, "
              << arith << " arithmetic\n";
    Tfield det_f = run<Tfield>("RealField", n, deg, reps);
    Tring det_r = run<Tring>("RealRing ", n, deg, reps);
    // Both runs must produce the same determinant.
    if (det_f != UniversalScalarConversion<Tfield, Tring>(det_r)) {
      std::cerr << "The field and the ring disagree on the determinant\n";
      std::cerr << "det_f=" << det_f << " det_r=" << det_r << "\n";
      throw TerminalException{1};
    }
    std::cout << "The two runs agree on the determinant\n";
    std::cout << "Normal termination of Bench_real_ring\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Bench_real_ring\n";
    exit(e.eVal);
  }
}
