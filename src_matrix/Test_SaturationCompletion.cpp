// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

// Tests of the repetitive solvers and the saturation machinery:
// --- SolutionMatRepetitive: agreement with SolutionMat on solvable and
//     unsolvable right hand sides.
// --- RecSolutionIntMat: agreement with SolutionIntMat on integral
//     combinations, scaled combinations and off-lattice vectors, and
//     the matrix interfaces get_solution_m / is_containing_m.
// --- IntegralSpaceSaturation: same rational span, and every integral
//     vector of the rational span is an integral combination.
// --- SubspaceCompletionInt: the completion of a saturated subspace to
//     a unimodular basis.

template <typename T> MyMatrix<T> RandomMatrix(int n_row, int n_col, int amp) {
  MyMatrix<T> A(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++)
      A(i, j) = T((random() % (2 * amp + 1)) - amp);
  return A;
}

template <typename T> MyVector<T> RandomVector(int n, int amp) {
  MyVector<T> V(n);
  for (int i = 0; i < n; i++)
    V(i) = T((random() % (2 * amp + 1)) - amp);
  return V;
}

// Over the field types: SolutionMatRepetitive against SolutionMat.
template <typename T> void process_repetitive(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    int n_col = n + 1 + (i % 2);
    MyMatrix<T> Basis(n, n_col);
    while (true) {
      Basis = RandomMatrix<T>(n, n_col, 5);
      if (RankMat(Basis) == n)
        break;
    }
    SolutionMatRepetitive<T> smr(Basis);
    for (int k = 0; k < 8; k++) {
      MyVector<T> V(n_col);
      if (k % 2 == 0) {
        MyVector<T> x = RandomVector<T>(n, 4);
        V = Basis.transpose() * x;
      } else {
        V = RandomVector<T>(n_col, 5);
      }
      std::optional<MyVector<T>> opt_rep = smr.GetSolution(V);
      std::optional<MyVector<T>> opt_dir = SolutionMat(Basis, V);
      if (opt_rep.has_value() != opt_dir.has_value()) {
        std::cerr << "SolutionMatRepetitive disagrees with SolutionMat\n";
        WriteMatrix(std::cerr, Basis);
        throw TerminalException{1};
      }
      if (opt_rep) {
        MyVector<T> V_check = Basis.transpose() * (*opt_rep);
        if (V_check != V) {
          std::cerr << "The solution of SolutionMatRepetitive is incorrect\n";
          WriteMatrix(std::cerr, Basis);
          throw TerminalException{1};
        }
      }
    }
  }
}

// Over the ring types: RecSolutionIntMat against SolutionIntMat.
template <typename T> void process_rec_int(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    int n_row = n + (i % 3) - 1;
    MyMatrix<T> A = RandomMatrix<T>(n_row, n, 4);
    RecSolutionIntMat<T> rec(A);
    for (int k = 0; k < 8; k++) {
      MyVector<T> V(n);
      if (k % 3 == 0) {
        MyVector<T> x = RandomVector<T>(n_row, 3);
        V = A.transpose() * x;
      } else {
        V = RandomVector<T>(n, 5);
      }
      std::optional<MyVector<T>> opt_dir = SolutionIntMat(A, V);
      bool has_rec = rec.has_solution_v(V);
      std::optional<MyVector<T>> opt_rec = rec.get_solution_v(V);
      if (opt_dir.has_value() != has_rec ||
          opt_dir.has_value() != opt_rec.has_value()) {
        std::cerr << "RecSolutionIntMat disagrees with SolutionIntMat\n";
        WriteMatrix(std::cerr, A);
        throw TerminalException{1};
      }
      if (opt_rec) {
        MyVector<T> V_check = A.transpose() * (*opt_rec);
        if (V_check != V) {
          std::cerr << "The solution of RecSolutionIntMat is incorrect\n";
          WriteMatrix(std::cerr, A);
          throw TerminalException{1};
        }
      }
    }
    // The matrix interfaces on a stack of integral combinations.
    MyMatrix<T> Xs = RandomMatrix<T>(3, n_row, 3);
    MyMatrix<T> Vs = Xs * A;
    if (!rec.is_containing_m(Vs)) {
      std::cerr << "is_containing_m rejected contained vectors\n";
      throw TerminalException{1};
    }
    std::optional<MyMatrix<T>> opt_m = rec.get_solution_m(Vs);
    if (!opt_m) {
      std::cerr << "get_solution_m failed on contained vectors\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Vs_check = (*opt_m) * A;
    if (!TestEqualityMatrix(Vs_check, Vs)) {
      std::cerr << "The solution of get_solution_m is incorrect\n";
      throw TerminalException{1};
    }
  }
}

// Over the ring types: saturation and completion.
template <typename T> void process_saturation(int n, int nb) {
  using Tfield = typename overlying_field<T>::field_type;
  for (int i = 0; i < nb; i++) {
    int r = 1 + (random() % (n - 1));
    MyMatrix<T> M(r, n);
    while (true) {
      M = RandomMatrix<T>(r, n, 4);
      if (RankMat(M) == r)
        break;
    }
    MyMatrix<T> S = IntegralSpaceSaturation(M);
    if (RankMat(S) != r || S.rows() != r) {
      std::cerr << "The saturation has the wrong rank\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    // The same rational span.
    MyMatrix<Tfield> Mf = UniversalMatrixConversion<Tfield, T>(M);
    MyMatrix<Tfield> Sf = UniversalMatrixConversion<Tfield, T>(S);
    if (!TestEqualitySpannedSpaces(Mf, Sf)) {
      std::cerr << "The saturation does not span the same space\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    // Saturated: an integral vector of the rational span is an integral
    // combination of the saturation basis.
    MyVector<Tfield> wf = ZeroVector<Tfield>(n);
    for (int k = 0; k < r; k++) {
      Tfield num = Tfield((random() % 9) - 4);
      Tfield den = Tfield(1 + (random() % 4));
      for (int j = 0; j < n; j++)
        wf(j) += (num / den) * Sf(k, j);
    }
    if (!IsZeroVector(wf)) {
      MyVector<Tfield> wf_int = RemoveFractionVector(wf);
      MyVector<T> w = UniversalVectorConversion<T, Tfield>(wf_int);
      if (!SolutionIntMat(S, w)) {
        std::cerr << "The saturation basis is not saturated\n";
        WriteMatrix(std::cerr, M);
        throw TerminalException{1};
      }
    }
    // The completion of the saturated subspace to a unimodular basis.
    MyMatrix<T> C = SubspaceCompletionInt(S, n);
    if (C.rows() != n - r || C.cols() != n) {
      std::cerr << "The completion has the wrong dimensions\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Full(n, n);
    for (int k = 0; k < r; k++)
      Full.row(k) = S.row(k);
    for (int k = 0; k < n - r; k++)
      Full.row(r + k) = C.row(k);
    T det = DeterminantMat(Full);
    if (det != 1 && det != -1) {
      std::cerr << "The completed basis is not unimodular, det=" << det
                << "\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
  }
}

template <typename T> void process(int n) {
  int nb = 30;
  if constexpr (is_ring_field<T>::value) {
    process_repetitive<T>(n, nb);
    std::cerr << "process_repetitive done\n";
  } else {
    process_rec_int<T>(n, nb);
    std::cerr << "process_rec_int done\n";
    process_saturation<T>(n, nb);
    std::cerr << "process_saturation done\n";
  }
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_SaturationCompletion [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(n);
      if (arith == "safe_integer")
        return process<SafeInt64>(n);
      if (arith == "boost_cpp_int")
        return process<boost::multiprecision::cpp_int>(n);
      if (arith == "mpq_class")
        return process<mpq_class>(n);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(n);
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_SaturationCompletion\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_SaturationCompletion\n";
    exit(e.eVal);
  }
  runtime(time);
}
