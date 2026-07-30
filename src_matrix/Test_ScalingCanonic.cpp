// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

// Tests of the fraction removal and the scalar canonicalization:
// --- RemoveFractionVectorPlusCoeff / RemoveFractionMatrixPlusCoeff: the
//     output is TheMult times the input, integral, of content one, with
//     a positive multiplier.
// --- NonUniqueScaleToIntegerVector: integral positive multiple.
// --- CanonicalizeVector: equal on positive multiples of a vector and
//     idempotent.
// --- ScalarCanonicalizationMatrix: equal on positive multiples.

template <typename T> MyVector<T> RandomFractionalVector(int n) {
  MyVector<T> V(n);
  for (int i = 0; i < n; i++) {
    T num = T((random() % 21) - 10);
    T den = T(1 + (random() % 8));
    V(i) = num / den;
  }
  return V;
}

template <typename T> MyMatrix<T> RandomFractionalMatrix(int n_row, int n_col) {
  MyMatrix<T> M(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++) {
      T num = T((random() % 21) - 10);
      T den = T(1 + (random() % 8));
      M(i, j) = num / den;
    }
  return M;
}

template <typename T> T RandomPositiveScalar() {
  T num = T(1 + (random() % 9));
  T den = T(1 + (random() % 9));
  return num / den;
}

// The content of an integral vector, as computed over the field type.
template <typename T> T VectorContent(MyVector<T> const &V) {
  T eGCD(0);
  for (int i = 0; i < V.size(); i++)
    eGCD = GcdPair(eGCD, V(i));
  return T_abs(eGCD);
}

template <typename T> void process_vector(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    MyVector<T> V = RandomFractionalVector<T>(n);
    if (IsZeroVector(V))
      continue;
    FractionVector<T> fr = RemoveFractionVectorPlusCoeff(V);
    MyVector<T> scaled = fr.TheMult * V;
    if (scaled != fr.TheVect || !IsIntegralVector(fr.TheVect)) {
      std::cerr << "RemoveFractionVectorPlusCoeff is inconsistent\n";
      WriteVector(std::cerr, V);
      throw TerminalException{1};
    }
    if (fr.TheMult <= 0 || VectorContent(fr.TheVect) != 1) {
      std::cerr << "RemoveFractionVectorPlusCoeff: wrong multiplier or "
                   "content\n";
      WriteVector(std::cerr, V);
      throw TerminalException{1};
    }
    MyVector<T> Vint = NonUniqueScaleToIntegerVector(V);
    if (!IsIntegralVector(Vint)) {
      std::cerr << "NonUniqueScaleToIntegerVector is not integral\n";
      throw TerminalException{1};
    }
    // The canonical vector does not depend on positive rescaling and is
    // idempotent.
    T lambda = RandomPositiveScalar<T>();
    MyVector<T> Vscal = lambda * V;
    MyVector<T> can1 = CanonicalizeVector(V);
    MyVector<T> can2 = CanonicalizeVector(Vscal);
    MyVector<T> can3 = CanonicalizeVector(can1);
    if (can1 != can2 || can1 != can3) {
      std::cerr << "CanonicalizeVector is not canonical\n";
      WriteVector(std::cerr, V);
      throw TerminalException{1};
    }
  }
}

template <typename T> void process_matrix(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    MyMatrix<T> M = RandomFractionalMatrix<T>(n + 1, n);
    FractionMatrix<T> fr = RemoveFractionMatrixPlusCoeff(M);
    MyMatrix<T> scaled = fr.TheMult * M;
    if (!TestEqualityMatrix(scaled, fr.TheMat) ||
        !IsIntegralMatrix(fr.TheMat)) {
      std::cerr << "RemoveFractionMatrixPlusCoeff is inconsistent\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    if (fr.TheMult <= 0) {
      std::cerr << "RemoveFractionMatrixPlusCoeff: negative multiplier\n";
      throw TerminalException{1};
    }
    T content(0);
    for (int i_row = 0; i_row <= n; i_row++)
      for (int i_col = 0; i_col < n; i_col++)
        content = GcdPair(content, fr.TheMat(i_row, i_col));
    if (T_abs(content) != 1) {
      std::cerr << "RemoveFractionMatrixPlusCoeff: wrong content\n";
      throw TerminalException{1};
    }
    // The scalar canonicalization does not depend on positive rescaling.
    T lambda = RandomPositiveScalar<T>();
    MyMatrix<T> Mscal = lambda * M;
    MyMatrix<T> can1 = ScalarCanonicalizationMatrix(M);
    MyMatrix<T> can2 = ScalarCanonicalizationMatrix(Mscal);
    if (!TestEqualityMatrix(can1, can2)) {
      std::cerr << "ScalarCanonicalizationMatrix is not canonical\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
  }
}

template <typename T> void process(int n) {
  int nb = 100;
  process_vector<T>(n, nb);
  std::cerr << "process_vector done\n";
  process_matrix<T>(n, nb);
  std::cerr << "process_matrix done\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_ScalingCanonic [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f = [&]() -> void {
      if (arith == "mpq_class")
        return process<mpq_class>(n);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(n);
      if (arith == "boost_cpp_rational")
        return process<boost::multiprecision::cpp_rational>(n);
      if (arith == "boost_mpq_rational")
        return process<boost::multiprecision::mpq_rational>(n);
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_ScalingCanonic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_ScalingCanonic\n";
    exit(e.eVal);
  }
  runtime(time);
}
