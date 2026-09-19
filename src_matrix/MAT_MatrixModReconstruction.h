// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXMODRECONSTRUCTION_H_
#define SRC_MATRIX_MAT_MATRIXMODRECONSTRUCTION_H_

/*
  Multiple prime reconstruction of integral computations.

  The scheme is always the same one:
  ---We compute an a priori upper bound on the size of the result
     (a Hadamard style bound on the determinant).
  ---We compute the result modulo a sequence of primes p1, p2, ....
  ---We combine the modular results (by the Chinese Remainder Theorem
     for the determinant, by accumulating the product of the primes
     giving the same answer for the row/column selection) until the
     product of the primes used is large enough for the answer to be
     determined without ambiguity.
 */

// clang-format off
#include "MAT_MatrixMod.h"
#include <limits>
#include <map>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_MOD_RECONSTRUCTION
#endif

struct SelectionRowColData {
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
  size_t TheRank;

  bool operator<(const SelectionRowColData &other) const {
    if (TheRank != other.TheRank)
      return TheRank < other.TheRank;
    if (ListColSelect != other.ListColSelect)
      return ListColSelect < other.ListColSelect;
    return ListRowSelect < other.ListRowSelect;
  }
};

// We must have n_row >= n_col.
template <typename T>
T HadamardUpperBoundRectangularKernel(MyMatrix<T> const &TheMat) {
  int n_row = TheMat.rows();
  int n_col = TheMat.cols();

  if (n_row == 0 || n_col == 0) {
    return T(1);
  }

  std::multimap<T, int> row_norms_squared; // Use multimap to handle duplicates

  for (int i = 0; i < n_row; i++) {
    T row_norm_squared(0);
    for (int j = 0; j < n_col; j++) {
      T val = TheMat(i, j);
      row_norm_squared += val * val;
    }
    // #ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    //     std::cerr << "DETHADAMARD: i=" << i << " row_norm_squared=" <<
    //     row_norm_squared << "\n";
    // #endif
    row_norms_squared.insert({row_norm_squared, i});
  }

  // Take product of the n_col largest row norms squared
  T bound(1);
  auto it = row_norms_squared.rbegin(); // Start from largest
  for (int count = 0; count < n_col && it != row_norms_squared.rend();
       ++count, ++it) {
    T val = it->first;
    // #ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    //     std::cerr << "count=" << count << " val=" << val << "\n";
    // #endif
    if (val != 0) {
      bound *= val;
    }
  }

  return bound;
}

template <typename T>
T HadamardUpperBoundRectangular(MyMatrix<T> const &TheMat) {
  int n_row = TheMat.rows();
  int n_col = TheMat.cols();

  if (n_row >= n_col) {
    return HadamardUpperBoundRectangularKernel(TheMat);
  } else {
    MyMatrix<T> TheMatTr = TransposedMat(TheMat);
    return HadamardUpperBoundRectangularKernel(TheMatTr);
  }
}

template <typename T> T SquareHadamardUpperBound(MyMatrix<T> const &TheMat) {
  int n = TheMat.rows();
  if (n != TheMat.cols()) {
    std::cerr << "SquareHadamardUpperBound: Matrix must be square\n";
    throw TerminalException{1};
  }

  if (n == 0) {
    return T(1);
  }

  T bound = T(1);

  for (int i = 0; i < n; i++) {
    T row_norm_squared(0);

    // Compute squared Euclidean norm of row i
    for (int j = 0; j < n; j++) {
      T val = TheMat(i, j);
      row_norm_squared += val * val;
    }

    // Multiply to the bound
    bound *= row_norm_squared;
  }

  return bound;
}

template <typename T> T DeterminantMatHadamard(MyMatrix<T> const &TheMat) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");

  int n = TheMat.rows();
  if (n != TheMat.cols()) {
    std::cerr << "DeterminantMatHadamard: Matrix must be square\n";
    throw TerminalException{1};
  }

  if (n == 0) {
    return T(1);
  }

  if (n == 1) {
    return TheMat(0, 0);
  }

  // Step 1: Compute Hadamard upper bound
  T bound = SquareHadamardUpperBound(TheMat);
  T target_product = T(9) * bound;

#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
  std::cerr << "DETHADAMARD: Hadamard bound = " << bound << "\n";
  std::cerr << "DETHADAMARD: Target product = " << target_product << "\n";
#endif

  // Step 2: Generate primes and compute determinants modulo primes
  PrimeGenerator<T> prime_gen;
  std::vector<T> primes;
  std::vector<T> det_mods;
  T product(1);

  while (product * product < target_product) {
    T prime = prime_gen.get_prime();
    primes.push_back(prime);

    T det_mod = DeterminantMatMod(TheMat, prime);
    det_mods.push_back(det_mod);

    product *= prime;

#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    std::cerr << "DETHADAMARD: Prime = " << prime << ", det mod = " << det_mod
              << ", product = " << product << "\n";
#endif
  }

#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
  std::cerr << "DETHADAMARD: Used " << primes.size() << " primes\n";
  std::cerr << "DETHADAMARD: Final product = " << product << "\n";
#endif

  // Step 3: Use Chinese Remainder Theorem to find determinant mod product
  T det_crt = chinese_remainder_theorem(det_mods, primes);

  // Step 4: Find the value nearest to 0
  // If det_crt > product/2, then det_crt - product is closer to 0
  if (2 * det_crt > product) {
    det_crt = det_crt - product;
  }

#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
  std::cerr << "DETHADAMARD: Final determinant = " << det_crt << "\n";
#endif

  return det_crt;
}

template <typename T>
SelectionRowColData SelectRowColDataMatMod_inner(MyMatrix<T> const &TheMat,
                                                 T const &TheMod) {
  SelectionRowCol<T> result = SelectRowColMatMod(TheMat, TheMod);

  SelectionRowColData data;
  data.ListColSelect = result.ListColSelect;
  data.ListRowSelect = result.ListRowSelect;
  data.TheRank = result.TheRank;

  return data;
}

template <typename T>
SelectionRowColData SelectRowColDataMatMod(MyMatrix<T> const &TheMat,
                                           T const &TheMod) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");

  // Compute TheMod * TheMod for comparison
  T mod_squared = TheMod * TheMod;

  // Check if we can use int16_t
  int16_t val_16 = std::numeric_limits<int16_t>::max();
  T max_int16 = UniversalScalarConversion<T, int16_t>(val_16);
  if (mod_squared < max_int16) {
#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    std::cerr << "SELECTMOD: Using int16_t optimization\n";
#endif
    int16_t mod_small = UniversalScalarConversion<int16_t, T>(TheMod);
    MyMatrix<int16_t> mat_small = UniversalMatrixConversion<int16_t, T>(TheMat);
    return SelectRowColDataMatMod_inner(mat_small, mod_small);
  }

  // Check if we can use int32_t
  int32_t val_32 = std::numeric_limits<int32_t>::max();
  T max_int32 = UniversalScalarConversion<T, int32_t>(val_32);
  if (mod_squared < max_int32) {
#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    std::cerr << "SELECTMOD: Using int32_t optimization\n";
#endif
    int32_t mod_small = UniversalScalarConversion<int32_t, T>(TheMod);
    MyMatrix<int32_t> mat_small = UniversalMatrixConversion<int32_t, T>(TheMat);
    return SelectRowColDataMatMod_inner(mat_small, mod_small);
  }

  // Check if we can use int64_t
  int64_t val_64 = std::numeric_limits<int64_t>::max();
  T max_int64 = UniversalScalarConversion<T, int64_t>(val_64);
  if (mod_squared < max_int64) {
#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    std::cerr << "SELECTMOD: Using int64_t optimization\n";
#endif
    int64_t mod_small = UniversalScalarConversion<int64_t, T>(TheMod);
    MyMatrix<int64_t> mat_small = UniversalMatrixConversion<int64_t, T>(TheMat);
    return SelectRowColDataMatMod_inner(mat_small, mod_small);
  }

  // No optimization possible, use original type
#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
  std::cerr << "SELECTMOD: No optimization, using original type\n";
#endif
  return SelectRowColDataMatMod_inner(TheMat, TheMod);
}

template <typename T>
SelectionRowColData SelectRowColMatModHadamard(MyMatrix<T> const &TheMat) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");

  int n = TheMat.rows();
  int m = TheMat.cols();

  std::map<SelectionRowColData, T> results;

  if (n == 0 || m == 0) {
    SelectionRowColData data;
    data.TheRank = 0;
    return data;
  }

  // Step 1: Compute Hadamard upper bound
  T bound = HadamardUpperBoundRectangular(TheMat);

  T target_product = T(9) * bound;

#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
  std::cerr << "SELECTHADAMARD: Hadamard bound = " << bound << "\n";
  std::cerr << "SELECTHADAMARD: Target product = " << target_product << "\n";
#endif

  // Step 2: Generate primes and compute SelectRowColMatMod results
  PrimeGenerator<T> prime_gen;

  while (true) {
    T prime = prime_gen.get_prime();

#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    std::cerr << "prime=" << prime << " before\n";
#endif
    SelectionRowColData data = SelectRowColDataMatMod(TheMat, prime);
#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    std::cerr << "prime=" << prime << " after\n";
#endif
    T &value = results[data];

    // Update the product for this data entry
    if (value == 0) {
      value = prime;
    } else {
      value *= prime;
    }

#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
    std::cerr << "SELECTHADAMARD: Prime = " << prime
              << ", rank = " << data.TheRank << ", product = " << value << "\n";
#endif

    // Step 4: Check if this entry has product larger than Hadamard bound
    if (value * value > target_product) {
#ifdef DEBUG_MATRIX_MOD_RECONSTRUCTION
      std::cerr << "SELECTHADAMARD: Product " << value << " exceeds target "
                << target_product << "\n";
      std::cerr << "SELECTHADAMARD: Returning first qualifying entry\n";
#endif
      return data;
    }
  }
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXMODRECONSTRUCTION_H_
// clang-format on
