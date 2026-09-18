// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXINTSPARSE_H_
#define SRC_MATRIX_MAT_MATRIXINTSPARSE_H_

// Integer normal forms taking their input as a MySparseMatrix.
//
// The motivation is the Smith normal form of the boundary matrices of a
// chain complex. Those are enormous and almost empty -- a 11568 x 12119
// case holds 81720 nonzero entries, 0.06% of its cells -- and the whole
// dense pipeline is proportional to the number of CELLS rather than to
// the number of entries: parsing the file, constructing that many number
// objects, scanning for the density, building the sparse image of the
// elimination. Measured on that matrix, the elimination itself takes
// 1.7 s out of a 7.9 s run, so about 78% of the time is spent handling a
// dense form of a matrix that is not dense.
//
// Feeding MySparseMatrix straight into the elimination removes all of it.
// What comes back out is a core that is normally far smaller than the
// input, and only that core is ever made dense for the backend.

// clang-format off
#include "MAT_MatrixInt.h"
#include "MAT_SparseMatrix.h"
#include <utility>
#include <vector>
// clang-format on

// The working structure of the unit-pivot elimination, filled directly
// from the sparse entries.
template <typename T>
SmithSparseWork<T> SmithSparseWorkFromSparse(MySparseMatrix<T> const &M) {
  SmithSparseWork<T> W(M.rows(), M.cols());
  for (int k = 0; k < M.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(M, k); it; ++it)
      if (it.value() != 0)
        W.set_entry(it.row(), it.col(), it.value());
  return W;
}

template <typename T>
MySparseMatrix<T> SmithCoreSparse(SmithSparseWork<T> const &W,
                                  SmithUnitPivotLayout<T> const &layout) {
  std::vector<int> col_pos(W.nbCol, -1);
  for (size_t idx = 0; idx < layout.keep_cols.size(); idx++)
    col_pos[layout.keep_cols[idx]] = idx;
  int n_row_core = layout.keep_rows.size();
  int n_col_core = layout.keep_cols.size();
  using Ttrip = Eigen::Triplet<T>;
  std::vector<Ttrip> tripletList;
  tripletList.reserve(W.nnz);
  for (int idx = 0; idx < n_row_core; idx++)
    for (auto &ent : W.rows[layout.keep_rows[idx]])
      tripletList.push_back(Ttrip(idx, col_pos[ent.first], ent.second));
  MySparseMatrix<T> core(n_row_core, n_col_core);
  core.setFromTriplets(tripletList.begin(), tripletList.end());
  return core;
}

template <typename T> struct SmithUnitPivotReductionSparse {
  std::vector<T> pivots;
  MySparseMatrix<T> core;
  size_t core_nnz;
};

template <typename T>
SmithUnitPivotReductionSparse<T>
SmithUnitPivotEliminate_sparse(MySparseMatrix<T> const &M, int nb_candidates,
                               double switch_density) {
  SmithSparseWork<T> W = SmithSparseWorkFromSparse(M);
  SmithUnitPivotLayout<T> layout =
      SmithUnitPivotEliminate_Work(W, nb_candidates, switch_density);
  size_t core_nnz = W.nnz;
  MySparseMatrix<T> core = SmithCoreSparse(W, layout);
  return {std::move(layout.pivots), std::move(core), core_nnz};
}

// The invariant factors of a sparse matrix.
//
// The unit-pivot pre-elimination runs on the sparse form, and only what
// it cannot consume is expanded for the backend. The result is the same
// vector as SmithNormalFormInvariant on the dense form of M: same
// algorithm, same invariant factors, the representation of the input
// being the only difference.
template <typename T>
MyVector<T> SmithNormalFormInvariant_sparse(MySparseMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  int min_dim = std::min(nbRow, nbCol);
  if (min_dim == 0)
    return MyVector<T>(0);
  if constexpr (use_unit_pivot_preelimination<T>::value) {
    SmithUnitPivotReductionSparse<T> red = SmithUnitPivotEliminate_sparse(
        M, smith_unit_pivot_nb_candidates, smith_unit_pivot_switch_density);
    if (!red.pivots.empty()) {
      // Only a core with something left in it reaches the backend, and
      // only then is it made dense.
      if (red.core_nnz == 0) {
        MyMatrix<T> empty(0, 0);
        return SmithInvariantFromCore(min_dim, red.pivots, empty, 0);
      }
      MyMatrix<T> core_dense = MyMatrixFromSparseMatrix(red.core);
      return SmithInvariantFromCore(min_dim, red.pivots, core_dense,
                                    red.core_nnz);
    }
  }
  // Nothing could be peeled off, so the backend gets the whole matrix.
  MyMatrix<T> M_dense = MyMatrixFromSparseMatrix(M);
  return SmithNormalFormInvariant_Kernel(M_dense);
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXINTSPARSE_H_
// clang-format on
