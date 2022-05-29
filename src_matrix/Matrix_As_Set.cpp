// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "MAT_Matrix.h"

int main() {
  size_t n_rows = 10;
  size_t n_cols = 5;
  using T = int;
  MyMatrix<T> M(n_rows, n_cols);
  for (size_t i_row = 0; i_row < n_rows; i_row++)
    for (size_t i_col = 0; i_col < n_cols; i_col++) {
      T val = random() % 100;
      M(i_row, i_col) = val;
    }

  MyMatrix<T> VectorContain(1, n_cols);
  for (size_t i_col = 0; i_col < n_cols; i_col++)
    VectorContain(0, i_col) = 1;
  ContainerMatrix<T> Cont(M, VectorContain);
  //  Cont.SetPtr(&VectorContain);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    std::pair<bool, size_t> pair = Cont.GetIdx();
    std::cerr << "i_row=" << i_row << " pair.first=" << pair.first << " / "
              << pair.second << "\n";
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      VectorContain(0, i_col) = M(i_row, i_col);
  }
}
