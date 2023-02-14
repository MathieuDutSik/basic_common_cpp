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
  MyVector<T> V(n_cols);
  for (size_t i_col = 0; i_col < n_cols; i_col++)
    V(i_col) = 1;
  ContainerMatrix<T> Cont(M);
  auto test_fct=[&](int case_call, MyVector<T> const& W) -> void {
    std::optional<size_t> opt = Cont.GetIdx_v(V);
    if (opt) {
      std::cerr << "case_call=" << case_call << " pos=" << *opt << "\n";
    } else {
      std::cerr << "case_call=" << case_call << " pos=missing\n";
    }
  };
  test_fct(-1, V);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      V(i_col) = M(i_row, i_col);
    test_fct(i_row, V);
  }
}
