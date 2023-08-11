// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "Fp.h"
#include "MAT_Matrix.h"
// clang-format on


template<typename T>
std::string full_process_type(MyMatrix<int> const& M) {
  int n_row = M.rows();
  int n_col = M.cols();
  MyMatrix<T> M_T(n_row,n_col);
  for (int i=0; i<n_row; i++) {
    for (int j=0; j<n_col; j++) {
      M_T(i,j) = M(i,j);
    }
  }
  //  MyMatrix<T> M_T = UniversalMatrixConversion<T,int>(M);
  MyMatrix<T> TheKer = NullspaceMat(M_T);
  //
  std::stringstream os;
  WriteMatrix(os, TheKer);
  std::string converted(os.str());
  return converted;
}


std::string process(std::string const& arith, MyMatrix<int> const& M) {
  /*
  if (arith == "Fp") {
    using T = Fp<long, 2147389441>;
    return full_process_type<T>(M);
  }
  */
  if (arith == "rational<SafeInt64>") {
    using T = Rational<SafeInt64>;
    return full_process_type<T>(M);
  }
  if (arith == "rational<long>") {
    using T = Rational<long>;
    return full_process_type<T>(M);
  }
  if (arith == "rational") {
    using T = mpq_class;
    return full_process_type<T>(M);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return full_process_type<T>(M);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return full_process_type<T>(M);
  }
  std::optional<std::string> opt_realalgebraic =
    get_postfix(arith, "RealAlgebraic=");
  if (opt_realalgebraic) {
    std::string const &FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    using T_rat = mpq_class;
    HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return full_process_type<T>(M);
  }
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}


void process_listm_listarith(std::vector<MyMatrix<int>> const& ListM, std::vector<std::string> const& ListArith) {
  for (auto & eM : ListM) {
    std::unordered_map<std::string, std::vector<std::string>> map;
    std::cerr << "|eM|=" << eM.rows() << " / " << eM.cols() << "\n";
    for (auto & arith : ListArith) {
      std::string e_str = process(arith, eM);
      map[e_str].push_back(arith);
    }
    if (map.size() != 1) {
      std::cerr << "inconsistency in the computation |map|=" << map.size() << "\n";
      for (auto & kv : map) {
        std::cerr << "Result=" << kv.first << " for arithmetics =";
        for (auto & arith : kv.second)
          std::cerr << " " << arith;
        std::cerr << "\n";
      }
    }
  }
  for (auto arith : ListArith) {
    HumanTime time;
    for (auto & eM : ListM)
      process(arith, eM);
    std::cerr << "Result for arithmetic arith=" << arith << " time=" << time << "\n";
  }
}


MyMatrix<int> get_random_matrix(int m, int n) {
  MyMatrix<int> M(m, n);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      int val = rand() % 11 - 5;
      M(i,j) = val;
    }
  }
  return M;
}



int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    //    std::vector<std::string> ListArith = {"Fp", "rational<SafeInt64>", "rational<long>", "rational", "Qsqrt5", "Qsqrt2"};
    std::vector<std::string> ListArith = {"rational<SafeInt64>", "rational<long>", "rational", "Qsqrt5", "Qsqrt2"};


    std::vector<MyMatrix<int>> ListM;
    auto insert=[&](int m, int n) -> void {
      ListM.push_back(get_random_matrix(m,n));
      ListM.push_back(get_random_matrix(n,m));
    };
    insert(4,3);
    insert(5,3);
    insert(6,2);
    insert(4,2);
    insert(6,7);
    /*
    insert(8,9);
    insert(10,11);
    insert(15,14);
    insert(20,21);
    insert(25,24);
    */
    //
    process_listm_listarith(ListM, ListArith);
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  }
  runtime(time1);
}
