// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIX_FP_H_
#define SRC_MATRIX_MAT_MATRIX_FP_H_

#include "MAT_Matrix.h"
#include "Boost_bitset.h"
#include "Fp.h"

MyVector<mpz_class> RescaleVec(MyVector<mpq_class> const &v) {
  int cols = v.size();
  std::vector<mpz_class> dens(cols, 1);
  MyVector<mpz_class> vret(cols);
  for (int iCol = 0; iCol < cols; iCol++) {
    dens[iCol] = v(iCol).get_den();
  }
  mpz_class scale = LCMlist(dens);
  for (int iCol = 0; iCol < cols; iCol++) {
    vret(iCol) = (scale / v(iCol).get_den()) * v(iCol).get_num();
  }
  return vret;
}

MyMatrix<mpz_class> RescaleRows(MyMatrix<mpq_class> const &M) {
  int rows = M.rows();
  int cols = M.cols();
  std::vector<mpz_class> dens(cols, 1);
  MyMatrix<mpz_class> Mret(rows, cols);
  for (int iRow = 0; iRow < rows; iRow++) {
    for (int iCol = 0; iCol < cols; iCol++) {
      dens[iCol] = M(iRow, iCol).get_den();
    }
    mpz_class scale = LCMlist(dens);
    for (int iCol = 0; iCol < cols; iCol++) {
      Mret(iRow, iCol) =
        (scale / M(iRow, iCol).get_den()) * M(iRow, iCol).get_num();
    }
  }
  return Mret;
}

template<typename T>
struct SubsetRankOneSolver {
public:
  using Tint = T;
  MyMatrix<T> const& EXT;
  int nbCol;

  SubsetRankOneSolver(MyMatrix<Tint> const& _EXT) : EXT(_EXT) {
    nbCol = EXT.cols();
  }
  MyVector<Tint> GetKernelVector(Face const& sInc) {
    int nb = sInc.count();
    boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
    auto f = [&](MyMatrix<T> &M, size_t eRank,
                 [[maybe_unused]] size_t iRow) -> void {
      for (int iCol=0; iCol<nbCol; iCol++)
        M(eRank,iCol) = EXT(jRow,iCol);
      jRow = sInc.find_next(jRow);
    };
    return NullspaceTrMatTargetOne_Kernel<T, decltype(f)>(nb, nbCol, f);
  }
};



template<>
struct SubsetRankOneSolver<mpq_class> {
public:
  using T = mpq_class;
  using Tint = mpz_class;
  using Tfast = Fp<long, 2147389441>;
  MyMatrix<Tint> const& EXT;
  MyMatrix<Tfast> EXT_fast;
  MyMatrix<long> EXT_long;
  bool try_int;
  size_t max_bits;
  int nbRow;
  int nbCol;

  size_t get_bit(mpz_class const& v) const {
    return mpz_sizeinbase(v.get_mpz_t(), 2);
  }

  SubsetRankOneSolver(MyMatrix<Tint> const& _EXT) : EXT(_EXT) {
    nbRow = EXT.rows();
    nbCol = EXT.cols();
    //
    // Faster modular version of EXT_red
    //
    max_bits = 0;
    EXT_fast = MyMatrix<Tfast>(nbRow, nbCol);
    EXT_long = MyMatrix<long>(nbRow, nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        Tint const& val = EXT(iRow, iCol);
        max_bits = std::max(get_bit(val), max_bits);
        EXT_long(iRow, iCol) = val.get_si();
        EXT_fast(iRow, iCol) = Tfast(EXT_long(iRow, iCol));
      }
    }
    try_int = (max_bits <= 30);
    max_bits += get_bit(mpz_class(nbCol));
  }

  MyVector<Tint> GetKernelVector(Face const& sInc) {
    size_t nb = sInc.count();
    MyVector<Tint> Vkernel(nbCol);
    bool failed_int = false;
    if (try_int) {
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<Tfast> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        M.row(eRank) = EXT_fast.row(jRow);
        jRow = sInc.find_next(jRow);
      };
      MyVector<Tfast> Vzero_Tfast =
        NullspaceTrMatTargetOne_Kernel<Tfast, decltype(f)>(nb, nbCol, f);
      // check result at full precision in case of overflows
      bool allzero = true;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (Vzero_Tfast(iCol) != 0) {
          allzero = false;
          break;
        }
      }
      if (allzero) {
        failed_int = true;
      } else {
        MyVector<long> VZ_long(nbCol);
        // reconstruct the vector
        size_t max_bits_NSP = 0;
        std::vector<long> nums(nbCol, 0);
        std::vector<long> dens(nbCol, 1);
        for (int iCol = 0; iCol < nbCol; iCol++) {
          Rational<long> val = Vzero_Tfast(0, iCol).rational_lift();
          nums[iCol] = val.get_num();
          dens[iCol] = val.get_den();
        }
        long lcm = LCMlist(dens);
        for (int iCol = 0; iCol < nbCol; iCol++) {
          VZ_long(iCol) = nums[iCol] * (lcm / dens[iCol]);
          Vkernel(iCol) = Tint(VZ_long(iCol));
          max_bits_NSP = std::max(max_bits_NSP, get_bit(Vkernel(iCol)));
        }
        // check if elements are small enough to do computation in
        if (max_bits + max_bits_NSP <= 60) {
          // check if part of kernel
          jRow = sInc.find_first();
          for (size_t iRow = 0; iRow < nb; iRow++) {
            auto row = EXT_long.row(jRow);
            jRow = sInc.find_next(jRow);
            long sm = 0;
            for (int iCol = 0; iCol < nbCol; iCol++) {
              sm += VZ_long(iCol) * row(iCol);
            }
            if (sm != 0) {
              failed_int = true;
              break;
            }
          }
        } else {
          failed_int = true;
        }
      }
    }

    if (failed_int || !try_int) {
      std::cerr << "Rational<long> strategy failed, retrying with mpq_class\n";
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<T> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        for (int iCol=0; iCol<nbCol; iCol++)
          M(eRank,iCol) = UniversalScalarConversion<T,Tint>(EXT(jRow,iCol));
        jRow = sInc.find_next(jRow);
      };
      Vkernel = RescaleVec(NullspaceTrMatTargetOne_Kernel<T, decltype(f)>(nb, nbCol, f));
    }
    return Vkernel;
  }
};






// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIX_FP_H_
// clang-format on
