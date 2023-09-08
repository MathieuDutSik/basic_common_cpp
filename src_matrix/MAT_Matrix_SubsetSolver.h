// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIX_FP_H_
#define SRC_MATRIX_MAT_MATRIX_FP_H_

#include "MAT_MatrixInt.h"
#include "Boost_bitset.h"
#include "Fp.h"

template<typename T>
struct SubsetRankOneSolver_Field {
public:
  using Tint = T;
  MyMatrix<T> const& EXT;
  int nbCol;

  SubsetRankOneSolver_Field(MyMatrix<Tint> const& _EXT) : EXT(_EXT), nbCol(EXT.cols()) {
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


// The acceleration scheme is using reduction to Fp.
// techniques for the computation of the Kernel.
//
// The type T should be an implementation of Q.
// Then Tint is some implementation of Z.
// We then do computation over a fast modulo ring Tfast.
//
// The scheme should be failsafe, that is not throw any
// more error than the type T.
template<typename T>
struct SubsetRankOneSolver_Acceleration {
public:
  using Tint = typename underlying_ring<T>::ring_type;
  using Tlift = int64_t;
  using Tfast = Fp<Tlift, 2147389441>;
  MyMatrix<Tint> const& EXT;
  int nbRow;
  int nbCol;
  std::vector<std::pair<Tlift,Tlift>> lifts;
  MyMatrix<Tfast> EXT_fast;
  MyMatrix<Tlift> EXT_lift;
  bool try_int;
  size_t max_bits;
  SubsetRankOneSolver_Acceleration(MyMatrix<Tint> const& _EXT) : EXT(_EXT), nbRow(EXT.rows()), nbCol(EXT.cols()), lifts(nbCol) {
    //
    // Faster modular version of EXT_red
    //
    max_bits = 0;
    EXT_fast = MyMatrix<Tfast>(nbRow, nbCol);
    EXT_lift = MyMatrix<Tlift>(nbRow, nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        Tint const& val = EXT(iRow, iCol);
        max_bits = std::max(get_bit(val), max_bits);
        EXT_lift(iRow, iCol) = UniversalScalarConversion<Tlift,Tint>(val);
        EXT_fast(iRow, iCol) = Tfast(EXT_lift(iRow, iCol));
      }
    }
    try_int = (max_bits <= 30);
    max_bits += get_bit(static_cast<int64_t>(nbCol));
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
        MyVector<Tlift> VZ_lift(nbCol);
        // reconstruct the vector
        size_t max_bits_NSP = 0;
        lifts[0] = Vzero_Tfast(0, 0).rational_lift();
        Tlift lcm = lifts[0].second;
        for (int iCol = 1; iCol < nbCol; iCol++) {
          lifts[iCol] = Vzero_Tfast(0, iCol).rational_lift();
          lcm = LCMpair(lcm, lifts[iCol].second);
        }
        for (int iCol = 0; iCol < nbCol; iCol++) {
          VZ_lift(iCol) = lifts[iCol].first * (lcm / lifts[iCol].second);
          Vkernel(iCol) = UniversalScalarConversion<Tint,Tlift>(VZ_lift(iCol));
          max_bits_NSP = std::max(max_bits_NSP, get_bit(VZ_lift(iCol)));
        }
        // check if elements are small enough to do computation in
        if (max_bits + max_bits_NSP <= 60) {
          // check if part of kernel
          jRow = sInc.find_first();
          for (size_t iRow = 0; iRow < nb; iRow++) {
            Tlift sm = 0;
            for (int iCol = 0; iCol < nbCol; iCol++) {
              sm += VZ_lift(iCol) * EXT_lift(jRow, iCol);
            }
            if (sm != 0) {
              failed_int = true;
              break;
            }
            jRow = sInc.find_next(jRow);
          }
        } else {
          failed_int = true;
        }
      }
    }

    if (failed_int || !try_int) {
      std::cerr << "Lifting strategy failed, retrying with mpq_class\n";
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<T> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        for (int iCol=0; iCol<nbCol; iCol++)
          M(eRank,iCol) = UniversalScalarConversion<T,Tint>(EXT(jRow,iCol));
        jRow = sInc.find_next(jRow);
      };
      Vkernel = NonUniqueRescaleVecRing(NullspaceTrMatTargetOne_Kernel<T, decltype(f)>(nb, nbCol, f));
    }
    return Vkernel;
  }
};


// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIX_FP_H_
// clang-format on
