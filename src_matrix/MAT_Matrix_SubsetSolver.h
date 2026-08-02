// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIX_SUBSETSOLVER_H_
#define SRC_MATRIX_MAT_MATRIX_SUBSETSOLVER_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "Boost_bitset.h"
#include "Fp.h"
#include <algorithm>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_SUBSET_SOLVER
#endif

#ifdef DEBUG
#define DEBUG_MATRIX_SUBSET_SOLVER
#endif

// The sign fix of a kernel vector: the first row outside the subset
// with a non-zero scalar product decides the orientation, the rows
// lying on the hyperplane being skipped.
template <typename Tint>
void SubsetRankOneSolver_SignFix(MyMatrix<Tint> const &EXT, Face const &sInc,
                                 MyVector<Tint> &V) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (sInc[iRow] == 0) {
      Tint scal(0);
      for (int iCol = 0; iCol < nbCol; iCol++)
        scal += EXT(iRow, iCol) * V(iCol);
      if (scal > 0)
        return;
      if (scal < 0) {
        V = -V;
        return;
      }
    }
  }
  std::cerr << "SubsetRankOneSolver: no row has a non-zero scalar product\n";
  throw TerminalException{1};
}

// The sign fix together with the incidence of the hyperplane: a single
// pass computes all the scalar products, the zero ones forming the
// incidence and the first non-zero one deciding the orientation.
template <typename Tint>
Face SubsetRankOneSolver_SignFixAndFace(MyMatrix<Tint> const &EXT,
                                        MyVector<Tint> &V) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  Face f(nbRow);
  int sign = 0;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    Tint scal(0);
    for (int iCol = 0; iCol < nbCol; iCol++)
      scal += EXT(iRow, iCol) * V(iCol);
    if (scal == 0) {
      f[iRow] = 1;
    } else {
      if (sign == 0)
        sign = (scal > 0) ? 1 : -1;
    }
  }
  if (sign == -1)
    V = -V;
  return f;
}

template <typename T> struct SubsetRankOneSolver_Field {
public:
  using Tint = T;
  MyMatrix<T> const &EXT;
  int nbRow;
  int nbCol;

  SubsetRankOneSolver_Field(MyMatrix<Tint> const &_EXT)
      : EXT(_EXT), nbRow(EXT.rows()), nbCol(EXT.cols()) {}
  MyVector<Tint> GetKernelVector(Face const &sInc) {
    int nb = sInc.count();
    boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
    auto f = [&](MyMatrix<T> &M, size_t eRank,
                 [[maybe_unused]] size_t iRow) -> void {
      for (int iCol = 0; iCol < nbCol; iCol++)
        M(eRank, iCol) = EXT(jRow, iCol);
      jRow = sInc.find_next(jRow);
    };
    return NullspaceTrMatTargetOne_Kernel<T, decltype(f)>(nb, nbCol, f);
  }
  MyVector<Tint> GetPositiveKernelVector(Face const &sInc) {
    MyVector<Tint> V = GetKernelVector(sInc);
    SubsetRankOneSolver_SignFix(EXT, sInc, V);
    return V;
  }
  std::pair<MyVector<Tint>, Face> GetPositiveKernelVectorAndFace(Face const &sInc) {
    MyVector<Tint> V = GetKernelVector(sInc);
    Face f = SubsetRankOneSolver_SignFixAndFace(EXT, V);
    return {std::move(V), std::move(f)};
  }
};

// The acceleration scheme is using reduction to Fp
// techniques for the computation of the Kernel.
//
// The type T is an implementation of Q or of Z; Tint is the
// corresponding implementation of Z. The computation runs over a fast
// modulo ring Tfast and is verified by lifting; on failure the exact
// computation is used (over the field for a Q type, over the ring for
// a Z type), so the scheme does not throw any more errors than the
// exact types themselves.
template <typename T> struct SubsetRankOneSolver_Acceleration {
public:
  using Tint = typename underlying_ring<T>::ring_type;
  using Tlift = int64_t;
  using Tfast = Fp<Tlift, 2147389441>;
  MyMatrix<Tint> const &EXT;
  int nbRow;
  int nbCol;
  std::vector<std::pair<Tlift, Tlift>> lifts;
  MyMatrix<Tfast> EXT_fast;
  MyMatrix<Tlift> EXT_lift;
  bool try_int;
  size_t max_bits;
  // The matrix determines the lifting regime: when the kernel entries
  // exceed the reconstruction range of the prime, every attempt fails.
  // After several consecutive failures the attempts stop, so only the
  // exact computation is paid.
  static const size_t max_consecutive_fails = 5;
  size_t n_consecutive_fail = 0;
  SubsetRankOneSolver_Acceleration(MyMatrix<Tint> const &_EXT)
      : EXT(_EXT), nbRow(EXT.rows()), nbCol(EXT.cols()), lifts(nbCol) {
    //
    // Faster modular version of EXT_red
    //
    max_bits = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 0; iCol < nbCol; iCol++)
        max_bits = std::max(get_bit(EXT(iRow, iCol)), max_bits);
    try_int = (max_bits <= 30);
    if (try_int) {
      EXT_fast = MyMatrix<Tfast>(nbRow, nbCol);
      EXT_lift = MyMatrix<Tlift>(nbRow, nbCol);
      for (int iRow = 0; iRow < nbRow; iRow++) {
        for (int iCol = 0; iCol < nbCol; iCol++) {
          EXT_lift(iRow, iCol) =
              UniversalScalarConversion<Tlift, Tint>(EXT(iRow, iCol));
          EXT_fast(iRow, iCol) = Tfast(EXT_lift(iRow, iCol));
        }
      }
    }
    max_bits += get_bit(static_cast<int64_t>(nbCol));
  }

private:
  // The kernel vector together with its lifted representation when the
  // fast integer path succeeded: the successful lift comes with the
  // bit-size guarantee max_bits + max_bits_NSP <= 60, so the sign and
  // incidence post-processing can run on int64 scalar products instead
  // of the expensive Tint ones. The exact arithmetic is used only in
  // the failing scenario.
  struct KernelVectorLift {
    MyVector<Tint> V;
    bool lift_valid;
    MyVector<Tlift> VZ_lift;
  };
  KernelVectorLift GetKernelVectorLift(Face const &sInc) {
    size_t nb = sInc.count();
    MyVector<Tint> Vkernel(nbCol);
    MyVector<Tlift> VZ_lift(nbCol);
    bool attempt_lift = try_int && n_consecutive_fail < max_consecutive_fails;
    bool failed_int = false;
    if (attempt_lift) {
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
        // reconstruct the vector
        size_t max_bits_NSP = 0;
        lifts[0] = Vzero_Tfast(0, 0).rational_lift();
        Tlift lcm = lifts[0].second;
        for (int iCol = 1; iCol < nbCol; iCol++) {
          lifts[iCol] = Vzero_Tfast(iCol).rational_lift();
          lcm = LCMpair(lcm, lifts[iCol].second);
        }
        for (int iCol = 0; iCol < nbCol; iCol++) {
          VZ_lift(iCol) = lifts[iCol].first * (lcm / lifts[iCol].second);
          Vkernel(iCol) = UniversalScalarConversion<Tint, Tlift>(VZ_lift(iCol));
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

    if (failed_int || !attempt_lift) {
      if (failed_int)
        n_consecutive_fail++;
#ifdef DEBUG_MATRIX_SUBSET_SOLVER
      std::cerr << "Lifting strategy failed, retrying with the exact "
                   "computation\n";
#endif
      using Tfield = typename overlying_field<Tint>::field_type;
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<Tfield> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        for (int iCol = 0; iCol < nbCol; iCol++)
          M(eRank, iCol) =
              UniversalScalarConversion<Tfield, Tint>(EXT(jRow, iCol));
        jRow = sInc.find_next(jRow);
      };
      Vkernel = NonUniqueRescaleVecRing(
          NullspaceTrMatTargetOne_Kernel<Tfield, decltype(f)>(nb, nbCol, f));
      return {std::move(Vkernel), false, std::move(VZ_lift)};
    }
    n_consecutive_fail = 0;
    return {std::move(Vkernel), true, std::move(VZ_lift)};
  }

public:
  MyVector<Tint> GetKernelVector(Face const &sInc) {
    return GetKernelVectorLift(sInc).V;
  }
  MyVector<Tint> GetPositiveKernelVector(Face const &sInc) {
    KernelVectorLift res = GetKernelVectorLift(sInc);
    if (res.lift_valid) {
      for (int iRow = 0; iRow < nbRow; iRow++) {
        if (sInc[iRow] == 0) {
          Tlift sm = 0;
          for (int iCol = 0; iCol < nbCol; iCol++)
            sm += res.VZ_lift(iCol) * EXT_lift(iRow, iCol);
          if (sm > 0)
            return res.V;
          if (sm < 0)
            return -res.V;
        }
      }
      std::cerr << "SubsetRankOneSolver: no row has a non-zero scalar "
                   "product\n";
      throw TerminalException{1};
    }
    SubsetRankOneSolver_SignFix(EXT, sInc, res.V);
    return res.V;
  }
  std::pair<MyVector<Tint>, Face> GetPositiveKernelVectorAndFace(Face const &sInc) {
    KernelVectorLift res = GetKernelVectorLift(sInc);
    if (res.lift_valid) {
      Face f(nbRow);
      int sign = 0;
      for (int iRow = 0; iRow < nbRow; iRow++) {
        Tlift sm = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          sm += res.VZ_lift(iCol) * EXT_lift(iRow, iCol);
        if (sm == 0) {
          f[iRow] = 1;
        } else {
          if (sign == 0)
            sign = (sm > 0) ? 1 : -1;
        }
      }
      if (sign == -1)
        res.V = -res.V;
      return {std::move(res.V), std::move(f)};
    }
    Face f = SubsetRankOneSolver_SignFixAndFace(EXT, res.V);
    return {std::move(res.V), std::move(f)};
  }
};

// The exact kernel vector of the corank one subset over a euclidean
// ring: the saturated integral kernel of the selected rows
// (NullspaceIntTrMat), with a content one output. Used by the ring
// variant.
template <typename Tint>
MyVector<Tint> SubsetRankOneSolver_KernelRing(MyMatrix<Tint> const &EXT,
                                              Face const &sInc) {
  int nb = sInc.count();
  int nbCol = EXT.cols();
  MyMatrix<Tint> Msel(nb, nbCol);
  boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
  for (int i = 0; i < nb; i++) {
    Msel.row(i) = EXT.row(jRow);
    jRow = sInc.find_next(jRow);
  }
  MyMatrix<Tint> NSP = NullspaceIntTrMat(Msel);
#ifdef SANITY_CHECK_MATRIX_SUBSET_SOLVER
  if (NSP.rows() != 1) {
    std::cerr << "The subset does not have corank one: |NSP|=" << NSP.rows()
              << "\n";
    throw TerminalException{1};
  }
#endif
  return GetMatrixRow(NSP, 0);
}

// The ring variant: everything over the ring itself, no conversion to
// the overlying field. It serves the euclidean rings that do not have
// the Fp acceleration.
template <typename T> struct SubsetRankOneSolver_Ring {
public:
  using Tint = T;
  MyMatrix<T> const &EXT;
  int nbRow;
  int nbCol;

  SubsetRankOneSolver_Ring(MyMatrix<Tint> const &_EXT)
      : EXT(_EXT), nbRow(EXT.rows()), nbCol(EXT.cols()) {}
  MyVector<Tint> GetKernelVector(Face const &sInc) {
    return SubsetRankOneSolver_KernelRing(EXT, sInc);
  }
  MyVector<Tint> GetPositiveKernelVector(Face const &sInc) {
    MyVector<Tint> V = GetKernelVector(sInc);
    SubsetRankOneSolver_SignFix(EXT, sInc, V);
    return V;
  }
  std::pair<MyVector<Tint>, Face> GetPositiveKernelVectorAndFace(Face const &sInc) {
    MyVector<Tint> V = GetKernelVector(sInc);
    Face f = SubsetRankOneSolver_SignFixAndFace(EXT, V);
    return {std::move(V), std::move(f)};
  }
};

template <typename T>
using subsetsolver_type = std::conditional_t<
    has_reduction_subset_solver<T>::value,
    SubsetRankOneSolver_Acceleration<T>,
    std::conditional_t<is_ring_field<T>::value, SubsetRankOneSolver_Field<T>,
                       SubsetRankOneSolver_Ring<T>>>;

// The public solver, for the repeated computation of corank one kernel
// vectors on a fixed matrix: the one-shot use case is FindFacetInequality.
// The constructor takes the matrix over T and does everything itself:
// the canonical per-row integer rescaling for the types with a
// reduction, so the Fp acceleration is available right away, and the
// dispatch to the fastest implementation (Fp accelerated for the
// rational types, exact ring for the euclidean rings, field otherwise).
template <typename T> class SubsetRankOneSolver {
  using T_solver = subsetsolver_type<T>;

public:
  using Tint = typename T_solver::Tint;

private:
  MyMatrix<Tint> EXT_int;
  T_solver subsetsolver;
  static MyMatrix<Tint> get_ext_int(MyMatrix<T> const &EXT) {
    if constexpr (has_reduction_subset_solver<T>::value) {
      return UniqueRescaleRowsRing(EXT);
    } else {
      return EXT;
    }
  }

  static MyVector<T> to_T(MyVector<Tint> const &V) {
    if constexpr (std::is_same_v<T, Tint>) {
      return V;
    } else {
      return UniversalVectorConversion<T, Tint>(V);
    }
  }

public:
  SubsetRankOneSolver(MyMatrix<T> const &EXT)
      : EXT_int(get_ext_int(EXT)), subsetsolver(EXT_int) {}
  MyMatrix<Tint> const &GetEXT_int() const { return EXT_int; }
  // The results are expressed over T, the internal Tint arithmetic being
  // an implementation matter of the variants.
  MyVector<T> GetKernelVector(Face const &sInc) {
    return to_T(subsetsolver.GetKernelVector(sInc));
  }
  MyVector<T> GetPositiveKernelVector(Face const &sInc) {
    return to_T(subsetsolver.GetPositiveKernelVector(sInc));
  }
  std::pair<MyVector<T>, Face> GetPositiveKernelVectorAndFace(Face const &sInc) {
    std::pair<MyVector<Tint>, Face> pair =
        subsetsolver.GetPositiveKernelVectorAndFace(sInc);
    return {to_T(pair.first), std::move(pair.second)};
  }
};

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIX_SUBSETSOLVER_H_
// clang-format on
