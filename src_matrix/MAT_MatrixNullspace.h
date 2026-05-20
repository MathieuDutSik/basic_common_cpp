// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXNULLSPACE_H_
#define SRC_MATRIX_MAT_MATRIXNULLSPACE_H_

// clang-format off
#include "MAT_MatrixFund.h"
#include "MAT_NonUniqueRescale.h"
#include <limits>
#include <utility>
#include <vector>
// clang-format on

template <typename T>
MyMatrix<T> NullspaceMatSingleVector(MyVector<T> const &V) {
  int n = V.size();
  for (int i = 0; i < n; i++) {
    if (V(i) != 0) {
      MyMatrix<T> NSP = ZeroMatrix<T>(n - 1, n);
      int pos = 0;
      for (int j = 0; j < n; j++) {
        if (i != j) {
          NSP(pos, i) = -V(j);
          NSP(pos, j) = V(i);
          pos++;
        }
      }
      return NSP;
    }
  }
  std::cerr << "Failed to find a non-zero index\n";
  throw TerminalException{1};
}

template <typename T, typename F>
MyMatrix<T> NullspaceTrMat_Kernel(size_t nbRow, size_t nbCol, F f) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in NullspaceTrMat_Kernel");
  size_t maxRank = nbRow;
  if (nbCol < maxRank)
    maxRank = nbCol;
  size_t sizMat = maxRank + 1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      size_t eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = eCol; iCol < nbCol; iCol++) {
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
        }
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        if (provMat(eRank, iCol) != 0)
          return iCol;
      }
      return std::numeric_limits<size_t>::max();
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != std::numeric_limits<size_t>::max()) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[FirstNonZeroCol] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          size_t StartCol = ListColSelect[iRank];
          for (size_t iCol = StartCol; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
    }
  }
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = -1;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        size_t eCol = ListColSelect[iRank];
        NSP(nbVect, eCol) = provMat(iRank, iCol);
      }
      nbVect++;
    }
  return NSP;
}

//
template <typename T, typename F>
MyMatrix<T> NullspaceTrMatTarget_Kernel(size_t nbRow, size_t nbCol,
                                        size_t target_zero, F f) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in NullspaceTrMat_Kernel");
  size_t target_rank = nbCol - target_zero;
  MyMatrix<T> provMat(target_rank, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      size_t eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = eCol; iCol < nbCol; iCol++) {
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
        }
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        if (provMat(eRank, iCol) != 0)
          return iCol;
      }
      return std::numeric_limits<size_t>::max();
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != std::numeric_limits<size_t>::max()) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[FirstNonZeroCol] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          size_t StartCol = ListColSelect[iRank];
          for (size_t iCol = StartCol; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
      if (eRank == target_rank) {
        // We reach the expected rank. So now returning.
        MyMatrix<T> NSP = ZeroMatrix<T>(target_zero, nbCol);
        size_t nbVect = 0;
        for (size_t iCol = 0; iCol < nbCol; iCol++) {
          if (ListColSelect01[iCol] == 0) {
            NSP(nbVect, iCol) = -1;
            for (size_t iRank = 0; iRank < eRank; iRank++) {
              size_t eCol = ListColSelect[iRank];
              NSP(nbVect, eCol) = provMat(iRank, iCol);
            }
            nbVect++;
          }
        }
        return NSP;
      }
    }
  }
  // The target was not achieved, we get a larger kernel than expected.
  // It could be a problem down the line, but such is life.
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = -1;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        size_t eCol = ListColSelect[iRank];
        NSP(nbVect, eCol) = provMat(iRank, iCol);
      }
      nbVect++;
    }
  return NSP;
}

template <typename T, typename F>
MyVector<T> NullspaceTrMatTargetOne_Kernel(size_t nbRow, size_t nbCol, F f) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in NullspaceTrMat_Kernel");
  size_t target_rank = nbCol - 1;
  MyMatrix<T> provMat(target_rank, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      size_t eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = eCol; iCol < nbCol; iCol++) {
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
        }
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        if (provMat(eRank, iCol) != 0)
          return iCol;
      }
      return std::numeric_limits<size_t>::max();
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != std::numeric_limits<size_t>::max()) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[FirstNonZeroCol] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          size_t StartCol = ListColSelect[iRank];
          for (size_t iCol = StartCol; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
      if (eRank == target_rank) {
        // We reach the expected rank. So now returning.
        MyVector<T> Vzero = ZeroVector<T>(nbCol);
        for (size_t iCol = 0; iCol < nbCol; iCol++) {
          if (ListColSelect01[iCol] == 0) {
            Vzero(iCol) = -1;
            for (size_t iRank = 0; iRank < eRank; iRank++) {
              size_t eCol = ListColSelect[iRank];
              Vzero(eCol) = provMat(iRank, iCol);
            }
            return Vzero;
          }
        }
      }
    }
  }
  // The target was not achieved, we get a larger kernel than expected.
  // We select one vector and it has to be processed down the line
  MyVector<T> Vzero = ZeroVector<T>(nbCol);
  auto set_vzero = [&]() -> void {
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      if (ListColSelect01[iCol] == 0) {
        Vzero(iCol) = -1;
        for (size_t iRank = 0; iRank < eRank; iRank++) {
          size_t eCol = ListColSelect[iRank];
          Vzero(eCol) = provMat(iRank, iCol);
        }
        return;
      }
    }
  };
  set_vzero();
  return Vzero;
}

template <typename T>
requires is_ring_field<T>::value
inline MyMatrix<T> NullspaceTrMat(MyMatrix<T> const &Input) {
  size_t nbRow = Input.rows();
  size_t nbCol = Input.cols();
  auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
    M.row(eRank) = Input.row(iRow);
  };
  return NullspaceTrMat_Kernel<T, decltype(f)>(nbRow, nbCol, f);
}

template <typename T>
requires (!is_ring_field<T>::value)
inline MyMatrix<T> NullspaceTrMat(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  size_t nbRow = Input.rows();
  int nbCol_i = Input.cols();
  size_t nbCol = nbCol_i;
  auto f = [&](MyMatrix<Tfield> &M, size_t eRank, size_t iRow) -> void {
    for (int iCol=0; iCol<nbCol_i; iCol++) {
      M(eRank, iCol) = UniversalScalarConversion<Tfield,T>(Input(iRow,iCol));
    }
  };
  MyMatrix<Tfield> NSP_field = NullspaceTrMat_Kernel<Tfield, decltype(f)>(nbRow, nbCol, f);
  int dim = NSP_field.rows();
  MyMatrix<T> NSP(dim, NSP_field.cols());
  for (int iNSP=0; iNSP<dim; iNSP++) {
    MyVector<Tfield> V1 = GetMatrixRow(NSP_field, iNSP);
    MyVector<Tfield> V2 = NonUniqueScaleToIntegerVector(V1);
    MyVector<T> V3 = UniversalVectorConversion<T, Tfield>(V2);
    AssignMatrixRow(NSP, iNSP, V3);
  }
  return NSP;
}

template <typename T> MyMatrix<T> NullspaceMat(MyMatrix<T> const &M) {
  return NullspaceTrMat(TransposedMat(M));
}


template <typename T>
MyMatrix<T> NullspaceTrMat_no_division(MyMatrix<T> const &Input) {
  // No division allowed. Maybe faster if we can allow for it using mpz_class.
  size_t nbRow = Input.rows();
  size_t nbCol = Input.cols();
  size_t maxRank = nbRow;
  if (nbCol < maxRank)
    maxRank = nbCol;
  size_t sizMat = maxRank + 1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  size_t miss_val = std::numeric_limits<size_t>::max();
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      provMat(eRank, iCol) = Input(iRow, iCol);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      int eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      T eVal2 = provMat(iRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = 0; iCol < nbCol; iCol++)
          provMat(eRank, iCol) =
              provMat(eRank, iCol) * eVal2 - provMat(iRank, iCol) * eVal1;
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        if (provMat(eRank, iCol) != 0)
          return iCol;
      return miss_val;
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != miss_val) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[size_t(FirstNonZeroCol)] = 1;
      T eVal2 = provMat(eRank, FirstNonZeroCol);
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          for (size_t iCol = 0; iCol < nbCol; iCol++)
            provMat(iRank, iCol) =
                provMat(iRank, iCol) * eVal2 - provMat(eRank, iCol) * eVal1;
        }
      }
      eRank++;
    }
  }
  // CODE INCOMPLETE BELOW. SOMETHING IS NEEDED.
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = 1;
      T prodVal(1);
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        size_t eCol = ListColSelect[iRank];
        T pivotVal = provMat(iRank, iCol);
        if (pivotVal != 0) {
          for (size_t jRank = 0; jRank < iRank; jRank++) {
            size_t fCol = ListColSelect[jRank];
            NSP(nbVect, fCol) *= provMat(iRank, eCol);
          }
          NSP(nbVect, iCol) *= provMat(iRank, eCol);
          //
          NSP(nbVect, eCol) = -prodVal * provMat(iRank, iCol);
          prodVal *= provMat(iRank, eCol);
        }
      }
      nbVect++;
    }
  return NSP;
}

template <typename T>
MyMatrix<T> SubspaceCompletionRational(MyMatrix<T> const &M, int const &n) {
  if (M.rows() == 0)
    return IdentityMat<T>(n);
  return NullspaceTrMat(M);
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXNULLSPACE_H_
// clang-format on
