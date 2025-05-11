// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXINT_H_
#define SRC_MATRIX_MAT_MATRIXINT_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>
// clang-format on

//  #undef TRACK_MAXIMUM_SIZE_COEFF

#ifdef DEBUG
#define DEBUG_MATRIX_INT
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_INT
#endif

// Now declarations of generic code.
// The code below generally requires the field T to be the ring (or fraction
// ring) of a factorial ring. Operations may work for fields and rings as well.
template <typename T> bool IsIntegralVector(MyVector<T> const &V) {
  int nbCol = V.size();
  for (int i = 0; i < nbCol; i++)
    if (!IsInteger(V(i)))
      return false;
  return true;
}

template <typename T> bool IsIntegralMatrix(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (!IsInteger(M(iRow, iCol)))
        return false;
  return true;
}

// T must be an integral domain with a norm function
// for computing GCD.
// eMat is a m x n matrix which defines a submodule L of T^n.
// The function computes the index i of L in T^n.
template <typename T> T Int_IndexLattice(MyMatrix<T> const &eMat) {
  static_assert(is_euclidean_domain<T>::value,
                "Requires T to be an Euclidean domain in Int_IndexLattice");
  using Treal = typename underlying_totally_ordered_ring<T>::real_type;
  size_t iRowF = 0, iColF = 0;
  MyMatrix<T> eMatW = eMat;
  size_t nbCol = eMat.cols();
  size_t nbRow = eMat.rows();
  std::vector<int> colStat(nbCol, 1);
  std::vector<int> rowStat(nbRow, 1);
  size_t nbDone = 0;
  T TheIndex = 1;
  while (true) {
    bool IsFirst = true;
    Treal MinPivot(0);
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      if (colStat[iCol] == 1) {
        for (size_t iRow = 0; iRow < nbRow; iRow++) {
          if (rowStat[iRow] == 1) {
            T eVal = eMatW(iRow, iCol);
            if (eVal != 0) {
              Treal eValA = T_NormGen(eVal);
              if (IsFirst) {
                iRowF = iRow;
                iColF = iCol;
                MinPivot = eValA;
              } else {
                if (eValA < MinPivot) {
                  iRowF = iRow;
                  iColF = iCol;
                  MinPivot = eValA;
                }
              }
              IsFirst = false;
            }
          }
        }
      }
    }
    if (IsFirst)
      return 0;
#ifdef SANITY_CHECK_MATRIX_INT
    if (MinPivot == 0) {
      std::cerr << "Clear error in the code of IndexLattice\n";
      throw TerminalException{1};
    }
#endif
    T ThePivot = eMatW(iRowF, iColF);
    bool IsFinished = true;
    for (size_t iRow = 0; iRow < nbRow; iRow++) {
      if (rowStat[iRow] == 1 && iRow != iRowF) {
        T eVal = eMatW(iRow, iColF);
        if (eVal != 0) {
          IsFinished = false;
          T TheQ = QuoInt(eVal, ThePivot);
          if (TheQ != 0)
            eMatW.row(iRow) -= TheQ * eMatW.row(iRowF);
        }
      }
    }
    if (IsFinished) {
      colStat[iColF] = 0;
      rowStat[iRowF] = 0;
      nbDone++;
      TheIndex = TheIndex * ThePivot;
    }
    if (nbDone == nbCol)
      return TheIndex;
  }
}

// This is the return type for GCD dot computation
// In input a list of entires x=(x_0, ...., x_m)
// In return we have:
// ---gcd: The greatest common divisor
// ---A vector v
//    x \dot v = gcd
template <typename T> struct GCD_dot {
  MyVector<T> V;
  T gcd;
};

template <typename T> GCD_dot<T> ComputeGcdDot(MyVector<T> const &x) {
  int siz = x.size();
  T gcd = x(0);
  MyVector<T> V(siz);
  V(0) = 1;
  for (int u = 1; u < siz; u++) {
    PairGCD_dot<T> res = ComputePairGcdDot(gcd, x(u));
    gcd = res.gcd;
    V(u) = res.b;
    for (int v = 0; v < u; v++) {
      V(v) = V(v) * res.a;
    }
  }
#ifdef SANITY_CHECK_MATRIX_INT
  T sum(0);
  for (int u = 0; u < siz; u++) {
    sum += V(u) * x(u);
  }
  if (sum != gcd) {
    std::cerr << "A: sum=" << sum << " gcd=" << gcd << "\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(V), gcd};
}

template <typename T>
GCD_dot<T> PositivityNormalizeGcdDot(GCD_dot<T> const &x) {
  if (x.gcd > 0) {
    return x;
  }
  return {-x.V, -x.gcd};
}

// This is the return type for the second extended GCD computations
// In input a list of entries x=(x_0, ...., x_m)
// In return we have:
// ---gcd: The greatest common divisor
// ---Pmat: An unimodular matrix P such that
//    x P = (gcd, 0, ....., 0)
template <typename T> struct GCD_int {
  MyMatrix<T> Pmat;
  T gcd;
};

template <typename T>
void WriteGCD_int(std::ostream &os, GCD_int<T> const &eGCD) {
  os << "gcd=" << eGCD.gcd << "\n";
  os << "Pmat=";
  WriteMatrix(os, eGCD.Pmat);
}

template <typename T>
inline typename std::enable_if<!is_mpz_class<T>::value, GCD_int<T>>::type
ComputePairGcd(T const &m, T const &n) {
  static_assert(is_euclidean_domain<T>::value,
                "Requires T to be an Euclidean domain in ComputePairGcd");
  T f, g, h, fm, gm, hm, q;
  if (n == 0 && m == 0) {
    f = 0;
    MyMatrix<T> Pmat = IdentityMat<T>(2);
    return {std::move(Pmat), f};
  }
  if (m >= 0) {
    f = m;
    fm = 1;
  } else {
    f = -m;
    fm = -1;
  }
  if (n >= 0) {
    g = n;
    gm = 0;
  } else {
    g = -n;
    gm = 0;
  }
  while (g != 0) {
    q = QuoInt(f, g);
    h = g;
    hm = gm;
    g = f - q * g;
    gm = fm - q * gm;
    f = h;
    fm = hm;
  }
  T eCoeff1, eCoeff2;
  if (n == 0) {
    eCoeff1 = fm;
    eCoeff2 = 0;
  } else {
    eCoeff1 = fm;
    eCoeff2 = (f - fm * m) / n;
  }
  MyMatrix<T> Pmat(2, 2);
  Pmat(0, 0) = eCoeff1;
  Pmat(1, 0) = eCoeff2;
  Pmat(0, 1) = -n / f;
  Pmat(1, 1) = m / f;
#ifdef SANITY_CHECK_MATRIX_INT
  T diff1 = f - Pmat(0, 0) * m - Pmat(1, 0) * n;
  T diff2 = Pmat(0, 1) * m + Pmat(1, 1) * n;
  if (diff1 != 0 || diff2 != 0) {
    std::cerr << "A: diff1=" << diff1 << " diff2=" << diff2 << "\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(Pmat), f};
}

#ifdef SRC_NUMBER_NUMBERTHEORYGMP_H_

template <typename T>
inline typename std::enable_if<is_mpz_class<T>::value, GCD_int<T>>::type
ComputePairGcd(T const &m, T const &n) {
  mpz_class eGCD;
  if (n == 0 && m == 0) {
    eGCD = 0;
    MyMatrix<T> Pmat = IdentityMat<T>(2);
    return {std::move(Pmat), eGCD};
  }
  mpz_class s, t;
  mpz_gcdext(eGCD.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), m.get_mpz_t(),
             n.get_mpz_t());
  MyMatrix<T> Pmat(2, 2);
  Pmat(0, 0) = s;
  Pmat(1, 0) = t;
  Pmat(0, 1) = -n / eGCD;
  Pmat(1, 1) = m / eGCD;
#ifdef SANITY_CHECK_MATRIX_INT
  T diff1 = eGCD - Pmat(0, 0) * m - Pmat(1, 0) * n;
  T diff2 = Pmat(0, 1) * m + Pmat(1, 1) * n;
  if (diff1 != 0 || diff2 != 0) {
    std::cerr << "m=" << m << " n=" << n << "\n";
    std::cerr << "s=" << s << " t=" << t << "  eGCD=" << eGCD << "\n";
    std::cerr << "B: diff1=" << diff1 << " diff2=" << diff2 << "\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(Pmat), eGCD};
}

#endif

template <typename T> T ComputeLCM(std::vector<T> const &eVect) {
  int siz = eVect.size();
  T eLCM(1);
  for (int i = 0; i < siz; i++)
    eLCM = LCMpair(eLCM, eVect[i]);
  return eLCM;
}

template <typename T, typename Tint>
MyVector<Tint> RescaleVec(MyVector<T> const &v) {
  static_assert(is_implementation_of_Q<T>::value,
                "Requires T to be an implementation of Q");
  static_assert(is_implementation_of_Q<Tint>::value ||
                    is_implementation_of_Z<Tint>::value,
                "Requires Tint to be an implementation of Z or Q");
  int cols = v.size();
  std::vector<Tint> dens(cols, 1);
  MyVector<Tint> vret = MyVector<Tint>(cols);
  for (int iCol = 0; iCol < cols; iCol++) {
    dens[iCol] = v(iCol).get_den();
  }
  Tint scale = LCMlist(dens);
  for (int iCol = 0; iCol < cols; iCol++) {
    vret(iCol) = (scale / v(iCol).get_den()) * v(iCol).get_num();
  }
  return vret;
}

template <typename T, typename Tint>
MyMatrix<Tint> RescaleRows(MyMatrix<T> const &M) {
  static_assert(is_implementation_of_Q<T>::value,
                "Requires T to be an implementation of Q");
  static_assert(is_implementation_of_Q<Tint>::value ||
                    is_implementation_of_Z<Tint>::value,
                "Requires Tint to be an implementation of Z or Q");
  int rows = M.rows();
  int cols = M.cols();
  std::vector<Tint> dens(cols, 1);
  MyMatrix<Tint> Mret = MyMatrix<Tint>(rows, cols);
  for (int iRow = 0; iRow < rows; iRow++) {
    for (int iCol = 0; iCol < cols; iCol++) {
      dens[iCol] = M(iRow, iCol).get_den();
    }
    Tint scale = LCMlist(dens);
    for (int iCol = 0; iCol < cols; iCol++) {
      Mret(iRow, iCol) =
          (scale / M(iRow, iCol).get_den()) * M(iRow, iCol).get_num();
    }
  }
  return Mret;
}

template <typename T> struct FractionMatrix {
  T TheMult;
  MyMatrix<T> TheMat;
};

template <typename T>
FractionMatrix<T> RemoveFractionMatrixPlusCoeff(MyMatrix<T> const &M) {
  if (IsZeroMatrix(M)) {
    T TheMult(1);
    return {TheMult, M};
  }
  int nbRow = M.rows();
  int nbCol = M.cols();
  using Tring = typename underlying_ring<T>::ring_type;
  Tring eLCM_ring(1);
  // iRow is inner loop because of cache locality
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      eLCM_ring = LCMpair(eLCM_ring, GetDenominator_z(M(iRow, iCol)));
  T eLCM = eLCM_ring;
  MyMatrix<T> M1 = eLCM * M;
  T eGCD = M1(0, 0);
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      eGCD = GcdPair(eGCD, M1(iRow, iCol));
  MyMatrix<T> M2 = M1 / eGCD;
  T TheMult = eLCM / eGCD;
  return {TheMult, std::move(M2)};
}

template<typename T>
T GetDenominatorMatrix(MyMatrix<T> const& M) {
  return RemoveFractionMatrixPlusCoeff(M).TheMult;
}

template <typename T> MyMatrix<T> RemoveFractionMatrix(MyMatrix<T> const &M) {
  return RemoveFractionMatrixPlusCoeff(M).TheMat;
}

template<typename T>
inline typename std::enable_if<is_ring_field<T>::value, MyMatrix<T>>::type
ScaledInverse(MyMatrix<T> const& M) {
  return Inverse(M);
}

template<typename T>
inline typename std::enable_if<!is_ring_field<T>::value, MyMatrix<T>>::type
ScaledInverse(MyMatrix<T> const& M) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> M_field = UniversalMatrixConversion<Tfield, T>(M);
  MyMatrix<Tfield> Minv_field = Inverse(M_field);
  MyMatrix<Tfield> MinvRescal_field = RemoveFractionMatrix(Minv_field);
  MyMatrix<T> MinvRescal = UniversalMatrixConversion<T,Tfield>(MinvRescal_field);
  return MinvRescal;
}

template <typename T>
FractionMatrix<T>
CanonicalizationSmallestCoefficientMatrixPlusCoeff(MyMatrix<T> const &M) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  T the_sma = GetSmallestMatrixCoefficient(M);
  MyMatrix<T> M2 = M / the_sma;
  T TheMult = 1 / the_sma;
  return {TheMult, std::move(M2)};
}

template <typename T> struct FractionVector {
  T TheMult;
  MyVector<T> TheVect;
};

template <typename T>
FractionVector<T>
CanonicalizationSmallestCoefficientVectorPlusCoeff(MyVector<T> const &V) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  T the_sma = GetSmallestVectorCoefficient(V);
  MyVector<T> V2 = V / the_sma;
  T TheMult = 1 / the_sma;
  return {TheMult, std::move(V2)};
}

template <typename T>
FractionMatrix<T> ScalarCanonicalizationMatrixPlusCoeff(MyMatrix<T> const &M) {
  using Tfield = typename overlying_field<T>::field_type;
  if constexpr (is_implementation_of_Q<Tfield>::value) {
    return RemoveFractionMatrixPlusCoeff(M);
  } else {
    return CanonicalizationSmallestCoefficientMatrixPlusCoeff(M);
  }
}

template <typename T>
MyMatrix<T> ScalarCanonicalizationMatrix(MyMatrix<T> const &M) {
  return ScalarCanonicalizationMatrixPlusCoeff(M).TheMat;
}

template <typename T>
FractionVector<T> RemoveFractionVectorPlusCoeff(MyVector<T> const &V) {
  if (IsZeroVector(V)) {
    T TheMult(1);
    return {TheMult, V};
  }
  int n = V.size();
  std::vector<T> eVect(n);
  T eLCM = GetDenominator(V(0));
  for (int i = 1; i < n; i++)
    eLCM = LCMpair(eLCM, GetDenominator(V(i)));
  MyVector<T> V1 = eLCM * V;
  T eGCD = V1(0);
  for (int i = 1; i < n; i++)
    eGCD = GcdPair(eGCD, V1(i));
  MyVector<T> Vret = V1 / eGCD;
  T TheMult = eLCM / eGCD;
  return {TheMult, std::move(Vret)};
}

template <typename T>
MyVector<typename underlying_ring<T>::ring_type>
NonUniqueRescaleVecRing(MyVector<T> const &V) {
  // It is non-unique because if we have V = (4, 6, 8) then we return (4, 6, 8)
  // while for RemoveFraction it returns (2,3,4).
  int n = V.size();
  using Tring = typename underlying_ring<T>::ring_type;
  std::vector<Tring> Lden(n);
  for (int i = 0; i < n; i++) {
    Lden[i] = GetDenominator_z(V(i));
  }
  Tring scale = LCMlist(Lden);
  MyVector<Tring> Vret(n);
  for (int i = 0; i < n; i++) {
    Vret(i) = (scale / Lden[i]) * GetNumerator_z(V(i));
  }
  return Vret;
}

template <typename T>
MyMatrix<typename underlying_ring<T>::ring_type>
NonUniqueRescaleRowsRing(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  using Tring = typename underlying_ring<T>::ring_type;
  MyMatrix<Tring> Mret(nbRow, nbCol);
  std::vector<Tring> dens(nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int iCol = 0; iCol < nbCol; iCol++) {
      dens[iCol] = GetDenominator_z(M(iRow, iCol));
    }
    Tring scale = LCMlist(dens);
    for (int iCol = 0; iCol < nbCol; iCol++) {
      Mret(iRow, iCol) = (scale / dens[iCol]) * GetNumerator_z(M(iRow, iCol));
    }
  }
  return Mret;
}

template <typename T>
MyMatrix<typename underlying_ring<T>::ring_type>
UniqueRescaleRowsRing(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  using Tring = typename underlying_ring<T>::ring_type;
  MyMatrix<Tring> Mret(nbRow, nbCol);
  std::vector<Tring> dens(nbCol), Vret(nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int iCol = 0; iCol < nbCol; iCol++) {
      dens[iCol] = GetDenominator_z(M(iRow, iCol));
    }
    Tring scale = LCMlist(dens);
    for (int iCol = 0; iCol < nbCol; iCol++) {
      Vret[iCol] = (scale / dens[iCol]) * GetNumerator_z(M(iRow, iCol));
    }
    Tring eGCD = Vret[0];
    for (int iCol = 1; iCol < nbCol; iCol++)
      eGCD = GcdPair(eGCD, Vret[iCol]);
    for (int iCol = 0; iCol < nbCol; iCol++) {
      Mret(iRow, iCol) = Vret[iCol] / eGCD;
    }
  }
  return Mret;
}

template <typename T>
FractionVector<T> ScalarCanonicalizationVectorPlusCoeff(MyVector<T> const &M) {
  using Tfield = typename overlying_field<T>::field_type;
  if constexpr (is_implementation_of_Q<Tfield>::value) {
    return RemoveFractionVectorPlusCoeff(M);
  } else {
    return CanonicalizationSmallestCoefficientVectorPlusCoeff(M);
  }
}

template <typename T>
MyVector<T> ScalarCanonicalizationVector(MyVector<T> const &M) {
  return ScalarCanonicalizationVectorPlusCoeff(M).TheVect;
}

template <typename T>
FractionVector<T>
NonUniqueScaleToIntegerVectorPlusCoeff_Kernel(MyVector<T> const &V) {
  using Tresidual = typename T::Tresidual;
  using Tring = typename underlying_ring<Tresidual>::ring_type;
  int siz = V.size();
  Tring eLCM_ring = ScalingInteger<Tring, T>(V(0));
  for (int i = 1; i < siz; i++)
    eLCM_ring = LCMpair(eLCM_ring, ScalingInteger<Tring, T>(V(i)));
  Tresidual eLCM_res = UniversalScalarConversion<Tresidual, Tring>(eLCM_ring);
  T eLCM(eLCM_res);
  MyVector<T> Vret = V * eLCM;
  return {eLCM, std::move(Vret)};
}

template <typename T>
FractionVector<T> NonUniqueScaleToIntegerVectorPlusCoeff(MyVector<T> const &V) {
  if constexpr (is_implementation_of_Q<T>::value) {
    return RemoveFractionVectorPlusCoeff(V);
  } else {
    return NonUniqueScaleToIntegerVectorPlusCoeff_Kernel(V);
  }
}

template <typename T>
MyVector<T> NonUniqueScaleToIntegerVector(MyVector<T> const &V) {
  return NonUniqueScaleToIntegerVectorPlusCoeff(V).TheVect;
}

template <typename T>
inline
    typename std::enable_if<is_implementation_of_Z<T>::value, MyVector<T>>::type
    CanonicalizeVector(MyVector<T> const &V) {
  return RemoveFractionVectorPlusCoeff(V).TheVect;
}

template <typename T>
inline
    typename std::enable_if<!is_float_arithmetic<T>::value, MyVector<T>>::type
    RemoveFractionVector(MyVector<T> const &V) {
  return RemoveFractionVectorPlusCoeff(V).TheVect;
}

// In this function we do not care about the invertible elements of the ring
// Below is the Z-case where we just have +1, -1.
template <typename T>
inline typename std::enable_if<is_totally_ordered<T>::value, MyVector<T>>::type
CanonicalizeVectorToInvertible(MyVector<T> const &V) {
  MyVector<T> eVect = CanonicalizeVector(V);
  int len = eVect.size();
  int FirstNZ = -1;
  for (int u = 0; u < len; u++)
    if (eVect(u) != 0 && FirstNZ == -1)
      FirstNZ = u;
  if (FirstNZ == -1)
    return eVect;
  if (eVect(FirstNZ) < 0)
    return -eVect;
  return eVect;
}

template <typename T> bool IsVectorPrimitive(MyVector<T> const &TheV) {
  size_t n = TheV.size();
  T TheGCD = TheV(0);
  for (size_t i = 1; i < n; i++) {
    T val = TheV(i);
    TheGCD = PairGcd(TheGCD, val);
  }
  return T_abs(TheGCD) == 1;
}

template <typename T>
void CheckGCD_information(GCD_int<T> const &gi, std::vector<T> const &ListX) {
  auto print_inf = [&]() -> void {
    WriteGCD_int(std::cerr, gi);
    std::cerr << "ListX =";
    for (auto &val : ListX)
      std::cerr << " " << val;
    std::cerr << "\n";
  };
  if (!IsIntegralMatrix(gi.Pmat)) {
    std::cerr << "Matrix should be integral\n";
    print_inf();
  }
  if (T_abs(DeterminantMat(gi.Pmat)) != 1) {
    std::cerr << "The determinant should be 1\n";
    print_inf();
  }
  size_t len = ListX.size();
  MyVector<T> V(len);
  for (size_t i = 0; i < len; i++)
    V(i) = ListX[i];
  MyVector<T> eP = gi.Pmat.transpose() * V;
  for (size_t i = 0; i < len; i++) {
    T val = 0;
    if (i == 0) {
      val = gi.gcd;
    }
    if (val != eP(i)) {
      std::cerr << "We should have val=" << val << " but we have eP=" << eP(i)
                << "\n";
      print_inf();
    }
  }
}

template <typename T>
GCD_int<T> ComputeGCD_information(std::vector<T> const &ListX) {
  static_assert(
      is_euclidean_domain<T>::value,
      "Requires T to be an Euclidean domain in ComputeGCD_information");
  size_t siz = ListX.size();
  if (siz == 1) {
    T gcd = ListX[0];
    MyMatrix<T> Pmat = IdentityMat<T>(1);
    GCD_int<T> eGCD_int{std::move(Pmat), gcd};
    return eGCD_int;
  }
  if (siz == 2)
    return ComputePairGcd(ListX[0], ListX[1]);
  std::vector<T> ListXred(ListX.begin(), ListX.end() - 1);
  GCD_int<T> eGCD_int = ComputeGCD_information(ListXred);
  GCD_int<T> eGCD2 = ComputePairGcd(eGCD_int.gcd, ListX[siz - 1]);
  MyMatrix<T> Pmat = MyMatrix<T>(siz, siz);
  // 1 : the column for the GCD
  for (size_t i = 0; i < siz - 1; i++)
    Pmat(i, 0) = eGCD2.Pmat(0, 0) * eGCD_int.Pmat(i, 0);
  Pmat(siz - 1, 0) = eGCD2.Pmat(1, 0);
  // 2 : The columns from the previous block
  for (size_t i = 0; i < siz - 2; i++) {
    for (size_t j = 0; j < siz - 1; j++)
      Pmat(j, i + 1) = eGCD_int.Pmat(j, i + 1);
    Pmat(siz - 1, i + 1) = 0;
  }
  // 3 : The zero columns
  for (size_t i = 0; i < siz - 1; i++)
    Pmat(i, siz - 1) = eGCD2.Pmat(0, 1) * eGCD_int.Pmat(i, 0);
  Pmat(siz - 1, siz - 1) = eGCD2.Pmat(1, 1);
  //
  return {std::move(Pmat), eGCD2.gcd};
}

// See https://en.wikipedia.org/wiki/Hermite_normal_form
// for the Row-style Hermite normal form
//
// The matrix M is rewritten as U A = H
// * H is upper triangular with H_{ij}=0 for i > j.
// * The leading coefficient of H of a row is strictly to the right of the above
// one.
// * The elements below pivots are zero and elements above pivots are
// nonnegative and strictly smaller than the pivot.
//
// This implementation is fairly naive. But its advantage is that it follows the
// same Template structure as other codes in this library.
// The output of the call of ComputeRowHermiteNormalForm ( M )
// The return (U,H) is U an unimodular matrix with U M = H
template <typename T, typename F>
void ComputeRowHermiteNormalForm_Kernel(MyMatrix<T> &H, F f) {
  size_t nbRow = H.rows();
  size_t nbCol = H.cols();
  Face StatusRow(nbRow);
  //
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  T MaxSizeCoeff = 0;
#endif

  size_t TopPosition = 0;
  MyMatrix<T> TheMat1 = IdentityMat<T>(nbRow);
  for (size_t iCol = 0; iCol < nbCol; iCol++) {
    std::vector<T> ListX;
    std::vector<size_t> ListIdx;
    bool HasNonZero = false;
    for (size_t iRow = 0; iRow < nbRow; iRow++)
      if (StatusRow[iRow] == 0) {
        ListIdx.push_back(iRow);
        T eVal = H(iRow, iCol);
        ListX.push_back(eVal);
        if (eVal != 0)
          HasNonZero = true;
      }
    if (HasNonZero) {
      //
      // Ensuring that the column has a pivot and that everything below is ZERO
      size_t siz = ListIdx.size();
      GCD_int<T> eGCD = ComputeGCD_information(ListX);
      for (size_t i = 0; i < siz; i++)
        for (size_t j = 0; j < siz; j++)
          TheMat1(ListIdx[i], ListIdx[j]) = eGCD.Pmat(j, i);
      auto fct1 = [&](MyMatrix<T> &m) -> void { m = TheMat1 * m; };
      f(fct1);
#ifdef TRACK_MAXIMUM_SIZE_COEFF
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(TheMat1));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(U));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(H));
#endif
      for (size_t i = 0; i < siz; i++)
        for (size_t j = 0; j < siz; j++) {
          if (i != j)
            TheMat1(ListIdx[i], ListIdx[j]) = 0;
          else
            TheMat1(ListIdx[i], ListIdx[j]) = 1;
        }
      //
      // Ensuring that the pivot is strictly positive
      // (in the case of integer. For other rings this is a different story)
      T eCanUnit = CanonicalizationUnit(H(TopPosition, iCol));
      if (eCanUnit != 1) {
        auto fct2 = [&](MyMatrix<T> &m) -> void {
          m.row(TopPosition) = eCanUnit * m.row(TopPosition);
        };
        f(fct2);
      }
      //
      // Putting the coefficients over the pivot
      T ThePivot = H(TopPosition, iCol);
      for (size_t iRow = 0; iRow < TopPosition; iRow++) {
        T eVal = H(iRow, iCol);
        T TheQ = QuoInt(eVal, ThePivot);
        if (TheQ != 0) {
          auto fct3 = [&](MyMatrix<T> &m) -> void {
            m.row(iRow) -= TheQ * m.row(TopPosition);
          };
          f(fct3);
        }
      }
#ifdef TRACK_MAXIMUM_SIZE_COEFF
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(U));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(H));
#endif
      //
      // Increasing the index
      StatusRow[TopPosition] = 1;
      TopPosition++;
    }
  }
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  std::cerr << "MaxSizeCoeff = " << MaxSizeCoeff << "\n";
#endif
}

template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
ComputeRowHermiteNormalForm(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  MyMatrix<T> H = M;
  MyMatrix<T> U = IdentityMat<T>(nbRow);
  auto f = [&](auto g) -> void {
    g(H);
    g(U);
  };
  ComputeRowHermiteNormalForm_Kernel(H, f);
  return {std::move(U), std::move(H)};
}

template <typename T>
MyMatrix<T> ComputeRowHermiteNormalForm_second(MyMatrix<T> const &M) {
  MyMatrix<T> H = M;
  auto f = [&](auto g) -> void { g(H); };
  ComputeRowHermiteNormalForm_Kernel(H, f);
  return H;
}

template <typename T, typename F>
void ComputeColHermiteNormalForm_Kernel(MyMatrix<T> &H, F f) {
  size_t nbRow = H.rows();
  size_t nbCol = H.cols();
  Face StatusRow(nbCol);
  //
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  T MaxSizeCoeff = 0;
#endif

  size_t TopPosition = 0;
  MyMatrix<T> TheMat1 = IdentityMat<T>(nbCol);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    std::vector<T> ListX;
    std::vector<size_t> ListIdx;
    bool HasNonZero = false;
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      if (StatusRow[iCol] == 0) {
        ListIdx.push_back(iCol);
        T eVal = H(iRow, iCol);
        ListX.push_back(eVal);
        if (eVal != 0)
          HasNonZero = true;
      }
    if (HasNonZero) {
      //
      // Ensuring that the column has a pivot and that everything below is ZERO
      size_t siz = ListIdx.size();
      GCD_int<T> eGCD = ComputeGCD_information(ListX);
      for (size_t i = 0; i < siz; i++)
        for (size_t j = 0; j < siz; j++)
          TheMat1(ListIdx[i], ListIdx[j]) = eGCD.Pmat(i, j);
      auto fct1 = [&](MyMatrix<T> &m) -> void { m = m * TheMat1; };
      f(fct1);
#ifdef TRACK_MAXIMUM_SIZE_COEFF
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(TheMat1));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(U));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(H));
#endif
      for (size_t i = 0; i < siz; i++)
        for (size_t j = 0; j < siz; j++) {
          if (i != j)
            TheMat1(ListIdx[i], ListIdx[j]) = 0;
          else
            TheMat1(ListIdx[i], ListIdx[j]) = 1;
        }
      //
      // Ensuring that the pivot is strictly positive
      // (in the case of integer. For other rings this is a different story)
      T eCanUnit = CanonicalizationUnit(H(iRow, TopPosition));
      if (eCanUnit != 1) {
        auto fct2 = [&](MyMatrix<T> &m) -> void {
          m.col(TopPosition) = eCanUnit * m.col(TopPosition);
        };
        f(fct2);
      }
      //
      // Putting the coefficients over the pivot
      T ThePivot = H(iRow, TopPosition);
      for (size_t iCol = 0; iCol < TopPosition; iCol++) {
        T eVal = H(iRow, iCol);
        T TheQ = QuoInt(eVal, ThePivot);
        if (TheQ != 0) {
          auto fct3 = [&](MyMatrix<T> &m) -> void {
            m.col(iCol) -= TheQ * m.col(TopPosition);
          };
          f(fct3);
        }
      }
#ifdef TRACK_MAXIMUM_SIZE_COEFF
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(U));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(H));
#endif
      //
      // Increasing the index
      StatusRow[TopPosition] = 1;
      TopPosition++;
    }
  }
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  std::cerr << "MaxSizeCoeff = " << MaxSizeCoeff << "\n";
#endif
}

template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
ComputeColHermiteNormalForm(MyMatrix<T> const &M) {
  int nbCol = M.cols();
  MyMatrix<T> H = M;
  MyMatrix<T> U = IdentityMat<T>(nbCol);
  auto f = [&](auto g) -> void {
    g(H);
    g(U);
  };
  ComputeColHermiteNormalForm_Kernel(H, f);
  return {std::move(U), std::move(H)};
}

template <typename T>
MyMatrix<T> ComputeColHermiteNormalForm_second(MyMatrix<T> const &M) {
  MyMatrix<T> H = M;
  auto f = [&](auto g) -> void { g(H); };
  ComputeColHermiteNormalForm_Kernel(H, f);
  return H;
}

/*
  Smith Normal form is needed for a number of applications.
  We apply here a fairly naive algorithm by acting on rows and columns.
  The result is a pair of matrices A and B such that A M B = Mred
  With Mred a reduced matrix.
  ---
  Action on rows correspond to multiplying on the left (by A)
  Action on columns correspond to multiplying on the right (by B)





 */
template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>> SmithNormalForm(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyMatrix<T> H = M;
  MyMatrix<T> ROW = IdentityMat<T>(nbRow);
  MyMatrix<T> COL = IdentityMat<T>(nbCol);
#ifdef SANITY_CHECK_MATRIX_INT
  auto check_consistency = [&](std::string const &mesg) -> void {
    MyMatrix<T> eProd = ROW * M * COL;
    if (eProd != H) {
      std::cerr << "ROW=\n";
      WriteMatrix(std::cerr, ROW);
      std::cerr << "COL=\n";
      WriteMatrix(std::cerr, COL);
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      std::cerr << "eProd=\n";
      WriteMatrix(std::cerr, eProd);
      std::cerr << "H=\n";
      WriteMatrix(std::cerr, H);
      std::cerr << "Error at stage mesg=" << mesg << "\n";
      throw TerminalException{1};
    }
  };
  std::string mesg;
#endif
  int posDone = 0;
  while (true) {
    struct choice {
      int iRow;
      int iCol;
      T norm;
    };
    std::optional<choice> opt;
    for (int iRow = posDone; iRow < nbRow; iRow++) {
      for (int iCol = posDone; iCol < nbCol; iCol++) {
        T eVal = T_abs(H(iRow, iCol));
        if (eVal != 0) {
          if (!opt) {
            opt = {iRow, iCol, eVal};
          } else {
            if (eVal < opt->norm) {
              opt = {iRow, iCol, eVal};
            }
          }
        }
      }
    }
    if (!opt)
      break;
    int iRowF = opt->iRow;
    int iColF = opt->iCol;
    T ThePivot = H(iRowF, iColF);
    bool NonZeroResidue = false;
    for (int iRow = posDone; iRow < nbRow; iRow++) {
      if (iRow != iRowF) {
        T eVal = H(iRow, iColF);
        if (eVal != 0) {
          T TheQ = QuoInt(eVal, ThePivot);
          H.row(iRow) -= TheQ * H.row(iRowF);
          ROW.row(iRow) -= TheQ * ROW.row(iRowF);
#ifdef SANITY_CHECK_MATRIX_INT
          mesg = "1 : Error_at iRowF=" + std::to_string(iRowF) +
                 " iColF=" + std::to_string(iColF) +
                 " iRpw=" + std::to_string(iRow) +
                 " TheQ=" + std::to_string(TheQ);
          check_consistency(mesg);
#endif
          if (H(iRow, iColF) != 0)
            NonZeroResidue = true;
        }
      }
    }
    for (int iCol = posDone; iCol < nbCol; iCol++) {
      if (iCol != iColF) {
        T eVal = H(iRowF, iCol);
        if (eVal != 0) {
          T TheQ = QuoInt(eVal, ThePivot);
          H.col(iCol) -= TheQ * H.col(iColF);
          COL.col(iCol) -= TheQ * COL.col(iColF);
#ifdef SANITY_CHECK_MATRIX_INT
          mesg = "2 : Error_at iRowF=" + std::to_string(iRowF) +
                 " iColF=" + std::to_string(iColF) +
                 " iCol=" + std::to_string(iCol) +
                 " TheQ=" + std::to_string(TheQ);
          check_consistency(mesg);
#endif
          if (H(iRowF, iCol) != 0)
            NonZeroResidue = true;
        }
      }
    }
    if (!NonZeroResidue) {
      if (iRowF != posDone) {
        MyMatrix<T> Trans = TranspositionMatrix<T>(nbRow, posDone, iRowF);
        ROW = Trans * ROW;
        H = Trans * H;
#ifdef SANITY_CHECK_MATRIX_INT
        mesg = "3 : Error_at iRowF=" + std::to_string(iRowF) +
               " posDone=" + std::to_string(posDone);
        check_consistency(mesg);
#endif
      }
      if (iColF != posDone) {
        MyMatrix<T> Trans = TranspositionMatrix<T>(nbCol, posDone, iColF);
        COL = COL * Trans;
        H = H * Trans;
#ifdef SANITY_CHECK_MATRIX_INT
        mesg = "4 : Error_at iColF=" + std::to_string(iColF) +
               " posDone=" + std::to_string(posDone);
        check_consistency(mesg);
#endif
      }
      T CanUnit = CanonicalizationUnit(H(posDone, posDone));
      if (CanUnit != 1) {
        ROW.row(posDone) = CanUnit * ROW.row(posDone);
        H.row(posDone) = CanUnit * H.row(posDone);
      }
      posDone++;
    }
  }
#ifdef DEBUG_MATRIX_INT
  MyMatrix<T> Test = ROW * M * COL;
  auto show_res = [&]() -> void {
    std::cerr << "Test=\n";
    WriteMatrix(std::cerr, Test);
    std::cerr << "ROW=\n";
    WriteMatrix(std::cerr, ROW);
    std::cerr << "COL=\n";
    WriteMatrix(std::cerr, COL);
    std::cerr << "M=\n";
    WriteMatrix(std::cerr, M);
    std::cerr << "Please debug\n";
    throw TerminalException{1};
  };
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      if (iRow != iCol && Test(iRow, iCol) != 0)
        show_res();
      if (iRow == iCol && Test(iRow, iCol) < 0)
        show_res();
    }
#endif
  return {ROW, COL};
}

template<typename T>
MyVector<T> SmithNormalFormInvariant(MyMatrix<T> const &M) {
  std::pair<MyMatrix<T>, MyMatrix<T>> pair = SmithNormalForm(M);
  MyMatrix<T> RedMat = pair.first * M * pair.second;
  int nbRow = M.rows();
  int nbCol = M.cols();
  int minDim = nbRow;
  if (nbCol < nbRow) {
    minDim = nbCol;
  }
  MyVector<T> VectInv(minDim);
  for (int i=0; i<minDim; i++) {
    VectInv(i) = RedMat(i,i);
  }
  return VectInv;
}

/*
  After thinking, it would seem that we need to use the
  SmithNormalForm in order to do those SubspaceCompletionInt
  operations.
 */
template <typename T>
MyMatrix<T> SubspaceCompletionInt(MyMatrix<T> const &M, int const &n) {
  int nbRow = M.rows();
#ifdef DEBUG_MATRIX_INT
  if (nbRow > 0) {
    if (M.cols() != n) {
      std::cerr << "MAT_INT: We should have M.cols() = n\n";
      throw TerminalException{1};
    }
  }
#endif
  //  int nbCol = M.cols();
  if (nbRow == 0)
    return IdentityMat<T>(n);
  std::pair<MyMatrix<T>, MyMatrix<T>> PairRed = SmithNormalForm(M);
  MyMatrix<T> const &A1 = PairRed.first;
  MyMatrix<T> const &A2 = PairRed.second;
  MyMatrix<T> TheProd = A1 * M * A2;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (TheProd(iRow, iRow) != 1) {
      std::cerr << "The Smith normal form indicates that the matrix is not "
                   "adequate\n";
      std::cerr << "The subspace spanned is not saturated\n";
      throw TerminalException{1};
    }
  }
  MyMatrix<T> A1bis = ZeroMatrix<T>(n, n);
  int d = nbRow;
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
      A1bis(i, j) = A1(i, j);
  for (int i = d; i < n; i++)
    A1bis(i, i) = 1;
  MyMatrix<T> FullBasis = Inverse(A1bis) * Inverse(A2);
  MyMatrix<T> TheCompletion(n - d, n);
  for (int i = d; i < n; i++)
    TheCompletion.row(i - d) = FullBasis.row(i);
  return TheCompletion;
}

template <typename T>
bool IsColumnNonEmpty(MyMatrix<T> const &eMat, int const &minAllowed,
                      int const &iCol) {
  int nbRow = eMat.rows();
  for (int iRow = minAllowed; iRow < nbRow; iRow++) {
    T eVal = eMat(iRow, iCol);
    if (eVal != 0)
      return true;
  }
  return false;
}

// T should be a ring of integers with factorization,
// e.g. Gaussian Integers, Eisenstein Integers, etc.
// but we need to provide QuoInt and T_Norm functions.
//
// The matrix in output is the matrix NullspaceIntMat(TransposedMat(M))
template <typename T> MyMatrix<T> NullspaceIntTrMat(MyMatrix<T> const &eMat) {
  static_assert(is_euclidean_domain<T>::value,
                "Requires T to be an Euclidean domain in NullspaceIntTrMat");
  auto INT_ClearColumn = [](MyMatrix<T> &eMat, size_t const &iCol,
                            size_t const &MinAllowedRow,
                            size_t &iRowFound) -> void {
    using Treal = typename underlying_totally_ordered_ring<T>::real_type;
    size_t nbRow = eMat.rows();
    while (true) {
      Treal MinVal(-1);
      size_t nbFound = 0;
      for (size_t iRow = MinAllowedRow; iRow < nbRow; iRow++) {
        T eVal = eMat(iRow, iCol);
        if (eVal != 0) {
          Treal AbsEVal = T_NormGen(eVal);
          if (nbFound == 0) {
            MinVal = AbsEVal;
            iRowFound = iRow;
          } else {
            if (AbsEVal < MinVal) {
              MinVal = AbsEVal;
              iRowFound = iRow;
            }
          }
          nbFound++;
        }
      }
#ifdef SANITY_CHECK_MATRIX_INT
      if (nbFound == 0) {
        std::cerr << "The column is zero. No work possible\n";
        throw TerminalException{1};
      }
#endif
      T ThePivot = eMat(iRowFound, iCol);
      for (size_t iRow = 0; iRow < nbRow; iRow++)
        if (iRow != iRowFound) {
          T eVal = eMat(iRow, iCol);
          T TheQ = QuoInt(eVal, ThePivot);
          if (TheQ != 0)
            eMat.row(iRow) -= TheQ * eMat.row(iRowFound);
        }
      if (nbFound == 1)
        return;
    }
  };
  auto SwitchRow = [](MyMatrix<T> &eMat, int const &iRow,
                      int const &jRow) -> void {
    int nbCol = eMat.cols();
    if (iRow == jRow)
      return;
    for (int iCol = 0; iCol < nbCol; iCol++) {
      T eVal1 = eMat(iRow, iCol);
      T eVal2 = eMat(jRow, iCol);
      eMat(iRow, iCol) = eVal2;
      eMat(jRow, iCol) = eVal1;
    }
  };
  MyMatrix<T> eMatW = eMat;
  size_t nbCol = eMat.cols();
  std::vector<size_t> ListIndex;
  std::vector<size_t> ListNonIndex;
  size_t eRank = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (IsColumnNonEmpty(eMatW, eRank, iCol)) {
      ListIndex.push_back(iCol);
      size_t iRowFound = 444;
      INT_ClearColumn(eMatW, iCol, eRank, iRowFound);
      SwitchRow(eMatW, eRank, iRowFound);
      eRank++;
    } else {
      ListNonIndex.push_back(iCol);
    }
  size_t dimSpace = ListNonIndex.size();
  std::vector<std::vector<T>> TheBasis;
  for (size_t i = 0; i < dimSpace; i++) {
    std::vector<T> eVect;
    for (size_t j = 0; j < dimSpace; j++) {
      if (i == j) {
        eVect.push_back(T(1));
      } else {
        eVect.push_back(T(0));
      }
    }
    TheBasis.push_back(eVect);
  }
  for (size_t iRank = 0; iRank < eRank; iRank++) {
    size_t iRow = eRank - 1 - iRank;
    size_t iCol = ListIndex[iRow];
    std::vector<T> ListX;
    T eVal = eMatW(iRow, iCol);
    ListX.push_back(eVal);
    std::vector<size_t> ListRelIndex;
    for (size_t jRow = iRow + 1; jRow < eRank; jRow++) {
      size_t jCol = ListIndex[jRow];
      ListRelIndex.push_back(jCol);
    }
    for (auto &jCol : ListNonIndex) {
      ListRelIndex.push_back(jCol);
    }
    size_t sizRelIndex = ListRelIndex.size();
    for (size_t iVect = 0; iVect < dimSpace; iVect++) {
      std::vector<T> eVect = TheBasis[iVect];
      T eSum(0);
      for (size_t iRel = 0; iRel < sizRelIndex; iRel++) {
        size_t jCol = ListRelIndex[iRel];
        T fVal = eMatW(iRow, jCol);
        eSum += eVect[iRel] * fVal;
      }
      ListX.push_back(eSum);
    }
    GCD_int<T> eGCD = ComputeGCD_information(ListX);
    //    CheckGCD_information(eGCD, ListX);
    std::vector<std::vector<T>> NewBasis;
    for (size_t iVect = 0; iVect < dimSpace; iVect++) {
      std::vector<T> eVectNew(sizRelIndex + 1, T(0));
      eVectNew[0] = eGCD.Pmat(0, iVect + 1);
      for (size_t i = 1; i <= dimSpace; i++) {
        T fVal = eGCD.Pmat(i, iVect + 1);
        std::vector<T> basVect = TheBasis[i - 1];
        for (size_t j = 0; j < sizRelIndex; j++)
          eVectNew[j + 1] += fVal * basVect[j];
      }
      NewBasis.push_back(eVectNew);
    }
    TheBasis = NewBasis;
  }
  MyMatrix<T> retNSP(dimSpace, nbCol);
  for (size_t iVect = 0; iVect < dimSpace; iVect++) {
    std::vector<T> eVect = TheBasis[iVect];
    size_t idx = 0;
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      size_t iCol = ListIndex[iRank];
      retNSP(iVect, iCol) = eVect[idx];
      idx++;
    }
    for (size_t iDim = 0; iDim < dimSpace; iDim++) {
      size_t iCol = ListNonIndex[iDim];
      retNSP(iVect, iCol) = eVect[idx];
      idx++;
    }
  }
#ifdef SANITY_CHECK_MATRIX_INT
  size_t nbRow = eMat.rows();
  for (size_t iVect = 0; iVect < dimSpace; iVect++)
    for (size_t iRow = 0; iRow < nbRow; iRow++) {
      T eSum(0);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        eSum += eMat(iRow, iCol) * retNSP(iVect, iCol);
      if (eSum != 0) {
        std::cerr << "There are remaining errors in NullspaceIntTrMat\n";
        throw TerminalException{1};
      }
    }
#endif
  return retNSP;
}

template <typename T> MyMatrix<T> NullspaceIntMat(MyMatrix<T> const &eMat) {
  return NullspaceIntTrMat(TransposedMat(eMat));
}

template<typename T>
MyMatrix<T> NullspaceIntVect(MyVector<T> const& V) {
  int n = V.size();
  MyMatrix<T> eMat(1,n);
  for (int i=0; i<n; i++) {
    eMat(0,i) = V(i);
  }
  return NullspaceIntTrMat(eMat);
}


// Given a vector v of T^n which is primitive
// the function returns a family of vector v1, ...., v(n-1)
// such that (v1, ...., v(n-1), v) is a T-basis of T^n
template <typename T> MyMatrix<T> ComplementToBasis(MyVector<T> const &TheV) {
  static_assert(is_euclidean_domain<T>::value,
                "Requires T to be an Euclidean domain in ComplementToBasis");
  using Treal = typename underlying_totally_ordered_ring<T>::real_type;
  std::vector<int> OperCol1;
  std::vector<int> OperCol2;
  std::vector<T> OperCoef;
  Treal AbsVal = -400;
  MyVector<T> TheVcopy = TheV;
  int n = TheV.size();
  int idxSelect;
  while (true) {
    bool IsFirst = true;
    int nbDiffZero = 0;
    idxSelect = -1;
    for (int i = 0; i < n; i++)
      if (TheVcopy(i) != 0) {
        nbDiffZero++;
        Treal hVal = T_NormGen(TheVcopy(i));
        if (IsFirst) {
          IsFirst = false;
          AbsVal = hVal;
          idxSelect = i;
        } else {
          if (hVal < AbsVal) {
            idxSelect = i;
            AbsVal = hVal;
          }
        }
      }
#ifdef SANITY_CHECK_MATRIX_INT
    if (idxSelect == -1) {
      std::cerr << "Inconsistency in computation of value\n";
      throw TerminalException{1};
    }
#endif
    if (nbDiffZero == 1) {
#ifdef SANITY_CHECK_MATRIX_INT
      if (AbsVal != 1) {
        std::cerr << "Wrong value for AbsVal\n";
        throw TerminalException{1};
      }
#endif
      break;
    }
    for (int j = 0; j < n; j++)
      if (TheVcopy(j) != 0 && idxSelect != j) {
        T eVal1 = TheVcopy(idxSelect);
        T eVal2 = TheVcopy(j);
        std::pair<T, T> ePair = ResQuoInt(eVal2, eVal1);
        TheVcopy(j) = ePair.first;
        OperCol1.push_back(idxSelect);
        OperCol2.push_back(j);
        OperCoef.push_back(ePair.second);
      }
  }
  size_t nbOper = OperCol1.size();
  MyMatrix<T> TheReturn = MyMatrix<T>(n - 1, n);
  int idx = 0;
  for (int i = 0; i < n; i++)
    if (i != idxSelect) {
      for (int j = 0; j < n; j++) {
        T eVal = 0;
        if (j == i)
          eVal = 1;
        TheReturn(idx, j) = eVal;
      }
      idx++;
    }
  for (size_t iOper = 0; iOper < nbOper; iOper++) {
    size_t jOper = nbOper - 1 - iOper;
    int i1 = OperCol1[jOper];
    int i2 = OperCol2[jOper];
    T eCoeff = OperCoef[jOper];
    for (int iRow = 0; iRow < n - 1; iRow++) {
      T eVal = TheReturn(iRow, i2) + eCoeff * TheReturn(iRow, i1);
      TheReturn(iRow, i2) = eVal;
    }
    T eVal = TheVcopy(i2) + eCoeff * TheVcopy(i1);
    TheVcopy(i2) = eVal;
  }
#ifdef SANITY_CHECK_MATRIX_INT
  if (!TestEquality(TheVcopy, TheV)) {
    std::cerr << "TheVcopy =";
    WriteVectorNoDim(std::cerr, TheVcopy);
    std::cerr << "TheV =";
    WriteVectorNoDim(std::cerr, TheV);
    std::cerr << "Bookkeeping error\n";
    throw TerminalException{1};
  }
#endif
  return TheReturn;
}

// We have two matrices M1 and M2 and we check if they define
// the same subspace of T^n
template <typename T>
bool TestEqualitySpaces(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
  using Treal = typename underlying_totally_ordered_ring<T>::real_type;
  size_t idxSelect = std::numeric_limits<size_t>::max();
  size_t k = M1.rows();
  size_t n = M1.cols();
  MyMatrix<T> M1copy = M1;
  MyMatrix<T> M2copy = M2;
  std::vector<int> StatusRow(k, 0);
  int idxSearch = 0;
  for (size_t iK = 0; iK < k; iK++) {
    while (true) {
      int nbDiff = 0;
      for (size_t j = 0; j < k; j++)
        if (StatusRow[j] == 0)
          if (M1copy(j, idxSearch) != 0)
            nbDiff++;
      if (nbDiff > 0)
        break;
      idxSearch++;
    }
    while (true) {
      int nbDiff = 0;
      bool IsFirst = true;
      Treal AbsVal = 0;
      for (size_t j = 0; j < k; j++) {
        T eVal = M1copy(j, idxSearch);
        if (eVal != 0 && StatusRow[j] == 0) {
          nbDiff++;
          Treal hVal = T_NormGen(eVal);
          if (IsFirst) {
            IsFirst = false;
            AbsVal = hVal;
            idxSelect = j;
          } else {
            if (hVal < AbsVal) {
              idxSelect = j;
              AbsVal = hVal;
            }
          }
        }
      }
      if (nbDiff == 1)
        break;
      for (size_t j = 0; j < k; j++)
        if (j != idxSelect) {
          T eVal1 = M1copy(idxSelect, idxSearch);
          T eVal2 = M1copy(j, idxSearch);
          T TheQ = QuoInt(eVal2, eVal1);
          if (TheQ != 0)
            M1copy.row(j) -= TheQ * M1copy.row(idxSelect);
        }
    }
    StatusRow[idxSelect] = 1;
    T eVal1 = M1copy(idxSelect, idxSearch);
    for (size_t j = 0; j < k; j++) {
      T eVal2 = M2copy(j, idxSearch);
      std::pair<T, T> ePair = ResQuoInt(eVal2, eVal1);
      if (ePair.first != 0)
        return false;
      if (ePair.second != 0)
        M2copy.row(j) = M2copy.row(j) - ePair.second * M1copy.row(idxSelect);
    }
    idxSearch++;
  }
  for (size_t j = 0; j < k; j++)
    for (size_t i = 0; i < n; i++)
      if (M2copy(j, i) != 0)
        return false;
  return true;
}

template <typename T>
int PositionSubspace(std::vector<MyMatrix<T>> const &ListSubspace,
                     MyMatrix<T> const &OneSubspace) {
  int nbSpace = ListSubspace.size();
  for (int iSub = 0; iSub < nbSpace; iSub++) {
    bool test = TestEqualitySpaces(ListSubspace[iSub], OneSubspace);
    if (test)
      return iSub;
  }
  return -1;
}

template <typename Ti, typename Td>
void FuncInsertSubspace(std::vector<MyMatrix<Ti>> &ListSubspace,
                        std::vector<Td> &ListDet,
                        MyMatrix<Ti> const &OneSubspace, Td const &OneDet) {
  if (PositionSubspace(ListSubspace, OneSubspace) == -1) {
    ListSubspace.push_back(OneSubspace);
    ListDet.push_back(OneDet);
  }
}

// T1 is integer type and T2 is real kind type
template <typename Ti, typename Td>
void ReorderSubspaceDet(std::vector<MyMatrix<Ti>> &ListSubspace,
                        std::vector<Td> &ListDet) {
  int nbSub = ListSubspace.size();
  for (int i = 0; i < nbSub - 1; i++)
    for (int j = i + 1; j < nbSub; j++)
      if (ListDet[i] > ListDet[j]) {
        MyMatrix<Ti> eSub = ListSubspace[i];
        Td eDet = ListDet[i];
        ListSubspace[i] = ListSubspace[j];
        ListDet[i] = ListDet[j];
        ListSubspace[j] = eSub;
        ListDet[j] = eDet;
      }
}

/* test if ListSub2 is a subset of ListSub1 */
template <typename T>
bool TestInclusionFamilySubspace(std::vector<MyMatrix<T>> const &ListSub1,
                                 std::vector<MyMatrix<T>> const &ListSub2) {
  int nbSub2 = ListSub2.size();
  for (int i = 0; i < nbSub2; i++)
    if (PositionSubspace(ListSub1, ListSub2[i]) == -1)
      return false;
  return true;
}

template <typename T>
bool TestEqualityFamilySubspace(std::vector<MyMatrix<T>> const &ListSub1,
                                std::vector<MyMatrix<T>> const &ListSub2) {
  int nbSub1 = ListSub1.size();
  int nbSub2 = ListSub2.size();
  if (nbSub1 != nbSub2)
    return false;
  return TestInclusionFamilySubspace(ListSub1, ListSub2);
}

template <typename T>
MyMatrix<T>
GetNoncontainedSubspace(std::vector<MyMatrix<T>> const &ListSubBig,
                        std::vector<MyMatrix<T>> const &ListSubSma) {
  int nbSubBig = ListSubBig.size();
  int nbSubSma = ListSubSma.size();
  for (int i = 0; i < nbSubBig; i++)
    if (PositionSubspace(ListSubSma, ListSubBig[i]) == -1)
      return ListSubBig[i];
  std::cerr << "We did not find any subspace\n";
  std::cerr << "nbSubBig=" << nbSubBig << " nbSumSma=" << nbSubSma << "\n";
  throw TerminalException{1};
}

template <typename T>
std::string ResultSolutionIntMat_to_GAP(const std::optional<MyVector<T>> &res) {
  if (res) {
    std::stringstream s;
    WriteVectorGAP(s, *res);
    return s.str();
  }
  return "fail";
}

// Find an integral solution of the equation Y = X A
// if it exists.
template <typename T>
std::optional<MyVector<T>> SolutionIntMat(MyMatrix<T> const &TheMat,
                                          MyVector<T> const &TheVect) {
  static_assert(is_euclidean_domain<T>::value,
                "Requires T to be an Euclidean domain in SolutionIntMat");
  using Treal = typename underlying_totally_ordered_ring<T>::real_type;
  int nbDiff;
  int len = TheVect.size();
  int nbVect = TheMat.rows();
  int nbCol = TheMat.cols();
  if (len != nbCol) {
    std::cerr << "Error in SolutionIntMat : The number of column of TheMat "
                 "should be equal to the size of TheVect\n";
    std::cerr << "nbCol(TheMat)=" << TheMat.cols()
              << " |TheVect|=" << TheVect.size() << "\n";
    throw TerminalException{1};
  }
  if (nbVect == 0) {
    MyVector<T> eSol;
    if (IsZeroVector(TheVect)) {
      return eSol;
    } else {
      return {};
    }
  }
  MyVector<T> eSol = ZeroVector<T>(nbVect);
  MyMatrix<T> eEquivMat = IdentityMat<T>(nbVect);
  MyMatrix<T> TheMatWork = TheMat;
  MyVector<T> TheVectWork = TheVect;
  std::vector<int> VectStatus(nbVect, 1);
  for (int i = 0; i < nbCol; i++) {
    int iVectFound = -1;
    while (true) {
      bool IsFirst = true;
      Treal MinValue(0);
      nbDiff = 0;
      for (int iVect = 0; iVect < nbVect; iVect++)
        if (VectStatus[iVect] == 1) {
          T prov1 = TheMatWork(iVect, i);
          Treal eNorm = T_NormGen(prov1);
          if (prov1 != 0) {
            nbDiff++;
            if (IsFirst) {
              IsFirst = false;
              MinValue = eNorm;
              iVectFound = iVect;
            } else {
              if (eNorm < MinValue) {
                MinValue = eNorm;
                iVectFound = iVect;
              }
            }
          }
        }
      if (nbDiff == 1 || nbDiff == 0) {
        break;
      }
#ifdef SANITY_CHECK_MATRIX_INT
      if (MinValue == 0) {
        std::cerr << "MinValue should not be zero\n";
        throw TerminalException{1};
      }
#endif
      for (int iVect = 0; iVect < nbVect; iVect++)
        if (VectStatus[iVect] == 1 && iVect != iVectFound) {
          T prov1b = TheMatWork(iVectFound, i);
          T prov2 = TheMatWork(iVect, i);
          T TheQ = QuoInt(prov2, prov1b);
          if (TheQ != 0) {
            TheMatWork.row(iVect) -= TheQ * TheMatWork.row(iVectFound);
            eEquivMat.row(iVect) -= TheQ * eEquivMat.row(iVectFound);
          }
        }
    }
    if (nbDiff == 1) {
#ifdef SANITY_CHECK_MATRIX_INT
      if (iVectFound == -1) {
        std::cerr << "Clear error in the program\n";
        throw TerminalException{1};
      }
#endif
      VectStatus[iVectFound] = 0;
      T prov1 = TheVectWork(i);
      T prov2 = TheMatWork(iVectFound, i);
      T TheQ = QuoInt(prov1, prov2);
      if (TheQ != 0) {
        for (int j = 0; j < nbCol; j++)
          TheVectWork(j) -= TheQ * TheMatWork(iVectFound, j);
        for (int iVect = 0; iVect < nbVect; iVect++)
          eSol(iVect) += TheQ * eEquivMat(iVectFound, iVect);
      }
    }
    if (TheVectWork(i) != 0)
      return {};
  }
  return eSol;
}

template <typename T>
struct RecSolutionIntMat {
private:
  int nbRow;
  int nbCol;
  std::vector<int> ListRow;
  MyMatrix<T> TheMatWork;
  MyMatrix<T> eEquivMat;
public:
  RecSolutionIntMat(MyMatrix<T> const &TheMat) {
    static_assert(is_euclidean_domain<T>::value,
                  "Requires T to be an Euclidean domain in "
                  "ComputeCanonicalFormFastReduction");
    using Treal = typename underlying_totally_ordered_ring<T>::real_type;
    int nbDiff;
    nbRow = TheMat.rows();
    nbCol = TheMat.cols();
    ListRow = std::vector<int>(nbCol);
    eEquivMat = IdentityMat<T>(nbRow);
    TheMatWork = TheMat;
#ifdef SANITY_CHECK_MATRIX_INT
    if (nbRow == 0) {
      std::cerr << "Need to write the code here\n";
      throw TerminalException{1};
    }
#endif
    std::vector<int> VectStatus(nbRow, 1);
    for (int i = 0; i < nbCol; i++) {
      int iVectFound = -1;
      while (true) {
        bool IsFirst = true;
        Treal MinValue(0);
        nbDiff = 0;
        for (int iVect = 0; iVect < nbRow; iVect++) {
          if (VectStatus[iVect] == 1) {
            T prov1 = TheMatWork(iVect, i);
            Treal eNorm = T_NormGen(prov1);
            if (prov1 != 0) {
              nbDiff++;
              if (IsFirst) {
                IsFirst = false;
                MinValue = eNorm;
                iVectFound = iVect;
              } else {
                if (eNorm < MinValue) {
                  MinValue = eNorm;
                  iVectFound = iVect;
                }
              }
            }
          }
        }
        if (nbDiff == 1 || nbDiff == 0) {
          break;
        }
#ifdef SANITY_CHECK_MATRIX_INT
        if (MinValue == 0) {
          std::cerr << "MinValue should not be zero\n";
          throw TerminalException{1};
        }
#endif
        for (int iVect = 0; iVect < nbRow; iVect++)
          if (VectStatus[iVect] == 1 && iVect != iVectFound) {
            T prov1b = TheMatWork(iVectFound, i);
            T prov2 = TheMatWork(iVect, i);
            T TheQ = QuoInt(prov2, prov1b);
            if (TheQ != 0) {
              TheMatWork.row(iVect) -= TheQ * TheMatWork.row(iVectFound);
              eEquivMat.row(iVect) -= TheQ * eEquivMat.row(iVectFound);
            }
          }
      }
      int eVal;
      if (nbDiff == 1) {
        eVal = iVectFound;
#ifdef SANITY_CHECK_MATRIX_INT
        if (iVectFound == -1) {
          std::cerr << "Clear error in the program\n";
          throw TerminalException{1};
        }
#endif
        VectStatus[iVectFound] = 0;
      } else {
        eVal = -1;
      }
      ListRow[i] = eVal;
    }
  }
  bool has_solution_v(MyVector<T> const &TheVect) const {
    MyVector<T> TheVectWork = TheVect;
    for (int i = 0; i < nbCol; i++) {
      int iRow = ListRow[i];
      if (iRow >= 0) {
        T const& prov2 = TheMatWork(iRow, i);
        T TheQ = QuoInt(TheVectWork(i), prov2);
        if (TheQ != 0) {
          for (int j = 0; j < nbCol; j++) {
            TheVectWork(j) -= TheQ * TheMatWork(iRow, j);
          }
        }
      }
      if (TheVectWork(i) != 0) {
        return false;
      }
    }
    return true;
  }
  std::optional<MyVector<T>> get_solution_v(MyVector<T> const &TheVect) const {
    int nbVect = TheMatWork.rows();
    int nbCol = TheMatWork.cols();
    MyVector<T> TheVectWork = TheVect;
    MyVector<T> eSol = ZeroVector<T>(nbVect);
    for (int i = 0; i < nbCol; i++) {
      int iRow = ListRow[i];
      if (iRow >= 0) {
        T const& prov2 = TheMatWork(iRow, i);
        T TheQ = QuoInt(TheVectWork(i), prov2);
        if (TheQ != 0) {
          for (int j = 0; j < nbCol; j++) {
            TheVectWork(j) -= TheQ * TheMatWork(iRow, j);
          }
          for (int iVect = 0; iVect < nbVect; iVect++) {
            eSol(iVect) += TheQ * eEquivMat(iRow, iVect);
          }
        }
      }
      if (TheVectWork(i) != 0) {
        return {};
      }
    }
    return eSol;
  }
  std::optional<MyMatrix<T>> get_solution_m(MyMatrix<T> const &TheMat) const {
    int n_row = TheMat.rows();
    MyMatrix<T> RetMat(n_row, nbRow);
    for (int i = 0; i < n_row; i++) {
      MyVector<T> V = GetMatrixRow(TheMat, i);
      std::optional<MyVector<T>> opt = get_solution_v(V);
      if (!opt)
        return {};
      AssignMatrixRow(RetMat, i, *opt);
    }
    return RetMat;
  }
  bool is_containing_m(MyMatrix<T> const &TheMat) const {
    int n_row = TheMat.rows();
    MyVector<T> V(nbCol);
    for (int i_row=0; i_row<n_row; i_row++) {
      for (int i_col=0; i_col<nbCol; i_col++) {
        V(i_col) = TheMat(i_row, i_col);
      }
      bool test = has_solution_v(V);
      if (!test) {
        return false;
      }
    }
    return true;
  }
};

template <typename T> struct BasisReduction {
  MyMatrix<T> TheBasisReduced;
  MyMatrix<T> Pmat;
  std::vector<int> IdxVector;
  MyMatrix<T> TheBasisReord;
};

// We need TheBasis to be of full rank.
// This code is for the Copositivity.
//
// Following operation is done on the matrix
// The matrix TheBasis is transformed by operations
// on the columns in order to diminish the size
// of the coefficients.
//
// Specifically, the matrix TheBasisReord is changed so that
// it is of the form
// (a11 0 ....... 0)
// ( x  a22 0 ... 0)
//     .
//     .
// ( x  .... x  ann)
// We have aII > 0
// This is done in two ways
// Applying an integral matrix to the column
// and then a permutation on the rows.
template <typename T>
BasisReduction<T> ComputeBasisReduction(MyMatrix<T> const &TheBasis) {
  static_assert(
      is_euclidean_domain<T>::value,
      "Requires T to be an Euclidean domain in ComputeBasisReduction");
  using Treal = typename underlying_totally_ordered_ring<T>::real_type;
  size_t nbCol = TheBasis.cols();
  size_t nbRow = TheBasis.rows();
  std::vector<int> colStat(nbCol, 1);
  std::vector<int> rowStat(nbRow, 1);
  MyMatrix<T> TheBasisReduced = TheBasis;
  MyMatrix<T> Pmat = IdentityMat<T>(nbCol);
  MyMatrix<T> TheBasisReord(nbRow, nbCol);
  std::vector<int> IdxVector;
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto FindMinGCDrow = [&](int const &iRank) -> size_t {
    size_t iRowSearch = miss_val;
    bool IsFirst = true;
    Treal AbsVal = 0;
    for (size_t iRow = 0; iRow < nbRow; iRow++)
      if (rowStat[iRow] == 1) {
        std::vector<T> eRowRed(nbCol - iRank);
        for (size_t iCol = 0; iCol < nbCol - iRank; iCol++) {
          T eVal = TheBasisReduced(iRow, iCol + iRank);
          eRowRed[iCol] = eVal;
        }
        GCD_int<T> eGCD = ComputeGCD_information(eRowRed);
        Treal eNorm = T_NormGen(eGCD.gcd);
        if (IsFirst) {
          IsFirst = false;
          AbsVal = eNorm;
          iRowSearch = iRow;
        } else {
          if (eNorm < AbsVal) {
            AbsVal = eNorm;
            iRowSearch = iRow;
          }
        }
      }
    return iRowSearch;
  };
  auto SingleMultiplicationUpdate = [&](MyMatrix<T> const &PartMat) -> void {
    Pmat = Pmat * PartMat;
    TheBasisReduced = TheBasisReduced * PartMat;
  };
  auto UpdateMatrices = [&](size_t const &iRank,
                            size_t const &iRowSearch) -> void {
    rowStat[iRowSearch] = 0;
    IdxVector.push_back(iRowSearch);
    std::vector<T> eRowRed(nbCol - iRank);
    for (size_t iCol = 0; iCol < nbCol - iRank; iCol++)
      eRowRed[iCol] = TheBasisReduced(iRowSearch, iCol + iRank);
    GCD_int<T> eGCD = ComputeGCD_information(eRowRed);
    MyMatrix<T> PartMat = IdentityMat<T>(nbCol);
    for (size_t iCol = iRank; iCol < nbCol; iCol++)
      for (size_t iRow = iRank; iRow < nbCol; iRow++)
        PartMat(iRow, iCol) = eGCD.Pmat(iRow - iRank, iCol - iRank);

    SingleMultiplicationUpdate(PartMat);
    for (size_t iCol = 0; iCol < iRank; iCol++) {
      T a = TheBasisReduced(iRowSearch, iCol);
      T b = TheBasisReduced(iRowSearch, iRank);
      T q = QuoInt(a, b);
      MyMatrix<T> PartMatB = IdentityMat<T>(nbCol);
      PartMatB(iRank, iCol) = -q;
      SingleMultiplicationUpdate(PartMatB);
    }
    if (TheBasisReduced(iRowSearch, iRank) < 0) {
      MyMatrix<T> PartMatC = IdentityMat<T>(nbCol);
      PartMatC(iRank, iRank) = -1;
      SingleMultiplicationUpdate(PartMatC);
    }
    TheBasisReord.row(iRank) = TheBasisReduced.row(iRowSearch);
  };
  size_t TheRank = std::min(nbCol, nbRow);
  for (size_t iRank = 0; iRank < TheRank; iRank++) {
    size_t iRowSearch = FindMinGCDrow(iRank);
    UpdateMatrices(iRank, iRowSearch);
  }
  BasisReduction<T> eRed{TheBasisReduced, Pmat, IdxVector, TheBasisReord};
  return eRed;
}

struct AffineBasisResult {
  bool result;
  std::vector<int> ListIdx;
};

// Given a family of points EXT, find a subset (v1, ...., vN) of EXT
// such that for every point of v of EXT there exist lambda1, ..., lambdaN in T
// with v=lambda1 v1 + .... + lambdaN vN
// This may not exist
template <typename T>
AffineBasisResult Kernel_ComputeAffineBasis(MyMatrix<T> const &EXT) {
  static_assert(is_ring_field<T>::value,
                "Requires T to have inverses in Kernel_ComputeAffineBasis");
  size_t nbRow = EXT.rows();
  size_t nbCol = EXT.cols();
  size_t n = nbCol;
  MyMatrix<T> ListExp = ZeroMatrix<T>(nbRow, nbCol);
  MyMatrix<T> EXTwork = EXT;
  std::vector<int> RowStatus(nbRow, 0);
  std::vector<int> ColumnStatus(nbCol, 1);
  std::cerr << "Starting Kernel_ComputeAffineBasis\n";
  MyVector<T> V1(nbCol);
  MyVector<T> V2(nbCol);
  MyVector<T> eExpr(nbCol);
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto fInsertValue = [&](size_t const &idx, int const &iVect) -> bool {
    size_t eCol = miss_val;
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      if (eCol == miss_val && EXTwork(iVect, iCol) != 0 &&
          ColumnStatus[iCol] == 1)
        eCol = iCol;
#ifdef SANITY_CHECK_MATRIX_INT
    if (eCol == miss_val) {
      std::cerr << "This should not be selected\n";
      std::cerr << "nbCol=" << nbCol << "\n";
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        std::cerr << " iCol=" << iCol << " stat=" << ColumnStatus[iCol]
                  << " val=" << EXTwork(iVect, iCol) << "\n";
      }
      throw TerminalException{1};
    }
#endif
    V1 = EXTwork.row(iVect);
    std::cerr << "V1 rows=" << V1.rows() << " cols=" << V1.cols() << "\n";
    ListExp(iVect, idx) = 1;
    for (size_t iRow = 0; iRow < nbRow; iRow++)
      if (RowStatus[iRow] == 0) {
        V2 = EXTwork.row(iRow);
        bool test = IsVectorMultiple(V1, V2);
        if (test) {
          T eQuot = V2(eCol) / V1(eCol);
          eExpr = ListExp.row(iRow) - eQuot * ListExp.row(iVect);
          for (size_t iCol = 0; iCol < nbCol; iCol++)
            if (!IsInteger(eExpr(iCol)))
              return false;
        }
      }
    ColumnStatus[eCol] = 0;
    RowStatus[iVect] = 1;
    for (size_t iRow = 0; iRow < nbRow; iRow++)
      if (RowStatus[iRow] == 0) {
        V2 = EXTwork.row(iRow);
        T eQuot = V2(eCol) / V1(eCol);
        ListExp.row(iRow) = ListExp.row(iRow) - eQuot * ListExp.row(iVect);
        EXTwork.row(iRow) = EXTwork.row(iRow) - eQuot * EXTwork.row(iVect);
        V2 = EXTwork.row(iRow);
        bool IsZero = IsZeroVector(V2);
        if (IsZero)
          RowStatus[iRow] = 1;
      }
    size_t nbFinished = 0;
    for (size_t iRow = 0; iRow < nbRow; iRow++)
      if (RowStatus[iRow] == 1)
        nbFinished++;
    std::cerr << "nbFinished=" << nbFinished << "\n";
    return true;
  };
  int nbIter = 1000;
  std::vector<int> ListIdx(n);
  std::vector<int> UsedNumber(nbRow, 0);
  auto GetRandomNumber = [&]() -> int {
    for (int iter = 0; iter < nbIter; iter++) {
      int eVal = random() % nbRow;
      if (UsedNumber[eVal] == 0 && RowStatus[eVal] == 0)
        return eVal;
    }
    return -1;
  };
  auto SetLocallyCorrectIndex = [&](size_t const &idx) -> int {
    for (int iter = 0; iter < nbIter; iter++) {
      int eVal = GetRandomNumber();
      if (eVal == -1)
        return -1;
      bool res = fInsertValue(idx, eVal);
      UsedNumber[eVal] = 1;
      if (res) {
        ListIdx[idx] = eVal;
        return 0;
      }
    }
    return -1;
  };
  for (size_t i = 0; i < n; i++) {
    std::cerr << "i=" << i << "\n";
    int eVal = SetLocallyCorrectIndex(i);
    if (eVal == -1)
      return {false, {}};
  }
  return {true, std::move(ListIdx)};
}

template <typename T> MyMatrix<T> RandomUnimodularMatrix(int const &n) {
  MyMatrix<T> RetMat = IdentityMat<T>(n);
  int n_iter = 3 * n;
  for (int iter = 0; iter < n_iter; iter++) {
    MyMatrix<T> eMat = IdentityMat<T>(n);
    int idx1 = random() % n;
    int idx2 = random() % n;
    if (idx1 != idx2) {
      int pivot = (random() % 21) - 10;
      eMat(idx1, idx2) = pivot;
    }
    RetMat = eMat * RetMat;
  }
  return RetMat;
}

template <typename T>
AffineBasisResult ComputeAffineBasis(MyMatrix<T> const &EXT) {
  int nbIter = 1000;
  for (int iter = 0; iter < nbIter; iter++) {
    AffineBasisResult eAffRes = Kernel_ComputeAffineBasis(EXT);
    if (eAffRes.result)
      return eAffRes;
  }
  return {false, {}};
}

// Compute the translation classes.
// Two classes eV and fV are equivalent if there exists a vector w integer such
// that eV - fV = w M that is we need to compute (eV - fV) M^(-1)
template <typename T, typename Tout>
std::vector<MyVector<Tout>>
Kernel_ComputeTranslationClasses(MyMatrix<T> const &M) {
  int n = M.rows();
  MyMatrix<T> eInv = Inverse(M);
  std::vector<MyVector<Tout>> ListClasses;
  std::vector<uint8_t> ListStatus;
  auto IsEquivalent = [&](MyVector<Tout> const &eV,
                          MyVector<Tout> const &fV) -> bool {
    MyVector<Tout> diff = eV - fV;
    for (int i = 0; i < n; i++) {
      T eVal = 0;
      for (int j = 0; j < n; j++)
        eVal += diff(j) * eInv(j, i);
      if (!IsInteger(eVal))
        return false;
    }
    return true;
  };
  auto FuncInsert = [&](MyVector<Tout> const &eV) -> void {
    for (auto &fV : ListClasses) {
      if (IsEquivalent(eV, fV))
        return;
    }
    ListClasses.push_back(eV);
    ListStatus.push_back(1);
  };
  MyVector<Tout> zerV = ZeroVector<Tout>(n);
  FuncInsert(zerV);
  while (true) {
    bool IsFinished = true;
    size_t nbClass = ListClasses.size();
    for (size_t iClass = 0; iClass < nbClass; iClass++) {
      if (ListStatus[iClass] == 1) {
        ListStatus[iClass] = 0;
        IsFinished = false;
        MyVector<Tout> eClass = ListClasses[iClass];
        for (int i = 0; i < n; i++) {
          MyVector<Tout> fClass = eClass;
          fClass[i]++;
          FuncInsert(fClass);
        }
      }
    }
    if (IsFinished)
      break;
  }
#ifdef SANITY_CHECK_MATRIX_INT
  T det = T_abs(DeterminantMat(M));
  T n_class = ListClasses.size();
  if (det != n_class) {
    std::cerr << "The determinant det=" << det
              << " does not coincide with the number of classes =" << n_class
              << "\n";
    std::cerr << "This ought to be considered a bug\n";
    throw TerminalException{1};
  }
#endif
  return ListClasses;
}

template <typename T, typename Tout>
inline typename std::enable_if<is_ring_field<T>::value,
                               std::vector<MyVector<Tout>>>::type
ComputeTranslationClasses(MyMatrix<T> const &Input) {
  return Kernel_ComputeTranslationClasses<T, Tout>(Input);
}

template <typename T, typename Tout>
inline typename std::enable_if<!is_ring_field<T>::value,
                               std::vector<MyVector<Tout>>>::type
ComputeTranslationClasses(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> Input_field = UniversalMatrixConversion<Tfield, T>(Input);
  return Kernel_ComputeTranslationClasses<Tfield, Tout>(Input_field);
}

template <typename T>
std::vector<size_t> GetActionOnClasses(std::vector<MyVector<T>> const &l_v,
                                       MyMatrix<T> const &Transform,
                                       MyMatrix<T> const &M) {
  MyMatrix<T> eInv = Inverse(M);
  int n = M.rows();
  auto IsEquivalent = [&](MyVector<T> const &eV,
                          MyVector<T> const &fV) -> bool {
    MyVector<T> diff = eV - fV;
    for (int i = 0; i < n; i++) {
      T eVal = 0;
      for (int j = 0; j < n; j++)
        eVal += diff(j) * eInv(j, i);
      if (!IsInteger(eVal))
        return false;
    }
    return true;
  };
  size_t len = l_v.size();
  auto get_position = [&](MyVector<T> const &eV) -> size_t {
    for (size_t i = 0; i < len; i++)
      if (IsEquivalent(eV, l_v[i]))
        return i;
    std::cerr << "Failed to find\n";
    throw TerminalException{1};
  };
  std::vector<size_t> V(len);
  for (size_t i = 0; i < len; i++) {
    MyVector<T> eV_img = Transform.transpose() * l_v[i];
    V[i] = get_position(eV_img);
  }
  return V;
}

// Given a family of vector v1, ....., vN
// Find an integral family of vector w1, ...., wR
// such that all vi are expressed integrally in term of wi
// and uniquely.
// This always exist as opposed to affine basis.
//
// As it happens, the algorithms works canonically.
// That is if you replace the family of vectors V by VP
// for some invertible P then we have
// Zbasis(VP) = Zbasis(V) P
// But of course this depends on the ordering of the operations.
//
template <typename T> MyMatrix<T> GetZbasis(MyMatrix<T> const &ListElement) {
  static_assert(is_euclidean_domain<T>::value,
                "Requires T to be an Euclidean domain in GetZbasis");
  using Treal = typename underlying_totally_ordered_ring<T>::real_type;
  int TheDim = ListElement.cols();
  MyMatrix<T> ListEqua;
  MyMatrix<T> InvMatrix;
  MyMatrix<T> InvMatrixTr;
  MyMatrix<T> TheBasis;
  std::vector<int> eSet;
  auto fGetOneBasis = [&](MyVector<T> const &eSol) -> MyMatrix<T> {
    int DimLoc = TheBasis.rows();
    MyMatrix<T> TheRedMat = ZeroMatrix<T>(DimLoc + 1, DimLoc);
    for (int i = 0; i < DimLoc; i++)
      TheRedMat(i, i) = 1;
    for (int i = 0; i < DimLoc; i++)
      TheRedMat(DimLoc, i) = eSol(i);
    MyMatrix<T> NSP = NullspaceIntMat(TheRedMat);
#ifdef SANITY_CHECK_MATRIX_INT
    if (NSP.rows() != 1) {
      std::cerr << "|NSP|=" << NSP.rows() << " when it should be 1\n";
      std::cerr << "TheRedMat:\n";
      WriteMatrix(std::cerr, TheRedMat);
      std::cerr << "TheRedMat:\n";
      WriteMatrix(std::cerr, NSP);
      std::cerr << "Inconsistency that needs to be corrected\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> eVect_pre = GetMatrixRow(NSP, 0);
    MyVector<T> eVect = CanonicalizeVectorToInvertible(eVect_pre);
    int n2 = DimLoc + 1;
    while (true) {
      std::vector<int> ListIdxZ;
      std::vector<int> ListIdxNZ;
      for (int i = 0; i < n2; i++) {
        if (eVect(i) == 0) {
          ListIdxZ.push_back(i);
        } else {
          ListIdxNZ.push_back(i);
        }
      }
      if (ListIdxNZ.size() == 1)
        return SelectRow(TheRedMat, ListIdxZ);
      std::vector<int> AbsList;
      bool IsFirst = true;
      int ThePivot = -1;
      Treal TheMin = -1;
      for (auto &eVal : ListIdxNZ) {
        T eVal_T = eVect(eVal);
        Treal eAbs = T_NormGen(eVal_T);
        if (IsFirst) {
          TheMin = eAbs;
          ThePivot = eVal;
        } else {
          if (eAbs < TheMin) {
            TheMin = eAbs;
            ThePivot = eVal;
          }
        }
        IsFirst = false;
      }
      for (int iCol = 0; iCol < n2; iCol++)
        if (iCol != ThePivot) {
          T TheQ = QuoInt(eVect(iCol), eVect(ThePivot));
          if (TheQ != 0) {
            TheRedMat.row(ThePivot) += TheQ * TheRedMat.row(iCol);
            eVect(iCol) -= TheQ * eVect(ThePivot);
          }
        }
    }
  };
  auto fComputeSpeed = [&]() -> void {
    int dimSpace = TheBasis.rows();
    if (dimSpace == 0) {
      ListEqua = IdentityMat<T>(TheDim);
      eSet = {};
    } else {
      ListEqua = NullspaceTrMat(TheBasis);
      eSet = ColumnReductionSet(TheBasis);
      InvMatrix = Inverse(SelectColumn(TheBasis, eSet));
      InvMatrixTr = InvMatrix.transpose();
    }
  };
  auto IsInSpace = [&](MyVector<T> const &eElt) -> bool {
    int nbEqua = ListEqua.rows();
    for (int iEqua = 0; iEqua < nbEqua; iEqua++) {
      T eSum = 0;
      for (int i = 0; i < TheDim; i++)
        eSum += ListEqua(iEqua, i) * eElt(i);
      if (eSum != 0)
        return false;
    }
    return true;
  };
  auto fInsert = [&](MyVector<T> const &eElt) -> void {
    bool test = IsInSpace(eElt);
    if (!test) {
      TheBasis = ConcatenateMatVec(TheBasis, eElt);
      fComputeSpeed();
    } else {
      if (TheBasis.rows() == 0)
        return;
      MyVector<T> eEltRed = SelectColumnVector(eElt, eSet);
      MyVector<T> eSol = InvMatrixTr * eEltRed;
      if (IsIntegralVector(eSol))
        return;
      MyMatrix<T> NewBasis = fGetOneBasis(eSol);
      TheBasis = NewBasis * TheBasis;
      fComputeSpeed();
    }
  };
  fComputeSpeed();
  int nbElt = ListElement.rows();
  for (int iElt = 0; iElt < nbElt; iElt++) {
    MyVector<T> eElt = GetMatrixRow(ListElement, iElt);
    fInsert(eElt);
  }

#ifdef SANITY_CHECK_MATRIX_INT
  int DimSpace = TheBasis.rows();
  for (int iBas = 0; iBas < DimSpace; iBas++) {
    MyVector<T> eLine = GetMatrixRow(TheBasis, iBas);
    std::optional<MyVector<T>> opt = SolutionIntMat(ListElement, eLine);
    if (!opt) {
      std::cerr << "Error in GetZbasis 1\n";
      throw TerminalException{1};
    }
  }
  for (int iElt = 0; iElt < nbElt; iElt++) {
    MyVector<T> eElt = GetMatrixRow(ListElement, iElt);
    std::optional<MyVector<T>> opt = SolutionIntMat(TheBasis, eElt);
    if (!opt) {
      std::cerr << "Error in GetZbasis 2\n";
      throw TerminalException{1};
    }
  }
#endif
  return TheBasis;
}

/*
  WRONG IDEA:
  M1 spans a lattice L1
  M2 spans a lattice L2
  We want to find a basis of the lattice L1 cap L2.
  ---
  We have the formula (L1 \cap L2)* = L1* + L2*
  This allows to apply the GetZbasis function

  CORRECT SOLUTION:
  write the equation system:
  sum_i lambda_i v^1_i = sum_j mu_j v^2_j
  with lambda and mu integer and deduce from there.
  ---
  This solution works even if the lattices are not full
  dimensional. But this requires a little bit of additional processing.
 */
template <typename T>
MyMatrix<T> IntersectionLattice(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
  int dim1 = M1.rows();
  int dim2 = M2.rows();
  int n = M1.cols();
  if (n != M2.cols()) {
    std::cerr << "The dimension of M1 and M2 should be the same\n";
    throw TerminalException{1};
  }
  MyMatrix<T> M1_M2 = Concatenate(M1, M2);
  MyMatrix<T> NSP = NullspaceIntMat(M1_M2);
  std::vector<int> L(dim1);
  for (int i = 0; i < dim1; i++)
    L[i] = i;
  MyMatrix<T> NSPred = SelectColumn(NSP, L);
  MyMatrix<T> SpannSet = NSPred * M1;
  if (dim1 == n && dim2 == n) {
    if (SpannSet.rows() != n) {
      std::cerr << "The dimension should be exactl n. Bug to be solved\n";
      throw TerminalException{1};
    }
  }
  return SpannSet;
}

template <typename T>
MyMatrix<T> IntersectionLattice_VectorSpace(MyMatrix<T> const &Latt,
                                            MyMatrix<T> const &Space) {
  int n = Latt.rows();
  int n_spa = Space.rows();
  MyMatrix<T> eBasis = ExtendToBasis(Space);
  MyMatrix<T> Latt2 = Latt * Inverse(eBasis);
  std::vector<int> V(n - n_spa);
  for (int i = 0; i < n - n_spa; i++)
    V[i] = i + n_spa;
  MyMatrix<T> Latt3 = SelectColumn(Latt2, V);
  MyMatrix<T> NSP = NullspaceIntMat(Latt3);
#ifdef SANITY_CHECK_MATRIX_INT
  if (!IsIntegralMatrix(NSP)) {
    std::cerr << "NSP should be integral\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> IntBasis = NSP * Latt;
#ifdef SANITY_CHECK_MATRIX_INT
  for (int i_s = 0; i_s < n_spa; i_s++) {
    MyVector<T> v = GetMatrixRow(IntBasis, i_s);
    std::optional<MyVector<T>> opt1 = SolutionIntMat(Latt, v);
    if (!opt1) {
      std::cerr << "The vector should be expressed integrally in terms of the "
                   "lattice\n";
      throw TerminalException{1};
    }
    std::optional<MyVector<T>> opt2 = SolutionMat(Space, v);
    if (!opt2) {
      std::cerr
          << "The vector should be expressed integrally in terms of the spacen";
      throw TerminalException{1};
    }
  }
#endif
  return IntBasis;
}

template <typename Tint>
MyMatrix<Tint> SYMPL_ComputeSymplecticBasis(MyMatrix<Tint> const &M) {
  int nb_row = M.rows();
  int n = M.cols() / 2;
  MyMatrix<Tint> Mwork = M;
  MyMatrix<Tint> SympFormMat = ZeroMatrix<Tint>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    SympFormMat(i, n + i) = 1;
    SympFormMat(n + i, i) = -1;
  }
  auto GetInitialVector = [&]() -> MyMatrix<Tint> {
    int pos = 0;
    std::vector<Tint> ListX(2 * n);
    while (true) {
      MyVector<Tint> V = GetMatrixRow(Mwork, pos);
      if (!IsZeroVector(V))
        return CanonicalizeVector(V);
      pos++;
      if (pos == 2 * n)
        break;
    }
    std::cerr << "Failed to find non-zero vector\n";
    throw TerminalException{1};
  };
  auto GetPairVector = [&](MyVector<Tint> const &w1) -> MyVector<Tint> {
    std::vector<Tint> ListScal(nb_row);
    for (int i_row = 0; i_row < nb_row; i_row++) {
      Tint eScal = 0;
      for (int i = 0; i < 2 * n; i++)
        eScal += w1(i) * Mwork(i_row, i);
      ListScal[i_row] = eScal;
    }
    GCD_int<Tint> eGCD = ComputeGCD_information(ListScal);
    if (T_abs(eGCD.gcd) != 1) {
      std::cerr << "The gcd should be equal to 1\n";
      throw TerminalException{1};
    }
    MyVector<Tint> SumVect = ZeroVector<Tint>(2 * n);
    for (int i_row = 0; i_row < nb_row; i_row++)
      SumVect += GetMatrixRow(Mwork, i_row) * eGCD.Pmat(i_row, 0);
    return SumVect;
  };
  MyMatrix<Tint> CompleteBasis(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    MyVector<Tint> w1 = GetInitialVector();
    MyVector<Tint> wN = GetPairVector(w1);
    CompleteBasis.row(i) = w1;
    CompleteBasis.row(n + i) = wN;
    MyVector<Tint> J_w1 = SympFormMat * w1;
    MyVector<Tint> J_wN = SympFormMat * wN;
    for (int i_row = 0; i_row < nb_row; i_row++) {
      MyVector<Tint> eRow = Mwork.row(i_row);
      Tint scal1 = eRow.dot(J_w1);
      Tint scalN = eRow.dot(J_wN);
      MyVector<Tint> NewRow = eRow - scalN * w1 + scal1 * wN;
      Mwork.row(i_row) = NewRow;
    }
  }
  return CompleteBasis;
}

template <typename T>
MyMatrix<T> CanonicalizeOrderedMatrix_Kernel(const MyMatrix<T> &M) {
  static_assert(is_ring_field<T>::value,
                "Requires T to have inverses in Kernel_ComputeAffineBasis");
  MyMatrix<T> Basis = RowReduction(M);
  MyMatrix<T> M1 = M * Inverse(Basis);
  return ScalarCanonicalizationMatrix(M1);
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, MyMatrix<T>>::type
CanonicalizeOrderedMatrix(MyMatrix<T> const &Input) {
  return CanonicalizeOrderedMatrix_Kernel(Input);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, MyMatrix<T>>::type
CanonicalizeOrderedMatrix(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  MyMatrix<Tfield> OutputF = CanonicalizeOrderedMatrix_Kernel(InputF);
  return UniversalMatrixConversion<T, Tfield>(OutputF);
}

template <typename Tint, typename Tfloat>
std::pair<Tfloat, MyVector<Tint>>
FindBestIntegerApproximation(MyVector<Tfloat> const &V, int const &N) {
  int n = V.size();
  Tfloat SumAbsVal = 0;
  for (int i = 0; i < n; i++)
    SumAbsVal += T_abs(V(i));
  Tfloat MinErr = SumAbsVal;
  MyVector<Tint> MinimalApprox(n);
  //
  for (int iter = 1; iter < N; iter++) {
    MyVector<Tint> CandApprox(n);
    Tfloat SumErr = 0;
    Tfloat TheMult = iter;
    for (int i = 0; i < n; i++) {
      Tfloat val = TheMult * V(i);
      Tint val_i = UniversalNearestScalarInteger<Tint, Tfloat>(val);
      CandApprox(i) = val_i;
      Tfloat valApprox_d =
          UniversalScalarConversion<Tfloat, Tint>(val_i) / TheMult;
      SumErr += T_abs(valApprox_d - V(i));
    }
    //
    if (SumErr < MinErr) {
      MinErr = SumErr;
      MinimalApprox = CandApprox;
    }
  }
  return {MinErr, std::move(MinimalApprox)};
}

template<typename T>
MyMatrix<T> GetZbasisColumn(MyMatrix<T> const& Amat) {
  MyMatrix<T> AmatTr = Amat.transpose();
  MyMatrix<T> Xmat = GetZbasis(AmatTr);
  MyMatrix<T> XmatTr = Xmat.transpose();
  return XmatTr;
}

// We want to consider the equation X A = B
// The equation is potentially underdefined and B in Z^*
// We look for the number d>0 such that for all B in Z^*
// the equation has a solution in Z^* / d
template<typename T>
T GetDenominatorQuotientSolution(MyMatrix<T> const& Amat) {
  MyMatrix<T> AmatRed1 = GetZbasis(Amat);
  MyMatrix<T> AmatRed2 = GetZbasisColumn(AmatRed1);
  return T_abs(DeterminantMat(AmatRed2));
}

template<typename T>
MyMatrix<T> IntegralSpaceSaturation(MyMatrix<T> const& TheSpace) {
  int nbRow = TheSpace.rows();
#ifdef SANITY_CHECK_MATRIX_INT
  if (nbRow != RankMat(TheSpace)) {
    std::cerr << "TheSpace should be full dimensional\n";
    throw TerminalException{1};
  }
#endif
  int nbCol = TheSpace.cols();
  if (nbRow == 0) {
    return MyMatrix<T>(0, nbCol);
  }
  if (nbRow == nbCol) {
    return IdentityMat<T>(nbRow);
  }
  MyMatrix<T> eOrth = NullspaceTrMat(TheSpace);
  MyMatrix<T> TheSpaceInt = NullspaceIntTrMat(eOrth);
  return TheSpaceInt;
}

template<typename T>
MyMatrix<T> RemoveFractionMatrixRows(MyMatrix<T> const& M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyMatrix<T> Mret(nbRow, nbCol);
  for (int u=0; u<nbRow; u++) {
    MyVector<T> eV = GetMatrixRow(M, u);
    MyVector<T> fV = RemoveFractionVector(eV);
    AssignMatrixRow(Mret, u, fV);
  }
  return Mret;
}


/*
  Given a vector eVect and a linear space W, we find
  a vector w = v + W * t  such that w has the minimal
  number of coefficients.
 */
template<typename T>
MyVector<T> EliminateSuperfluousPrimeDenominators(MyVector<T> const& eVect, MyMatrix<T> const& ListVect) {
  int dim = eVect.size();
  if (ListVect.rows() == 0) {
    // Nothing possible really
    return eVect;
  }
  if (RankMat(ListVect) == dim) {
    return ZeroVector<T>(dim);
  }
  MyMatrix<T> ListVect1 = RemoveFractionMatrixRows(ListVect);
  MyMatrix<T> ListVect2 = NullspaceTrMat(ListVect1);
  MyMatrix<T> ListVect3 = RemoveFractionMatrixRows(ListVect2);
  MyMatrix<T> ListVect4 = NullspaceIntTrMat(ListVect3);
  //
  MyMatrix<T> TheCompl = SubspaceCompletionInt(ListVect4, dim);
  MyMatrix<T> TheBasis = Concatenate(TheCompl, ListVect4);
  MyMatrix<T> TheBasisInv = Inverse(TheBasis);
  MyVector<T> eSol = TheBasisInv.transpose() * eVect;
  int dimCompl = TheCompl.rows();
  MyVector<T> eSolRed(dimCompl);
  for (int u=0; u<dimCompl; u++) {
    eSolRed(u) = eSol(u);
  }
  MyVector<T> TheRet = TheCompl.transpose() * eSolRed;
  return TheRet;
}

template<typename T>
MyMatrix<T> EliminateSuperfluousPrimeDenominators_Matrix(MyMatrix<T> const& eMat, std::vector<MyMatrix<T>> const& ListMat) {
  if (ListMat.size() == 0) {
    return eMat;
  }
  int dim = eMat.rows();
  MyVector<T> eVect = MatrixToVector(eMat);
  std::vector<MyVector<T>> ListVect;
  for (auto & eMat : ListMat) {
    MyVector<T> eV = MatrixToVector(eMat);
    ListVect.push_back(eV);
  }
  MyMatrix<T> MatVect = MatrixFromVectorFamily(ListVect);
  MyVector<T> TheVect = EliminateSuperfluousPrimeDenominators(eVect, MatVect);
  MyMatrix<T> TheMat = VectorToMatrix(TheVect, dim);
  return TheMat;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXINT_H_
// clang-format on
