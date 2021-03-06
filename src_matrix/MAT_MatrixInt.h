#ifndef TEMP_MATRIX_INTEGRAL_RING
#define TEMP_MATRIX_INTEGRAL_RING

#include "MAT_Matrix.h"
#include "Boost_bitset.h"
#include "NumberTheory.h"

#undef TRACK_MAXIMUM_SIZE_COEFF
#undef DEBUG

// Now declarations of generic code.
// The code below generally requires the field T to be the ring (or fraction ring) of
// a factorial ring.
// Operations may work for fields and rings as well.
template<typename T>
bool IsIntegralVector(MyVector<T> const& V)
{
  int nbCol=V.size();
  for (int i=0; i<nbCol; i++)
    if (!IsInteger(V(i)))
      return false;
  return true;
}


template<typename T>
bool IsIntegralMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  for (int iCol=0; iCol<nbCol; iCol++)
    for (int iRow=0; iRow<nbRow; iRow++)
      if (!IsInteger(M(iRow, iCol)))
	return false;
  return true;
}


// T must be an integral domain with a norm function
// for computing GCD.
// eMat is a m x n matrix which defines a submodule L of T^n.
// The function computes the index i of L in T^n.
template<typename T>
T Int_IndexLattice(MyMatrix<T> const& eMat)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in Int_IndexLattice");
  size_t iRowF=0, iColF=0;
  MyMatrix<T> eMatW=eMat;
  size_t nbCol=eMat.cols();
  size_t nbRow=eMat.rows();
  std::vector<int> colStat(nbCol,1);
  std::vector<int> rowStat(nbRow,1);
  size_t nbDone=0;
  T TheIndex=1;
  while(true) {
    bool IsFirst=true;
    int MinPivot=0;
    for (size_t iCol=0; iCol<nbCol; iCol++)
      if (colStat[iCol] == 1)
	for (size_t iRow=0; iRow<nbRow; iRow++)
	  if (rowStat[iRow] == 1) {
	    T eVal=eMatW(iRow, iCol);
	    if (eVal != 0) {
	      int eValA=T_Norm(eVal);
	      if (IsFirst) {
		iRowF = iRow;
		iColF = iCol;
		MinPivot=eValA;
	      } else {
		if (eValA < MinPivot) {
		  iRowF = iRow;
		  iColF = iCol;
		  MinPivot=eValA;
		}
	      }
	      IsFirst=false;
	    }
	  }
    if (IsFirst)
      return 0;
#ifdef DEBUG
    if (MinPivot == 0) {
      std::cerr << "Clear error in the code of IndexLattice\n";
      throw TerminalException{1};
    }
#endif
    //    std::cerr << "Before row operations\n";
    T ThePivot=eMatW(iRowF, iColF);
    bool IsFinished=true;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      if (rowStat[iRow] == 1 && iRow != iRowF) {
	T eVal=eMatW(iRow, iColF);
	if (eVal != 0) {
	  IsFinished=false;
	  T TheQ=QuoInt(eVal, ThePivot);
          if (TheQ != 0)
            eMatW.row(iRow) -= TheQ*eMatW.row(iRowF);
	}
      }
    }
    if (IsFinished) {
      colStat[iColF]=0;
      rowStat[iRowF]=0;
      nbDone++;
      TheIndex=TheIndex*ThePivot;
    }
    if (nbDone == nbCol)
      return TheIndex;
  }
}


// This is the return type for the GCD computations
// In input a list of entries x=(x_0, ...., x_m)
// In return we have
// ---gcd: The greatest common dovisor
// ---Pmat: A matrix P unimodulaire such that
//    V P = (gcd, 0, ....., 0)
template<typename T>
struct GCD_int {
  MyMatrix<T> Pmat;
  T gcd;
};


template<typename T>
void WriteGCD_int(std::ostream & os, GCD_int<T> const& eGCD)
{
  os << "gcd=" << eGCD.gcd << "\n";
  os << "Pmat=";
  WriteMatrix(os, eGCD.Pmat);
}


template<typename T>
inline typename std::enable_if<(not is_mpz_class<T>::value),GCD_int<T>>::type ComputePairGcd(T const& m, T const& n)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in ComputePairGcd");
  T f, g, h, fm, gm, hm, q;
  if (n == 0 && m == 0) {
    f=0;
    MyMatrix<T> Pmat=IdentityMat<T>(2);
    return {std::move(Pmat), f};
  }
  if (0 <= m) {
    f=m; fm=1;
  } else {
    f=-m; fm=-1;
  }
  if (0 <= n) {
    g=n; gm=0;
  } else {
    g=-n; gm=0;
  }
  while (g != 0) {
    q = QuoInt( f, g );
    h = g;          hm = gm;
    g = f - q * g;  gm = fm - q * gm;
    f = h;          fm = hm;
  }
  T eCoeff1, eCoeff2;
  if (n == 0) {
    eCoeff1=fm;
    eCoeff2=0;
  } else {
    eCoeff1=fm;
    eCoeff2=(f - fm * m) / n;
  }
  MyMatrix<T> Pmat(2,2);
  Pmat(0,0) = eCoeff1;
  Pmat(1,0) = eCoeff2;
  Pmat(0,1) = -n/f;
  Pmat(1,1) = m/f;
#ifdef DEBUG
  T diff1 = f - Pmat(0,0) * m - Pmat(1,0) * n;
  T diff2 = Pmat(0,1) * m + Pmat(1,1) * n;
  if (diff1 != 0 || diff2 != 0) {
    std::cerr << "A: diff1=" << diff1 << " diff2=" << diff2 << "\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(Pmat), f};
}

template<typename T>
inline typename std::enable_if<is_mpz_class<T>::value,GCD_int<T>>::type ComputePairGcd(T const& m, T const& n)
{
  mpz_class eGCD;
  if (n == 0 && m == 0) {
    eGCD=0;
    MyMatrix<T> Pmat=IdentityMat<T>(2);
    return {std::move(Pmat), eGCD};
  }
  mpz_class s, t;
  mpz_gcdext(eGCD.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), m.get_mpz_t(), n.get_mpz_t());
  //  std::cerr << "n=" << n << " m=" << m << " eGCD=" << eGCD << "\n";
  MyMatrix<T> Pmat(2,2);
  Pmat(0,0) = s;
  Pmat(1,0) = t;
  Pmat(0,1) = -n / eGCD;
  Pmat(1,1) =  m / eGCD;
#ifdef DEBUG
  T diff1 = eGCD - Pmat(0,0) * m - Pmat(1,0) * n;
  T diff2 = Pmat(0,1) * m + Pmat(1,1) * n;
  if (diff1 != 0 || diff2 != 0) {
    std::cerr << "m=" << m << " n=" << n << "\n";
    std::cerr << "s=" << s << " t=" << t << "  eGCD=" << eGCD << "\n";
    std::cerr << "B: diff1=" << diff1 << " diff2=" << diff2 << "\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(Pmat), eGCD};
}



template<typename T>
inline typename std::enable_if<is_mpz_class<T>::value,T>::type KernelGcdPair(T const& a, T const& b)
{
  mpz_class eGCD;
  mpz_gcd(eGCD.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return eGCD;
}

template<typename T>
inline typename std::enable_if<(not is_mpz_class<T>::value),T>::type KernelGcdPair(T const& a, T const& b)
{
  GCD_int<T> eGCD=ComputePairGcd(a, b);
  return eGCD.gcd;
}



template<typename T>
inline typename std::enable_if<is_totally_ordered<T>::value,T>::type GcdPair(T const& a, T const& b)
{
  T eGCD=KernelGcdPair(a,b);
  if (eGCD > 0)
    return eGCD;
  return -eGCD;
}

template<typename T>
inline typename std::enable_if<(not is_totally_ordered<T>::value),T>::type GcdPair(T const& a, T const& b)
{
  return KernelGcdPair(a,b);
}


template<typename T>
inline typename std::enable_if<(not is_mpz_class<T>::value),T>::type KernelLCMpair(T const& a, T const& b)
{
  if (a == 0)
    return b;
  if (b == 0)
    return a;
  return a * b / KernelGcdPair(a,b);
}

template<typename T>
inline typename std::enable_if<is_mpz_class<T>::value,T>::type KernelLCMpair(T const& a, T const& b)
{
  mpz_class eLCM;
  mpz_lcm(eLCM.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return eLCM;
}


template<typename T>
inline typename std::enable_if<(not is_totally_ordered<T>::value),T>::type LCMpair(T const& a, T const& b)
{
  return KernelLCMpair(a, b);
}

template<typename T>
inline typename std::enable_if<is_totally_ordered<T>::value,T>::type LCMpair(T const& a, T const& b)
{
  T eLCM = KernelLCMpair(a, b);
  if (eLCM > 0)
    return eLCM;
  return -eLCM;
}



template<typename T>
T ComputeLCM(std::vector<T> const& eVect)
{
  int siz=eVect.size();
  T eLCM=1;
  for (int i=0; i<siz; i++)
    eLCM = LCMpair(eLCM, eVect[i]);
  return eLCM;
}






template<typename T>
struct FractionMatrix {
  T TheMult;
  MyMatrix<T> TheMat;
};


template<typename T>
FractionMatrix<T> RemoveFractionMatrixPlusCoeff(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  using Tring = typename underlying_ring<T>::ring_type;
  Tring eLCM_ring = 1;
  // iRow is inner loop because of cache locality
  for (int iCol=0; iCol<nbCol; iCol++)
    for (int iRow=0; iRow<nbRow; iRow++)
      eLCM_ring = LCMpair(eLCM_ring, GetDenominator_z(M(iRow,iCol)));
  T eLCM = eLCM_ring;
  MyMatrix<T> Mret = eLCM * M;
  return {eLCM, std::move(Mret)};
}

template<typename T>
MyMatrix<T> RemoveFractionMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int nbCol=M.cols();
  using Tring = typename underlying_ring<T>::ring_type;
  Tring eLCM_ring = 1;
  // iRow is inner loop because of cache locality
  for (int iCol=0; iCol<nbCol; iCol++)
    for (int iRow=0; iRow<nbRow; iRow++)
      eLCM_ring = LCMpair(eLCM_ring, GetDenominator_z(M(iRow,iCol)));
  T eLCM = eLCM_ring;
  return eLCM * M;
}

template<typename T>
struct FractionVector {
  T TheMult;
  MyVector<T> TheVect;
};


template<typename T>
FractionVector<T> RemoveFractionVectorPlusCoeff(MyVector<T> const& V)
{
  int n=V.size();
  std::vector<T> eVect(n);
  T eLCM=1;
  for (int i=0; i<n; i++)
    eLCM = LCMpair(eLCM, GetDenominator(V(i)));
  MyVector<T> Vret = eLCM * V;
  return {eLCM, std::move(Vret)};
}

template<typename T>
inline typename std::enable_if<is_implementation_of_Z<T>::value,MyVector<T>>::type CanonicalizeVector(MyVector<T> const& V)
{
  return RemoveFractionVectorPlusCoeff(V).TheVect;
}


template<typename T>
inline typename std::enable_if<(not is_float_arithmetic<T>::value),MyVector<T>>::type RemoveFractionVector(MyVector<T> const& V)
{
  return RemoveFractionVectorPlusCoeff(V).TheVect;
}


// In this function we do not care about the invertible elements of the ring
// Below is the Z-case where we just have +1, -1.
template<typename T>
inline typename std::enable_if<is_totally_ordered<T>::value,MyVector<T>>::type CanonicalizeVectorToInvertible(MyVector<T> const& V)
{
  MyVector<T> eVect = CanonicalizeVector(V);
  int len=eVect.size();
  int FirstNZ = -1;
  for (int u=0; u<len; u++)
    if (eVect(u) != 0 && FirstNZ == -1)
      FirstNZ = u;
  if (FirstNZ == -1)
    return eVect;
  if (eVect(FirstNZ) < 0)
    return -eVect;
  return eVect;
}









int IsVectorPrimitive(MyVector<int> const& TheV)
{
  size_t n=TheV.size();
  int TheGCD=TheV(0);
  for (size_t i=1; i<n; i++) {
    int eValI=TheV(i);
    GCD_int<int> eRec=ComputePairGcd(TheGCD, eValI);
    TheGCD=eRec.gcd;
  }
  if (abs(TheGCD) == 1)
    return 1;
  return 0;
}






template<typename T>
GCD_int<T> ComputeGCD_information(std::vector<T> const& ListX)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in ComputeGCD_information");
  size_t siz=ListX.size();
  if (siz == 1) {
    T gcd=ListX[0];
    MyMatrix<T> Pmat=IdentityMat<T>(1);
    GCD_int<T> eGCD_int{std::move(Pmat), gcd};
    return eGCD_int;
  }
  if (siz == 2)
    return ComputePairGcd(ListX[0], ListX[1]);
  std::vector<T> ListXred(ListX.begin(), ListX.end()-1);
  GCD_int<T> eGCD_int=ComputeGCD_information(ListXred);
  GCD_int<T> eGCD2=ComputePairGcd(eGCD_int.gcd, ListX[siz-1]);
  MyMatrix<T> Pmat=MyMatrix<T>(siz,siz);
  // 1 : the column for the GCD
  for (size_t i=0; i<siz-1; i++)
    Pmat(i,0) = eGCD2.Pmat(0,0) * eGCD_int.Pmat(i,0);
  Pmat(siz-1, 0) = eGCD2.Pmat(1,0);
  // 2 : The columns from the previous block
  for (size_t i=0; i<siz-2; i++) {
    for (size_t j=0; j<siz-1; j++)
      Pmat(j, i+1) = eGCD_int.Pmat(j, i+1);
    Pmat(siz-1, i+1) = 0;
  }
  // 3 : The zero columns
  for (size_t i=0; i<siz-1; i++)
    Pmat(i, siz-1) = eGCD2.Pmat(0,1) * eGCD_int.Pmat(i,0);
  Pmat(siz-1, siz-1) = eGCD2.Pmat(1,1);
  //
  return {std::move(Pmat), eGCD2.gcd};
}


// See https://en.wikipedia.org/wiki/Hermite_normal_form
// for the Row-style Hermite normal form
//
// The matrix M is rewritten as U A = H
// * H is upper triangular with H_{ij}=0 for i > j.
// * The leading coefficient of H of a row is strictly to the right of the above one.
// * The elements below pivots are zero and elements above pivots are nonnegative and strictly smaller than the pivot.
//
// This implementation is fairly naive. But its advantage is that it follows the
// same Template structure as other codes in this library.
// The output of the call of ComputeRowHermiteNormalForm ( M )
// The return (U,H) is U an unimodular matrix with U M = H
template<typename T, typename F>
void ComputeRowHermiteNormalForm_Kernel(MyMatrix<T> & H, F f)
{
  size_t nbRow=H.rows();
  size_t nbCol=H.cols();
  Face StatusRow(nbRow);
  //
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  T MaxSizeCoeff = 0;
#endif

  size_t TopPosition=0;
  MyMatrix<T> TheMat1 = IdentityMat<T>(nbRow);
  for (size_t iCol=0; iCol<nbCol; iCol++) {
    std::vector<T> ListX;
    std::vector<size_t> ListIdx;
    bool HasNonZero=false;
    for (size_t iRow=0; iRow<nbRow; iRow++)
      if (StatusRow[iRow] == 0) {
        ListIdx.push_back(iRow);
        T eVal = H(iRow,iCol);
        ListX.push_back(eVal);
        if (eVal != 0)
          HasNonZero=true;
      }
    if (HasNonZero) {
      //
      // Ensuring that the column has a pivot and that everything below is ZERO
      size_t siz = ListIdx.size();
      GCD_int<T> eGCD = ComputeGCD_information(ListX);
      for (size_t i=0; i<siz; i++)
        for (size_t j=0; j<siz; j++)
          TheMat1(ListIdx[i],ListIdx[j]) = eGCD.Pmat(j,i);
      auto fct1=[&](MyMatrix<T> & m) -> void {
        m = TheMat1 * m;
      };
      f(fct1);
#ifdef TRACK_MAXIMUM_SIZE_COEFF
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(TheMat1));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(U));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(H));
#endif
      for (size_t i=0; i<siz; i++)
        for (size_t j=0; j<siz; j++) {
          if (i != j)
            TheMat1(ListIdx[i],ListIdx[j]) = 0;
          else
            TheMat1(ListIdx[i],ListIdx[j]) = 1;
        }
      //
      // Ensuring that the pivot is strictly positive
      // (in the case of integer. For other rings this is a different story)
      T eCanUnit = CanonicalizationUnit(H(TopPosition, iCol));
      if (eCanUnit != 1) {
        auto fct2=[&](MyMatrix<T> & m) -> void {
          m.row(TopPosition) = eCanUnit * m.row(TopPosition);
        };
        f(fct2);
      }
      //
      // Putting the coefficients over the pivot
      T ThePivot = H(TopPosition, iCol);
      for (size_t iRow=0; iRow<TopPosition; iRow++) {
        T eVal = H(iRow, iCol);
        T TheQ = QuoInt(eVal, ThePivot);
        if (TheQ != 0) {
          auto fct3=[&](MyMatrix<T> & m) -> void {
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
      StatusRow[TopPosition]=1;
      TopPosition++;
     }
  }
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  std::cerr << "MaxSizeCoeff = " << MaxSizeCoeff << "\n";
#endif
}

template<typename T>
std::pair<MyMatrix<T>, MyMatrix<T>> ComputeRowHermiteNormalForm(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  MyMatrix<T> H = M;
  MyMatrix<T> U = IdentityMat<T>(nbRow);
  auto f=[&](auto g) -> void {
    g(H);
    g(U);
  };
  ComputeRowHermiteNormalForm_Kernel(H, f);
  return {std::move(U), std::move(H)};
}



template<typename T, typename F>
void ComputeColHermiteNormalForm_Kernel(MyMatrix<T> & H, F f)
{
  size_t nbRow=H.rows();
  size_t nbCol=H.cols();
  Face StatusRow(nbCol);
  //
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  T MaxSizeCoeff = 0;
#endif

  size_t TopPosition=0;
  MyMatrix<T> TheMat1 = IdentityMat<T>(nbCol);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    std::vector<T> ListX;
    std::vector<size_t> ListIdx;
    bool HasNonZero=false;
    for (size_t iCol=0; iCol<nbCol; iCol++)
      if (StatusRow[iCol] == 0) {
        ListIdx.push_back(iCol);
        T eVal = H(iRow,iCol);
        ListX.push_back(eVal);
        if (eVal != 0)
          HasNonZero=true;
      }
    if (HasNonZero) {
      //
      // Ensuring that the column has a pivot and that everything below is ZERO
      size_t siz = ListIdx.size();
      GCD_int<T> eGCD = ComputeGCD_information(ListX);
      for (size_t i=0; i<siz; i++)
        for (size_t j=0; j<siz; j++)
          TheMat1(ListIdx[i],ListIdx[j]) = eGCD.Pmat(i,j);
      auto fct1=[&](MyMatrix<T> & m) -> void {
        m = m * TheMat1;
      };
      f(fct1);
#ifdef TRACK_MAXIMUM_SIZE_COEFF
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(TheMat1));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(U));
      MaxSizeCoeff = T_max(MaxSizeCoeff, Linfinity_norm_mat(H));
#endif
      for (size_t i=0; i<siz; i++)
        for (size_t j=0; j<siz; j++) {
          if (i != j)
            TheMat1(ListIdx[i],ListIdx[j]) = 0;
          else
            TheMat1(ListIdx[i],ListIdx[j]) = 1;
        }
      //
      // Ensuring that the pivot is strictly positive
      // (in the case of integer. For other rings this is a different story)
      T eCanUnit = CanonicalizationUnit(H(iRow, TopPosition));
      if (eCanUnit != 1) {
        auto fct2=[&](MyMatrix<T> & m) -> void {
          m.col(TopPosition) = eCanUnit * m.col(TopPosition);
        };
        f(fct2);
      }
      //
      // Putting the coefficients over the pivot
      T ThePivot = H(iRow, TopPosition);
      for (size_t iCol=0; iCol<TopPosition; iCol++) {
        T eVal = H(iRow, iCol);
        T TheQ = QuoInt(eVal, ThePivot);
        if (TheQ != 0) {
          auto fct3=[&](MyMatrix<T> & m) -> void {
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
      StatusRow[TopPosition]=1;
      TopPosition++;
     }
  }
#ifdef TRACK_MAXIMUM_SIZE_COEFF
  std::cerr << "MaxSizeCoeff = " << MaxSizeCoeff << "\n";
#endif
}



template<typename T>
std::pair<MyMatrix<T>, MyMatrix<T>> ComputeColHermiteNormalForm(MyMatrix<T> const& M)
{
  int nbCol=M.cols();
  MyMatrix<T> H = M;
  MyMatrix<T> U = IdentityMat<T>(nbCol);
  auto f=[&](auto g) -> void {
    g(H);
    g(U);
  };
  ComputeColHermiteNormalForm_Kernel(H, f);
  return {std::move(U), std::move(H)};
}

template<typename T>
MyMatrix<T> ComputeColHermiteNormalForm_second(MyMatrix<T> const& M)
{
  MyMatrix<T> H = M;
  auto f=[&](auto g) -> void {
    g(H);
  };
  ComputeColHermiteNormalForm_Kernel(H, f);
  return H;
}



template<typename T>
void SwitchRow(MyMatrix<T> & eMat, int const& iRow, int const& jRow)
{
  int nbCol=eMat.cols();
  if (iRow == jRow)
    return;
  for (int iCol=0; iCol<nbCol; iCol++) {
    T eVal1=eMat(iRow, iCol);
    T eVal2=eMat(jRow, iCol);
    eMat(iRow, iCol)=eVal2;
    eMat(jRow, iCol)=eVal1;
  }
}

template<typename T>
void INT_ClearColumn(MyMatrix<T> & eMat, size_t const& iCol, size_t const& MinAllowedRow, size_t & iRowFound)
{
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  size_t nbRow=eMat.rows();
  while(true) {
    Treal MinVal=-1;
    size_t nbFound=0;
    for (size_t iRow=MinAllowedRow; iRow<nbRow; iRow++) {
      T eVal=eMat(iRow, iCol);
      if (eVal != 0) {
	Treal AbsEVal=T_NormGen(eVal);
	if (nbFound == 0) {
	  MinVal=AbsEVal;
	  iRowFound=iRow;
	} else {
	  if (AbsEVal < MinVal) {
	    MinVal=AbsEVal;
	    iRowFound=iRow;
	  }
	}
	nbFound++;
      }
    }
#ifdef DEBUG
    if (nbFound == 0) {
      std::cerr << "The column is zero. No work possible\n";
      throw TerminalException{1};
    }
#endif
    T ThePivot=eMat(iRowFound, iCol);
    for (size_t iRow=0; iRow<nbRow; iRow++)
      if (iRow != iRowFound) {
	T eVal=eMat(iRow, iCol);
	T TheQ=QuoInt(eVal, ThePivot);
        if (TheQ != 0)
          eMat.row(iRow) -= TheQ*eMat.row(iRowFound);
      }
    if (nbFound == 1)
      return;
  }
}


template<typename T>
bool IsColumnNonEmpty(MyMatrix<T> const& eMat, int const& minAllowed, int const& iCol)
{
  int nbRow=eMat.rows();
  for (int iRow=minAllowed; iRow<nbRow; iRow++) {
    T eVal=eMat(iRow, iCol);
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
template<typename T>
MyMatrix<T> NullspaceIntTrMat(MyMatrix<T> const& eMat)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in NullspaceIntTrMat");
  MyMatrix<T> eMatW=eMat;
  size_t nbCol=eMat.cols();
  std::vector<size_t> ListIndex;
  std::vector<size_t> ListNonIndex;
  size_t eRank=0;
  for (size_t iCol=0; iCol<nbCol; iCol++)
    if (IsColumnNonEmpty(eMatW, eRank, iCol)) {
      ListIndex.push_back(iCol);
      size_t iRowFound=444;
      INT_ClearColumn(eMatW, iCol, eRank, iRowFound);
      SwitchRow(eMatW, eRank, iRowFound);
      eRank++;
    } else {
      ListNonIndex.push_back(iCol);
    }
  size_t dimSpace=ListNonIndex.size();
  std::vector<std::vector<T>> TheBasis;
  for (size_t i=0; i<dimSpace; i++) {
    std::vector<T> eVect;
    for (size_t j=0; j<dimSpace; j++) {
      if (i == j) {
	eVect.push_back(1);
      } else {
	eVect.push_back(0);
      }
    }
    TheBasis.push_back(eVect);
  }
  for (size_t iRank=0; iRank<eRank; iRank++) {
    size_t iRow=eRank-1-iRank;
    size_t iCol=ListIndex[iRow];
    std::vector<T> ListX;
    T eVal=eMatW(iRow, iCol);
    ListX.push_back(eVal);
    std::vector<size_t> ListRelIndex;
    for (size_t jRow=iRow+1; jRow<eRank; jRow++) {
      size_t jCol=ListIndex[jRow];
      ListRelIndex.push_back(jCol);
    }
    for (auto & jCol : ListNonIndex) {
      ListRelIndex.push_back(jCol);
    }
    size_t sizRelIndex=ListRelIndex.size();
    for (size_t iVect=0; iVect<dimSpace; iVect++) {
      std::vector<T> eVect=TheBasis[iVect];
      T eSum=0;
      for (size_t iRel=0; iRel<sizRelIndex; iRel++) {
	size_t jCol=ListRelIndex[iRel];
	T fVal=eMatW(iRow, jCol);
	eSum += eVect[iRel]*fVal;
      }
      ListX.push_back(eSum);
    }
    GCD_int<T> eGCD=ComputeGCD_information(ListX);
    std::vector<std::vector<T>> NewBasis;
    for (size_t iVect=0; iVect<dimSpace; iVect++) {
      std::vector<T> eVectNew(sizRelIndex+1,0);
      eVectNew[0]=eGCD.Pmat(0,iVect+1);
      for (size_t i=1; i<=dimSpace; i++) {
	T fVal=eGCD.Pmat(i,iVect+1);
	std::vector<T> basVect=TheBasis[i-1];
	for (size_t j=0; j<sizRelIndex; j++)
	  eVectNew[j+1] += fVal*basVect[j];
      }
      NewBasis.push_back(eVectNew);
    }
    TheBasis=NewBasis;
  }
  MyMatrix<T> retNSP(dimSpace,nbCol);
  for (size_t iVect=0; iVect<dimSpace; iVect++) {
    std::vector<T> eVect=TheBasis[iVect];
    size_t idx=0;
    for (size_t iRank=0; iRank<eRank; iRank++) {
      size_t iCol=ListIndex[iRank];
      retNSP(iVect, iCol)=eVect[idx];
      idx++;
    }
    for (size_t iDim=0; iDim<dimSpace; iDim++) {
      size_t iCol=ListNonIndex[iDim];
      retNSP(iVect, iCol)=eVect[idx];
      idx++;
    }
  }
#ifdef DEBUG
  size_t nbRow=eMat.rows();
  for (size_t iVect=0; iVect<dimSpace; iVect++)
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      T eSum=0;
      for (size_t iCol=0; iCol<nbCol; iCol++)
	eSum += eMat(iRow, iCol) * retNSP(iVect, iCol);
      if (eSum != 0) {
	std::cerr << "There are remaining errors in NullspaceIntTrMat\n";
	throw TerminalException{1};
      }
    }
#endif
  return retNSP;
}

template<typename T>
MyMatrix<T> NullspaceIntMat(MyMatrix<T> const& eMat)
{
  return NullspaceIntTrMat(TransposedMat(eMat));
}



// Given a vector v of T^n which is primitive
// the function returns a family of vector v1, ...., v(n-1)
// such that (v1, ...., v(n-1), v) is a T-basis of T^n
template<typename T>
MyMatrix<T> ComplementToBasis(MyVector<T> const& TheV)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in ComplementToBasis");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  std::vector<int> OperCol1;
  std::vector<int> OperCol2;
  std::vector<T> OperCoef;
  Treal AbsVal=-400;
  MyVector<T> TheVcopy=TheV;
  int n=TheV.size();
  int idxSelect;
  while(true) {
    bool IsFirst=true;
    int nbDiffZero=0;
    idxSelect=-1;
    for (int i=0; i<n; i++)
      if (TheVcopy(i) != 0) {
	nbDiffZero++;
	Treal hVal=T_NormGen(TheVcopy(i));
	if (IsFirst) {
	  IsFirst=false;
	  AbsVal=hVal;
	  idxSelect=i;
	} else {
	  if (hVal < AbsVal)
	    {
	      idxSelect=i;
	      AbsVal=hVal;
	    }
	}
      }
#ifdef DEBUG
    if (idxSelect == -1) {
      std::cerr << "Inconsistency in computation of value\n";
      throw TerminalException{1};
    }
#endif
    if (nbDiffZero == 1) {
#ifdef DEBUG
      if (AbsVal != 1) {
	std::cerr << "Wrong value for AbsVal\n";
	throw TerminalException{1};
      }
#endif
      break;
    }
    for (int j=0; j<n; j++)
      if (TheVcopy(j) != 0 && idxSelect !=j) {
        T eVal1=TheVcopy(idxSelect);
        T eVal2=TheVcopy(j);
        std::pair<T,T> ePair = ResQuoInt(eVal2, eVal1);
        TheVcopy(j)=ePair.first;
        OperCol1.push_back(idxSelect);
        OperCol2.push_back(j);
        OperCoef.push_back(ePair.second);
      }
  }
  size_t nbOper=OperCol1.size();
  MyMatrix<T> TheReturn=MyMatrix<T>(n-1, n);
  int idx=0;
  for (int i=0; i<n; i++)
    if (i != idxSelect) {
      for (int j=0; j<n; j++) {
	T eVal=0;
	if (j == i)
	  eVal=1;
	TheReturn(idx, j)=eVal;
      }
      idx++;
    }
  for (size_t iOper=0; iOper<nbOper; iOper++) {
    size_t jOper=nbOper-1-iOper;
    int i1=OperCol1[jOper];
    int i2=OperCol2[jOper];
    T eCoeff=OperCoef[jOper];
    for (int iRow=0; iRow<n-1; iRow++) {
      T eVal=TheReturn(iRow, i2) + eCoeff*TheReturn(iRow, i1);
      TheReturn(iRow, i2)=eVal;
    }
    T eVal=TheVcopy(i2) + eCoeff*TheVcopy(i1);
    TheVcopy(i2)=eVal;
  }
#ifdef DEBUG
  if (!TestEquality(TheVcopy, TheV)) {
    std::cerr << "TheVcopy =";
    WriteVector(std::cerr, TheVcopy);
    std::cerr << "TheV =";
    WriteVector(std::cerr, TheV);
    std::cerr << "Bookkeeping error\n";
    throw TerminalException{1};
  }
#endif
  return TheReturn;
}


// We have two matrices M1 and M2 and we check if they define
// the same subspace of T^n
template<typename T>
bool TestEqualitySpaces(MyMatrix<T> const& M1, MyMatrix<T> const& M2)
{
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int idxSelect=-1;
  size_t k=M1.rows();
  size_t n=M1.cols();
  MyMatrix<T> M1copy=M1;
  MyMatrix<T> M2copy=M2;
  std::vector<int> StatusRow(k, 0);
  int idxSearch=0;
  for (size_t iK=0; iK<k; iK++) {
    while(true) {
      int nbDiff=0;
      for (size_t j=0; j<k; j++)
	if (StatusRow[j] == 0)
	  if (M1copy(j, idxSearch) != 0)
	    nbDiff++;
      if (nbDiff > 0)
	break;
      idxSearch++;
    }
    while(true) {
      int nbDiff=0;
      bool IsFirst=true;
      Treal AbsVal=0;
      for (size_t j=0; j<k; j++) {
	T eVal=M1copy(j,idxSearch);
	if (eVal != 0 && StatusRow[j] == 0) {
	  nbDiff++;
	  Treal hVal=T_NormGen(eVal);
	  if (IsFirst) {
	    IsFirst=false;
	    AbsVal=hVal;
	    idxSelect=j;
	  } else {
	    if (hVal < AbsVal) {
	      idxSelect=j;
	      AbsVal=hVal;
	    }
	  }
	}
      }
      if (nbDiff == 1)
	break;
      for (size_t j=0; j<k; j++)
	if (j != idxSelect) {
	  T eVal1=M1copy(idxSelect, idxSearch);
	  T eVal2=M1copy(j, idxSearch);
	  T TheQ=QuoInt(eVal2, eVal1);
          if (TheQ != 0)
            M1copy.row(j) -= TheQ*M1copy.row(idxSelect);
	}
    }
    StatusRow[idxSelect]=1;
    T eVal1=M1copy(idxSelect, idxSearch);
    for (size_t j=0; j<k; j++) {
      T eVal2=M2copy(j, idxSearch);
      std::pair<T,T> ePair = ResQuoInt(eVal2, eVal1);
      if (ePair.first != 0)
	return false;
      if (ePair.second != 0)
        M2copy.row(j)=M2copy.row(j) - ePair.second*M1copy.row(idxSelect);
    }
    idxSearch++;
  }
  for (size_t j=0; j<k; j++)
    for (size_t i=0; i<n; i++)
      if (M2copy(j, i) != 0)
	return false;
  return true;
}


template<typename T>
int PositionSubspace(std::vector<MyMatrix<T>> const& ListSubspace, MyMatrix<T> const& OneSubspace)
{
  int nbSpace=ListSubspace.size();
  for (int iSub=0; iSub<nbSpace; iSub++) {
    bool test=TestEqualitySpaces(ListSubspace[iSub], OneSubspace);
    if (test)
      return iSub;
  }
  return -1;
}

template<typename Ti, typename Td>
void FuncInsertSubspace(std::vector<MyMatrix<Ti>> & ListSubspace, std::vector<Td> & ListDet, MyMatrix<Ti> const& OneSubspace, Td const& OneDet)
{
  if (PositionSubspace(ListSubspace, OneSubspace) == -1) {
    ListSubspace.push_back(OneSubspace);
    ListDet.push_back(OneDet);
  }
}


// T1 is integer type and T2 is real kind type
template<typename Ti, typename Td>
void ReorderSubspaceDet(std::vector<MyMatrix<Ti>> &ListSubspace, std::vector<Td> &ListDet)
{
  int nbSub=ListSubspace.size();
  for (int i=0; i<nbSub-1; i++)
    for (int j=i+1; j<nbSub; j++)
      if (ListDet[i] > ListDet[j])
	{
	  MyMatrix<Ti> eSub=ListSubspace[i];
	  Td eDet=ListDet[i];
	  ListSubspace[i]=ListSubspace[j];
	  ListDet[i]=ListDet[j];
	  ListSubspace[j]=eSub;
	  ListDet[j]=eDet;
	}
}

/* test if ListSub2 is a subset of ListSub1 */
template<typename T>
bool TestInclusionFamilySubspace(std::vector<MyMatrix<T>> const& ListSub1, std::vector<MyMatrix<T>> const& ListSub2)
{
  int nbSub2=ListSub2.size();
  for (int i=0; i<nbSub2; i++)
    if (PositionSubspace(ListSub1, ListSub2[i]) == -1)
      return false;
  return true;
}


template<typename T>
bool TestEqualityFamilySubspace(std::vector<MyMatrix<T>> const& ListSub1, std::vector<MyMatrix<T>> const& ListSub2)
{
  int nbSub1=ListSub1.size();
  int nbSub2=ListSub2.size();
  if (nbSub1 != nbSub2)
    return false;
  return TestInclusionFamilySubspace(ListSub1, ListSub2);
}



template<typename T>
MyMatrix<T> GetNoncontainedSubspace(std::vector<MyMatrix<T>> const& ListSubBig, std::vector<MyMatrix<T>> const& ListSubSma)
{
  int nbSubBig=ListSubBig.size();
  int nbSubSma=ListSubSma.size();
  for (int i=0; i<nbSubBig; i++)
    if (PositionSubspace(ListSubSma, ListSubBig[i]) == -1)
      return ListSubBig[i];
  std::cerr << "We did not find any subspace\n";
  std::cerr << "nbSubBig=" << nbSubBig << " nbSumSma=" << nbSubSma << "\n";
  throw TerminalException{1};
}


template<typename T>
struct ResultSolutionIntMat {
  bool TheRes;
  MyVector<T> eSol;
};



// Find an integral solution of the equation Y = X A
// if it exists.
template<typename T>
ResultSolutionIntMat<T> SolutionIntMat(MyMatrix<T> const& TheMat, MyVector<T> const& TheVect)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in SolutionIntMat");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int nbDiff;
  int nbVect=TheMat.rows();
  int nbCol=TheMat.cols();
  if (nbVect == 0) {
    MyVector<T> eSol;
    if (IsZeroVector(TheVect)) {
      return {true, std::move(eSol)};
    } else {
      return {false, {}};
    }
  }
  MyVector<T> eSol=ZeroVector<T>(nbVect);
  MyMatrix<T> eEquivMat=IdentityMat<T>(nbVect);
  MyMatrix<T> TheMatWork=TheMat;
  MyVector<T> TheVectWork=TheVect;
  std::vector<int> VectStatus(nbVect,1);
  for (int i=0; i<nbCol; i++) {
    int iVectFound=-1;
    while(true) {
      bool IsFirst=true;
      Treal MinValue=0;
      nbDiff=0;
      for (int iVect=0; iVect<nbVect; iVect++)
        if (VectStatus[iVect] == 1) {
          T prov1=TheMatWork(iVect, i);
          Treal eNorm=T_NormGen(prov1);
	  if (prov1 != 0) {
	    nbDiff++;
	    if (IsFirst) {
	      IsFirst=false;
	      MinValue=eNorm;
	      iVectFound=iVect;
	    } else {
	      if (eNorm < MinValue) {
		MinValue=eNorm;
		iVectFound=iVect;
	      }
	    }
	  }
	}
      if (nbDiff == 1 || nbDiff == 0)
	break;
#ifdef DEBUG
      if (MinValue == 0) {
	std::cerr << "MinValue should not be zero\n";
	throw TerminalException{1};
      }
#endif
      for (int iVect=0; iVect<nbVect; iVect++)
	if (VectStatus[iVect] == 1 && iVect != iVectFound) {
	  T prov1b=TheMatWork(iVectFound, i);
	  T prov2=TheMatWork(iVect, i);
	  T TheQ=QuoInt(prov2, prov1b);
          if (TheQ != 0) {
            TheMatWork.row(iVect) -= TheQ*TheMatWork.row(iVectFound);
            eEquivMat.row(iVect) -= TheQ*eEquivMat.row(iVectFound);
          }
	}
    }
    if (nbDiff == 1) {
#ifdef DEBUG
      if (iVectFound == -1) {
        std::cerr << "Clear error in the program\n";
	throw TerminalException{1};
      }
#endif
      VectStatus[iVectFound]=0;
      T prov1=TheVectWork(i);
      T prov2=TheMatWork(iVectFound, i);
      T TheQ=QuoInt(prov1, prov2);
      if (TheQ != 0) {
        for (int j=0; j<nbCol; j++)
          TheVectWork(j) -= TheQ*TheMatWork(iVectFound, j);
        for (int iVect=0; iVect<nbVect; iVect++)
          eSol(iVect) += TheQ*eEquivMat(iVectFound,iVect);
      }
    }
    if (TheVectWork(i) != 0)
      return {false, {}};
  }
  return {true, std::move(eSol)};
}


template<typename T>
struct CanSolIntMat {
  std::vector<int> ListRow;
  MyMatrix<T> TheMatWork;
  MyMatrix<T> eEquivMat;
};




template<typename T>
CanSolIntMat<T> ComputeCanonicalFormFastReduction(MyMatrix<T> const& TheMat)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in ComputeCanonicalFormFastReduction");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int nbDiff;
  int nbVect=TheMat.rows();
  int nbCol=TheMat.cols();
#ifdef DEBUG
  if (nbVect == 0) {
    std::cerr << "Need to write the code here\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> eEquivMat=IdentityMat<T>(nbVect);
  MyMatrix<T> TheMatWork=TheMat;
  std::vector<int> VectStatus(nbVect,1);
  std::vector<int> ListRow(nbCol);
  for (int i=0; i<nbCol; i++) {
    int iVectFound=-1;
    while(true) {
      bool IsFirst=true;
      Treal MinValue=0;
      nbDiff=0;
      for (int iVect=0; iVect<nbVect; iVect++)
        if (VectStatus[iVect] == 1) {
          T prov1=TheMatWork(iVect, i);
          Treal eNorm=T_NormGen(prov1);
	  if (prov1 != 0) {
	    nbDiff++;
	    if (IsFirst) {
	      IsFirst=false;
	      MinValue=eNorm;
	      iVectFound=iVect;
	    } else {
	      if (eNorm < MinValue) {
		MinValue=eNorm;
		iVectFound=iVect;
	      }
	    }
	  }
	}
      if (nbDiff == 1 || nbDiff == 0)
	break;
#ifdef DEBUG
      if (MinValue == 0) {
	std::cerr << "MinValue should not be zero\n";
	throw TerminalException{1};
      }
#endif
      for (int iVect=0; iVect<nbVect; iVect++)
	if (VectStatus[iVect] == 1 && iVect != iVectFound) {
	  T prov1b=TheMatWork(iVectFound, i);
	  T prov2=TheMatWork(iVect, i);
	  T TheQ=QuoInt(prov2, prov1b);
          if (TheQ != 0) {
            TheMatWork.row(iVect) -= TheQ*TheMatWork.row(iVectFound);
            eEquivMat.row(iVect)  -= TheQ*eEquivMat.row(iVectFound);
          }
	}
    }
    int eVal;
    if (nbDiff == 1) {
      eVal=iVectFound;
#ifdef DEBUG
      if (iVectFound == -1) {
        std::cerr << "Clear error in the program\n";
	throw TerminalException{1};
      }
#endif
      VectStatus[iVectFound]=0;
    } else {
      eVal=-1;
    }
    ListRow[i]=eVal;
  }
  return {std::move(ListRow), std::move(TheMatWork), std::move(eEquivMat)};
}

template<typename T>
bool CanTestSolutionIntMat(CanSolIntMat<T> const& eCan, MyVector<T> const& TheVect)
{
  int nbCol=eCan.TheMatWork.cols();
  MyVector<T> TheVectWork=TheVect;
  for (int i=0; i<nbCol; i++) {
    int iRow=eCan.ListRow[i];
    if (iRow >= 0) {
      T prov1=TheVectWork(i);
      T prov2=eCan.TheMatWork(iRow, i);
      T TheQ=QuoInt(prov1, prov2);
      if (TheQ != 0) {
        for (int j=0; j<nbCol; j++)
          TheVectWork(j) -= TheQ*eCan.TheMatWork(iRow, j);
      }
    }
    if (TheVectWork(i) != 0)
      return false;
  }
  return true;
}

template<typename T>
ResultSolutionIntMat<T> CanSolutionIntMat(CanSolIntMat<T> const& eCan, MyVector<T> const& TheVect)
{
  int nbVect=eCan.TheMatWork.rows();
  int nbCol=eCan.TheMatWork.cols();
  MyVector<T> TheVectWork=TheVect;
  MyVector<T> eSol=ZeroVector<T>(nbVect);
  for (int i=0; i<nbCol; i++) {
    int iRow=eCan.ListRow[i];
    if (iRow >= 0) {
      T prov1=TheVectWork(i);
      T prov2=eCan.TheMatWork(iRow, i);
      T TheQ=QuoInt(prov1, prov2);
      if (TheQ != 0) {
        for (int j=0; j<nbCol; j++)
          TheVectWork(j) -= TheQ*eCan.TheMatWork(iRow, j);
        for (int iVect=0; iVect<nbVect; iVect++)
          eSol(iVect) += TheQ*eCan.eEquivMat(iRow,iVect);
      }
    }
    if (TheVectWork(i) != 0)
      return {false, {}};
  }
  return {true, std::move(eSol)};
}



template<typename T>
struct BasisReduction {
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
template<typename T>
BasisReduction<T> ComputeBasisReduction(MyMatrix<T> const& TheBasis)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in ComputeBasisReduction");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  size_t nbCol=TheBasis.cols();
  size_t nbRow=TheBasis.rows();
  std::vector<int> colStat(nbCol, 1);
  std::vector<int> rowStat(nbRow, 1);
  MyMatrix<T> TheBasisReduced = TheBasis;
  MyMatrix<T> Pmat=IdentityMat<T>(nbCol);
  MyMatrix<T> TheBasisReord(nbRow, nbCol);
  std::vector<int> IdxVector;
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto FindMinGCDrow=[&](int const& iRank) -> size_t {
    size_t iRowSearch = miss_val;
    bool IsFirst=true;
    Treal AbsVal=0;
    for (size_t iRow=0; iRow<nbRow; iRow++)
      if (rowStat[iRow] == 1) {
	std::vector<T> eRowRed(nbCol - iRank);
	for (size_t iCol=0; iCol<nbCol - iRank; iCol++) {
	  T eVal=TheBasisReduced(iRow, iCol+iRank);
	  eRowRed[iCol]=eVal;
	}
	GCD_int<T> eGCD=ComputeGCD_information(eRowRed);
	Treal eNorm=T_NormGen(eGCD.gcd);
	if (IsFirst) {
	  IsFirst=false;
	  AbsVal=eNorm;
	  iRowSearch=iRow;
	} else {
	  if (eNorm < AbsVal) {
	    AbsVal=eNorm;
	    iRowSearch=iRow;
	  }
	}
      }
    return iRowSearch;
  };
  auto SingleMultiplicationUpdate=[&](MyMatrix<T> const& PartMat) -> void {
    Pmat = Pmat*PartMat;
    TheBasisReduced = TheBasisReduced*PartMat;
  };
  auto UpdateMatrices=[&](size_t const& iRank, size_t const& iRowSearch) -> void {
    rowStat[iRowSearch]=0;
    IdxVector.push_back(iRowSearch);
    std::vector<T> eRowRed(nbCol - iRank);
    for (size_t iCol=0; iCol<nbCol - iRank; iCol++)
      eRowRed[iCol]=TheBasisReduced(iRowSearch, iCol+iRank);
    GCD_int<T> eGCD=ComputeGCD_information(eRowRed);
    MyMatrix<T> PartMat=IdentityMat<T>(nbCol);
    for (size_t iCol=iRank; iCol<nbCol; iCol++)
      for (size_t iRow=iRank; iRow<nbCol; iRow++)
	PartMat(iRow, iCol) = eGCD.Pmat(iRow-iRank, iCol-iRank);

    SingleMultiplicationUpdate(PartMat);
    for (size_t iCol=0; iCol<iRank; iCol++) {
      T a=TheBasisReduced(iRowSearch, iCol);
      T b=TheBasisReduced(iRowSearch, iRank);
      T q=QuoInt(a, b);
      MyMatrix<T> PartMatB=IdentityMat<T>(nbCol);
      PartMatB(iRank, iCol)=-q;
      SingleMultiplicationUpdate(PartMatB);
    }
    if (TheBasisReduced(iRowSearch, iRank) < 0) {
      MyMatrix<T> PartMatC=IdentityMat<T>(nbCol);
      PartMatC(iRank, iRank)=-1;
      SingleMultiplicationUpdate(PartMatC);
    }
    TheBasisReord.row(iRank)=TheBasisReduced.row(iRowSearch);
  };
  size_t TheRank=std::min(nbCol, nbRow);
  for (size_t iRank=0; iRank<TheRank; iRank++) {
    size_t iRowSearch=FindMinGCDrow(iRank);
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
template<typename T>
AffineBasisResult Kernel_ComputeAffineBasis(MyMatrix<T> const& EXT)
{
  static_assert(is_ring_field<T>::value, "Requires T to have inverses in Kernel_ComputeAffineBasis");
  size_t nbRow=EXT.rows();
  size_t nbCol=EXT.cols();
  size_t n=nbCol;
  MyMatrix<T> ListExp=ZeroMatrix<T>(nbRow, nbCol);
  MyMatrix<T> EXTwork=EXT;
  std::vector<int> RowStatus(nbRow, 0);
  std::vector<int> ColumnStatus(nbCol,1);
  std::cerr << "Starting Kernel_ComputeAffineBasis\n";
  MyVector<T> V1(nbCol);
  MyVector<T> V2(nbCol);
  MyVector<T> eExpr(nbCol);
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto fInsertValue=[&](int const& idx, int const& iVect) -> bool {
    size_t eCol=miss_val;
    for (size_t iCol=0; iCol<nbCol; iCol++)
      if (eCol == miss_val && EXTwork(iVect, iCol) != 0 && ColumnStatus[iCol] == 1)
	eCol = iCol;
#ifdef DEBUG
    std::cerr << "eCol=" << eCol << "\n";
    if (eCol == miss_val) {
      std::cerr << "This should not be selected\n";
      std::cerr << "nbCol=" << nbCol << "\n";
      for (int iCol=0; iCol<nbCol; iCol++) {
	std::cerr << " iCol=" << iCol << " stat=" << ColumnStatus[iCol] << " val=" << EXTwork(iVect,iCol) << "\n";
      }
      throw TerminalException{1};
    }
#endif
    V1=EXTwork.row(iVect);
    std::cerr << "V1 rows=" << V1.rows() << " cols=" << V1.cols() << "\n";
    ListExp(iVect, idx)=1;
    for (size_t iRow=0; iRow<nbRow; iRow++)
      if (RowStatus[iRow] == 0) {
	V2=EXTwork.row(iRow);
	bool test=IsVectorMultiple(V1, V2);
	if (test) {
	  T eQuot=V2(eCol)/V1(eCol);
	  eExpr=ListExp.row(iRow) - eQuot*ListExp.row(iVect);
	  for (size_t iCol=0; iCol<nbCol; iCol++)
	    if (!IsInteger(eExpr(iCol)))
	      return false;
	}
      }
    ColumnStatus[eCol]=0;
    RowStatus[iVect]=1;
    for (size_t iRow=0; iRow<nbRow; iRow++)
      if (RowStatus[iRow] == 0) {
	V2=EXTwork.row(iRow);
	T eQuot=V2(eCol)/V1(eCol);
	ListExp.row(iRow)=ListExp.row(iRow) - eQuot*ListExp.row(iVect);
	EXTwork.row(iRow)=EXTwork.row(iRow) - eQuot*EXTwork.row(iVect);
	V2=EXTwork.row(iRow);
	bool IsZero=IsZeroVector(V2);
	if (IsZero)
	  RowStatus[iRow]=1;
      }
    size_t nbFinished=0;
    for (size_t iRow=0; iRow<nbRow; iRow++)
      if (RowStatus[iRow] == 1)
	nbFinished++;
    std::cerr << "nbFinished=" << nbFinished << "\n";
    return true;
  };
  int nbIter=1000;
  std::vector<int> ListIdx(n);
  std::vector<int> UsedNumber(nbRow, 0);
  auto GetRandomNumber=[&]() -> int {
    for (int iter=0; iter<nbIter; iter++) {
      int eVal=rand() % nbRow;
      if (UsedNumber[eVal] == 0 && RowStatus[eVal] == 0)
	return eVal;
    }
    return -1;
  };
  auto SetLocallyCorrectIndex=[&](int const& idx) -> int {
    for (int iter=0; iter<nbIter; iter++) {
      int eVal=GetRandomNumber();
      if (eVal == -1)
	return -1;
      bool res=fInsertValue(idx, eVal);
      UsedNumber[eVal]=1;
      if (res) {
	ListIdx[idx]=eVal;
	return 0;
      }
    }
    return -1;
  };
  for (int i=0; i<n; i++) {
    std::cerr << "i=" << i << "\n";
    int eVal=SetLocallyCorrectIndex(i);
    if (eVal == -1)
      return {false, {}};
  }
  return {true, std::move(ListIdx)};
}


template<typename T>
MyMatrix<T> RandomUnimodularMatrix(int const& n)
{
  MyMatrix<T> RetMat = IdentityMat<T>(n);
  for (int iter=0; iter<10; iter++) {
    MyMatrix<T> eMat = IdentityMat<T>(n);
    int idx1 = rand() % n;
    int idx2 = rand() % n;
    if (idx1 != idx2) {
      int pivot = (rand() % 21) - 10;
      eMat(idx1, idx2) = pivot;
    }
    RetMat = eMat * RetMat;
  }
  return RetMat;
}


template<typename T>
AffineBasisResult ComputeAffineBasis(MyMatrix<T> const& EXT)
{
  int nbIter=1000;
  for (int iter=0; iter<nbIter; iter++) {
    AffineBasisResult eAffRes=Kernel_ComputeAffineBasis(EXT);
    if (eAffRes.result)
      return eAffRes;
  }
  return {false, {}};
}


// Compute the translation classes.
// Two classes eV and fV are equivalent if there exists a vector w integer such that
// eV - fV = w M
// that is we need to compute (eV - fV) M^(-1)
template<typename T>
std::vector<MyVector<int>> ComputeTranslationClasses(MyMatrix<T> const& M)
{
  int n=M.rows();
  MyMatrix<T> eInv=Inverse(M);
  std::vector<MyVector<int>> ListClasses;
  std::vector<uint8_t> ListStatus;
  auto IsEquivalent=[&](MyVector<int> const& eV, MyVector<int> const& fV) -> bool {
    MyVector<int> diff=eV - fV;
    for (int i=0; i<n; i++) {
      T eVal=0;
      for (int j=0; j<n; j++)
	eVal += diff(j) * eInv(j,i);
      if (!IsInteger(eVal))
	return false;
    }
    return true;
  };
  auto FuncInsert=[&](MyVector<int> const& eV) -> void {
    for (auto & fV : ListClasses) {
      if (IsEquivalent(eV, fV))
	return;
    }
    ListClasses.push_back(eV);
    ListStatus.push_back(1);
  };
  MyVector<int> zerV = ZeroVector<int>(n);
  FuncInsert(zerV);
  while(true) {
    bool IsFinished=true;
    size_t nbClass=ListClasses.size();
    for (size_t iClass=0; iClass<nbClass; iClass++) {
      if (ListStatus[iClass] == 1) {
	ListStatus[iClass]=0;
	IsFinished=false;
	MyVector<int> eClass=ListClasses[iClass];
	for (int i=0; i<n; i++) {
	  MyVector<int> fClass=eClass;
	  fClass[i]++;
	  FuncInsert(fClass);
	}
      }
    }
    if (IsFinished)
      break;
  }
  return ListClasses;
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
template<typename T>
MyMatrix<T> GetZbasis(MyMatrix<T> const& ListElement)
{
  static_assert(is_euclidean_domain<T>::value, "Requires T to be an Euclidean domain in GetZbasis");
  using Treal=typename underlying_totally_ordered_ring<T>::real_type;
  int TheDim=ListElement.cols();
  MyMatrix<T> ListEqua;
  MyMatrix<T> InvMatrix;
  MyMatrix<T> InvMatrixTr;
  MyMatrix<T> TheBasis;
  std::vector<int> eSet;
  auto fGetOneBasis=[&](MyVector<T> const& eSol) -> MyMatrix<T> {
    int DimLoc=TheBasis.rows();
    //    std::cerr << "DimLoc=" << DimLoc << " |eSol|=" << eSol.size() << "\n";
    MyMatrix<T> TheRedMat=ZeroMatrix<T>(DimLoc+1, DimLoc);
    for (int i=0; i<DimLoc; i++)
      TheRedMat(i,i)=1;
    for (int i=0; i<DimLoc; i++)
      TheRedMat(DimLoc,i)=eSol(i);
    //    std::cerr << "After TheRedMat construction\n";
    MyMatrix<T> NSP=NullspaceIntMat(TheRedMat);
    //    std::cerr << "We have NSP\n";
#ifdef DEBUG
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
    MyVector<T> eVect_pre=GetMatrixRow(NSP,0);
    MyVector<T> eVect = CanonicalizeVectorToInvertible(eVect_pre);
    /*    std::cerr << "eVect=";
	  WriteVector(std::cerr, eVect);*/
    int n2=DimLoc+1;
    while(true) {
      std::vector<int> ListIdxZ;
      std::vector<int> ListIdxNZ;
      for (int i=0; i<n2; i++) {
	if (eVect(i) == 0) {
	  ListIdxZ.push_back(i);
	} else {
	  ListIdxNZ.push_back(i);
	}
      }
      if (ListIdxNZ.size() == 1)
	return SelectRow(TheRedMat, ListIdxZ);
      std::vector<int> AbsList;
      bool IsFirst=true;
      int ThePivot=-1;
      Treal TheMin=-1;
      for (auto & eVal : ListIdxNZ) {
	T eVal_T=eVect(eVal);
	Treal eAbs=T_NormGen(eVal_T);
	if (IsFirst) {
	  TheMin=eAbs;
	  ThePivot=eVal;
	} else {
	  if (eAbs < TheMin) {
	    TheMin=eAbs;
	    ThePivot=eVal;
	  }
	}
	IsFirst=false;
      }
      for (int iCol=0; iCol<n2; iCol++)
	if (iCol != ThePivot) {
	  T TheQ=QuoInt(eVect(iCol), eVect(ThePivot));
          if (TheQ != 0) {
            TheRedMat.row(ThePivot) += TheQ*TheRedMat.row(iCol);
            eVect(iCol) -= TheQ*eVect(ThePivot);
          }
	}
    }
  };
  auto fComputeSpeed=[&]() -> void {
    int dimSpace=TheBasis.rows();
    if (dimSpace == 0) {
      ListEqua=IdentityMat<T>(TheDim);
      eSet={};
    } else {
      /*      std::cerr << "TheBasis=\n";
	      WriteMatrix(std::cerr, TheBasis);*/
      ListEqua=NullspaceTrMat(TheBasis);
      eSet=ColumnReductionSet(TheBasis);
      //      std::cerr << "|eSet|=" << eSet.size() << "\n";
      InvMatrix=Inverse(SelectColumn(TheBasis, eSet));
      InvMatrixTr=InvMatrix.transpose();
    }
  };
  auto IsInSpace=[&](MyVector<T> const& eElt) -> bool {
    int nbEqua=ListEqua.rows();
    //    std::cerr << "nbEqua=" << nbEqua << "\n";
    for (int iEqua=0; iEqua<nbEqua; iEqua++) {
      T eSum=0;
      for (int i=0; i<TheDim; i++)
	eSum += ListEqua(iEqua,i)*eElt(i);
      if (eSum != 0)
	return false;
    }
    return true;
  };
  auto fInsert=[&](MyVector<T> const& eElt) -> void {
    bool test=IsInSpace(eElt);
    //    std::cerr << "test=" << test << "\n";
    if (!test) {
      TheBasis=ConcatenateMatVec(TheBasis, eElt);
      fComputeSpeed();
    } else {
      if (TheBasis.rows() == 0)
	return;
      MyVector<T> eEltRed=SelectColumnVector(eElt, eSet);
      MyVector<T> eSol=InvMatrixTr*eEltRed;
      if (IsIntegralVector(eSol))
	return;
      MyMatrix<T> NewBasis=fGetOneBasis(eSol);
      TheBasis=NewBasis*TheBasis;
      fComputeSpeed();
    }
  };
  fComputeSpeed();
  int nbElt=ListElement.rows();
  //  std::cerr << "nbElt=" << nbElt << "\n";
  for (int iElt=0; iElt<nbElt; iElt++) {
    //    std::cerr << "iElt=" << iElt << "\n";
    MyVector<T> eElt=GetMatrixRow(ListElement, iElt);
    fInsert(eElt);
    //    std::cerr << "After fInsert\n";
  }

#ifdef DEBUG
  int DimSpace=TheBasis.rows();
  for (int iBas=0; iBas<DimSpace; iBas++) {
    MyVector<T> eLine=GetMatrixRow(TheBasis, iBas);
    //      std::cerr << "Before SolutionIntMat, iBas=" << iBas << "\n";
    ResultSolutionIntMat<T> eResIntMat=SolutionIntMat(ListElement, eLine);
    /*      std::cerr << "ListElement=\n";
	    WriteMatrixGAP(std::cerr, ListElement);
	    std::cerr << "eLine=\n";
	    WriteVectorGAP(std::cerr, eLine);
	    std::cerr << "After SolutionIntMat 1\n";*/
    if (!eResIntMat.TheRes) {
      std::cerr << "Error in GetZbasis 1\n";
      throw TerminalException{1};
    }
  }
  for (int iElt=0; iElt<nbElt; iElt++) {
    MyVector<T> eElt=GetMatrixRow(ListElement, iElt);
    //      std::cerr << "Before SolutionIntMat, iElt=" << iElt << "\n";
    ResultSolutionIntMat<T> eResIntMat=SolutionIntMat(TheBasis, eElt);
    /*      std::cerr << "TheBasis=\n";
	    WriteMatrixGAP(std::cerr, TheBasis);
	    std::cerr << "eElt=\n";
	    WriteVectorGAP(std::cerr, eElt);
	    std::cerr << "After SolutionIntMat 2 eResIntMat.TheRes=" << eResIntMat.TheRes << "\n";*/
    if (!eResIntMat.TheRes) {
      std::cerr << "Error in GetZbasis 2\n";
      throw TerminalException{1};
    }
  }
#endif
  return TheBasis;
}



template<typename Tint>
MyMatrix<Tint> SYMPL_ComputeSymplecticBasis(MyMatrix<Tint> const& M)
{
  int nb_row = M.rows();
  int n = M.cols() / 2;
  MyMatrix<Tint> Mwork = M;
  MyMatrix<Tint> SympFormMat = ZeroMatrix<Tint>(2*n, 2*n);
  for (int i=0; i<n; i++) {
    SympFormMat(i,n+i) = 1;
    SympFormMat(n+i,i) = -1;
  }
  auto GetInitialVector=[&]() -> MyMatrix<Tint> {
    int pos=0;
    std::vector<Tint> ListX(2*n);
    while(true) {
      MyVector<Tint> V = GetMatrixRow(Mwork, pos);
      if (!IsZeroVector(V))
        return CanonicalizeVector(V);
      pos++;
      if (pos == 2*n)
        break;
    }
    std::cerr << "Failed to find non-zero vector\n";
    throw TerminalException{1};
  };
  auto GetPairVector=[&](MyVector<Tint> const& w1) -> MyVector<Tint> {
    std::vector<Tint> ListScal(nb_row);
    for (int i_row=0; i_row<nb_row; i_row++) {
      Tint eScal=0;
      for (int i=0; i<2*n; i++)
        eScal += w1(i) * Mwork(i_row, i);
      ListScal[i_row] = eScal;
    }
    GCD_int<Tint> eGCD = ComputeGCD_information(ListScal);
    if (T_abs(eGCD.gcd) != 1) {
      std::cerr << "The gcd should be equal to 1\n";
      throw TerminalException{1};
    }
    MyVector<Tint> SumVect = ZeroVector<Tint>(2*n);
    for (int i_row=0; i_row<nb_row; i_row++)
      SumVect += GetMatrixRow(Mwork, i_row) * eGCD.Pmat(i_row, 0);
    return SumVect;
  };
  MyMatrix<Tint> CompleteBasis(2*n, 2*n);
  for (int i=0; i<n; i++) {
    MyVector<Tint> w1 = GetInitialVector();
    MyVector<Tint> wN = GetPairVector(w1);
    CompleteBasis.row(i) = w1;
    CompleteBasis.row(n + i) = wN;
    MyVector<Tint> J_w1 = SympFormMat * w1;
    MyVector<Tint> J_wN = SympFormMat * wN;
    for (int i_row=0; i_row<nb_row; i_row++) {
      MyVector<Tint> eRow = Mwork.row(i_row);
      Tint scal1 = eRow.dot(J_w1);
      Tint scalN = eRow.dot(J_wN);
      MyVector<Tint> NewRow = eRow - scalN * w1 + scal1 * wN;
      Mwork.row(i_row) = NewRow;
    }
  }
  return CompleteBasis;
}


template<typename T>
MyMatrix<T> CanonicalizeOrderedMatrix_Kernel(const MyMatrix<T>& M)
{
  static_assert(is_ring_field<T>::value, "Requires T to have inverses in Kernel_ComputeAffineBasis");
  MyMatrix<T> Basis = RowReduction(M);
  MyMatrix<T> M1 = M * Inverse(Basis);
  return RemoveFractionMatrix(M1);
}

template<typename T>
inline typename std::enable_if<is_ring_field<T>::value,MyMatrix<T>>::type CanonicalizeOrderedMatrix(MyMatrix<T> const& Input)
{
  return CanonicalizeOrderedMatrix_Kernel(Input);
}

template<typename T>
inline typename std::enable_if<(not is_ring_field<T>::value),MyMatrix<T>>::type CanonicalizeOrderedMatrix(MyMatrix<T> const& Input)
{
  using Tfield=typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield,T>(Input);
  MyMatrix<Tfield> OutputF = CanonicalizeOrderedMatrix_Kernel(InputF);
  return UniversalMatrixConversion<T,Tfield>(OutputF);
}




#endif
