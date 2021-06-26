#ifndef COMBINATORICS_FUNCTIONALITIES_INCLUDE
#define COMBINATORICS_FUNCTIONALITIES_INCLUDE

#include "COMB_Combinatorics_elem.h"
#include "COMB_Vectors.h"
#include "MAT_Matrix.h"
#include "Boost_bitset.h"
#include "COMB_Stor.h"





struct BlockIteration {
public:
  // no copy
  BlockIteration(const BlockIteration&) = delete;

  // no assign
  BlockIteration& operator=(const BlockIteration&) = delete;

  // no move
  BlockIteration(BlockIteration&&) = delete;

  // no default constructor
  BlockIteration() = delete;

  BlockIteration(int const& eDim, int const& eSize) : dim(eDim), size(eSize), eVect(dim,0)
  {
  }
  int IncrementShow()
  {
    for (int i=0; i<dim; i++)
      if (eVect[i] < size-1) {
	eVect[i]++;
	for (int j=0; j<i; j++)
	  eVect[j]=0;
	return i;
      }
    return -1;
  }
  void IncrementSilent()
  {
    for (int i=0; i<dim; i++)
      if (eVect[i] < size-1) {
	eVect[i]++;
	for (int j=0; j<i; j++)
	  eVect[j]=0;
	return;
      }
    std::cerr << "Should not reach that stage\n";
    throw TerminalException{1};
  }
  std::vector<int> GetVect() const
  {
    return eVect;
  }
  int GetNbPoss() const
  {
    int eRet=1;
    for (int i=0; i<dim; i++)
      eRet *= size;
    return eRet;
  }
private:
  int dim;
  int size;
  std::vector<int> eVect;
};



struct BlockCppIterator {
private:
  int dim;
  int size;
  struct IteratorContain {
  private:
    int dim_iter;
    int size_iter;
    std::vector<int> U;
    void single_increase()
    {
      for (int i=0; i<dim_iter; i++)
        if (U[i] < size_iter-1) {
          U[i]++;
          for (int j=0; j<i; j++)
            U[j]=0;
          return;
        }
      U = {};
    }
  public:
    IteratorContain(int const& eDim, int const& eSize, std::vector<int> const& eU) : dim_iter(eDim), size_iter(eSize), U(eU)
    {
    }
    std::vector<int> const& operator*()
    {
      return U;
    }
    IteratorContain & operator++()
    {
      single_increase();
      return *this;
    }
    IteratorContain operator++(int)
    {
      IteratorContain tmp = *this;
      single_increase();
      return tmp;
    }
    bool operator!=(IteratorContain const& iter)
    {
      if (iter.dim_iter != dim_iter)
        return true;
      if (iter.size_iter != size_iter)
        return true;
      if (iter.U.size() != U.size())
        return true;
      for (size_t i=0; i<U.size(); i++)
        if (iter.U[i] != U[1])
          return true;
      return false;
    }
    bool operator==(IteratorContain const& iter)
    {
      if (iter.dim_iter != dim_iter)
        return false;
      if (iter.size_iter != size_iter)
        return false;
      if (iter.U.size() != U.size())
        return false;
      for (size_t i=0; i<U.size(); i++)
        if (iter.U[i] != U[1])
          return false;
      return true;
    }
  };
public:
  // no copy
  BlockCppIterator(const BlockCppIterator&) = delete;

  // no assign
  BlockCppIterator& operator=(const BlockCppIterator&) = delete;

  // no move
  BlockCppIterator(BlockCppIterator&&) = delete;

  // no default constructor
  BlockCppIterator() = delete;

  BlockCppIterator(int const& eDim, int const& eSize) : dim(eDim), size(eSize)
  {
  }
  // The iterator business
  using iterator=IteratorContain;
  using const_iterator=IteratorContain;
  const_iterator cbegin() const
  {
    return { dim, size, std::vector<int>(dim,0) };
  }
  const_iterator cend() const
  {
    return {dim, size, {}};
  }
  const_iterator begin() const
  {
    return { dim, size, std::vector<int>(dim,0) };
  }
  const_iterator end() const
  {
    return {dim, size, {}};
  }
};






struct BlockIterationMultiple {
public:
  // no copy
  BlockIterationMultiple(const BlockIterationMultiple&) = delete;

  // no assign
  BlockIterationMultiple& operator=(const BlockIterationMultiple&) = delete;

  // no move
  BlockIterationMultiple(BlockIterationMultiple&&) = delete;

  // no default constructor
  BlockIterationMultiple() = delete;

  BlockIterationMultiple(std::vector<int> const& VectSiz) : dim(VectSiz.size()), VectSiz(VectSiz)
  {
  }
  //
  // hand written loop functionality
  //
  size_t size() const
  {
    size_t eRet=1;
    for (size_t i=0; i<dim; i++)
      eRet *= VectSiz[i];
    return eRet;
  }
  //
  // The iterator business
  //
  struct IteratorContain {
  private:
    size_t dim_iter;
    std::vector<int> VectSiz;
    std::vector<int> U;
    void single_increase()
    {
      for (size_t i=0; i<dim_iter; i++)
        if (U[i] < VectSiz[i] - 1) {
          U[i]++;
          for (size_t j=0; j<i; j++)
            U[j]=0;
          return;
        }
      U = {};
    }
  public:
    IteratorContain(std::vector<int> const& VectSiz, std::vector<int> const& U) : dim_iter(VectSiz.size()), VectSiz(VectSiz), U(U)
    {
    }
    std::vector<int> const& operator*()
    {
      return U;
    }
    IteratorContain & operator++()
    {
      single_increase();
      return *this;
    }
    IteratorContain operator++(int)
    {
      IteratorContain tmp = *this;
      single_increase();
      return tmp;
    }
    bool operator!=(IteratorContain const& iter)
    {
      if (dim_iter != iter.dim_iter)
        return true;
      for (size_t i=0; i<dim_iter; i++)
        if (VectSiz[i] != iter.VectSiz[i])
          return true;
      if (U.size() != iter.U.size())
        return true;
      for (size_t i=0; i<U.size(); i++)
        if (U[i] != iter.U[i])
          return true;
      return false;
    }
    bool operator==(IteratorContain const& iter)
    {
      if (dim_iter != iter.dim_iter)
        return false;
      for (size_t i=0; i<dim_iter; i++)
        if (VectSiz[i] != iter.VectSiz[i])
          return false;
      if (U.size() != iter.U.size())
        return false;
      for (size_t i=0; i<dim_iter; i++)
        if (U[i] != iter.U[i])
          return false;
      return true;
    }
  };
public:


  // The iterator business
  using iterator=IteratorContain;
  using const_iterator=IteratorContain;
  const_iterator cbegin() const
  {
    return IteratorContain(VectSiz, std::vector<int>(dim,0) );
  }
  const_iterator cend() const
  {
    return IteratorContain(VectSiz, {});
  }
  const_iterator begin() const
  {
    return IteratorContain(VectSiz, std::vector<int>(dim,0));
  }
  const_iterator end() const
  {
    return IteratorContain(VectSiz, {});
  }
private:
  size_t dim;
  std::vector<int> VectSiz;
};



template<typename Tint>
Tint ConvertVectorToNumber(std::vector<int> const& V, int const& N)
{
  Tint val=0;
  Tint expo=1;
  int len=V.size();
  for (int i=0; i<len; i++) {
    val += expo*V[i];
    expo *= N;
  }
  return val;
}


std::vector<int> RandomPermutation(int const& n)
{
  std::vector<int> RetList(n,-1);
  for (int i=0; i<n; i++) {
    int rnd_pos = rand() % (n - i);
    int idx=0;
    bool IsAssigned = false;
    for (int j=0; j<n; j++) {
      if (RetList[j] == -1 && !IsAssigned) {
        if (rnd_pos == idx) {
          RetList[j] = i;
          IsAssigned = true;
        }
        idx++;
      }
    }
  }
  return RetList;
}




template<typename Tint>
std::vector<int> ConvertNumberToVector(Tint const& val, int const& N, int const& len)
{
  std::vector<int> V(len);
  Tint valWork=val;
  for (int i=0; i<len; i++) {
    int res=valWork % N;
    V[i]=res;
    valWork = (valWork - res)/N;
  }
  return V;
}

template<typename Tint>
std::vector<int> ConvertNumberToPrimitiveVector(Tint const& val, int const& N, int const& len)
{
  Tint valWork=val;
  std::vector<int> V(len,0);
  Tint expo=1;
  for (int i=0; i<len; i++) {
    if (valWork < expo) {
      std::vector<int> V2=ConvertNumberToVector(valWork, N, i);
      for (int j=0; j<i; j++)
	V[j]=V2[j];
      V[i]=1;
      return V;
    }
    else {
      valWork -= expo;
      expo *= N;
    }
  }
  std::cerr << "val=" << val << " N=" << N << " len=" << len << "\n";
  std::cerr << "Failed to find an index in ConvertNumberToPrimitiveVector\n";
  throw TerminalException{1};
}

/* Primitive vector is a vector of the form
   (a1, ...., aj, 1, 0, ...., 0) */
template<typename Tint>
Tint ConvertPrimitiveVectorToNumber(std::vector<int> const& V, int const& N)
{
  int len=V.size();
  int ifound=-1;
  for (int i=0; i<len; i++) {
    if (V[i] != 0)
      ifound=i;
  }
  if (V[ifound] != 1) {
    std::cerr << "V=";
    for (int i=0; i<len; i++)
      std::cerr << " " << V[i];
    std::cerr << "\n";
    std::cerr << "ConvertPrimitiveVectorToNumber Inconsistent input\n";
    throw TerminalException{1};
  }
  Tint val=0;
  Tint expo=1;
  std::vector<int> Vred(ifound);
  for (int i=0; i<ifound; i++) {
    val += expo;
    expo *= N;
    Vred[i]=V[i];
  }
  return val + ConvertVectorToNumber<Tint>(Vred, N);
}










struct BlockIteration01 {
public:
  // no copy
  BlockIteration01(const BlockIteration01&) = delete;

  // no assign
  BlockIteration01& operator=(const BlockIteration01&) = delete;

  // no move
  BlockIteration01(BlockIteration01&&) = delete;

  BlockIteration01() = delete;
  
  BlockIteration01(int const& eDim) : dim(eDim), eFace(dim)
  {
  }
  int IncrementShow()
  {
    for (int i=0; i<dim; i++)
      if (eFace[i] == 0) {
	eFace[i]=1;
	for (int j=0; j<i; j++)
	  eFace[j]=0;
	return i;
      }
    return -1;
  }
  void IncrementSilent()
  {
    for (int i=0; i<dim; i++)
      if (eFace[i] == 0) {
	eFace[i]=1;
	for (int j=0; j<i; j++)
	  eFace[j]=0;
	return;
      }
    std::cerr << "Should not reach that stage\n";
    throw TerminalException{1};
  }
  Face GetFace() const
  {
    return eFace;
  }
  int GetNbPoss() const
  {
    int eRet=1;
    for (int i=0; i<dim; i++)
      eRet *= 2;
    return eRet;
  }
private:
  int dim;
  Face eFace;
};









template<typename Tint, typename TstorAll, typename TstorActive>
std::vector<Tint> ListRepresentativeConnectedComponent(Tint const& len, std::function<std::vector<Tint>(Tint const&)> const& GetAdjacent)
{
  TstorAll Vlist(len);
  for (Tint i=0; i<len; i++)
    Vlist.set(i, true);
  std::vector<Tint> ListOrbRepresentative;
  auto ComputeFullOrbit=[&](Tint const& ePoint) -> void {
    TstorActive Vactive(len);
    Vactive.set(ePoint,true);
    while(true) {
      if (Vactive.empty())
	break;
      Tint TheFirst=Vactive.find_first();
      Vactive.set(TheFirst,false);
      Vlist.set(TheFirst,false);
      std::vector<Tint> ListAdj=GetAdjacent(TheFirst);
      for (auto & eAdj : ListAdj)
	if (!Vactive.get(eAdj) && Vlist.get(eAdj))
	  Vactive.set(eAdj,true);
    }
  };
  while(true) {
    if (Vlist.empty())
      break;
    Tint TheFirst=Vlist.find_first();
    ListOrbRepresentative.push_back(TheFirst);
    ComputeFullOrbit(TheFirst);
  }
  return ListOrbRepresentative;
}



// Relation of binomial coefficients is
// C_n^k = (n/k) C_{n-1}^{k-1}
//
template<typename T, typename Tint>
T BinomialCoefficient(Tint const& n, Tint const& k)
{
  T eRet=1;
  for (Tint i=1; i<=k; i++) {
    T eNum = n - (k - i);
    T eDen = i;
    eRet *= eNum;
    eRet /= eDen;
  }
  return eRet;
}



// functionality for directly accessing to binomial coefficients
// by their number.
// Very little memory overhead.
//
// Used formula is
// C_n^k = C_{n-1}^k + C_{n-1}^{k-1}
// 
// So, in the structure we have 
// (1, ... ,1, 0, ... , 0) which is vector of position 0.
template<typename T>
struct IteratorBinomial {
public:
  IteratorBinomial() = delete;

  IteratorBinomial(int const& eN, int const& eK) : n(eN), k(eK)
  {
    ListBinom.setZero(n+1,n+1);
    ListBinom(0,0)=1;
    for (int eS=1; eS<=n; eS++)
      for (int fK=0; fK<=eS; fK++)
	ListBinom(eS, fK)=BinomialCoefficient<T,int>(eS, fK);
  }
  T size() const
  {
    return ListBinom(n,k);
  }
  Face first_face() const
  {
    Face eFace(n);
    for (int i=0; i<k; i++)
      eFace[i]=1;
    return eFace;
  }
  std::vector<int> first_stdvect() const
  {
    std::vector<int> eVect(k);
    for (int i=0; i<k; i++)
      eVect[i]=i;
    return eVect;
  }
  std::vector<Face> ListAllFace() const
  {
    int eBin=BinomialCoefficient<int,int>(n,k);
    std::vector<Face> ListFace(eBin);
    Face eFace=(*this).first_face();
    int idx=0;
    for (int i=0; i<eBin; i++) {
      ListFace[idx]=eFace;
      idx++;
      if (i<eBin-1)
	FaceIncrement(eFace);
    }
    return ListFace;
  }
  T GetIndexFromFace(Face const& eFace) const
  {
    T pos=0;
    int sizRemain=k;
    for (int i=0; i<n; i++) {
      int j=n-1-i;
      if (eFace[j] == 1) {
	T eBin=ListBinom[n-1-i][sizRemain];
	pos += eBin;
	sizRemain--;
      }
    }
    return pos;
  }
  Face GetFaceFromIndex(T const& pos) const
  {
    T posWork=pos;
    Face eFace(n);
    int sizRemain=k;
    for (int i=0; i<n; i++) {
      int j=n-1-i;
      T eBin=ListBinom(j,sizRemain);
      if (posWork >= eBin) {
	eFace[j]=1;
	posWork -= eBin;
	sizRemain--;
      }
    }
    return eFace;
  }
  bool StdvectIncrement(std::vector<int> & Tvect) const
  {
    if (k == 0) // Set is empty, no increment possible
      return false;
    Tvect[0]++;
    int xy2=1;
    while ((xy2 < k) && (Tvect[xy2-1] >= Tvect[xy2])) {
      Tvect[xy2]++;
      xy2++;
    }
    if (xy2 != 1) {
      for (int xy1=0; xy1<xy2-1; xy1++)
	Tvect[xy1]=xy1;
    }
    if (Tvect[k-1] == n)
      return false;
    return true;
  }
  bool IsCorrectStdVect(std::vector<int> const& eVect) const
  {
    for (int i=0; i<k; i++) {
      if (eVect[i] < 0)
	return false;
      if (eVect[i] >= n)
	return false;
    }
    for (int i=1; i<k; i++) {
      if (eVect[i-1] >= eVect[i])
	return false;
    }
    return true;
  }
  bool FaceIncrement(Face & eFace) const
  {
    std::vector<int> Tvect(k);
    int idx=0;
    for (int i=0; i<n; i++)
      if (eFace[i] == 1) {
	Tvect[idx]=i;
	idx++;
	eFace[i]=0;
      }
    bool res=StdvectIncrement(Tvect);
    if (!res)
      return false;
    for (int i=0; i<k; i++)
      eFace[Tvect[i]]=1;
    return true;
  }
private:
  int n;
  int k;
  MyMatrix<T> ListBinom;
};


template<typename T>
std::vector<T> MinimumDihedralOrbit(std::vector<T> const& eL)
{
  int len=eL.size();
  T eMin=eL[0];
  int posmin=0;
  for (int pos=1; pos<len; pos++)
    if (eL[pos] < eMin) {
      posmin=pos;
      eMin = eL[pos];
    }
  std::vector<T> V1(len);
  std::vector<T> V2(len);
  int pos1=posmin;
  int pos2=posmin;
  for (int u=0; u<len; u++) {
    V1[u] = eL[pos1];
    V2[u] = eL[pos2];
    pos1 = NextIdx(len, pos1);
    pos2 = PrevIdx(len, pos2);
  }
  for (int u=0; u<len; u++) {
    if (V1[u] < V2[u])
      return V1;
    if (V2[u] < V1[u])
      return V2;
  }
  return V1;
};

template<typename T1, typename T2>
void SortParallel(std::vector<T1> & list1, std::vector<T2> & list2)
{
  struct PairData {
    T1 val1;
    T2 val2;
  };
  int len=list1.size();
  std::vector<PairData> ListPair(len);
  for (int i=0; i<len; i++)
    ListPair[i] = {list1[i], list2[i]};
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const & a, PairData const& b) -> bool {
	 return a.val1 < b.val1;
       });
  for (int i=0; i<len; i++) {
    list1[i] = ListPair[i].val1;
    list2[i] = ListPair[i].val2;
  }
}



#endif
