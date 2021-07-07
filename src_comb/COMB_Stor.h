#ifndef POLYTOPE_STOR
#define POLYTOPE_STOR

#include "Temp_common.h"
#include "Boost_bitset.h"

//
// This is a set of functionality for storing bits
// with different complexity theories and memory requirements.
// Those classes are supposed to be used as template arguments.
//
// Each class implements following subset of dynamic_subset:
// a[i] = val with val a boolean, i.e. 0/1 or true/false.
// val  = a[i]
// find_first()  : first true element of the list
// empty()       : true if empty, false otherwise
// count()       : number of true entries
//
// Possible available classes:
// dynamic_bitset (aliases to Face)
//  ---find_first  : O(n)
//  ---data access : O(1)
//  ---count       : O(n)
// Memory expenses : 2^n + overhead
//
// DoubleList:
//  ---find_first  : O(1)
//  ---data access : O(1)
//  ---count       : O(n)
// Memory expenses : 2n * 2^n + overhead
//

struct IntegerSubsetStorage {
  size_t MaxElement;
  std::vector<size_t> ListNext;
  std::vector<size_t> ListPrev;
};


void VSLT_ZeroAssignment(IntegerSubsetStorage & VSLT)
{
  size_t MaxElement=VSLT.MaxElement;
  size_t Maxp2=MaxElement+2;
  size_t miss_val = std::numeric_limits<size_t>::max();
  for (size_t iVert=0; iVert<Maxp2; iVert++) {
    VSLT.ListNext[iVert] = miss_val;
    VSLT.ListPrev[iVert] = miss_val;
  }
  VSLT.ListNext[MaxElement  ] = MaxElement+1;
  VSLT.ListPrev[MaxElement+1] = MaxElement;
}


IntegerSubsetStorage VSLT_InitializeStorage(size_t const& MaxElement)
{
  IntegerSubsetStorage VSLT;
  size_t Maxp2 = MaxElement+2;
  std::vector<size_t> ListNext(Maxp2);
  std::vector<size_t> ListPrev(Maxp2);
  VSLT.ListNext = ListNext;
  VSLT.ListPrev = ListPrev;
  VSLT.MaxElement = MaxElement;
  VSLT_ZeroAssignment(VSLT);
  return VSLT;
}


int VSLT_NrElement(IntegerSubsetStorage const& VSLT)
{
  size_t MaxElt=VSLT.MaxElement+1;
  size_t pos=MaxElt-1;
  size_t NbElt=0;
  while(true) {
    size_t posNext=VSLT.ListNext[pos];
    if (posNext == MaxElt)
      return NbElt;
    pos = posNext;
    NbElt++;
  }
}

int VSLT_TheFirstPosition(IntegerSubsetStorage const& VSLT)
{
  return VSLT.ListNext[VSLT.MaxElement];
}

bool VSLT_IsItInSubset(IntegerSubsetStorage const& VSLT, size_t const& pos)
{
  size_t miss_val = std::numeric_limits<size_t>::max();
  return VSLT.ListNext[pos] != miss_val;
}

void VSLT_StoreValue(IntegerSubsetStorage & VSLT, size_t const& pos)
{
  size_t posAfter=VSLT.ListNext[VSLT.MaxElement];
  VSLT.ListNext[VSLT.MaxElement] = pos;
  VSLT.ListNext[pos] = posAfter;
  VSLT.ListPrev[posAfter] = pos;
  VSLT.ListPrev[pos] = VSLT.MaxElement;
}



void VSLT_RemoveValue(IntegerSubsetStorage & VSLT, size_t const& pos)
{
  size_t miss_val = std::numeric_limits<size_t>::max();
  size_t posNext=VSLT.ListNext[pos];
  size_t posPrev=VSLT.ListPrev[pos];
  VSLT.ListNext[posPrev] = posNext;
  VSLT.ListPrev[posNext] = posPrev;
  VSLT.ListNext[pos] = miss_val;
  VSLT.ListPrev[pos] = miss_val;
}


bool VSLT_IsEmpty(IntegerSubsetStorage const& VSLT)
{
  return VSLT.ListNext[VSLT.MaxElement] == VSLT.MaxElement+1;
}


// This is an artificial class for allowing the operation
// a[i]=a
template<typename Tint>
struct DoubleList {
public:
  DoubleList(const DoubleList&) = delete;
  DoubleList& operator=(const DoubleList&) = delete;
  DoubleList(DoubleList&&) = delete;
  DoubleList() = delete;
  DoubleList(Tint const& eMax) : MaxElement(eMax)
  {
    Tint Maxp2=eMax + 2;
    ListNext = std::vector<Tint>(Maxp2, Maxp2);
    ListPrev = std::vector<Tint>(Maxp2, Maxp2);
    ListNext[MaxElement]=MaxElement+1;
    ListPrev[MaxElement+1]=MaxElement;
  }
  ~DoubleList()
  {
  }
  Tint count() const
  {
    Tint Maxp1=MaxElement+1;
    Tint pos=MaxElement;
    Tint NbElt=0;
    while(true) {
      pos=ListNext[pos];
      if (pos == Maxp1)
	return NbElt;
      NbElt++;
    }
    std::cerr << "We should not reach that stage\n";
    throw TerminalException{1};
  }
  bool empty() const
  {
    Tint Maxp1=MaxElement+1;
    //    std::cerr << "MaxElement=" << MaxElement << " Maxp1=" << Maxp1 << " ListNext[]=" << ListNext[MaxElement] << "\n";
    if (ListNext[MaxElement] == Maxp1)
      return true;
    return false;
  }
  Tint find_first() const
  {
    return ListNext[MaxElement];
  }
  void set(Tint const& pos, bool const& eVal)
  {
    Tint Maxp2=MaxElement+2;
    if (eVal) {
      if (ListNext[pos] == Maxp2) {
	Tint posAfter=ListNext[MaxElement];
	ListNext[MaxElement]=pos;
	ListNext[pos]=posAfter;
	ListPrev[posAfter]=pos;
	ListPrev[pos]=MaxElement;
      }
    } else {
      if (ListNext[pos] != Maxp2) {
	Tint posNext=ListNext[pos];
	Tint posPrev=ListPrev[pos];
	ListNext[posPrev]=posNext;
	ListPrev[posNext]=posPrev;
	ListNext[pos]=Maxp2;
	ListPrev[pos]=Maxp2;
      }
    }
  }
  bool get(Tint const& pos) const
  {
    Tint Maxp2=MaxElement+2;
    if (ListNext[pos] == Maxp2)
      return false;
    return true;
  }
private:
  Tint MaxElement;
  std::vector<Tint> ListNext;
  std::vector<Tint> ListPrev;
};


struct EncapsDynamicBitset {
public:
  EncapsDynamicBitset(const EncapsDynamicBitset&) = delete;
  EncapsDynamicBitset& operator=(const EncapsDynamicBitset&) = delete;
  EncapsDynamicBitset(EncapsDynamicBitset&&) = delete;
  EncapsDynamicBitset() = delete;
  EncapsDynamicBitset(std::size_t const& eMax) : eVect(eMax)
  {
    TotalCount=0;
  }
  ~EncapsDynamicBitset()
  {
  }
  std::size_t count() const
  {
    return TotalCount;
  }
  std::size_t find_first() const
  {
    return eVect.find_first();
  }
  void set(std::size_t const& pos, bool const& eVal)
  {
    if (eVect[pos] && !eVal)
      TotalCount--;
    if (!eVect[pos] && eVal)
      TotalCount++;
    eVect[pos]=eVal;
  }
  bool get(std::size_t const& pos) const
  {
    return eVect[pos];
  }
  bool empty() const
  {
    if (TotalCount == 0)
      return true;
    return false;
  }
private:
  std::size_t TotalCount;
  boost::dynamic_bitset<> eVect;
};



template<typename T>
struct StackStorage {
public:
  StackStorage()
  {
    TheSize=0;
  }
  void push_back(T const& val)
  {
    if (TheSize == ListElt.size()) {
      TheSize++;
      ListElt.push_back(val);
    }
    else {
      ListElt[TheSize]=val;
      TheSize++;
    }
  }
  T pop()
  {
    T val=ListElt[TheSize-1];
    TheSize--;
    return val;
  }
  size_t size()
  {
    return TheSize;
  }
private:
  std::vector<T> ListElt;
  size_t TheSize;
};


#endif
