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
  int MaxElement;
  std::vector<int> ListNext;
  std::vector<int> ListPrev;
};


void VSLT_ZeroAssignment(IntegerSubsetStorage & VSLT)
{
  int MaxElement=VSLT.MaxElement;
  int Maxp2=MaxElement+2;
  for (int iVert=0; iVert<Maxp2; iVert++) {
    VSLT.ListNext[iVert]=-1;
    VSLT.ListPrev[iVert]=-1;
  }
  VSLT.ListNext[MaxElement  ] = MaxElement+1;
  VSLT.ListPrev[MaxElement+1] = MaxElement;
}


IntegerSubsetStorage VSLT_InitializeStorage(int const& MaxElement)
{
  IntegerSubsetStorage VSLT;
  int Maxp2 = MaxElement+2;
  std::vector<int> ListNext(Maxp2);
  std::vector<int> ListPrev(Maxp2);
  VSLT.ListNext = ListNext;
  VSLT.ListPrev = ListPrev;
  VSLT.MaxElement = MaxElement;
  VSLT_ZeroAssignment(VSLT);
  return VSLT;
}


int VSLT_NrElement(IntegerSubsetStorage const& VSLT)
{
  int MaxElt=VSLT.MaxElement+1;
  int pos=MaxElt-1;
  int NbElt=0;
  while(true) {
    int posNext=VSLT.ListNext[pos];
    if (posNext == MaxElt)
      return NbElt;
    pos=posNext;
    NbElt++;
  }
}

int VSLT_TheFirstPosition(IntegerSubsetStorage const& VSLT)
{
  return VSLT.ListNext[VSLT.MaxElement];
}

bool VSLT_IsItInSubset(IntegerSubsetStorage const& VSLT, int const& pos)
{
  return VSLT.ListNext[pos] != -1;
}

void VSLT_StoreValue(IntegerSubsetStorage & VSLT, int const& pos)
{
  int posAfter=VSLT.ListNext[VSLT.MaxElement];
  VSLT.ListNext[VSLT.MaxElement] = pos;
  VSLT.ListNext[pos] = posAfter;
  VSLT.ListPrev[posAfter] = pos;
  VSLT.ListPrev[pos] = VSLT.MaxElement;
}



void VSLT_RemoveValue(IntegerSubsetStorage & VSLT, int const& pos)
{
  int posNext=VSLT.ListNext[pos];
  int posPrev=VSLT.ListPrev[pos];
  VSLT.ListNext[posPrev] = posNext;
  VSLT.ListPrev[posNext] = posPrev;
  VSLT.ListNext[pos] = -1;
  VSLT.ListPrev[pos] = -1;
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
    //    std::cerr << "MaxElement=" << MaxElement << " Maxp1=" << Maxp1 << "\n";
    while(true) {
      //      std::cerr <<  "   pos=" << pos << " ListNext[pos]=" << ListNext[pos] << "\n";
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
    }
    else {
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


/* Basic bit operations */

static constexpr uint8_t kBitmask[] = {1, 2, 4, 8, 16, 32, 64, 128};

inline bool getbit(std::vector<uint8_t> const& V, size_t const& pos)
{
  return (V[pos >> 3] >> (pos & 0x07)) & 1;
}

inline void setbit(std::vector<uint8_t> & V, size_t const& pos, bool val) {
  V[pos / 8] ^= static_cast<uint8_t>(-static_cast<uint8_t>(val) ^ V[pos / 8]) & kBitmask[pos % 8];
}

/* Container of vector of faces */


struct vectface {
private:
  size_t n;
  size_t n_face;
  std::vector<uint8_t> V;
public:
  vectface() = delete;

  vectface(size_t const& _n) : n(_n), n_face(0)
  {}

  // vectface API similar to std::vector<Face>
  void push_back(Face f)
  {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    for (size_t i=curr_len; i<needed_len; i++)
      V.push_back(0);
    //
    size_t pos = n_face * n;
    for (size_t i=0; i<n; i++) {
      bool val = f[i];
      setbit(V, pos, val);
      pos++;
    }
    n_face++;
  }

  Face operator[](size_t i_orb) const
  {
    Face f(n);
    size_t pos = i_orb * n;
    for (size_t i=0; i<n; i++) {
      f[i] = getbit(V, pos);
      pos++;
    }
    return f;
  }

  size_t size() const
  {
    return n_face;
  }

  void pop_back()
  {
    n_face--;
  }

  Face pop()
  {
    n_face--;
    Face f(n);
    size_t pos = n_face * n;
    for (size_t i=0; i<n; i++) {
      f[i] = getbit(V, pos);
      pos++;
    }
    return f;
  }

  // non standard API
  template<typename F>
  void InsertFace(F fct)
  {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    for (size_t i=curr_len; i<needed_len; i++)
      V.push_back(0);
    //
    size_t pos = n_face * n;
    for (size_t i=0; i<n; i++) {
      bool val = fct(i);
      setbit(V, pos, val);
      pos++;
    }
    n_face++;
  }

  void SetFace(Face & f, size_t i_orb) const
  {
    size_t pos = i_orb * n;
    for (size_t i=0; i<n; i++) {
      f[i] = getbit(V, pos);
      pos++;
    }
  }

  void append(vectface const& w)
  {
    size_t curr_len = V.size();
    size_t n_bits = ( n_face + w.n_face ) * n;
    size_t needed_len = (n_bits + 7) / 8;
    for (size_t i=curr_len; i<needed_len; i++)
      V.push_back(0);
    // Now appending
    size_t pos = n_face * n;
    size_t depl = w.n_face * n;
    for (size_t i=0; i<depl; i++) {
      bool val = getbit(w.V, i);
      setbit(V, pos, val);
      pos++;
    }
  }

  // Iterating stuff
private:
  struct IteratorContain {
  private:
    const vectface & v;
    size_t pos;
    Face f;
  public:
    IteratorContain(vectface const& _v, size_t const& _pos) : v(_v), pos(_pos), f(_v.n)
    {}
    Face const& operator*()
    {
      v.SetFace(f, pos);
      return f;
    }
    IteratorContain& operator++()
    {
      pos++;
      return *this;
    }
    IteratorContain operator++(int)
    {
      IteratorContain tmp = *this;
      pos++;
      return tmp;
    }
    bool operator!=(IteratorContain const& iter)
    {
      return pos != iter.pos;
    }
    bool operator==(IteratorContain const& iter)
    {
      return pos == iter.pos;
    }
  };
public:
  using iterator = IteratorContain;
  using const_iterator = IteratorContain;
  const_iterator cbegin() const
  {
    return IteratorContain(*this, 0);
  }
  const_iterator cend() const
  {
    return IteratorContain(*this, n_face);
  }
  const_iterator begin() const
  {
    return IteratorContain(*this, 0);
  }
  const_iterator end() const
  {
    return IteratorContain(*this, n_face);
  }
};






#endif
