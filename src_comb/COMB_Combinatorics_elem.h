#ifndef INCLUDE_COMBINATORICS_ELEM
#define INCLUDE_COMBINATORICS_ELEM

#include "Temp_common.h"

std::vector<int> BinomialStdvect_First(int const& k)
{
  std::vector<int> eVect(k);
  for (int i=0; i<k; i++)
    eVect[i]=i;
  return eVect;
}


bool BinomialStdvect_Increment(int const&n, int const&k, std::vector<int> & Tvect)
{
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



struct SetCppIterator {
  using Tidx = uint16_t;
private:
  Tidx dim;
  Tidx size;
  struct IteratorContain {
  private:
    Tidx dim_iter;
    Tidx size_iter;
    std::vector<Tidx> V;
    void single_increase()
    {
      V[0]++;
      Tidx xy2=1;
      while ((xy2 < size_iter) && (V[xy2-1] >= V[xy2])) {
        V[xy2]++;
        xy2++;
      }
      if (xy2 != 1) {
        for (Tidx xy1=0; xy1<xy2-1; xy1++)
          V[xy1]=xy1;
      }
      if (V[size_iter-1] == dim_iter)
        V.clear();
    }
  public:
    IteratorContain(Tidx const& eDim, Tidx const& eSize, std::vector<Tidx> const& eV) : dim_iter(eDim), size_iter(eSize), V(eV)
    {
    }
    std::vector<Tidx> const& operator*()
    {
      return V;
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
      if (iter.V.size() != V.size())
        return true;
      for (size_t i=0; i<V.size(); i++)
        if (iter.V[i] != V[i])
          return true;
      return false;
    }
    bool operator==(IteratorContain const& iter)
    {
      if (iter.dim_iter != dim_iter)
        return false;
      if (iter.size_iter != size_iter)
        return false;
      if (iter.V.size() != V.size())
        return false;
      for (size_t i=0; i<V.size(); i++)
        if (iter.V[i] != V[i])
          return false;
      return true;
    }
  };
public:
  // no copy
  SetCppIterator(const SetCppIterator&) = delete;

  // no assign
  SetCppIterator& operator=(const SetCppIterator&) = delete;

  // no move
  SetCppIterator(SetCppIterator&&) = delete;

  // no default constructor
  SetCppIterator() = delete;

  SetCppIterator(const Tidx& eDim, const Tidx& eSize) : dim(eDim), size(eSize)
  {
    if (size == 0) {
      std::cerr << "We should have a non-zero size\n";
      throw TerminalException{1};
    }
  }

  // The iterator business
  using iterator=IteratorContain;
  using const_iterator=IteratorContain;
  const_iterator cbegin() const
  {
    if (size > dim)
      return {dim, size, {}};
    std::vector<Tidx> V(size);
    for (Tidx i=0; i<size; i++)
      V[i] = i;
    return {dim, size, V};
  }
  const_iterator cend() const
  {
    return {dim, size, {}};
  }
  const_iterator begin() const
  {
    if (size > dim)
      return {dim, size, {}};
    std::vector<Tidx> V(size);
    for (Tidx i=0; i<size; i++)
      V[i] = i;
    return {dim, size, V};
  }
  const_iterator end() const
  {
    return {dim, size, {}};
  }
};





#endif


