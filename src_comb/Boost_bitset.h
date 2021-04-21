#ifndef INCLUDE_FACE_BITSET
#define INCLUDE_FACE_BITSET

// Boost libraries

#include "Temp_common.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Face;




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





std::vector<int> FaceToVector(Face const& eSet)
{
  int nbVert=eSet.count();
  std::vector<int> eList(nbVert);
  int aRow=eSet.find_first();
  for (int i=0; i<nbVert; i++) {
    eList[i]=aRow;
    aRow=eSet.find_next(aRow);
  }
  return eList;
}


std::vector<int> FaceTo01vector(Face const& eSet)
{
  int nbVert=eSet.size();
  int siz=eSet.count();
  std::vector<int> eList(nbVert,0);
  int aRow=eSet.find_first();
  for (int i=0; i<siz; i++) {
    eList[aRow]=1;
    aRow=eSet.find_next(aRow);
  }
  return eList;
}







void WriteFace(std::ostream & os, Face const& eList)
{
  int len;
  len=eList.size();
  os << len;
  for (int i=0; i<len; i++) {
    int eVal=eList[i];
    os << " " << eVal;
  }
  os << "\n";
}

Face ReadFace(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "ReadFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  int len, eVal;
  is >> len;
  Face eFace(len);
  for (int i=0; i<len; i++) {
    is >> eVal;
    eFace[i]=eVal;
  }
  return eFace;
}


vectface ReadListFace(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "ReadListFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  int nbFace;
  is >> nbFace;
  if (nbFace == 0) {
    std::cerr << "We cannot handle that case because we need the base length\n";
    throw TerminalException{1};
  }
  Face f = ReadFace(is);
  vectface ListFace(f.size());
  ListFace.push_back(f);
  for (int iFace=0; iFace<nbFace; iFace++) {
    Face f2 = ReadFace(is);
    ListFace.push_back(f2);
  }
  return ListFace;
}

void WriteListFace(std::ostream & os, vectface const& ListFace)
{
  int nbFace=ListFace.size();
  os << nbFace << "\n";
  for (int iFace=0; iFace<nbFace; iFace++)
    WriteFace(os, ListFace[iFace]);
}


void WriteFaceGAP(std::ostream &os, Face const& f)
{
  int nb=f.count();
  //  int siz=f.size();
  os << "[";
  int aPos=f.find_first();
  for (int i=0; i<nb; i++) {
    if (i>0)
      os << ",";
    int eVal=aPos+1;
    os << eVal;
    aPos=f.find_next(aPos);
  }
  os << "]";
}


void WriteListFaceGAP(std::ostream & os, vectface const& ListFace)
{
  os << "[";
  bool IsFirst=true;
  for (auto & eFace : ListFace) {
    if (!IsFirst)
      os << ",";
    IsFirst=false;
    WriteFaceGAP(os, eFace);
  }
  os << "]";
}

void WriteListFaceGAPfile(std::string const& eFile, vectface const& ListFace)
{
  std::ofstream os(eFile);
  os << "return ";
  WriteListFaceGAP(os, ListFace);
  os << ";\n";
}




// We require x and y to be of the same size
bool operator<(Face const& x, Face const& y)
{
  int len=x.size();
  for (int i=0; i<len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}




void PrintVectInt(std::ostream &os, Face const& eList)
{
  int len, i;
  len=eList.size();
  for (i=0; i<len; i++)
    if (eList[i] == 1)
      os << " " << i;
  os << "\n";
}







Face FullFace(int const& len)
{
  Face eFace(len);
  for (int u=0; u<len; u++)
    eFace[u]=1;
  return eFace;
}



ulong FaceToUnsignedLong(Face const& f)
{
  int len=f.size();
  if (len > 32) {
    std::cerr << "Too large value, conversion impossible";
    throw TerminalException{1};
  }
  ulong pos=0;
  ulong pow=1;
  for (int i=0; i<len; i++) {
    pos += pow*f[i];
    pow *= 2;
  }
  return pos;
}

Face UnsignedLongToFace(int const& len, ulong const& eVal)
{
  if (len > 32) {
    std::cerr << "length error\n";
    throw TerminalException{1};
  }
  ulong eWork = eVal;
  Face eFace(len);
  ulong pow=1;
  for (int i=0; i<len; i++) {
    ulong Pow2=pow * 2;
    ulong res=eWork % Pow2;
    if (res == pow) {
      eFace[i]=1;
      eWork -= pow;
    }
    pow=Pow2;
  }
  return eFace;
}

void VectVectInt_Magma_Print(std::ostream &os, vectface const&ListOrbit)
{
  int nbOrbit=ListOrbit.size();
  os << "[";
  for (int iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
    if (iOrbit>0)
      os << ",\n";
    Face eRepr=ListOrbit[iOrbit];
    int siz=eRepr.count();
    os << "[";
    int eVal=eRepr.find_first();
    for (int i=0; i<siz; i++) {
      if (i>0)
	os << ",";
      os << eVal;
      eVal=eRepr.find_next(eVal);
    }
    os << "]";
  }
  os << "]\n";
}



#endif
