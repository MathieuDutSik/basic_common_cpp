#ifndef INCLUDE_FACE_BITSET_KERNEL
#define INCLUDE_FACE_BITSET_KERNEL

// Boost libraries

#include <bitset>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Face;

// Those are needed for the tsl::sparse_map
#define IMPLEMENT_COPY_OPERATOR



/* Basic bit operations */

static constexpr uint8_t kBitmask[] = {1, 2, 4, 8, 16, 32, 64, 128};

inline bool getbit(std::vector<uint8_t> const& V, size_t const& pos)
{
  return (V[pos >> 3] >> (pos & 0x07)) & 1;
}

inline void setbit(std::vector<uint8_t> & V, size_t const& pos, bool val) {
  V[pos / 8] ^= static_cast<uint8_t>(-static_cast<uint8_t>(val) ^ V[pos / 8]) & kBitmask[pos % 8];
}

inline void setbit_ptr(uint8_t* arr, size_t const& pos, bool val) {
  arr[pos / 8] ^= static_cast<uint8_t>(-static_cast<uint8_t>(val) ^ arr[pos / 8]) & kBitmask[pos % 8];
}

/* Container of vector of faces */

struct vectface {
public:
  size_t n;
  size_t n_face;
  std::vector<uint8_t> V;
  // Constructors, move operators and the like
  vectface() : n(0), n_face(0)
  {}

  vectface(size_t const& _n) : n(_n), n_face(0)
  {}

  vectface(vectface&& vf) : n(vf.n), n_face(vf.n_face), V(std::move(vf.V))
  {
  }

#ifdef IMPLEMENT_COPY_OPERATOR
  vectface(const vectface& vf) : n(vf.n), n_face(vf.n_face), V(vf.V)
  {
  }
  vectface& operator=(const vectface& vf) : n(vf.n), n_face(vf.n_face), V(vf.V)
  {
    return *this;
  }
#else
  vectface(const vectface&) = delete;
  vectface& operator=(const vectface&) = delete;
#endif
  
  

  // The actual API

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
    n_face += w.n_face;
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






template<typename T>
T getsetasint(const Face& face)
{
  size_t len = face.size();
  T eSum = 0;
  T pow = 1;
  for (size_t i=0; i<len; i++) {
    if (face[i] == 1)
      eSum += pow;
    pow *= 2;

  }
  return eSum;
}


#endif
