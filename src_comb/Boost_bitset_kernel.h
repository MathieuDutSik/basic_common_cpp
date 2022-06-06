// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COMB_BOOST_BITSET_KERNEL_H_
#define SRC_COMB_BOOST_BITSET_KERNEL_H_

// Boost libraries

#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <boost/serialization/nvp.hpp>
#include <utility>
#include <vector>

typedef boost::dynamic_bitset<> Face;

// Those are needed for the tsl::sparse_map

/* Basic bit operations */

static constexpr uint8_t kBitmask[] = {1, 2, 4, 8, 16, 32, 64, 128};

inline bool getbit(std::vector<uint8_t> const &V, size_t const &pos) {
  return (V[pos >> 3] >> (pos & 0x07)) & 1;
}

inline void setbit(std::vector<uint8_t> &V, size_t const &pos, bool val) {
  V[pos / 8] ^= static_cast<uint8_t>(-static_cast<uint8_t>(val) ^ V[pos / 8]) &
                kBitmask[pos % 8];
}

inline void setbit_ptr(uint8_t *arr, size_t const &pos, bool val) {
  arr[pos / 8] ^=
      static_cast<uint8_t>(-static_cast<uint8_t>(val) ^ arr[pos / 8]) &
      kBitmask[pos % 8];
}

/* Container of vector of faces */

struct vectface {
public:
  size_t n;
  size_t n_face;
  std::vector<uint8_t> V;
  size_t append_len;
  std::vector<uint8_t> Vappend;
  // Constructors, move operators and the like
  vectface() : n(0), n_face(0) {}

  vectface(size_t const &_n)
      : n(_n), n_face(0), append_len((n + 7) / 8),
        Vappend(std::vector<uint8_t>(append_len, 0)) {}

  vectface(vectface &&vf)
      : n(vf.n), n_face(vf.n_face), V(std::move(vf.V)), append_len((n + 7) / 8),
        Vappend(std::vector<uint8_t>(append_len, 0)) {}

  // Serialization related stuff
  const std::vector<uint8_t> &serial_get_std_vector_uint8_t() const {
    return V;
  }
  void build_vectface(const size_t &_n, const size_t &_n_face,
                      std::vector<uint8_t> &&_V) {
    n = _n;
    n_face = _n_face;
    V = std::move(_V);
    append_len = (n + 7) / 8;
    Vappend = std::vector<uint8_t>(append_len, 0);
  }

  vectface &operator=(const vectface &&vf) {
    n = vf.n;
    n_face = vf.n_face;
    V = std::move(vf.V);
    return *this;
  }

  vectface(const vectface &) = delete;
  vectface &operator=(const vectface &) = delete;

  // The actual API

  // vectface API similar to std::vector<Face>
  void push_back(const Face &f) {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      size_t delta = needed_len - curr_len;
      V.insert(V.end(), Vappend.begin(), Vappend.begin() + delta);
    }
    //
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      bool val = f[i];
      setbit(V, pos, val);
      pos++;
    }
    n_face++;
  }

  Face operator[](size_t i_orb) const {
    Face f(n);
    size_t pos = i_orb * n;
    for (size_t i = 0; i < n; i++) {
      f[i] = getbit(V, pos);
      pos++;
    }
    return f;
  }

  size_t size() const { return n_face; }

  size_t get_n() const { return n; }

  void pop_back() { n_face--; }

  Face pop() {
    n_face--;
    Face f(n);
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      f[i] = getbit(V, pos);
      pos++;
    }
    return f;
  }

  bool operator==(vectface const &vf) {
    if (n != vf.n)
      return false;
    if (n_face != vf.n_face)
      return false;
    size_t n_elt = n * n_face;
    size_t q   = n_elt / 8;
    size_t q8  = q * 8;
    for (size_t u = 0; u < q; u++)
      if (V[u] != vf.V[u])
        return false;
    for (size_t i = q8; i<n_elt; i++)
      if (getbit(V, i) != getbit(vf.V, i))
        return false;
    return true;
  }

  bool operator!=(vectface const &vf) {
    if (n != vf.n)
      return true;
    if (n_face != vf.n_face)
      return true;
    size_t n_elt = n * n_face;
    size_t q   = n_elt / 8;
    size_t q8  = q * 8;
    for (size_t u = 0; u < q; u++)
      if (V[u] != vf.V[u])
        return true;
    for (size_t i = q8; i<n_elt; i++)
      if (getbit(V, i) != getbit(vf.V, i))
        return true;
    return false;
  }

  // non standard API
  template <typename F> void InsertFace(F fct) {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      size_t delta = needed_len - curr_len;
      V.insert(V.end(), Vappend.begin(), Vappend.begin() + delta);
    }
    //
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      bool val = fct(i);
      setbit(V, pos, val);
      pos++;
    }
    n_face++;
  }

  template <typename F> void InsertFaceRef(F &fct) {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      size_t delta = needed_len - curr_len;
      V.insert(V.end(), Vappend.begin(), Vappend.begin() + delta);
    }
    //
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      bool val = fct(i);
      setbit(V, pos, val);
      pos++;
    }
    n_face++;
  }

  void SetFace(Face &f, size_t i_orb) const {
    size_t pos = i_orb * n;
    for (size_t i = 0; i < n; i++) {
      f[i] = getbit(V, pos);
      pos++;
    }
  }

  void append(vectface const &w) {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + w.n_face) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      size_t delta = needed_len - curr_len;
      size_t n_iter = delta / append_len;
      for (size_t i_iter = 0; i_iter < n_iter; i_iter++)
        V.insert(V.end(), Vappend.begin(), Vappend.begin() + append_len);
      size_t res = delta % append_len;
      if (res > 0)
        V.insert(V.end(), Vappend.begin(), Vappend.begin() + res);
    }
    // Now appending
    size_t pos = n_face * n;
    size_t depl = w.n_face * n;
    for (size_t i = 0; i < depl; i++) {
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
    const vectface &v;
    size_t pos;
    Face f;

  public:
    IteratorContain(vectface const &_v, size_t const &_pos)
        : v(_v), pos(_pos), f(_v.n) {}
    Face const &operator*() {
      v.SetFace(f, pos);
      return f;
    }
    IteratorContain &operator++() {
      pos++;
      return *this;
    }
    IteratorContain operator++(int) {
      IteratorContain tmp = *this;
      pos++;
      return tmp;
    }
    IteratorContain &operator--() {
      pos--;
      return *this;
    }
    IteratorContain operator--(int) {
      IteratorContain tmp = *this;
      pos--;
      return tmp;
    }
    bool operator!=(IteratorContain const &iter) { return pos != iter.pos; }
    bool operator==(IteratorContain const &iter) { return pos == iter.pos; }
    friend std::ptrdiff_t operator-(IteratorContain const &x,
                                    IteratorContain const &y) {
      return x.pos - y.pos;
    }
  };

public:
  using iterator = IteratorContain;
  using const_iterator = IteratorContain;
  const_iterator cbegin() const { return IteratorContain(*this, 0); }
  const_iterator cend() const { return IteratorContain(*this, n_face); }
  const_iterator begin() const { return IteratorContain(*this, 0); }
  const_iterator end() const { return IteratorContain(*this, n_face); }
};

template <> struct std::iterator_traits<typename vectface::iterator> {
  using value_type = Face;
  using difference_type = std::ptrdiff_t;
};

/*
template<>
void std::swap(Face & x, Face & y)
{
  size_t siz = x.size();
  if (siz != y.size()) {
    Face tmp = x;
    x = y;
    y = tmp;
  } else {
    for (size_t i=0; i<siz; i++) {
      int val = x[i];
      x[i] = y[i];
      y[i] = val;
    }
  }
}
*/

namespace boost::serialization {

template <class Archive>
inline void load(Archive &ar, vectface &val,
                 [[maybe_unused]] const unsigned int version) {
  size_t n, n_face, len_vect;
  ar &make_nvp("n", n);
  ar &make_nvp("n_face", n_face);
  ar &make_nvp("len_vect", len_vect);
  std::vector<uint8_t> V(len_vect);
  for (size_t u = 0; u < len_vect; u++)
    ar &make_nvp("Vu", V[u]);
  val.build_vectface(n, n_face, std::move(V));
}

template <class Archive>
inline void save(Archive &ar, vectface const &val,
                 [[maybe_unused]] const unsigned int version) {
  size_t n = val.get_n();
  size_t n_face = val.size();
  ar &make_nvp("n", n);
  ar &make_nvp("n_face", n_face);
  const std::vector<uint8_t> &V = val.serial_get_std_vector_uint8_t();
  size_t len_vect = V.size();
  ar &make_nvp("len_vect", len_vect);
  for (size_t u = 0; u < len_vect; u++)
    ar &make_nvp("Vu", V[u]);
}

template <class Archive>
inline void serialize(Archive &ar, vectface &val, const unsigned int version) {
  split_free(ar, val, version);
}

} // namespace boost::serialization

template <typename T> T getsetasint(const Face &face) {
  size_t len = face.size();
  T eSum = 0;
  T pow = 1;
  for (size_t i = 0; i < len; i++) {
    if (face[i] == 1)
      eSum += pow;
    pow *= 2;
  }
  return eSum;
}

// clang-format off
#endif  // SRC_COMB_BOOST_BITSET_KERNEL_H_
// clang-format on
