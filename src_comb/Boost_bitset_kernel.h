// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COMB_BOOST_BITSET_KERNEL_H_
#define SRC_COMB_BOOST_BITSET_KERNEL_H_

// Boost libraries

#include "boost_serialization.h"
#include "hash_functions.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <utility>
#include <vector>

typedef boost::dynamic_bitset<> Face;

namespace boost::serialization {

template <class Archive>
inline void load(Archive &ar, Face &val,
                 [[maybe_unused]] const unsigned int version) {
  size_t n;
  ar &make_nvp("n", n);
  val = Face(n);
  for (size_t u = 0; u < n; u++) {
    int scal;
    ar &make_nvp("Vu", scal);
    val[u] = scal;
  }
}

template <class Archive>
inline void save(Archive &ar, Face const &val,
                 [[maybe_unused]] const unsigned int version) {
  size_t n = val.size();
  ar &make_nvp("n", n);
  for (size_t u = 0; u < n; u++) {
    int scal = val[u];
    ar &make_nvp("Vu", scal);
  }
}

template <class Archive>
inline void serialize(Archive &ar, Face &val, const unsigned int version) {
  split_free(ar, val, version);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// Those are needed for the tsl::sparse_map

/* Basic bit operations */

static constexpr uint8_t kBitmask[] = {1, 2, 4, 8, 16, 32, 64, 128};

inline bool getbit_vector(std::vector<uint8_t> const &V, size_t const &pos) {
  return (V[pos >> 3] >> (pos & 0x07)) & 1;
}

inline void setbit_vector(std::vector<uint8_t> &V, size_t const &pos, bool val) {
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

  vectface(const size_t &_n, const std::vector<size_t> &l_n_face,
           std::vector<std::vector<uint8_t>> const &l_V) {
    n = _n;
    n_face = 0;
    for (auto &e_n : l_n_face)
      n_face += e_n;
    //
    size_t n_tot = (n_face * n + 7) / 8;
    V = std::vector<uint8_t>(n_tot);
    size_t pos = 0;
    for (size_t u = 0; u < l_n_face.size(); u++) {
      size_t const &e_n_face = l_n_face[u];
      std::vector<uint8_t> const &eV = l_V[u];
      for (size_t v = 0; v < e_n_face * n; v++) {
        bool val = getbit_vector(eV, v);
        setbit_vector(V, pos, val);
        pos++;
      }
    }
    append_len = (n + 7) / 8;
    Vappend = std::vector<uint8_t>(append_len, 0);
  }
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
  size_t hash() {
    // We can have some insertions followed by removals, this has to be handled
    // The relevant length according to the number of faces
    int rel_len = (n * n_face + 7) / 8;
    // Setting the bits to 0 (could have been inserted and then pop_back)
    int n_bit = n * n_face;
    int rel_bit = rel_len * 8;
    for (int i_bit=n_bit; i_bit<rel_bit; i_bit++)
      setbit_vector(V, i_bit, false);
    // Now computing with a standard seed.
    uint32_t seed = 0x1b853560;
    return robin_hood_hash_bytes(V.data(), rel_len, seed);
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
    // We have to handle the fact that vectface was built with default
    // constructor. And so we assign n
    n = f.size();
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      V.resize(needed_len, uint8_t(0));
    }
    //
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      bool val = f[i];
      setbit_vector(V, pos, val);
      pos++;
    }
    n_face++;
  }

  Face operator[](size_t i_orb) const {
    Face f(n);
    size_t pos = i_orb * n;
    for (size_t i = 0; i < n; i++) {
      f[i] = getbit_vector(V, pos);
      pos++;
    }
    return f;
  }

  // This exist because we did not manage to code a operator[] that allow modification.
  void AssignEntry(Face const& f, size_t i_orb) {
    size_t pos = i_orb * n;
    for (size_t i = 0; i < n; i++) {
      bool val = f[i];
      setbit_vector(V, pos, val);
      pos++;
    }
  }

  size_t size() const { return n_face; }

  size_t get_n() const { return n; }

  void pop_back() { n_face--; }

  void clear() { n_face = 0; }

  Face pop() {
    n_face--;
    Face f(n);
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      f[i] = getbit_vector(V, pos);
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
    size_t q = n_elt / 8;
    size_t q8 = q * 8;
    for (size_t u = 0; u < q; u++)
      if (V[u] != vf.V[u])
        return false;
    for (size_t i = q8; i < n_elt; i++)
      if (getbit_vector(V, i) != getbit_vector(vf.V, i))
        return false;
    return true;
  }

  bool operator!=(vectface const &vf) {
    return !(*this == vf);
  }

  // non standard API
  template <typename F> void InsertFace(F fct) {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      V.resize(needed_len, uint8_t(0));
    }
    //
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      bool val = fct(i);
      setbit_vector(V, pos, val);
      pos++;
    }
    n_face++;
  }

  template <typename F> void InsertFaceRef(F &fct) {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + 1) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      V.resize(needed_len, uint8_t(0));
    }
    //
    size_t pos = n_face * n;
    for (size_t i = 0; i < n; i++) {
      bool val = fct(i);
      setbit_vector(V, pos, val);
      pos++;
    }
    n_face++;
  }

  void SetFace(Face &f, size_t i_orb) const {
    size_t pos = i_orb * n;
    for (size_t i = 0; i < n; i++) {
      f[i] = getbit_vector(V, pos);
      pos++;
    }
  }

  void append(vectface const &w) {
    size_t curr_len = V.size();
    size_t n_bits = (n_face + w.n_face) * n;
    size_t needed_len = (n_bits + 7) / 8;
    if (curr_len < needed_len) {
      V.resize(needed_len, uint8_t(0));
    }
    // Now appending
    size_t pos = n_face * n;
    size_t depl = w.n_face * n;
    for (size_t i = 0; i < depl; i++) {
      bool val = getbit_vector(w.V, i);
      setbit_vector(V, pos, val);
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

vectface select_minimum_count(vectface const& vf) {
  size_t n = vf.n;
  size_t min_incd = std::numeric_limits<size_t>::max();
  for (auto & eFace : vf) {
    size_t incd = eFace.count();
    if (incd < min_incd)
      min_incd = incd;
  }
  vectface vf_ret(n);
  for (auto & eFace : vf) {
    size_t incd = eFace.count();
    if (incd == min_incd)
      vf_ret.push_back(eFace);
  }
  return vf_ret;
}

vectface sort_vectface(vectface const& vf) {
  size_t n = vf.get_n();
  std::set<Face> set;
  for (auto & f : vf) {
    set.insert(f);
  }
  vectface vf_ret(n);
  for (auto & f : set) {
    vf_ret.push_back(f);
  }
  return vf_ret;
}

vectface unicize_vectface(vectface const& vf) {
  return sort_vectface(vf);
}



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

// clang-format off
}  // namespace boost::serialization
// clang-format on

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
