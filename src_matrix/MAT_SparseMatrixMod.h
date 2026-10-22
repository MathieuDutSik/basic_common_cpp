// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_SPARSEMATRIXMOD_H_
#define SRC_MATRIX_MAT_SPARSEMATRIXMOD_H_

// clang-format off
#include "ExceptionsFunc.h"
#include <cstdint>
#include <iostream>
#include <functional>
#include <map>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_SPARSE_MATRIX_MOD
#endif

/*
  Sparse linear algebra modulo a prime p < 2^31, so that the products of
  two residues fit in an int64_t. This is the approach of Ash, Gunnells and
  McConnell for the cohomology of the congruence subgroups of SL_4(Z):
  the matrices of the chain complexes are very sparse, the dense rational
  elimination is out of reach for them, and working modulo a large prime
  gives the ranks and the Hecke polynomials modulo p.

  A sparse vector is the list of its nonzero entries (index, value) sorted
  by index, with the values in [1, p-1].
 */

using SparseVectorMod = std::vector<std::pair<int, int64_t>>;

// The inverse of a modulo the prime p, by the extended Euclidean algorithm.
inline int64_t InverseModPrime(int64_t const &a, int64_t const &p) {
  int64_t r0 = p, r1 = a % p;
  if (r1 < 0) {
    r1 += p;
  }
  int64_t s0 = 0, s1 = 1;
  while (r1 != 0) {
    int64_t q = r0 / r1;
    int64_t r2 = r0 - q * r1;
    r0 = r1;
    r1 = r2;
    int64_t s2 = s0 - q * s1;
    s0 = s1;
    s1 = s2;
  }
  int64_t inv = s0 % p;
  if (inv < 0) {
    inv += p;
  }
  return inv;
}

// The residue of a in [0, p-1].
inline int64_t ResidueModPrime(int64_t const &a, int64_t const &p) {
  int64_t r = a % p;
  if (r < 0) {
    r += p;
  }
  return r;
}

// c * v modulo p.
inline SparseVectorMod SparseVectorScaleMod(SparseVectorMod const &v, int64_t const &c,
                                            int64_t const &p) {
  SparseVectorMod w;
  int64_t c_red = ResidueModPrime(c, p);
  if (c_red == 0) {
    return w;
  }
  w.reserve(v.size());
  for (auto &e : v) {
    w.push_back({e.first, (e.second * c_red) % p});
  }
  return w;
}

// v1 + c * v2 modulo p.
inline SparseVectorMod SparseVectorAddMultipleMod(SparseVectorMod const &v1, int64_t const &c,
                                                  SparseVectorMod const &v2, int64_t const &p) {
  SparseVectorMod w;
  int64_t c_red = ResidueModPrime(c, p);
  if (c_red == 0) {
    return v1;
  }
  w.reserve(v1.size() + v2.size());
  size_t i1 = 0, i2 = 0;
  while (i1 < v1.size() || i2 < v2.size()) {
    if (i2 == v2.size() || (i1 < v1.size() && v1[i1].first < v2[i2].first)) {
      w.push_back(v1[i1]);
      i1 += 1;
    } else if (i1 == v1.size() || v2[i2].first < v1[i1].first) {
      w.push_back({v2[i2].first, (c_red * v2[i2].second) % p});
      i2 += 1;
    } else {
      int64_t val = (v1[i1].second + c_red * v2[i2].second) % p;
      if (val != 0) {
        w.push_back({v1[i1].first, val});
      }
      i1 += 1;
      i2 += 1;
    }
  }
  return w;
}

/*
  A row echelon form modulo p built incrementally. The rows have their
  leading entry equal to 1, at distinct pivot columns. The reduction of a
  vector subtracts the rows in the increasing order of the columns: since a
  row has no entry before its pivot column, the result has a zero at every
  pivot column. It is the normal form of the vector modulo the span of the
  rows: two vectors have the same normal form iff their difference is in
  the span.
 */
class SparseEchelonMod {
public:
  explicit SparseEchelonMod(int64_t const &_p) : p(_p) {}
  int64_t prime() const { return p; }
  size_t rank() const { return rows.size(); }
  bool is_pivot(int const &col) const { return pivot_row.count(col) > 0; }
  SparseVectorMod const &get_row(size_t const &i) const { return rows[i]; }
  // The pivot column of the row i.
  int get_pivot(size_t const &i) const { return rows[i][0].first; }
  /*
    The normal form of v. If coeff is not null, the coefficients c_i with
    v = normal form + sum_i c_i rows[i] are stored there (indexed by the
    rows).
   */
  SparseVectorMod reduce(SparseVectorMod const &v,
                         std::vector<std::pair<size_t, int64_t>> *coeff = nullptr) const {
    // A dense buffer of the residues with the indices to process in a
    // min-heap: a row only adds entries after its pivot, so each index is
    // processed once and in increasing order, which gives a sorted result.
    std::priority_queue<int, std::vector<int>, std::greater<int>> heap;
    auto f_add = [&](int const &idx, int64_t const &val) -> void {
      if (static_cast<size_t>(idx) >= buf.size()) {
        buf.resize(2 * idx + 1, 0);
        queued.resize(2 * idx + 1, 0);
      }
      buf[idx] = ResidueModPrime(buf[idx] + val, p);
      if (queued[idx] == 0) {
        queued[idx] = 1;
        heap.push(idx);
      }
    };
    for (auto &e : v) {
      f_add(e.first, e.second);
    }
    SparseVectorMod w;
    while (!heap.empty()) {
      int col = heap.top();
      heap.pop();
      queued[col] = 0;
      int64_t c = buf[col];
      buf[col] = 0;
      if (c == 0) {
        continue;
      }
      auto iter_piv = pivot_row.find(col);
      if (iter_piv == pivot_row.end()) {
        w.push_back({col, c});
        continue;
      }
      size_t i_row = iter_piv->second;
      if (coeff) {
        coeff->push_back({i_row, c});
      }
      for (auto &e : rows[i_row]) {
        if (e.first != col) {
          f_add(e.first, p - (c * e.second) % p);
        }
      }
    }
    return w;
  }
  // Inserts v. Returns whether v is independent of the rows.
  bool insert(SparseVectorMod const &v) {
    SparseVectorMod w = reduce(v);
    if (w.empty()) {
      return false;
    }
    int64_t inv = InverseModPrime(w[0].second, p);
    SparseVectorMod w_norm = SparseVectorScaleMod(w, inv, p);
    pivot_row[w_norm[0].first] = rows.size();
    rows.emplace_back(std::move(w_norm));
    return true;
  }

private:
  int64_t p;
  std::unordered_map<int, size_t> pivot_row;
  std::vector<SparseVectorMod> rows;
  // Work space of reduce, zero between the calls.
  mutable std::vector<int64_t> buf;
  mutable std::vector<uint8_t> queued;
};

/*
  A basis of the left kernel of the matrix whose rows are l_row (vectors
  with indices in [0, n_col)), modulo p: the combinations x with
  sum_i x_i l_row[i] = 0. Each row is extended by the unit vector e_i at
  the indices n_col + i and the extended rows are reduced on their first
  n_col coordinates; a row whose first part vanishes gives the kernel
  vector held in its second part.
 */
inline std::vector<SparseVectorMod> SparseLeftKernelMod(std::vector<SparseVectorMod> const &l_row,
                                                        int const &n_col, int64_t const &p) {
  SparseEchelonMod echelon(p);
  std::vector<SparseVectorMod> l_kernel;
  int n_row = l_row.size();
  for (int i = 0; i < n_row; i++) {
    SparseVectorMod v = l_row[i];
    v.push_back({n_col + i, 1});
    SparseVectorMod w = echelon.reduce(v);
    if (w.empty() || w[0].first >= n_col) {
      // The first part vanishes: the second part is a kernel vector.
      SparseVectorMod x;
      for (auto &e : w) {
        x.push_back({e.first - n_col, e.second});
      }
      l_kernel.emplace_back(std::move(x));
    } else {
      echelon.insert(w);
    }
  }
#ifdef SANITY_CHECK_SPARSE_MATRIX_MOD
  for (auto &x : l_kernel) {
    std::map<int, int64_t> sum;
    for (auto &e : x) {
      for (auto &f : l_row[e.first]) {
        int64_t &val = sum[f.first];
        val = ResidueModPrime(val + e.second * f.second, p);
      }
    }
    for (auto &kv : sum) {
      if (kv.second != 0) {
        std::cerr << "SPARSEMOD: a vector of the left kernel is not in the kernel\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return l_kernel;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_SPARSEMATRIXMOD_H_
// clang-format on
