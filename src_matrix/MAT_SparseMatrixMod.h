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
#include <set>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
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
  A reduced row echelon form modulo p built incrementally: every row has
  its pivot entry equal to 1 and no entry at the pivot columns of the
  other rows. The reduction of a vector is then a single pass (the rows
  added have no entry at a pivot column) and gives the normal form of the
  vector modulo the span of the rows, with a zero entry at every pivot
  column.

  This is efficient when the space spanned has a small codimension in the
  space of the entries met: the rows then have few entries (their pivot
  and the columns that are not pivots). The pivot of a new row is the
  column of the row met by the fewest other rows (the Markowitz criterion),
  which keeps the elimination of that column from the other rows cheap.
 */
class SparseRrefMod {
public:
  explicit SparseRrefMod(int64_t const &_p) : p(_p) {}
  int64_t prime() const { return p; }
  size_t rank() const { return rows.size(); }
  bool is_pivot(int const &col) const { return pivot_row.count(col) > 0; }
  SparseVectorMod reduce(SparseVectorMod const &v) const {
    std::vector<int> touched;
    auto f_add = [&](int const &idx, int64_t const &val) -> void {
      if (static_cast<size_t>(idx) >= buf.size()) {
        buf.resize(2 * idx + 1, 0);
        flag.resize(2 * idx + 1, 0);
      }
      if (flag[idx] == 0) {
        flag[idx] = 1;
        touched.push_back(idx);
      }
      buf[idx] = ResidueModPrime(buf[idx] + val, p);
    };
    for (auto &e : v) {
      auto iter = pivot_row.find(e.first);
      if (iter == pivot_row.end()) {
        f_add(e.first, e.second);
        continue;
      }
      int64_t c = ResidueModPrime(e.second, p);
      for (auto &f : rows[iter->second]) {
        if (f.first != e.first) {
          f_add(f.first, p - (c * f.second) % p);
        }
      }
    }
    std::sort(touched.begin(), touched.end());
    SparseVectorMod w;
    for (auto &idx : touched) {
      if (buf[idx] != 0) {
        w.push_back({idx, buf[idx]});
      }
      buf[idx] = 0;
      flag[idx] = 0;
    }
    return w;
  }
  // Inserts v. Returns whether v is independent of the rows.
  bool insert(SparseVectorMod const &v) {
    SparseVectorMod w = reduce(v);
    if (w.empty()) {
      return false;
    }
    // The pivot: the column of w contained in the fewest rows.
    int col = -1;
    size_t best = 0;
    for (auto &e : w) {
      auto iter = col_rows.find(e.first);
      size_t cnt = (iter == col_rows.end()) ? 0 : iter->second.size();
      if (col == -1 || cnt < best) {
        col = e.first;
        best = cnt;
      }
    }
    int64_t val_col = 0;
    for (auto &e : w) {
      if (e.first == col) {
        val_col = e.second;
      }
    }
    SparseVectorMod w_norm = SparseVectorScaleMod(w, InverseModPrime(val_col, p), p);
    size_t i_new = rows.size();
    // The column col is eliminated from the rows containing it.
    auto iter_col = col_rows.find(col);
    if (iter_col != col_rows.end()) {
      std::vector<size_t> l_row(iter_col->second.begin(), iter_col->second.end());
      for (auto &i_row : l_row) {
        SparseVectorMod const &r = rows[i_row];
        int64_t coef = 0;
        for (auto &e : r) {
          if (e.first == col) {
            coef = e.second;
          }
        }
        SparseVectorMod r_new = SparseVectorAddMultipleMod(r, p - coef, w_norm, p);
        // The changes of support of the row.
        size_t i1 = 0, i2 = 0;
        while (i1 < r.size() || i2 < r_new.size()) {
          if (i2 == r_new.size() || (i1 < r.size() && r[i1].first < r_new[i2].first)) {
            if (!is_pivot(r[i1].first)) {
              col_rows[r[i1].first].erase(i_row);
            }
            i1 += 1;
          } else if (i1 == r.size() || r_new[i2].first < r[i1].first) {
            col_rows[r_new[i2].first].insert(i_row);
            i2 += 1;
          } else {
            i1 += 1;
            i2 += 1;
          }
        }
        rows[i_row] = std::move(r_new);
      }
    }
    col_rows.erase(col);
    for (auto &e : w_norm) {
      if (e.first != col) {
        col_rows[e.first].insert(i_new);
      }
    }
    pivot_row[col] = i_new;
    rows.emplace_back(std::move(w_norm));
    return true;
  }

private:
  int64_t p;
  std::unordered_map<int, size_t> pivot_row;
  std::vector<SparseVectorMod> rows;
  // For the columns that are not pivots, the rows containing them.
  std::unordered_map<int, std::unordered_set<size_t>> col_rows;
  // Work space of reduce, zero between the calls.
  mutable std::vector<int64_t> buf;
  mutable std::vector<uint8_t> flag;
};

/*
  The quotient F_p^{n_col} / span(rows) of a set of sparse rows: its
  dimension, a basis made of unit vectors e_c (the basis columns) and the
  coordinates of the class of any vector in that basis. This is what the
  homology needs, the rows being the boundaries.

  The rows are eliminated in two phases.
  - The sparse phase: the pivots are chosen globally while eliminating, a
    row with a single entry if there is one (no fill-in), otherwise the
    column met by the fewest rows, with its shortest row. The pivot row,
    normalized to 1 at its pivot, is subtracted from the other rows
    containing the pivot column and removed. The pivot rows are kept in
    the order of the elimination: the row of the pivot k has no entry at
    the pivots before k, so subtracting them in that order reduces a
    vector to one supported on the columns that are not pivots.
  - The dense phase: for a matrix of a small corank with many more rows
    than columns (the boundaries of a cell complex), the remaining rows,
    mostly dependent, fill in on the last columns. Once at most
    max_dense_col columns remain, the annihilator K of the span of the
    remaining rows in the space of those columns is computed densely,
    starting from the whole space and imposing the rows one by one. The
    class of a vector is then given by its pairings with a basis of K.

  The basis columns are the columns that are neither pivots nor met by
  the remaining rows (their coordinate is the entry of the reduced
  vector) and the columns J of the dense phase for which the pairings of
  the e_c, c in J, with the basis of K form an invertible matrix P (the
  coordinates are P^{-1} times the pairings).
 */
class SparseQuotientMod {
public:
  SparseQuotientMod(std::vector<SparseVectorMod> const &l_row, int const &_n_col, int64_t const &_p,
                    size_t const &max_dense_col = 8000)
      : n_col(_n_col), p(_p) {
    int n_row = l_row.size();
    std::vector<SparseVectorMod> rows = l_row;
    std::vector<std::unordered_set<int>> col_rows(n_col);
    std::set<std::pair<size_t, int>> set_row;
    std::set<std::pair<size_t, int>> set_col;
    for (int i = 0; i < n_row; i++) {
      if (!rows[i].empty()) {
        set_row.insert({rows[i].size(), i});
      }
      for (auto &e : rows[i]) {
        if (e.first < 0 || e.first >= n_col) {
          std::cerr << "SPARSEMOD: SparseQuotientMod, an index is out of [0, n_col)\n";
          throw TerminalException{1};
        }
        col_rows[e.first].insert(i);
      }
    }
    for (int c = 0; c < n_col; c++) {
      if (!col_rows[c].empty()) {
        set_col.insert({col_rows[c].size(), c});
      }
    }
    auto f_col_change = [&](int const &c, bool const &add, int const &i) -> void {
      set_col.erase({col_rows[c].size(), c});
      if (add) {
        col_rows[c].insert(i);
      } else {
        col_rows[c].erase(i);
      }
      if (!col_rows[c].empty()) {
        set_col.insert({col_rows[c].size(), c});
      }
    };
    //
    // The sparse phase
    //
    while (!set_row.empty() && set_col.size() > max_dense_col) {
      int i_piv = -1;
      int c_piv = -1;
      if (set_row.begin()->first == 1) {
        i_piv = set_row.begin()->second;
        c_piv = rows[i_piv][0].first;
      } else {
        c_piv = set_col.begin()->second;
        size_t best = 0;
        for (auto &i : col_rows[c_piv]) {
          if (i_piv == -1 || rows[i].size() < best || (rows[i].size() == best && i < i_piv)) {
            i_piv = i;
            best = rows[i].size();
          }
        }
      }
      int64_t val_piv = 0;
      for (auto &e : rows[i_piv]) {
        if (e.first == c_piv) {
          val_piv = e.second;
        }
      }
      SparseVectorMod r_piv = SparseVectorScaleMod(rows[i_piv], InverseModPrime(val_piv, p), p);
      set_row.erase({rows[i_piv].size(), i_piv});
      for (auto &e : rows[i_piv]) {
        f_col_change(e.first, false, i_piv);
      }
      rows[i_piv].clear();
      std::vector<int> l_other(col_rows[c_piv].begin(), col_rows[c_piv].end());
      std::sort(l_other.begin(), l_other.end());
      for (auto &i : l_other) {
        SparseVectorMod const &r = rows[i];
        int64_t coef = 0;
        for (auto &e : r) {
          if (e.first == c_piv) {
            coef = e.second;
          }
        }
        SparseVectorMod r_new = SparseVectorAddMultipleMod(r, p - coef, r_piv, p);
        size_t i1 = 0, i2 = 0;
        while (i1 < r.size() || i2 < r_new.size()) {
          if (i2 == r_new.size() || (i1 < r.size() && r[i1].first < r_new[i2].first)) {
            f_col_change(r[i1].first, false, i);
            i1 += 1;
          } else if (i1 == r.size() || r_new[i2].first < r[i1].first) {
            f_col_change(r_new[i2].first, true, i);
            i2 += 1;
          } else {
            i1 += 1;
            i2 += 1;
          }
        }
        set_row.erase({r.size(), i});
        if (!r_new.empty()) {
          set_row.insert({r_new.size(), i});
        }
        rows[i] = std::move(r_new);
      }
      pivot_index[c_piv] = l_pivot.size();
      l_pivot.push_back({c_piv, std::move(r_piv)});
    }
    //
    // The dense phase, on the columns met by the remaining rows
    //
    std::vector<int> dense_pos(n_col, -1);
    for (auto &pair : set_col) {
      dense_pos[pair.second] = l_dense_col.size();
      l_dense_col.push_back(pair.second);
    }
    std::sort(l_dense_col.begin(), l_dense_col.end());
    for (size_t u = 0; u < l_dense_col.size(); u++) {
      dense_pos[l_dense_col[u]] = u;
    }
    size_t m = l_dense_col.size();
    std::vector<std::vector<int64_t>> l_kernel(m, std::vector<int64_t>(m, 0));
    for (size_t j = 0; j < m; j++) {
      l_kernel[j][j] = 1;
    }
    std::vector<int64_t> vals;
    for (auto &pair : set_row) {
      SparseVectorMod const &r = rows[pair.second];
      size_t dim_k = l_kernel.size();
      vals.assign(dim_k, 0);
      int j0 = -1;
      for (size_t j = 0; j < dim_k; j++) {
        int64_t sum = 0;
        for (auto &e : r) {
          sum = (sum + e.second * l_kernel[j][dense_pos[e.first]]) % p;
        }
        vals[j] = sum;
        if (sum != 0 && j0 == -1) {
          j0 = j;
        }
      }
      if (j0 == -1) {
        continue;
      }
      n_rank_dense += 1;
      int64_t inv_j0 = InverseModPrime(vals[j0], p);
      std::vector<int64_t> const &k0 = l_kernel[j0];
      for (size_t j = 0; j < dim_k; j++) {
        if (static_cast<int>(j) == j0 || vals[j] == 0) {
          continue;
        }
        int64_t c = (p - (vals[j] * inv_j0) % p) % p;
        std::vector<int64_t> &kj = l_kernel[j];
        for (size_t u = 0; u < m; u++) {
          if (k0[u] != 0) {
            kj[u] = (kj[u] + c * k0[u]) % p;
          }
        }
      }
      l_kernel.erase(l_kernel.begin() + j0);
    }
    //
    // The basis columns
    //
    std::vector<uint8_t> is_dense(n_col, 0);
    for (auto &c : l_dense_col) {
      is_dense[c] = 1;
    }
    for (int c = 0; c < n_col; c++) {
      if (is_dense[c] == 0 && pivot_index.count(c) == 0) {
        l_basis_col.push_back(c);
      }
    }
    n_free_sparse = l_basis_col.size();
    // The columns J of the dense phase: the pivots of an echelon form of
    // the basis of K, and the inverse of P(i, j) = K_i(J_j).
    size_t d = l_kernel.size();
    std::vector<std::vector<int64_t>> E = l_kernel;
    std::vector<size_t> l_j;
    for (size_t i = 0; i < d; i++) {
      size_t u0 = m;
      for (size_t u = 0; u < m; u++) {
        if (E[i][u] != 0) {
          u0 = u;
          break;
        }
      }
      if (u0 == m) {
        std::cerr << "SPARSEMOD: SparseQuotientMod, the basis of the annihilator is not independent\n";
        throw TerminalException{1};
      }
      l_j.push_back(u0);
      int64_t inv_u0 = InverseModPrime(E[i][u0], p);
      for (size_t i2 = i + 1; i2 < d; i2++) {
        if (E[i2][u0] != 0) {
          int64_t c = (p - (E[i2][u0] * inv_u0) % p) % p;
          for (size_t u = 0; u < m; u++) {
            E[i2][u] = (E[i2][u] + c * E[i][u]) % p;
          }
        }
      }
    }
    for (auto &u : l_j) {
      l_basis_col.push_back(l_dense_col[u]);
    }
    // P^{-1} by the Gauss-Jordan elimination of (P | I).
    std::vector<std::vector<int64_t>> A(d, std::vector<int64_t>(2 * d, 0));
    for (size_t i = 0; i < d; i++) {
      for (size_t j = 0; j < d; j++) {
        A[i][j] = l_kernel[i][l_j[j]];
      }
      A[i][d + i] = 1;
    }
    for (size_t col = 0; col < d; col++) {
      size_t piv = col;
      while (piv < d && A[piv][col] == 0) {
        piv += 1;
      }
      if (piv == d) {
        std::cerr << "SPARSEMOD: SparseQuotientMod, the matrix P is not invertible\n";
        throw TerminalException{1};
      }
      std::swap(A[piv], A[col]);
      int64_t inv_c = InverseModPrime(A[col][col], p);
      for (auto &x : A[col]) {
        x = (x * inv_c) % p;
      }
      for (size_t i = 0; i < d; i++) {
        if (i != col && A[i][col] != 0) {
          int64_t c = (p - A[i][col]) % p;
          for (size_t j = 0; j < 2 * d; j++) {
            A[i][j] = (A[i][j] + c * A[col][j]) % p;
          }
        }
      }
    }
    P_inv = std::vector<std::vector<int64_t>>(d, std::vector<int64_t>(d));
    for (size_t i = 0; i < d; i++) {
      for (size_t j = 0; j < d; j++) {
        P_inv[i][j] = A[i][d + j];
      }
    }
    kernel = std::move(l_kernel);
  }
  int64_t prime() const { return p; }
  size_t rank() const { return l_pivot.size() + n_rank_dense; }
  // The dimension of the quotient.
  size_t dim() const { return l_basis_col.size(); }
  std::vector<int> const &basis_columns() const { return l_basis_col; }
  size_t n_sparse_pivot() const { return l_pivot.size(); }
  size_t n_dense_col() const { return l_dense_col.size(); }
  // The number of entries of the pivot rows, for the fill-in.
  size_t n_entry() const {
    size_t n = 0;
    for (auto &pair : l_pivot) {
      n += pair.second.size();
    }
    return n;
  }
  // The coordinates of the class of v in the basis of the e_c for the
  // basis columns c.
  std::vector<int64_t> coordinates(SparseVectorMod const &v) const {
    std::vector<int64_t> buf(n_col, 0);
    std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> queue;
    std::unordered_set<size_t> set_queued;
    auto f_queue = [&](int const &idx) -> void {
      auto iter = pivot_index.find(idx);
      if (iter != pivot_index.end() && set_queued.insert(iter->second).second) {
        queue.push(iter->second);
      }
    };
    for (auto &e : v) {
      buf[e.first] = ResidueModPrime(buf[e.first] + e.second, p);
      f_queue(e.first);
    }
    while (!queue.empty()) {
      size_t k = queue.top();
      queue.pop();
      int64_t val = buf[l_pivot[k].first];
      if (val == 0) {
        continue;
      }
      for (auto &e : l_pivot[k].second) {
        buf[e.first] = ResidueModPrime(buf[e.first] - (val * e.second) % p, p);
        f_queue(e.first);
      }
    }
    size_t d = kernel.size();
    std::vector<int64_t> coord(n_free_sparse + d, 0);
    for (size_t i = 0; i < n_free_sparse; i++) {
      coord[i] = buf[l_basis_col[i]];
    }
    std::vector<int64_t> pairing(d, 0);
    for (size_t i = 0; i < d; i++) {
      int64_t sum = 0;
      for (size_t u = 0; u < l_dense_col.size(); u++) {
        int64_t val = buf[l_dense_col[u]];
        if (val != 0) {
          sum = (sum + val * kernel[i][u]) % p;
        }
      }
      pairing[i] = sum;
    }
    for (size_t i = 0; i < d; i++) {
      int64_t sum = 0;
      for (size_t j = 0; j < d; j++) {
        sum = (sum + P_inv[i][j] * pairing[j]) % p;
      }
      coord[n_free_sparse + i] = sum;
    }
    return coord;
  }

private:
  int n_col;
  int64_t p;
  // The pivots of the sparse phase (column, row normalized to 1 at the
  // column) in the order of the elimination.
  std::vector<std::pair<int, SparseVectorMod>> l_pivot;
  std::unordered_map<int, size_t> pivot_index;
  // The dense phase: its columns, the basis of the annihilator of the
  // remaining rows on them, its rank and the inverse of P.
  std::vector<int> l_dense_col;
  std::vector<std::vector<int64_t>> kernel;
  size_t n_rank_dense = 0;
  std::vector<std::vector<int64_t>> P_inv;
  // The basis columns: those of the sparse phase, then those of the dense
  // phase.
  std::vector<int> l_basis_col;
  size_t n_free_sparse = 0;
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
