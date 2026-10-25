// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "MAT_SparseMatrixMod.h"
#include "Timings.h"
#include <fstream>
#include <string>
#include <vector>
// clang-format on

/*
  Benchmark of the cycle lifting job modulo p: given sparse rows R_0, ...,
  R_{n-1} in F_p^m (in the homology of a cell complex, R_c is the cycle
  condition of the cell c: its boundary, or coboundary), find
  - the rank of the rows,
  - a set of free indices F with |F| = n - rank such that the vectors x
    with sum_c x_c R_c = 0 (the cycles) are determined by their values on
    F, and
  - for some f in F the cycle lift(f), equal to 1 at f and 0 on F - {f}.
  The homology needs the lifts of its basis cells only (a few of them).

  The instance file is
    p n_row n_col
    then one line per row: its number of entries, then the pairs
    (column, value) with the values in [0, p).

  The methods:
  - echelon: the incremental row echelon form of the rows (SparseEchelonMod,
    the columns ordered by increasing number of rows containing them), the
    rows independent of the previous ones being the set S = complement of
    F, then the echelon form of the rows of S marked by unit vectors, whose
    reduction of R_f expresses it in the rows of S. This is the method used
    up to October 2026 for the Hecke operators of polyhedral_common.
  - kernel: SparseKernelMod on the constraints, the columns of the matrix
    (R_c)_c, eliminated with global pivots.
  Both check their lifts against the rows.
 */

struct Instance {
  int64_t p;
  int n_row;
  int n_col;
  std::vector<SparseVectorMod> l_row;
};

Instance read_instance(std::string const &file) {
  std::ifstream is(file);
  if (!is.good()) {
    std::cerr << "Bench_SparseKernelMod: cannot open " << file << "\n";
    throw TerminalException{1};
  }
  Instance inst;
  is >> inst.p >> inst.n_row >> inst.n_col;
  inst.l_row.resize(inst.n_row);
  for (int i = 0; i < inst.n_row; i++) {
    int len;
    is >> len;
    inst.l_row[i].resize(len);
    for (int u = 0; u < len; u++) {
      is >> inst.l_row[i][u].first >> inst.l_row[i][u].second;
    }
  }
  if (!is.good()) {
    std::cerr << "Bench_SparseKernelMod: the file " << file << " is truncated\n";
    throw TerminalException{1};
  }
  return inst;
}

// Whether sum_c x_c R_c = 0.
bool is_cycle(Instance const &inst, SparseVectorMod const &x) {
  std::unordered_map<int, int64_t> sum;
  for (auto &e : x) {
    for (auto &f : inst.l_row[e.first]) {
      int64_t &val = sum[f.first];
      val = (val + e.second * f.second) % inst.p;
    }
  }
  for (auto &kv : sum) {
    if (kv.second != 0) {
      return false;
    }
  }
  return true;
}

struct LiftResult {
  size_t rank;
  std::vector<int> l_free;
  std::vector<SparseVectorMod> l_lift;
};

LiftResult method_echelon(Instance const &inst, size_t const &n_lift) {
  int64_t p = inst.p;
  int n_row = inst.n_row;
  int n_col = inst.n_col;
  MicrosecondTime time;
  // The columns by increasing number of rows containing them.
  std::vector<int> l_count(n_col, 0);
  for (auto &v : inst.l_row) {
    for (auto &e : v) {
      l_count[e.first] += 1;
    }
  }
  std::vector<int> l_order(n_col);
  for (int c = 0; c < n_col; c++) {
    l_order[c] = c;
  }
  std::stable_sort(l_order.begin(), l_order.end(),
                   [&](int const &a, int const &b) -> bool { return l_count[a] < l_count[b]; });
  std::vector<int> l_label(n_col);
  for (int c = 0; c < n_col; c++) {
    l_label[l_order[c]] = c;
  }
  auto f_relabel = [&](SparseVectorMod const &v) -> SparseVectorMod {
    SparseVectorMod w;
    w.reserve(v.size());
    for (auto &e : v) {
      w.push_back({l_label[e.first], e.second});
    }
    std::sort(w.begin(), w.end());
    return w;
  };
  SparseEchelonMod echelon(p);
  std::vector<uint8_t> in_S(n_row, 0);
  for (int c = 0; c < n_row; c++) {
    if (!inst.l_row[c].empty() && echelon.insert(f_relabel(inst.l_row[c]))) {
      in_S[c] = 1;
    }
  }
  std::cerr << "|echelon: echelon of the rows, rank=" << echelon.rank() << "|=" << time << "\n";
  SparseEchelonMod echelon_mark(p);
  size_t n_entry = 0;
  for (int c = 0; c < n_row; c++) {
    if (in_S[c] == 1) {
      SparseVectorMod v = f_relabel(inst.l_row[c]);
      v.push_back({n_col + c, 1});
      echelon_mark.insert(v);
    }
  }
  for (size_t i = 0; i < echelon_mark.rank(); i++) {
    n_entry += echelon_mark.get_row(i).size();
  }
  std::cerr << "|echelon: marked echelon, n_entry=" << n_entry << "|=" << time << "\n";
  LiftResult res;
  res.rank = echelon.rank();
  for (int c = 0; c < n_row; c++) {
    if (in_S[c] == 0) {
      res.l_free.push_back(c);
    }
  }
  for (size_t i = 0; i < std::min(n_lift, res.l_free.size()); i++) {
    int c = res.l_free[i];
    SparseVectorMod v{{c, 1}};
    if (!inst.l_row[c].empty()) {
      SparseVectorMod red = echelon_mark.reduce(f_relabel(inst.l_row[c]));
      for (auto &e : red) {
        if (e.first < n_col) {
          std::cerr << "Bench_SparseKernelMod: a row is not in the span of S\n";
          throw TerminalException{1};
        }
        v.push_back({e.first - n_col, e.second});
      }
      std::sort(v.begin(), v.end());
    }
    res.l_lift.emplace_back(std::move(v));
  }
  std::cerr << "|echelon: " << res.l_lift.size() << " lifts|=" << time << "\n";
  return res;
}

LiftResult method_kernel(Instance const &inst, size_t const &n_lift) {
  MicrosecondTime time;
  // The constraints: the columns of the matrix whose rows are the R_c.
  std::vector<SparseVectorMod> l_constraint(inst.n_col);
  for (int c = 0; c < inst.n_row; c++) {
    for (auto &e : inst.l_row[c]) {
      l_constraint[e.first].push_back({c, e.second});
    }
  }
  std::cerr << "|kernel: transposition|=" << time << "\n";
  SparseKernelMod ker(l_constraint, inst.n_row, inst.p);
  std::cerr << "|kernel: elimination, rank=" << ker.rank() << " n_sparse_pivot=" << ker.n_sparse_pivot()
            << " n_dense_col=" << ker.n_dense_col() << " n_entry=" << ker.n_entry() << "|=" << time << "\n";
  LiftResult res;
  res.rank = ker.rank();
  res.l_free = ker.free_columns();
  for (size_t i = 0; i < std::min(n_lift, res.l_free.size()); i++) {
    res.l_lift.emplace_back(ker.lift(i));
  }
  std::cerr << "|kernel: " << res.l_lift.size() << " lifts|=" << time << "\n";
  return res;
}

int main(int argc, char *argv[]) {
  HumanTime time_total;
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "Bench_SparseKernelMod [instance] [method] [n_lift]\n";
      std::cerr << "  method: echelon or kernel\n";
      std::cerr << "  n_lift: the number of lifts computed and checked (default 20)\n";
      return -1;
    }
    std::string file = argv[1];
    std::string method = argv[2];
    size_t n_lift = 20;
    if (argc == 4) {
      n_lift = std::stoul(argv[3]);
    }
    MicrosecondTime time;
    Instance inst = read_instance(file);
    size_t n_entry = 0;
    for (auto &v : inst.l_row) {
      n_entry += v.size();
    }
    std::cerr << "p=" << inst.p << " n_row=" << inst.n_row << " n_col=" << inst.n_col << " n_entry=" << n_entry
              << "\n";
    std::cerr << "|read|=" << time << "\n";
    LiftResult res;
    if (method == "echelon") {
      res = method_echelon(inst, n_lift);
    } else if (method == "kernel") {
      res = method_kernel(inst, n_lift);
    } else {
      std::cerr << "Bench_SparseKernelMod: unknown method " << method << "\n";
      return -1;
    }
    if (res.rank + res.l_free.size() != static_cast<size_t>(inst.n_row)) {
      std::cerr << "Bench_SparseKernelMod: rank + |F| = " << res.rank + res.l_free.size() << " and not n_row\n";
      return 1;
    }
    std::set<int> set_free(res.l_free.begin(), res.l_free.end());
    for (size_t i = 0; i < res.l_lift.size(); i++) {
      for (auto &e : res.l_lift[i]) {
        if (set_free.count(e.first) > 0 && (e.first != res.l_free[i] || e.second != 1)) {
          std::cerr << "Bench_SparseKernelMod: the lift " << i << " is not a unit vector on F\n";
          return 1;
        }
      }
      if (!is_cycle(inst, res.l_lift[i])) {
        std::cerr << "Bench_SparseKernelMod: the lift " << i << " is not a cycle\n";
        return 1;
      }
    }
    std::cout << "method=" << method << " rank=" << res.rank << " dim_kernel=" << res.l_free.size()
              << " n_lift_checked=" << res.l_lift.size() << "\n";
    std::cerr << "Normal termination of Bench_SparseKernelMod\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Bench_SparseKernelMod\n";
    exit(e.eVal);
  }
  runtime(time_total);
}
