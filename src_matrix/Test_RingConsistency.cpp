// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Consistency of a field and of its underlying ring, for the algebraic number
// fields whose underlying_ring is a genuine ring rather than the field itself.
//
// The generic code that moves a computation to underlying_ring<T> then runs it
// over a type that is neither a field nor a euclidean domain, a combination
// the matrix dispatches had not seen before, and the point of this test is
// that every such path returns the same answer as the field path it replaces.
//
// Two pairs are covered, and the suite is written once for both:
//   --- RealField over the cubic field of discriminant 49, whose generator is
//       2*cos(2*pi/7) of monic minimal polynomial X^3 + X^2 - 2X - 1, against
//       the order Z[x] on the powers of the generator;
//   --- QuadField<mpq_class, d> against QuadField<mpz_class, d>, that is
//       Q(sqrt(d)) against Z[sqrt(d)].

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryRealField.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "MAT_Matrix_SubsetSolver.h"
// clang-format on
#include <string>
#include <vector>

int const idx_field = 1;
using Tz = Tint_real_field;

// The dispatch branch this test is about: the ring is the third kind of type,
// neither a field nor a euclidean domain, and the field is its field of
// fractions.
template <typename Tfield> struct RingPair {
  using Tring = typename underlying_ring<Tfield>::ring_type;
  static_assert(!is_ring_field<Tring>::value, "the ring is not a field");
  static_assert(!is_euclidean_domain<Tring>::value,
                "the ring is not a euclidean domain");
  static_assert(std::is_same_v<typename overlying_field<Tring>::field_type,
                               Tfield>,
                "the overlying field of the ring is the field");
};

static int n_error = 0;

static void check(bool test, std::string const &name) {
  if (test) {
    std::cerr << "PASS: " << name << "\n";
  } else {
    std::cerr << "FAIL: " << name << "\n";
    n_error++;
  }
}

// A deterministic pseudo-random sequence, so a CI failure is reproducible.
struct SmallRandom {
  uint64_t state;
  explicit SmallRandom(uint64_t seed) : state(seed) {}
  int next(int modulo) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<int>((state >> 33) % modulo);
  }
  int centered(int amp) { return next(2 * amp + 1) - amp; }
};

template <typename Tfield, typename Tring>
static Tfield to_field(Tring const &x) {
  return UniversalScalarConversion<Tfield, Tring>(x);
}

template <typename Tfield, typename Tring>
static MyVector<Tfield> to_field(MyVector<Tring> const &V) {
  return UniversalVectorConversion<Tfield, Tring>(V);
}

template <typename Tfield, typename Tring>
static MyMatrix<Tfield> to_field(MyMatrix<Tring> const &M) {
  return UniversalMatrixConversion<Tfield, Tring>(M);
}

// The size of a canonical representative, used to check that the
// canonicalization really does bound the coordinates.
template <typename Tfield>
static Tfield SumAbsoluteEntries(MyVector<Tfield> const &V) {
  Tfield sum(0);
  for (int i = 0; i < V.size(); i++)
    sum += T_abs(V(i));
  return sum;
}

// A random ring element, one overload per shape of ring.
template <int i_field>
static RealRing<i_field> RandomRingElement(SmallRandom &rnd, int deg, int amp,
                                           RealRing<i_field> const &) {
  std::vector<Tz> V(deg);
  for (int u = 0; u < deg; u++)
    V[u] = rnd.centered(amp);
  return RealRing<i_field>(V);
}

template <typename T, int d>
static QuadField<T, d> RandomRingElement(SmallRandom &rnd,
                                         [[maybe_unused]] int deg, int amp,
                                         QuadField<T, d> const &) {
  return QuadField<T, d>(T(rnd.centered(amp)), T(rnd.centered(amp)));
}

template <typename Tring>
static Tring RandomRing(SmallRandom &rnd, int deg, int amp) {
  return RandomRingElement(rnd, deg, amp, Tring());
}

template <typename Tring>
static MyMatrix<Tring> RandomRingMatrix(SmallRandom &rnd, int n_row, int n_col,
                                        int deg, int amp) {
  MyMatrix<Tring> M(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++)
      M(i, j) = RandomRing<Tring>(rnd, deg, amp);
  return M;
}

// The two vectors are on the same line: u = c v with c non-zero, and with
// c > 0 when require_positive. The canonicalizations and the kernel vectors
// are only defined up to a scalar, so this is the comparison that applies to
// them: positive where the orientation is restored (the ring canonicalization
// scales the field one by a positive factor, and GetPositiveKernelVector fixes
// the sign), of either sign where it is not -- multiplying by a ring element
// of unknown sign reverses the canonical representative, which is why the
// callers of canonicalize_normal restore the orientation themselves.
template <typename Tfield>
static bool Proportional(MyVector<Tfield> const &u, MyVector<Tfield> const &v,
                         bool require_positive) {
  int n = u.size();
  if (v.size() != n)
    return false;
  int i_pivot = -1;
  for (int i = 0; i < n; i++)
    if (v(i) != 0) {
      i_pivot = i;
      break;
    }
  if (i_pivot == -1)
    return IsZeroVector(u);
  Tfield c = u(i_pivot) / v(i_pivot);
  if (c == 0)
    return false;
  if (require_positive && c < 0)
    return false;
  for (int i = 0; i < n; i++)
    if (u(i) != c * v(i))
      return false;
  return true;
}

template <typename Tfield>
static bool PositivelyProportional(MyVector<Tfield> const &u,
                                   MyVector<Tfield> const &v) {
  return Proportional(u, v, true);
}

// The square matrix operations: what the ring computes must be what the field
// computes, entry for entry.
template <typename Tfield>
static void process_square(int n, int deg, int nb) {
  using Tring = typename RingPair<Tfield>::Tring;
  SmallRandom rnd(20260906);
  bool det_ok = true, prod_ok = true, rank_ok = true, adj_ok = true;
  for (int i_test = 0; i_test < nb; i_test++) {
    MyMatrix<Tring> Mr = RandomRingMatrix<Tring>(rnd, n, n, deg, 3);
    MyMatrix<Tfield> Mf = to_field<Tfield>(Mr);
    // The determinant, which over the ring goes through Bareiss.
    if (to_field<Tfield>(DeterminantMat(Mr)) != DeterminantMat(Mf))
      det_ok = false;
    // The matrix product.
    MyMatrix<Tring> Pr = Mr * Mr;
    MyMatrix<Tfield> Pf = Mf * Mf;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (to_field<Tfield>(Pr(i, j)) != Pf(i, j))
          prod_ok = false;
    // The rank.
    if (RankMat(Mr) != RankMat(Mf))
      rank_ok = false;
    // The adjugate, the fraction-free core of the ring computation.
    if (DeterminantMat(Mf) != 0) {
      std::pair<MyMatrix<Tring>, Tring> pr = AdjugateDeterminant(Mr);
      std::pair<MyMatrix<Tfield>, Tfield> pf = AdjugateDeterminant(Mf);
      if (to_field<Tfield>(pr.second) != pf.second)
        adj_ok = false;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          if (to_field<Tfield>(pr.first(i, j)) != pf.first(i, j))
            adj_ok = false;
    }
  }
  check(det_ok, "DeterminantMat agrees between the ring and the field");
  check(prod_ok, "the matrix product agrees between the ring and the field");
  check(rank_ok, "RankMat agrees between the ring and the field");
  check(adj_ok, "AdjugateDeterminant agrees between the ring and the field");
}

// ScalarCanonicalizationVector / ScalarCanonicalizationMatrix: the ring has
// neither a gcd to reduce a content with nor a division to normalize with, so
// it normalizes through the overlying field and scales back. The result is the
// same direction as the field normalization.
template <typename Tfield>
static void process_canonicalization(int n, int deg, int nb) {
  using Tring = typename RingPair<Tfield>::Tring;
  SmallRandom rnd(777001);
  bool vec_ok = true, mat_ok = true, bounded_ok = true;
  for (int i_test = 0; i_test < nb; i_test++) {
    MyVector<Tring> Vr(n);
    for (int i = 0; i < n; i++)
      Vr(i) = RandomRing<Tring>(rnd, deg, 4);
    if (IsZeroVector(Vr))
      continue;
    MyVector<Tfield> Vf = to_field<Tfield>(Vr);
    MyVector<Tring> Cr = ScalarCanonicalizationVector(Vr);
    MyVector<Tfield> Cf = ScalarCanonicalizationVector(Vf);
    if (!PositivelyProportional(to_field<Tfield>(Cr), Cf))
      vec_ok = false;
    // The canonicalization is what keeps the coordinates from growing over the
    // repeated combinations the callers perform: re-canonicalizing a scaled up
    // vector comes back to the same representative, up to the sign that the
    // sign of the scaling factor carries.
    Tring big = RandomRing<Tring>(rnd, deg, 3);
    if (big != Tring(0)) {
      MyVector<Tring> Sr(n);
      for (int i = 0; i < n; i++)
        Sr(i) = Vr(i) * big;
      MyVector<Tring> Cr2 = ScalarCanonicalizationVector(Sr);
      if (!Proportional(to_field<Tfield>(Cr2), Cf, false))
        bounded_ok = false;
      // And the representative really is bounded: its entries do not carry
      // the scaling factor.
      if (!IsZeroVector(Cr2)) {
        Tfield n_scaled = SumAbsoluteEntries(to_field<Tfield>(Cr2));
        Tfield n_plain = SumAbsoluteEntries(to_field<Tfield>(Cr));
        if (n_scaled != n_plain)
          bounded_ok = false;
      }
    }
    MyMatrix<Tring> Mr = RandomRingMatrix<Tring>(rnd, n, n, deg, 4);
    MyMatrix<Tfield> Mf = to_field<Tfield>(Mr);
    MyMatrix<Tring> Kr = ScalarCanonicalizationMatrix(Mr);
    MyMatrix<Tfield> Kf = ScalarCanonicalizationMatrix(Mf);
    MyVector<Tfield> flat_r(n * n), flat_f(n * n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        flat_r(i * n + j) = to_field<Tfield>(Kr(i, j));
        flat_f(i * n + j) = Kf(i, j);
      }
    if (!PositivelyProportional(flat_r, flat_f))
      mat_ok = false;
  }
  check(vec_ok, "ScalarCanonicalizationVector agrees up to a positive scalar");
  check(mat_ok, "ScalarCanonicalizationMatrix agrees up to a positive scalar");
  check(bounded_ok,
        "the canonicalization brings a scaled up vector back to the same "
        "representative");
}

// CanonicalizeOrderedMatrix and, through it, SelectIndependentRowsRing. That
// selection reduces by cross multiplication, which doubles the size of the
// entries at every echelon level, and relies on NormalizeVectorContent to pull
// them back; the reduction used to be gated on is_euclidean_domain and so did
// nothing here. The canonical form must be the one the overlying field gives.
template <typename Tfield>
static void process_canonicalize_ordered(int n_row, int n_col, int deg,
                                         int nb) {
  using Tring = typename RingPair<Tfield>::Tring;
  SmallRandom rnd(5150);
  bool ok = true, sel_ok = true;
  int n_done = 0;
  for (int i_test = 0; i_test < nb; i_test++) {
    MyMatrix<Tring> Mr = RandomRingMatrix<Tring>(rnd, n_row, n_col, deg, 3);
    MyMatrix<Tfield> Mf = to_field<Tfield>(Mr);
    if (RankMat(Mf) < n_col)
      continue;
    // The independent rows selected over the ring are those selected over the
    // field: the selection is by rank, which no scaling changes.
    std::vector<int> sel_r = SelectIndependentRowsRing(Mr);
    std::vector<int> sel_f = SelectIndependentRowsRing(Mf);
    if (sel_r != sel_f)
      sel_ok = false;
    // And the canonical form is the one computed over the field, up to the
    // positive scalar the ring cannot divide by: the field representative
    // carries a rational denominator, the ring one is the integral multiple
    // of it. Being canonical on the ray is the whole contract, so this is the
    // comparison that applies; requiring equality would be requiring the ring
    // to leave itself.
    MyMatrix<Tring> Cr = CanonicalizeOrderedMatrix(Mr);
    MyMatrix<Tfield> Cf = CanonicalizeOrderedMatrix(Mf);
    if (Cr.rows() != Cf.rows() || Cr.cols() != Cf.cols()) {
      ok = false;
      continue;
    }
    MyVector<Tfield> flat_r(Cr.rows() * Cr.cols());
    MyVector<Tfield> flat_f(Cr.rows() * Cr.cols());
    for (int i = 0; i < Cr.rows(); i++)
      for (int j = 0; j < Cr.cols(); j++) {
        flat_r(i * Cr.cols() + j) = to_field<Tfield>(Cr(i, j));
        flat_f(i * Cr.cols() + j) = Cf(i, j);
      }
    if (!PositivelyProportional(flat_r, flat_f))
      ok = false;
    n_done++;
  }
  check(n_done > 0, "the ordered canonicalization test found usable matrices");
  check(sel_ok, "SelectIndependentRowsRing agrees between the ring and the field");
  check(ok, "CanonicalizeOrderedMatrix agrees with the field up to a positive scalar");
}

// SubsetRankOneSolver: the ring uses the SubsetRankOneSolver_RingOverField
// variant since the euclidean one needs a gcd. The kernel vector it returns
// must be the same ray as the one the field variant returns.
template <typename Tfield>
static void process_subset_solver(int n_row, int n_col, int deg, int nb) {
  using Tring = typename RingPair<Tfield>::Tring;
  SmallRandom rnd(424242);
  bool ok = true;
  int n_done = 0;
  for (int i_test = 0; i_test < nb; i_test++) {
    MyMatrix<Tring> EXTr = RandomRingMatrix<Tring>(rnd, n_row, n_col, deg, 3);
    MyMatrix<Tfield> EXTf = to_field<Tfield>(EXTr);
    if (RankMat(EXTf) < n_col)
      continue;
    // A subset of n_col - 1 rows of rank n_col - 1, that is of corank one.
    Face sInc(n_row);
    std::vector<int> chosen;
    for (int i = 0; i < n_row && static_cast<int>(chosen.size()) < n_col - 1;
         i++) {
      std::vector<int> cand = chosen;
      cand.push_back(i);
      MyMatrix<Tfield> Sub(cand.size(), n_col);
      for (size_t u = 0; u < cand.size(); u++)
        Sub.row(u) = EXTf.row(cand[u]);
      if (RankMat(Sub) == static_cast<int>(cand.size()))
        chosen = cand;
    }
    if (static_cast<int>(chosen.size()) != n_col - 1)
      continue;
    for (auto &i : chosen)
      sInc[i] = 1;
    SubsetRankOneSolver<Tring> solver_r(EXTr);
    SubsetRankOneSolver<Tfield> solver_f(EXTf);
    MyVector<Tring> Vr = solver_r.GetPositiveKernelVector(sInc);
    MyVector<Tfield> Vf = solver_f.GetPositiveKernelVector(sInc);
    if (!PositivelyProportional(to_field<Tfield>(Vr), Vf))
      ok = false;
    // It really is a kernel vector of the selected rows.
    for (auto &i : chosen) {
      Tring sum(0);
      for (int j = 0; j < n_col; j++)
        sum += EXTr(i, j) * Vr(j);
      if (sum != Tring(0))
        ok = false;
    }
    n_done++;
  }
  check(n_done > 0, "the subset solver test found usable configurations");
  check(ok, "SubsetRankOneSolver agrees between the ring and the field");
}

// The scaling of a field vector into the ring, which is how the callers enter
// the ring computation in the first place.
template <typename Tfield>
static void process_scaling(int n, int deg, int nb) {
  using Tring = typename RingPair<Tfield>::Tring;
  SmallRandom rnd(31337);
  bool ok = true, integral_ok = true;
  for (int i_test = 0; i_test < nb; i_test++) {
    MyVector<Tfield> Vf(n);
    for (int i = 0; i < n; i++) {
      Tring num = RandomRing<Tring>(rnd, deg, 4);
      int den = 1 + rnd.next(6);
      Vf(i) = to_field<Tfield>(num) / Tfield(den);
    }
    if (IsZeroVector(Vf))
      continue;
    FractionVectorRing<Tfield> fr = RemoveFractionVectorPlusCoeffRing(Vf);
    MyVector<Tring> const &Vr = fr.TheVect;
    if (!PositivelyProportional(to_field<Tfield>(Vr), Vf))
      ok = false;
    // The scaled vector is genuinely in the ring: converting it back to the
    // field and into the ring again is the identity.
    MyVector<Tring> Vr2 =
        UniversalVectorConversion<Tring, Tfield>(to_field<Tfield>(Vr));
    if (Vr2 != Vr)
      integral_ok = false;
  }
  check(ok, "RemoveFractionVectorPlusCoeffRing preserves the direction");
  check(integral_ok, "the scaled vector round trips through the field");
}

// The ordering and the sign determination, which the ring inherits from the
// same approximant ladder as the field.
template <typename Tfield>
static void process_ordering(int deg, int nb) {
  using Tring = typename RingPair<Tfield>::Tring;
  SmallRandom rnd(90210);
  bool cmp_ok = true, sign_ok = true, div_ok = true;
  for (int i_test = 0; i_test < nb; i_test++) {
    Tring a = RandomRing<Tring>(rnd, deg, 5);
    Tring b = RandomRing<Tring>(rnd, deg, 5);
    Tfield af = to_field<Tfield>(a);
    Tfield bf = to_field<Tfield>(b);
    if ((a < b) != (af < bf) || (a > b) != (af > bf) ||
        (a == b) != (af == bf) || (a <= b) != (af <= bf) ||
        (a >= b) != (af >= bf))
      cmp_ok = false;
    if (IsNonNegative(a) != IsNonNegative(af))
      sign_ok = false;
    if ((a > 0) != (af > 0) || (a < 0) != (af < 0))
      sign_ok = false;
    // The exact division: (a*b)/b is a, and it is the field quotient.
    if (b != Tring(0)) {
      Tring prod = a * b;
      Tring quot = prod / b;
      if (quot != a || to_field<Tfield>(quot) != to_field<Tfield>(prod) / bf)
        div_ok = false;
    }
  }
  check(cmp_ok, "the comparisons agree between the ring and the field");
  check(sign_ok, "the sign determination agrees between the ring and the field");
  check(div_ok, "the exact division agrees with the field quotient");
}

template <typename Tfield>
static void run_suite(std::string const &name, int deg, int nb) {
  std::cerr << "===== " << name << " =====\n";
  process_square<Tfield>(4, deg, nb);
  std::cerr << "STEP 1: the square matrix operations\n";
  process_canonicalization<Tfield>(4, deg, nb);
  std::cerr << "STEP 2: the canonicalizations\n";
  process_canonicalize_ordered<Tfield>(7, 4, deg, nb);
  std::cerr << "STEP 3: the ordered canonicalization and the row selection\n";
  process_subset_solver<Tfield>(6, 4, deg, nb);
  std::cerr << "STEP 4: the subset solver\n";
  process_scaling<Tfield>(5, deg, nb);
  std::cerr << "STEP 5: the scaling into the ring\n";
  process_ordering<Tfield>(deg, 10 * nb);
  std::cerr << "STEP 6: the ordering, the sign and the division\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T_rat = mpq_class;
    std::string eFile = "CI_tests/RealAlgebraicField/CubicFieldDisc_49";
    bool found = false;
    for (int level = 0; level <= 10; level++) {
      if (FILE_IsExistingFile(eFile)) {
        found = true;
        break;
      }
      eFile = "../" + eFile;
    }
    if (!found) {
      std::cerr << "Failed to find RealAlgebraicField test data after checking "
                   "paths from CI_tests/ up to 10 parent levels\n";
      throw TerminalException{1};
    }
    HelperClassRealField<T_rat> hcrf(eFile);
    insert_helper_real_algebraic_field(idx_field, hcrf);
    check(hcrf.is_monic(), "the minimal polynomial is monic");
    //
    int nb = 20;
    if (argc == 2)
      nb = ParseScalar<int>(std::string(argv[1]));
    std::cerr << "Running " << nb << " trials per section\n";
    run_suite<RealField<idx_field>>("RealField over the cubic field of "
                                    "discriminant 49, against Z[x]",
                                    hcrf.deg, nb);
    // The quadratic fields need no registration: sqrt(d) is a root of the
    // monic X^2 - d, so Z[sqrt(d)] is a ring for every d.
    run_suite<QuadField<mpq_class, 2>>("Q(sqrt(2)) against Z[sqrt(2)]", 2, nb);
    run_suite<QuadField<mpq_class, 3>>("Q(sqrt(3)) against Z[sqrt(3)]", 2, nb);
    run_suite<QuadField<mpq_class, 5>>("Q(sqrt(5)) against Z[sqrt(5)]", 2, nb);
    //
    std::cerr << "n_error=" << n_error << "\n";
    if (n_error > 0) {
      std::cerr << "Erroneous termination of Test_RingConsistency\n";
      return 1;
    }
    std::cerr << "Normal termination of Test_RingConsistency\n";
    std::cerr << "runtime = " << time << "\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
}
