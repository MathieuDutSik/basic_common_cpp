// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "ImagQuadRing.h"
// clang-format on

/*
  The checks of ImagQuadRing, every one against an independent computation:
  - nearest_integers against a brute force search in a large window, on a
    grid of points of the fundamental domain, and the covering radius
    against the largest distance met and the value at the deep hole;
  - the ring axioms, the norm and the division on pseudo-random elements;
  - in the norm-Euclidean rings, the Euclidean division, the extended gcd
    (u x + v y = g, g divides x and y) and gcd(x z, y z) = z gcd(x, y)
    up to a unit on many pairs;
  - in the other rings, the failure of the division of the deep hole.
 */
template <typename T> int process(int const &D, bool const &expected_euclidean) {
  using Ring = ImagQuadRing<T>;
  using Elt = typename Ring::Elt;
  Ring ring = Ring::from_discriminant(D);
  int n_error = 0;
  auto fail = [&](std::string const &msg) -> void {
    std::cerr << "D=" << D << " ERROR: " << msg << "\n";
    n_error += 1;
  };
  auto str = [](Elt const &x) -> std::string {
    std::stringstream ss;
    ss << "(" << x.first << "," << x.second << ")";
    return ss.str();
  };
  if (ring.discriminant() != T(D)) {
    fail("wrong discriminant");
  }
  if (ring.is_norm_euclidean() != expected_euclidean) {
    fail("wrong norm-Euclidean status");
  }
  size_t n_unit = ring.units().size();
  size_t n_unit_expected = (D == -4) ? 4 : ((D == -3) ? 6 : 2);
  if (n_unit != n_unit_expected) {
    fail("wrong number of units " + std::to_string(n_unit));
  }
  // nearest_integers against the brute force.
  auto brute_nearest = [&](Elt const &x) -> std::pair<T, std::vector<Elt>> {
    std::optional<T> best;
    std::vector<Elt> l_best;
    T fa = UniversalFloorScalarInteger<T, T>(x.first);
    T fb = UniversalFloorScalarInteger<T, T>(x.second);
    for (int a = -8; a <= 8; a++) {
      for (int b = -3; b <= 3; b++) {
        Elt q{fa + T(a), fb + T(b)};
        T val = ring.norm(Ring::sub(x, q));
        if (!best || val < *best) {
          best = val;
          l_best.clear();
        }
        if (val == *best) {
          l_best.push_back(q);
        }
      }
    }
    std::sort(l_best.begin(), l_best.end());
    return {*best, l_best};
  };
  T R2 = ring.covering_radius_squared();
  T max_dist(0);
  int K = 24;
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < K; j++) {
      // Also a shift by an integer, which must not matter.
      Elt x{T(i) / T(K) + T(3), T(j) / T(K) - T(2)};
      std::vector<Elt> l_near = ring.nearest_integers(x);
      std::sort(l_near.begin(), l_near.end());
      std::pair<T, std::vector<Elt>> pair = brute_nearest(x);
      if (l_near != pair.second) {
        fail("nearest_integers differs from the brute force at " + str(x));
      }
      T dist = pair.first;
      if (dist > max_dist) {
        max_dist = dist;
      }
    }
  }
  if (max_dist > R2) {
    fail("a point is farther than the covering radius");
  }
  Elt hole = ring.deep_hole();
  if (brute_nearest(hole).first != R2 || ring.distance_squared(hole) != R2) {
    fail("the deep hole is not at the covering radius");
  }
  // Pseudo-random elements.
  uint64_t seed = 12345 + static_cast<uint64_t>(-D);
  auto next_int = [&](int const &range) -> int {
    seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<int>((seed >> 33) % (2 * range + 1)) - range;
  };
  auto random_elt = [&](int const &range) -> Elt {
    return {T(next_int(range)), T(next_int(range))};
  };
  for (int iter = 0; iter < 500; iter++) {
    Elt x = random_elt(20), y = random_elt(20), z = random_elt(20);
    if (ring.mul(x, y) != ring.mul(y, x)) {
      fail("mul is not commutative");
    }
    if (ring.mul(ring.mul(x, y), z) != ring.mul(x, ring.mul(y, z))) {
      fail("mul is not associative");
    }
    if (ring.mul(x, Ring::add(y, z)) != Ring::add(ring.mul(x, y), ring.mul(x, z))) {
      fail("mul is not distributive");
    }
    if (ring.norm(ring.mul(x, y)) != ring.norm(x) * ring.norm(y)) {
      fail("the norm is not multiplicative");
    }
    if (ring.mul(x, ring.conj(x)) != Elt{ring.norm(x), T(0)}) {
      fail("x conj(x) is not N(x)");
    }
    if (!Ring::is_zero(y) && ring.mul(ring.div(x, y), y) != x) {
      fail("(x / y) y is not x");
    }
  }
  // The Euclidean division and the gcd.
  size_t n_division_fail = 0;
  for (int iter = 0; iter < 2000; iter++) {
    Elt x = random_elt(50), y = random_elt(50), z = random_elt(10);
    if (Ring::is_zero(y) || Ring::is_zero(z)) {
      continue;
    }
    std::optional<std::pair<Elt, Elt>> opt_div = ring.euclidean_division(x, y);
    if (!opt_div) {
      n_division_fail += 1;
      continue;
    }
    Elt const &q = opt_div->first;
    Elt const &r = opt_div->second;
    if (!Ring::is_integral(q) || Ring::add(ring.mul(q, y), r) != x ||
        ring.norm(r) >= ring.norm(y)) {
      fail("wrong Euclidean division of " + str(x) + " by " + str(y));
    }
    if (!expected_euclidean) {
      continue;
    }
    std::optional<std::vector<Elt>> opt_gcd = ring.extended_gcd(x, y);
    if (!opt_gcd) {
      fail("extended_gcd failed in a norm-Euclidean ring");
      continue;
    }
    Elt const &g = (*opt_gcd)[0];
    Elt const &u = (*opt_gcd)[1];
    Elt const &v = (*opt_gcd)[2];
    if (Ring::add(ring.mul(u, x), ring.mul(v, y)) != g) {
      fail("u x + v y is not the gcd for " + str(x) + ", " + str(y));
    }
    if (!ring.divides(g, x) || !ring.divides(g, y)) {
      fail("the gcd does not divide " + str(x) + ", " + str(y));
    }
    std::optional<Elt> opt_gz = ring.gcd(ring.mul(x, z), ring.mul(y, z));
    Elt zg = ring.mul(z, g);
    if (!opt_gz || !ring.divides(*opt_gz, zg) || !ring.divides(zg, *opt_gz)) {
      fail("gcd(x z, y z) is not z gcd(x, y) up to a unit");
    }
    if (ring.canonical_associate(g) != ring.canonical_associate(ring.mul(ring.units().back(), g))) {
      fail("canonical_associate is not invariant");
    }
  }
  if (expected_euclidean && n_division_fail > 0) {
    fail("a Euclidean division failed in a norm-Euclidean ring");
  }
  // The ideals and their generators: for the class number 1 the gcd always
  // exists, divides x and y and generates the same ideal; otherwise a
  // generator found generates the ideal, and the non principal ideals
  // give none.
  bool class_number_one = (D == -3 || D == -4 || D == -7 || D == -8 || D == -11 || D == -19 ||
                           D == -43 || D == -67 || D == -163);
  size_t n_non_principal = 0;
  for (int iter = 0; iter < 300; iter++) {
    Elt x = random_elt(30), y = random_elt(30);
    if (Ring::is_zero(x) && Ring::is_zero(y)) {
      continue;
    }
    std::optional<Elt> opt_g = ring.gcd(std::vector<Elt>{x, y});
    if (!opt_g) {
      n_non_principal += 1;
      if (class_number_one) {
        fail("no gcd in a ring of class number 1 for " + str(x) + ", " + str(y));
      }
      continue;
    }
    Elt const &g = *opt_g;
    if (!ring.divides(g, x) || !ring.divides(g, y)) {
      fail("the gcd does not divide " + str(x) + ", " + str(y));
    }
    if (ring.ideal_lattice(std::vector<Elt>{g}) != ring.ideal_lattice(std::vector<Elt>{x, y})) {
      fail("the gcd does not generate the ideal of " + str(x) + ", " + str(y));
    }
  }
  if (!class_number_one && n_non_principal == 0) {
    fail("no non principal ideal met in a ring of class number > 1");
  }
  // The deep hole: x / y = hole has no remainder of norm < 1 when R2 >= 1.
  T den = UniversalFloorScalarInteger<T, T>(T(-4 * D));
  Elt y{den, T(0)};
  Elt x = ring.mul(hole, y);
  if (!Ring::is_integral(x)) {
    fail("the multiple of the deep hole should be integral");
  }
  bool division_ok = ring.euclidean_division(x, y).has_value();
  if (division_ok != expected_euclidean) {
    fail("the division of the deep hole does not match the norm-Euclidean status");
  }
  std::cerr << "D=" << D << " s=" << ring.s << " p=" << ring.p << " R^2=" << R2
            << " euclidean=" << ring.is_norm_euclidean() << " n_unit=" << n_unit
            << " n_division_fail=" << n_division_fail << " n_non_principal=" << n_non_principal
            << " n_error=" << n_error << "\n";
  return n_error;
}

int main() {
  try {
    int n_error = 0;
    for (int D : {-3, -4, -7, -8, -11}) {
      n_error += process<mpq_class>(D, true);
#ifdef ENABLE_FLINT_SUPPORT
      n_error += process<fmpq_class>(D, true);
#endif
    }
    for (int D : {-15, -19, -20, -23, -24, -43, -67, -163}) {
      n_error += process<mpq_class>(D, false);
#ifdef ENABLE_FLINT_SUPPORT
      n_error += process<fmpq_class>(D, false);
#endif
    }
    if (n_error > 0) {
      std::cerr << "Test_ImagQuadRing: " << n_error << " errors\n";
      return 1;
    }
    std::cerr << "Normal termination of Test_ImagQuadRing\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_ImagQuadRing\n";
    exit(e.eVal);
  }
}
