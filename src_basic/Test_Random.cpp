// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Temp_common.h"
// clang-format on

// The generator of Basic_random.h gives the same numbers on every platform:
// the reference values below were produced with libc++ and libstdc++ alike,
// and a run on any system has to reproduce them.

void check(bool test, std::string const &context) {
  if (!test) {
    std::cerr << "Test_Random: failure in " << context << "\n";
    throw TerminalException{1};
  }
}

void test_standard_engine() {
  // The C++ standard fixes the 10000th output of a default constructed
  // std::mt19937_64.
  std::mt19937_64 engine;
  engine.discard(9999);
  check(engine() == 9981545732273789042ULL, "std::mt19937_64 10000th value");
}

void test_reference_values() {
  set_random_seed(12345);
  std::vector<uint64_t> ref_u64{970638883550249018ULL, 3765130718776470004ULL,
                                10451094180116843272ULL,
                                4752344318214878128ULL};
  for (uint64_t ref : ref_u64)
    check(random_u64() == ref, "random_u64 reference values");
  std::vector<uint64_t> ref_index{500, 321, 590, 648, 853,
                                  879, 717, 401, 263, 485};
  for (uint64_t ref : ref_index)
    check(random_index(1000) == ref, "random_index reference values");
  std::vector<int> V{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  random_shuffle_vector(V);
  check(V == std::vector<int>{4, 7, 6, 1, 3, 2, 8, 9, 0, 5},
        "random_shuffle_vector reference permutation");
}

void test_random_int() {
  set_random_seed(777);
  for (int ref : {1282271111, 287147645, 530360574, 1674523095})
    check(random_int() == ref, "random_int() reference values");
  for (int ref : {-5, -4, 8, -5, 4, -8, 9, -4, 9, 0})
    check(random_int(-10, 10) == ref, "random_int(-10, 10) reference values");
  for (int ref : {2002814411, 950268971, 1571890902})
    check(random_int(INT32_MIN, INT32_MAX) == ref,
          "random_int(INT32_MIN, INT32_MAX) reference values");
  // The bounds are reached and never passed.
  set_random_seed(31);
  bool seen_lo = false, seen_hi = false;
  for (int i = 0; i < 10000; i++) {
    int x = random_int(-3, 3);
    check(-3 <= x && x <= 3, "random_int(-3, 3) range");
    if (x == -3)
      seen_lo = true;
    if (x == 3)
      seen_hi = true;
  }
  check(seen_lo && seen_hi, "random_int(-3, 3) reaches its bounds");
  for (int i = 0; i < 1000; i++) {
    check(random_int(5, 5) == 5, "random_int(5, 5)");
    int x = random_int();
    check(0 <= x && x <= 2147483647, "random_int() range");
  }
  // The same draws as random_index over the same width.
  set_random_seed(4);
  std::vector<int> V1;
  for (int i = 0; i < 20; i++)
    V1.push_back(random_int(0, 999));
  set_random_seed(4);
  std::vector<int> V2;
  for (int i = 0; i < 20; i++)
    V2.push_back(static_cast<int>(random_index(1000)));
  check(V1 == V2, "random_int(0, n-1) agrees with random_index(n)");
  // random_int64 draws as random_int over a range of int, and covers ranges
  // beyond int, the full one included.
  set_random_seed(5);
  std::vector<int64_t> W1;
  for (int i = 0; i < 20; i++)
    W1.push_back(random_int(-1000, 1000));
  set_random_seed(5);
  std::vector<int64_t> W2;
  for (int i = 0; i < 20; i++)
    W2.push_back(random_int64(-1000, 1000));
  check(W1 == W2, "random_int64 agrees with random_int");
  int64_t big = 1000000000000000;
  for (int i = 0; i < 1000; i++) {
    int64_t x = random_int64(-big, big);
    check(-big <= x && x <= big, "random_int64(-10^15, 10^15) range");
    random_int64(INT64_MIN, INT64_MAX);
    check(random_int64(INT64_MAX, INT64_MAX) == INT64_MAX,
          "random_int64(INT64_MAX, INT64_MAX)");
  }
}

void test_local_engine() {
  // A stream of its own: the same values for the same engine seed, and no
  // effect on the draws of the thread.
  std::mt19937_64 eng(2026);
  for (int ref : {-20, 92, -52, -23, 70, 50})
    check(random_int(eng, -100, 100) == ref, "random_int(engine) references");
  for (double ref :
       {-0.54591112459201607, 2.6242259009096074, 1.8166949843926745})
    check(random_real(eng, -2.0, 3.0) == ref, "random_real references");
  // std::log and std::cos may round differently in the last bit.
  for (double ref : {-0.3058648000446515, -0.26468486319168938,
                     0.88951356600948661, -0.84270929941688111})
    check(std::abs(random_normal(eng) - ref) < 1e-12,
          "random_normal references");
  set_random_seed(12345);
  std::mt19937_64 eng2(1);
  for (int i = 0; i < 100; i++)
    random_int(eng2, 0, 9);
  check(random_u64() == 970638883550249018ULL,
        "a local engine leaves the thread draws alone");
  // Ranges and moments.
  set_random_seed(6);
  int n_draw = 400000;
  double sum = 0, sum2 = 0;
  for (int i = 0; i < n_draw; i++) {
    double x = random_real(-2.0, 3.0);
    check(-2.0 <= x && x < 3.0, "random_real range");
    double z = random_normal();
    check(std::isfinite(z), "random_normal is finite");
    sum += z;
    sum2 += z * z;
  }
  double mean = sum / n_draw;
  double var = sum2 / n_draw - mean * mean;
  check(std::abs(mean) < 0.01 && std::abs(var - 1) < 0.01,
        "random_normal mean 0 and variance 1");
}

void test_reseed() {
  // The same seed gives the same sequence, a different seed another one.
  auto draw = [](uint64_t seed) -> std::vector<uint64_t> {
    set_random_seed(seed);
    std::vector<uint64_t> V;
    for (int i = 0; i < 20; i++)
      V.push_back(random_u64());
    return V;
  };
  check(draw(7) == draw(7), "reseed reproduces the sequence");
  check(draw(7) != draw(8), "different seeds differ");
}

void test_ranges() {
  set_random_seed(2024);
  for (int i = 0; i < 1000; i++)
    check(random_index(1) == 0, "random_index(1)");
  // Moduli where the rejection threshold 2^64 mod n is large, on the 64-bit
  // kernel under random_index and random_int.
  for (uint64_t n : {uint64_t(3), (uint64_t(1) << 63) + 1, UINT64_MAX}) {
    for (int i = 0; i < 1000; i++)
      check(basic_random_detail::random_below(n) < n, "random_below(n) < n");
  }
  for (int i = 0; i < 1000; i++) {
    double x = random_unit();
    check(0 <= x && x < 1, "random_unit in [0, 1)");
  }
  // A die: each face within 2% of its expected count. The seed is fixed, so
  // the test is not flaky.
  int n_draw = 600000;
  std::vector<int> counts(6, 0);
  for (int i = 0; i < n_draw; i++)
    counts[random_index(6)]++;
  for (int c : counts)
    check(98000 < c && c < 102000, "random_index(6) uniformity");
  int n_true = 0;
  for (int i = 0; i < n_draw; i++)
    if (random_bool())
      n_true++;
  check(294000 < n_true && n_true < 306000, "random_bool balance");
  // The shuffle is a permutation.
  std::vector<int> W(100);
  for (int i = 0; i < 100; i++)
    W[i] = i;
  random_shuffle_vector(W);
  std::vector<int> Wsort = W;
  std::sort(Wsort.begin(), Wsort.end());
  for (int i = 0; i < 100; i++)
    check(Wsort[i] == i, "random_shuffle_vector is a permutation");
}

void test_threads() {
  // A new thread draws from its own engine, seeded from the main seed and its
  // thread index, and the reseed reaches it.
  set_random_seed(99);
  uint64_t main_first = random_u64();
  uint64_t thread_first = 0;
  std::thread thr([&]() -> void { thread_first = random_u64(); });
  thr.join();
  check(main_first != thread_first, "threads draw different sequences");
  set_random_seed(99);
  check(random_u64() == main_first, "reseed reaches the main thread");
}

int main() {
  try {
    test_standard_engine();
    test_reference_values();
    test_random_int();
    test_local_engine();
    test_reseed();
    test_ranges();
    test_threads();
    std::cerr << "Test_Random: all the checks passed\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
