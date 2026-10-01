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
  // Moduli where the rejection threshold 2^64 mod n is large.
  for (uint64_t n : {uint64_t(3), (uint64_t(1) << 63) + 1, UINT64_MAX}) {
    for (int i = 0; i < 1000; i++)
      check(random_index(n) < n, "random_index(n) < n");
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
    test_reseed();
    test_ranges();
    test_threads();
    std::cerr << "Test_Random: all the checks passed\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
