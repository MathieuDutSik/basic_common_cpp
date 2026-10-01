// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_RANDOM_H_
#define SRC_BASIC_BASIC_RANDOM_H_

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <utility>
#include <vector>
#ifndef WASM_PLATFORM
#include <thread>
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_BASIC_RANDOM
#endif

#ifdef SANITY_CHECK_BASIC_RANDOM
#include "ExceptionsFunc.h"
#include <iostream>
#endif

inline unsigned get_random_time_seed() {
#ifdef USE_NANOSECOND_RAND
  std::timespec ts;
  std::timespec_get(&ts, TIME_UTC);
  // A seed: the truncation to unsigned is intended.
  unsigned val = static_cast<unsigned>(ts.tv_nsec);
#else
  unsigned val = static_cast<unsigned>(time(nullptr));
#endif
  return val;
}

#ifndef WASM_PLATFORM

inline unsigned get_random_pid_seed() {
  // There seems to be no way of converting std::thread::id to size_t
  // even though the pid is going to be a normal integer. So, instead
  // we use the hash.
  std::thread::id this_id = std::this_thread::get_id();
  size_t hash = std::hash<std::thread::id>()(this_id);
  return static_cast<unsigned>(hash);
}

#endif

inline unsigned get_random_seed() {
  unsigned seed = get_random_time_seed();
#ifndef WASM_PLATFORM
  seed += get_random_pid_seed();
#endif
  return seed;
}

inline void srand_random_set() {
  unsigned val = get_random_seed();
  srand(val);
}

/*
  The random numbers of the code, the same on every platform.

  rand() and random() are two generators that POSIX keeps apart and whose
  sequences it does not specify: glibc shares one state between them, macOS
  does not, so srand() does not seed random() there and a seed gives
  different numbers on different systems. Instead every draw goes to a
  std::mt19937_64, whose output the C++ standard fixes exactly, and the
  ranges are formed here rather than by the std distributions, whose
  algorithms are left to each library. A seed therefore gives the same
  numbers with libstdc++, libc++ and MSVC.

  Each thread has its own engine, seeded from the main seed and the index of
  the thread (the order in which threads first draw), so threads draw
  without sharing state. The default seed is fixed: a run that does not
  call set_random_seed is deterministic, as random() is when unseeded.
  set_random_seed(seed) reseeds every thread at its next draw, and
  set_random_seed_nondeterministic() seeds from the time and the thread id.
 */
namespace basic_random_detail {

inline constexpr uint64_t default_seed = 0x5eed5eed5eed5eedULL;

inline std::atomic<uint64_t> &main_seed() {
  static std::atomic<uint64_t> seed{default_seed};
  return seed;
}

// Incremented by every set_random_seed, so that each thread engine sees that
// it has to be reseeded.
inline std::atomic<uint64_t> &seed_generation() {
  static std::atomic<uint64_t> generation{0};
  return generation;
}

inline std::atomic<uint64_t> &thread_counter() {
  static std::atomic<uint64_t> counter{0};
  return counter;
}

// The splitmix64 finalizer: spreads the bits of the seed and the thread index
// over the 64-bit engine seed.
inline uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

struct ThreadEngine {
  std::mt19937_64 engine;
  uint64_t index;
  // A generation that never occurs, so that the first access seeds.
  uint64_t generation = UINT64_MAX;
};

inline std::mt19937_64 &thread_engine() {
  thread_local ThreadEngine te{std::mt19937_64(),
                               thread_counter().fetch_add(1)};
  uint64_t generation = seed_generation().load(std::memory_order_acquire);
  if (te.generation != generation) {
    uint64_t seed = main_seed().load(std::memory_order_relaxed);
    te.engine.seed(splitmix64(seed ^ splitmix64(te.index)));
    te.generation = generation;
  }
  return te.engine;
}

} // namespace basic_random_detail

inline void set_random_seed(uint64_t seed) {
  basic_random_detail::main_seed().store(seed, std::memory_order_relaxed);
  basic_random_detail::seed_generation().fetch_add(1,
                                                   std::memory_order_release);
  // The calling thread takes its index now, so the thread that seeds is
  // thread 0 when it seeds before any other draw.
  basic_random_detail::thread_engine();
}

inline void set_random_seed_nondeterministic() {
  set_random_seed(get_random_seed());
}

// The engine of the calling thread, for the code that needs an engine. The
// std distributions applied to it are not portable: prefer the functions
// below.
inline std::mt19937_64 &get_random_engine() {
  return basic_random_detail::thread_engine();
}

inline uint64_t random_u64() { return basic_random_detail::thread_engine()(); }

namespace basic_random_detail {

// Uniform in [0, n), n > 0, without bias: the draws below 2^64 mod n are
// rejected, so that the accepted ones cover every residue equally often (the
// method of OpenBSD's arc4random_uniform, in 64-bit arithmetic only). All the
// bounded draws below go through it, with their conversions written once.
inline uint64_t random_below(uint64_t n) {
  uint64_t threshold = (0 - n) % n;
  while (true) {
    uint64_t x = random_u64();
    if (x >= threshold)
      return x % n;
  }
}

} // namespace basic_random_detail

static_assert(sizeof(int) >= 4, "random_int needs a 32-bit int");

// Uniform in the closed range [lo, hi], lo <= hi. The width hi - lo + 1 is at
// most 2^32 and is formed in 64 bits, so that the full int range is valid.
inline int random_int(int lo, int hi) {
#ifdef SANITY_CHECK_BASIC_RANDOM
  if (lo > hi) {
    std::cerr << "RANDOM: random_int with lo=" << lo << " > hi=" << hi << "\n";
    throw TerminalException{1};
  }
#endif
  int64_t lo64 = lo;
  int64_t hi64 = hi;
  uint64_t width = static_cast<uint64_t>(hi64 - lo64) + 1;
  int64_t offset =
      static_cast<int64_t>(basic_random_detail::random_below(width));
  return static_cast<int>(lo64 + offset);
}

// Uniform in [0, 2^31 - 1], the range of random(), on every platform.
inline int random_int() { return random_int(0, 2147483647); }

// Uniform in [0, n), n > 0: an index into a container of size n.
inline size_t random_index(size_t n) {
#ifdef SANITY_CHECK_BASIC_RANDOM
  if (n == 0) {
    std::cerr << "RANDOM: random_index with n=0\n";
    throw TerminalException{1};
  }
#endif
  return static_cast<size_t>(
      basic_random_detail::random_below(static_cast<uint64_t>(n)));
}

inline bool random_bool() { return (random_u64() >> 63) != 0; }

// Uniform in [0, 1), from the top 53 bits.
inline double random_unit() {
  return static_cast<double>(random_u64() >> 11) * 0x1.0p-53;
}

// Fisher-Yates, with random_index: the same permutation on every platform,
// unlike std::shuffle whose algorithm is left to the library.
template <typename T> void random_shuffle_vector(std::vector<T> &V) {
  size_t len = V.size();
  for (size_t i = len; i > 1; i--) {
    size_t j = random_index(i);
    using std::swap;
    swap(V[i - 1], V[j]);
  }
}

#ifdef _WIN32
// POSIX random() is not provided by the MinGW / MSVC C runtimes. Supply a
// portable shim with the same signature (long in [0, 2^31 - 1]) so call sites
// can use random() uniformly across platforms.
inline long random() {
  thread_local std::mt19937 gen{get_random_seed()};
  return static_cast<long>(gen() & 0x7FFFFFFFL);
}
#endif

// clang-format off
#endif  // SRC_BASIC_BASIC_RANDOM_H_
// clang-format on
