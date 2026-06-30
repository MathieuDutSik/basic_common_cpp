// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_RANDOM_H_
#define SRC_BASIC_BASIC_RANDOM_H_

#include <cstdlib>
#ifdef _WIN32
#include <random>
#endif
#ifndef WASM_PLATFORM
#include <thread>
#endif

inline unsigned get_random_time_seed() {
#ifdef USE_NANOSECOND_RAND
  std::timespec ts;
  std::timespec_get(&ts, TIME_UTC);
  unsigned val = ts.tv_nsec;
#else
  unsigned val = time(nullptr);
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
