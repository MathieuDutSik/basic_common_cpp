// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_RANDOM_H_
#define SRC_BASIC_BASIC_RANDOM_H_

#ifndef WASM_PLATFORM
#include <thread>
#endif

unsigned get_random_time_seed() {
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

unsigned get_random_pid_seed() {
  // There seems to be no way of converting std::thread::id to size_t
  // even though the pid is going to be a normal integer. So, instead
  // we use the hash.
  std::thread::id this_id = std::this_thread::get_id();
  size_t hash = std::hash<std::thread::id>()(this_id);
  return static_cast<unsigned>(hash);
}

#endif

unsigned get_random_seed() {
  unsigned seed = get_random_time_seed();
#ifndef WASM_PLATFORM
  seed += get_random_pid_seed();
#endif
  return seed;
}

void srand_random_set() {
  unsigned val = get_random_seed();
  srand(val);
}

// clang-format off
#endif  // SRC_BASIC_BASIC_RANDOM_H_
// clang-format on
