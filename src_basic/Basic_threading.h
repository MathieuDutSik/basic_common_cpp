// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_THREADING_H_
#define SRC_BASIC_BASIC_THREADING_H_

// Thread-dependent helpers split out from Temp_common.h and Timings.h so the
// rest of src_basic compiles on platforms where std::thread is unavailable
// (e.g. WebAssembly without pthread support).

#include "Temp_common.h"
#include <chrono>
#include <thread>

unsigned get_random_pid_seed() {
  // There seems to be no way of converting std::thread::id to size_t
  // even though the pid is going to be a normal integer. So, instead
  // we use the hash.
  std::thread::id this_id = std::this_thread::get_id();
  size_t hash = std::hash<std::thread::id>()(this_id);
  return static_cast<unsigned>(hash);
}

unsigned get_random_seed() {
  unsigned seed1 = get_random_time_seed();
  unsigned seed2 = get_random_pid_seed();
  return seed1 + seed2;
}

void srand_random_set() {
  unsigned val = get_random_seed();
  srand(val);
}

void my_sleep(size_t n_millisecond) {
  std::chrono::milliseconds timespan(n_millisecond);
  std::this_thread::sleep_for(timespan);
}

// clang-format off
#endif  // SRC_BASIC_BASIC_THREADING_H_
// clang-format on
