// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_THREADING_H_
#define SRC_BASIC_BASIC_THREADING_H_

// Thread-dependent helpers split out from Temp_common.h and Timings.h so the
// rest of src_basic compiles on platforms where std::thread is unavailable
// (e.g. WebAssembly without pthread support).

#include <chrono>
#include <thread>

void my_sleep(size_t n_millisecond) {
  std::chrono::milliseconds timespan(n_millisecond);
  std::this_thread::sleep_for(timespan);
}

// clang-format off
#endif  // SRC_BASIC_BASIC_THREADING_H_
// clang-format on
