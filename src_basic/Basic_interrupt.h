// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_INTERRUPT_H_
#define SRC_BASIC_BASIC_INTERRUPT_H_

// signal for early termination
#include <signal.h>

std::atomic<bool> ExitEvent;

void signal_callback_handler(int signum) {
  std::cout << "Caught signal " << signum << "\n";
  std::cout << "We are going to exit hopefully\n";
  ExitEvent = true;
}

void enable_sig_int() {
  ExitEvent = false;
  signal(SIGINT, signal_callback_handler);
}

// clang-format off
#endif  // SRC_BASIC_BASIC_INTERRUPT_H_
// clang-format on
