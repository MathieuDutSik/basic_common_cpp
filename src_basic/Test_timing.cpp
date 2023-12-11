// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Timings.h"

int main() {

  HumanTime ht;
  SecondTime st;
  MillisecondTime millist;
  MicrosecondTime microst;

  for (int i = 0; i < 100; i++) {
    size_t delta = rand() % 10000;
    my_sleep(delta);
    std::cerr << "delta=" << delta << " dateandtime=" << timeanddate() << "\n";
    std::cerr << "  second=" << st << " millisecond=" << millist
              << " microsecond=" << microst << " human time=" << ht << "\n";
  }
}
