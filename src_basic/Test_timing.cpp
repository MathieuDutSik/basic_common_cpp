// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Timings.h"

int main() {
  std::chrono::time_point<std::chrono::system_clock> time0 =
      get_cpp_time(1, 1, 1974);
  std::chrono::time_point<std::chrono::system_clock> time1 =
      get_cpp_time(2, 1, 2022);
  size_t dur =
      std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0)
          .count();
  std::cerr << "dur=" << dur << "\n";

  SecondTime st;
  MillisecondTime millist;
  MicrosecondTime microst;
  my_sleep(1000);
  std::cerr << "duration 1 s=" << st << " millisecond=" << millist
            << " microsecond=" << microst << " dateandtime=" << timeanddate()
            << "\n";
  my_sleep(1000);
  std::cerr << "duration 1 s=" << st << " millisecond=" << millist
            << " microsecond=" << microst << " dateandtime=" << timeanddate()
            << "\n";
}
