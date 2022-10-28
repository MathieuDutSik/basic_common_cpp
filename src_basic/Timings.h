// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_TIMINGS_H_
#define SRC_BASIC_TIMINGS_H_

#include <chrono>
#include <iostream>
#include <string>


std::chrono::time_point<std::chrono::system_clock> get_cpp_time(int day, int month, int year) {
  std::tm start{};
  start.tm_sec = 0;
  start.tm_min = 0;
  start.tm_hour = 0;
  start.tm_mday = day;
  start.tm_mon = month - 1;
  start.tm_year = year - 1900;
  return std::chrono::system_clock::from_time_t(std::mktime(&start));
}

// The SingletonTime

struct SingletonTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  SingletonTime() { time = std::chrono::system_clock::now(); }
};

std::string ms(SingletonTime const &s1, SingletonTime const &s2) {
  return std::to_string(
      std::chrono::duration_cast<std::chrono::microseconds>(s2.time - s1.time)
          .count());
}
std::string ms(SingletonTime const &s1) {
  return ms(s1, SingletonTime());
}

std::string s(SingletonTime const &s1, SingletonTime const &s2) {
  return std::to_string(
      std::chrono::duration_cast<std::chrono::seconds>(s2.time - s1.time)
          .count());
}
std::string s(SingletonTime const &s1) {
  return s(s1, SingletonTime());
}

int si(SingletonTime const &s1, SingletonTime const &s2) {
  return std::chrono::duration_cast<std::chrono::seconds>(s2.time - s1.time)
      .count();
}
int si(SingletonTime const &s1) {
  return si(s1, SingletonTime());
}

int msi(SingletonTime const &s1, SingletonTime const &s2) {
  return std::chrono::duration_cast<std::chrono::microseconds>(s2.time -
                                                               s1.time)
      .count();
}
int msi(SingletonTime const &s1) {
  return msi(s1, SingletonTime());
}

double sd(SingletonTime const &s1, SingletonTime const &s2) {
  int n_microsecs = std::chrono::duration_cast<std::chrono::microseconds>(s2.time - s1.time).count();
  int n_in_second = 1000 * 1000;
  return double(n_microsecs) / double(n_in_second);
}
double sd(SingletonTime const &s1) {
  return sd(s1, SingletonTime());
}

void runtime(SingletonTime const &start) {
  SingletonTime end;
  std::cerr << "runtime = " << s(start, end) << "\n";
}

void timing(std::string const &strTime, SingletonTime const &s1,
            SingletonTime const &s2) {
  std::cerr << "Timing |" << strTime << "|=" << s(s1, s2) << "\n";
}


struct SingletonTimeInc {
  SingletonTime time;
  SingletonTimeInc() : time() {}
  int get_ms() {
    SingletonTime timeNext;
    int val = msi(time, timeNext);
    time = timeNext;
    return val;
  }
};

// The MsTime

struct MilisecondTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  MilisecondTime() { time = std::chrono::system_clock::now(); }
};

std::ostream& operator<<(std::ostream& os, MilisecondTime & x) {
  std::chrono::time_point<std::chrono::system_clock> timeNew = std::chrono::system_clock::now();
  int n_ms = std::chrono::duration_cast<std::chrono::microseconds>(timeNew - x.time).count();
  x.time = timeNew;
  os << n_ms;
  return os;
}






// clang-format off
#endif  // SRC_BASIC_TIMINGS_H_
// clang-format on
