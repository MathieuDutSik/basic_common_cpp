// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_TIMINGS_H_
#define SRC_BASIC_TIMINGS_H_

#include <chrono>
#include <iostream>
#include <string>
#include <thread>

std::chrono::time_point<std::chrono::system_clock>
get_cpp_time(int day, int month, int year) {
  std::tm start{};
  start.tm_sec = 0;
  start.tm_min = 0;
  start.tm_hour = 0;
  start.tm_mday = day;
  start.tm_mon = month - 1;
  start.tm_year = year - 1900;
  return std::chrono::system_clock::from_time_t(std::mktime(&start));
}

// The sleep function

void my_sleep(size_t n_millisecond) {
  std::chrono::milliseconds timespan(n_millisecond);
  std::this_thread::sleep_for(timespan);
}

// The nice time printout

std::string timeanddate() {
  std::time_t result = std::time(nullptr);
  std::string estr = std::asctime(std::localtime(&result));
  return estr.substr(0, estr.size() - 1);
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
std::string ms(SingletonTime const &s1) { return ms(s1, SingletonTime()); }

std::string s(SingletonTime const &s1, SingletonTime const &s2) {
  return std::to_string(
      std::chrono::duration_cast<std::chrono::seconds>(s2.time - s1.time)
          .count());
}
std::string s(SingletonTime const &s1) { return s(s1, SingletonTime()); }

int si(SingletonTime const &s1, SingletonTime const &s2) {
  return std::chrono::duration_cast<std::chrono::seconds>(s2.time - s1.time)
      .count();
}
int si(SingletonTime const &s1) { return si(s1, SingletonTime()); }

int msi(SingletonTime const &s1, SingletonTime const &s2) {
  return std::chrono::duration_cast<std::chrono::microseconds>(s2.time -
                                                               s1.time)
      .count();
}
int msi(SingletonTime const &s1) { return msi(s1, SingletonTime()); }

double sd(SingletonTime const &s1, SingletonTime const &s2) {
  int n_microsecs =
      std::chrono::duration_cast<std::chrono::microseconds>(s2.time - s1.time)
          .count();
  int n_in_second = 1000 * 1000;
  return double(n_microsecs) / double(n_in_second);
}
double sd(SingletonTime const &s1) { return sd(s1, SingletonTime()); }

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

// The MicrosecondTime

struct MicrosecondTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  MicrosecondTime() { time = std::chrono::system_clock::now(); }
  int eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int delta =
        std::chrono::duration_cast<std::chrono::microseconds>(timeNew - time)
            .count();
    time = timeNew;
    return delta;
  }
};

std::ostream &operator<<(std::ostream &os, MicrosecondTime &x) {
  os << x.eval();
  return os;
}

// The SecondTime

struct SecondTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  SecondTime() { time = std::chrono::system_clock::now(); }
  int eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int delta = std::chrono::duration_cast<std::chrono::seconds>(timeNew - time)
                    .count();
    time = timeNew;
    return delta;
  }
};

std::ostream &operator<<(std::ostream &os, SecondTime &x) {
  os << x.eval();
  return os;
}

// The MillisecondTime

struct MillisecondTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  MillisecondTime() { time = std::chrono::system_clock::now(); }
  int eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int delta =
        std::chrono::duration_cast<std::chrono::milliseconds>(timeNew - time)
            .count();
    time = timeNew;
    return delta;
  }
};

std::ostream &operator<<(std::ostream &os, MillisecondTime &x) {
  os << x.eval();
  return os;
}

// The NanosecondTime

struct NanosecondTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  NanosecondTime() { time = std::chrono::system_clock::now(); }
  int eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int delta =
        std::chrono::duration_cast<std::chrono::nanoseconds>(timeNew - time)
            .count();
    time = timeNew;
    return delta;
  }
};

std::ostream &operator<<(std::ostream &os, NanosecondTime &x) {
  os << x.eval();
  return os;
}

// clang-format off
#endif  // SRC_BASIC_TIMINGS_H_
// clang-format on
