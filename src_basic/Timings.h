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
  return static_cast<double>(n_microsecs) / static_cast<double>(n_in_second);
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

// The TimeEval

template<typename T>
struct TimeEval {
  std::chrono::time_point<std::chrono::system_clock> time;
  TimeEval() { time = std::chrono::system_clock::now(); }
  int eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int delta = std::chrono::duration_cast<T>(timeNew - time).count();
    time = timeNew;
    return delta;
  }
};

template<typename T>
std::ostream &operator<<(std::ostream &os, TimeEval<T> &x) {
  os << x.eval();
  return os;
}

// The instantiations

using SecondTime = TimeEval<std::chrono::seconds>;
using MillisecondTime = TimeEval<std::chrono::milliseconds>;
using MicrosecondTime = TimeEval<std::chrono::microseconds>;
using NanosecondTime = TimeEval<std::chrono::nanoseconds>;

// The HumanTime

struct HumanTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  HumanTime() { time = std::chrono::system_clock::now(); }
  std::string eval() {
    std::chrono::time_point<std::chrono::system_clock> timeNew =
        std::chrono::system_clock::now();
    int64_t delta = std::chrono::duration_cast<std::chrono::nanoseconds>(timeNew - time).count();
    time = timeNew;
    int64_t one_microsecond = 1000;
    int64_t one_millisecond = 1000 * one_microsecond;
    int64_t one_second = 1000 * one_millisecond;
    int64_t one_minute = 60 * one_second;
    int64_t one_hour = 60 * one_minute;
    int64_t one_day = 24 * one_hour;
    if (delta < one_microsecond) {
      std::string reply = std::to_string(delta) + "ns";
      return reply;
    }
    if (delta < one_millisecond) {
      int64_t res = delta / one_microsecond;
      std::string reply = std::to_string(res) + "\u039Cs";
      int64_t res2 = delta - one_microsecond * res;
      if (res2 > 0) {
        reply += " " + std::to_string(res2) + "ns";
      }
      return reply;
    }
    if (delta < one_second) {
      int64_t res = delta / one_millisecond;
      std::string reply = std::to_string(res) + "ms";
      int64_t res2 = (delta - one_millisecond * res) / one_microsecond;
      if (res2 > 0) {
        reply += " " + std::to_string(res2) + "\u039Cs";
      }
      return reply;
    }
    if (delta < one_minute) {
      int64_t res = delta / one_second;
      std::string reply = std::to_string(res) + "s";
      int64_t res2 = (delta - one_second * res) / one_millisecond;
      if (res2 > 0) {
        reply += " " + std::to_string(res2) + "ms";
      }
      return reply;
    }
    if (delta < one_day) {
      int64_t res = delta / one_hour;
      std::string reply = std::to_string(res) + "h";
      int64_t res2 = (delta - one_day * res) / one_minute;
      if (res2 > 0) {
        reply += " " + std::to_string(res2) + "min";
      }
      return reply;
    }
    int64_t res = delta / one_day;
    std::string reply = std::to_string(res) + "day";
    int64_t res2 = (delta - one_day * res) / one_hour;
    if (res2 > 0) {
      reply += " " + std::to_string(res2) + "h";
    }
    return reply;
  }
};

std::ostream &operator<<(std::ostream &os, HumanTime &x) {
  os << x.eval();
  return os;
}

// clang-format off
#endif  // SRC_BASIC_TIMINGS_H_
// clang-format on
