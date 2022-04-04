#ifndef INCLUDE_TIMINGS_H_
#define INCLUDE_TIMINGS_H_

#include <chrono>

struct SingletonTime {
  std::chrono::time_point<std::chrono::system_clock> time;
  SingletonTime() {
    time = std::chrono::system_clock::now();
  }
};


std::string ms(SingletonTime const& s1, SingletonTime const& s2)
{
  return std::to_string(std::chrono::duration_cast<std::chrono::microseconds>(s2.time - s1.time).count());
}

std::string s(SingletonTime const& s1, SingletonTime const& s2)
{
  return std::to_string(std::chrono::duration_cast<std::chrono::seconds>(s2.time - s1.time).count());
}

int si(SingletonTime const& s1, SingletonTime const& s2)
{
  return std::chrono::duration_cast<std::chrono::microseconds>(s2.time - s1.time).count();
}

void runtime(SingletonTime const& start)
{
  SingletonTime end;
  std::cerr << "runtime = " << s(start,end) << "\n";
}





#endif