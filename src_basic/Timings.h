#ifndef INCLUDE_TIMINGS_H_
#define INCLUDE_TIMINGS_H_

#include <chrono>

struct SingletonTime {
  std::chrono::time_point<std::chrono::system_clock> time;
};


std::string ms(SingletonTime const& ms1, SingletonTime const& ms2)
{
  return std::to_string(std::chrono::duration_cast<std::chrono::microseconds>(ms2.time - ms1.time).count());
}







#endif
