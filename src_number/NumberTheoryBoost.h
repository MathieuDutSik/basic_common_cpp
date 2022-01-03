#ifndef INCLUDE_NUMBER_THEORY_BOOST
#define INCLUDE_NUMBER_THEORY_BOOST



#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/eigen.hpp>

#include "NumberTheory.h"





template <>
struct is_euclidean_domain<boost::multiprecision::cpp_int> {
  static const bool value = true;
};



template<>
struct is_exact_arithmetic<boost::multiprecision::cpp_int> {
  static const bool value = true;
};


inline boost::multiprecision::cpp_int CanonicalizationUnit(boost::multiprecision::cpp_int const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}




inline boost::multiprecision::cpp_int ResInt(boost::multiprecision::cpp_int const& a, boost::multiprecision::cpp_int const& b)
{
  using T = boost::multiprecision::cpp_int;
  T q = a / b;
  T res = a - q * b;
  return res;
}


inline boost::multiprecision::cpp_int QuoInt(boost::multiprecision::cpp_int const& a, boost::multiprecision::cpp_int const& b)
{
  using T = boost::multiprecision::cpp_int;
  T q = a / b;
  return q;
}



inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, double & a2)
{
  std::cerr << "Missing code, write here\n";
  throw TerminalException{1};
}


inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, int & a2)
{
  a2 = a1.template convert_to<int>();
}

inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, long & a2)
{
  a2 = a1.template convert_to<long>();
}

inline void TYPE_CONVERSION(int const& a1, boost::multiprecision::cpp_int & a2)
{
  a2 = a1;
}







#endif
