// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_TYPECONVERSIONFINAL_H_
#define SRC_NUMBER_TYPECONVERSIONFINAL_H_

#include <string>

// This is a block of type conversion after the types have been defined
// It can be coded only if the preceding include statements have all been
// clarified. Hence the moniker "final".

template <typename Ti, typename To>
inline void TYPE_CONVERSION_STRING(Ti const &a1, To &a2) {
  std::stringstream s;
  s << a1;
  std::string converted(s.str());
  //
  std::istringstream is(converted);
  is >> a2;
}

#if defined INCLUDE_NUMBER_THEORY_GMP &&                                       \
    defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT

// Nothing clever for this combination
inline void TYPE_CONVERSION(stc<mpq_class> const &a1,
                            boost::multiprecision::cpp_rational &a2) {
  TYPE_CONVERSION_STRING(a1.val, a2);
}

template<typename T1, typename T2>
void TYPE_CONVERSION_IsInteger(stc<T1> const& a1, T2 & a2) {
  if (!IsInteger(a1.val)) {
    std::string str_ret = "a1=" + std::to_string(a1.val) + " is not an integer";
    throw ConversionException{str_ret};
  }
  TYPE_CONVERSION_STRING(a1.val, a2);
}


inline void TYPE_CONVERSION(stc<mpq_class> const &a1,
                            boost::multiprecision::cpp_int &a2) {
  TYPE_CONVERSION_IsInteger(a1, a2);
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1,
                            boost::multiprecision::cpp_rational &a2) {
  TYPE_CONVERSION_STRING(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1,
                            boost::multiprecision::cpp_int &a2) {
  TYPE_CONVERSION_STRING(a1.val, a2);
}

// Now the other direction

inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_rational> const &a1,
                            mpq_class &a2) {
  TYPE_CONVERSION_STRING(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_rational> const &a1,
                            mpz_class &a2) {
  TYPE_CONVERSION_IsInteger(a1, a2);
}

inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_int> const &a1,
                            mpq_class &a2) {
  TYPE_CONVERSION_STRING(a1.val, a2);
}

inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_int> const &a1,
                            mpz_class &a2) {
  TYPE_CONVERSION_STRING(a1.val, a2);
}

#endif

// clang-format off
#endif  // SRC_NUMBER_TYPECONVERSIONFINAL_H_
// clang-format on
