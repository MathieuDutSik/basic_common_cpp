#ifndef INCLUDE_TYPE_CONVERSION_FINAL
#define INCLUDE_TYPE_CONVERSION_FINAL

// This is a block of type conversion after the types have been defined
// It can be coded only if the preceding include statements have all been
// clarified. Hence the moniker "final".

template<typename Ti, typename To>
inline void TYPE_CONVERSION_STRING(Ti const& a1, To & a2)
{
  std::stringstream s;
  s << a1;
  std::string converted(s.str());
  //
  std::istringstream is(converted);
  is >> a2;
}




# if defined INCLUDE_NUMBER_THEORY_GMP && defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT

// Nothing clever for this combination
inline void TYPE_CONVERSION(mpq_class const& a1, boost::multiprecision::cpp_rational & a2)
{
  TYPE_CONVERSION_STRING(a1, a2);
}

inline void TYPE_CONVERSION(mpq_class const& a1, boost::multiprecision::cpp_int & a2)
{
  if (!IsInteger(a1)) {
    std::string str_ret = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str_ret};
  }
  TYPE_CONVERSION_STRING(a1, a2);
}

inline void TYPE_CONVERSION(mpz_class const& a1, boost::multiprecision::cpp_rational & a2)
{
  TYPE_CONVERSION_STRING(a1, a2);
}

inline void TYPE_CONVERSION(mpz_class const& a1, boost::multiprecision::cpp_int & a2)
{
  TYPE_CONVERSION_STRING(a1, a2);
}

// Now the other direction

inline void TYPE_CONVERSION(boost::multiprecision::cpp_rational const& a1, mpq_class const& a2)
{
  TYPE_CONVERSION_STRING(a1, a2);
}

inline void TYPE_CONVERSION(boost::multiprecision::cpp_rational const& a1, mpz_class const& a2)
{
  if (!IsInteger(a1)) {
    std::string str_ret = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str_ret};
  }
  TYPE_CONVERSION_STRING(a1, a2);
}

inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, mpq_class const& a2)
{
  TYPE_CONVERSION_STRING(a1, a2);
}

inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, mpz_class const& a2)
{
  TYPE_CONVERSION_STRING(a1, a2);
}




# endif





#endif
