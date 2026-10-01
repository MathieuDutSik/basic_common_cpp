// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYBOOSTGMPINT_H_
#define SRC_NUMBER_NUMBERTHEORYBOOSTGMPINT_H_
#define INCLUDE_NUMBER_THEORY_BOOST_GMP_INT
// clang-format off
#include "BasicNumberTypes.h"
#include "NumberTheoryBoostFormat.h"
#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
#include "TypeConversion.h"
#include "boost_serialization.h"
#include <boost/multiprecision/gmp.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_free.hpp>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
// clang-format on

// The FMA form: mpz_int fuses `acc += a*b` via its expression templates, so
// the direct form is best (measured).
// clang-format off
#define NT_BOOST_INT boost::multiprecision::mpz_int
#define NT_BOOST_RAT boost::multiprecision::mpq_rational
#define NT_BOOST_IS_INT is_boost_mpz_int
#define NT_BOOST_IS_RAT is_boost_mpq_rational
#define NT_BOOST_INT_FMA_PREFERED true
#define NT_BOOST_INT_TO_SMALL mpz_int_to_small_integer
#define NT_BOOST_CEIL_RAT Ceil_mpq_rational
#define NT_BOOST_FLOOR_RAT Floor_mpq_rational
#include "NumberTheoryBoost_impl.h"
// clang-format on

// to_string


// boost serialization

namespace boost::serialization {

// boost::multiprecision::mpq_rational

template <class Archive>
inline void load(Archive &ar, boost::multiprecision::mpq_rational &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("mpq_rational", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, boost::multiprecision::mpq_rational const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::ostringstream os;
  os << val;
  std::string str = os.str();
  ar &make_nvp("mpq_rational", str);
}

template <class Archive>
inline void serialize(Archive &ar, boost::multiprecision::mpq_rational &val,
                      [[maybe_unused]] const unsigned int version) {
  split_free(ar, val, version);
}

// boost::multiprecision::mpz_int

template <class Archive>
inline void load(Archive &ar, boost::multiprecision::mpz_int &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("mpz_int", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, boost::multiprecision::mpz_int const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::ostringstream os;
  os << val;
  std::string str = os.str();
  ar &make_nvp("mpz_int", str);
}

template <class Archive>
inline void serialize(Archive &ar, boost::multiprecision::mpz_int &val,
                      const unsigned int version) {
  split_free(ar, val, version);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// See T_frexp in TypeConversion.h. The backend is a GMP integer.
inline double T_frexp(boost::multiprecision::mpz_int const &x, long &e) {
  return mpz_get_d_2exp(&e, x.backend().data());
}

// See T_frexp in TypeConversion.h. The backend is a GMP rational.
inline double T_frexp(boost::multiprecision::mpq_rational const &x, long &e) {
  long e_num, e_den;
  double m_num = mpz_get_d_2exp(&e_num, mpq_numref(x.backend().data()));
  double m_den = mpz_get_d_2exp(&e_den, mpq_denref(x.backend().data()));
  return T_frexp_quotient(m_num, e_num, m_den, e_den, e);
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYBOOSTGMPINT_H_
// clang-format on
