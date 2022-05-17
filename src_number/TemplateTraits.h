// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_TEMPLATETRAITS_H_
#define SRC_NUMBER_TEMPLATETRAITS_H_

// Trait for specific types.
// the std::is_same<T,mpz_class> is not adequate because it requires the type
// mpz_class to be known in the scope.

template <typename T> struct is_mpq_class { static const bool value = false; };

template <typename T> struct is_mpz_class { static const bool value = false; };

template <typename T> struct is_boost_cpp_int {
  static const bool value = false;
};

template <typename T> struct is_boost_cpp_rational {
  static const bool value = false;
};

template <typename T> struct is_boost_mpz_int {
  static const bool value = false;
};

template <typename T> struct is_boost_mpq_rational {
  static const bool value = false;
};

// is_euclidean_domain

template <typename T> struct is_euclidean_domain {
  static const bool value = false;
};

template <> struct is_euclidean_domain<int64_t> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<int32_t> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<int16_t> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<int8_t> {
  static const bool value = true;
};

// overlying_field block

template <typename T> struct overlying_field {};

template <> struct overlying_field<double> { typedef double field_type; };

template <> struct overlying_field<float> { typedef float field_type; };

// underlying_ring

template <typename T> struct underlying_ring {};

template <> struct underlying_ring<int64_t> { typedef int64_t ring_type; };

template <> struct underlying_ring<int32_t> { typedef int32_t ring_type; };

template <> struct underlying_ring<int16_t> { typedef int16_t ring_type; };

template <> struct underlying_ring<int8_t> { typedef int8_t ring_type; };

// Trait definition for subset of integers

template <typename T> struct is_implementation_of_Z {};

template <> struct is_implementation_of_Z<double> {
  static const bool value = false;
};

template <> struct is_implementation_of_Z<float> {
  static const bool value = false;
};

template <> struct is_implementation_of_Z<int64_t> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<int32_t> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<int16_t> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<int8_t> {
  static const bool value = true;
};

// Trait definition for exactness

template <typename T> struct is_exact_arithmetic {};

template <> struct is_exact_arithmetic<double> {
  static const bool value = false;
};

template <> struct is_exact_arithmetic<float> {
  static const bool value = false;
};

// Trait definition for fields

template <typename T> struct is_ring_field {};

template <> struct is_ring_field<int64_t> { static const bool value = false; };

template <> struct is_ring_field<int32_t> { static const bool value = false; };

template <> struct is_ring_field<int16_t> { static const bool value = false; };

template <> struct is_ring_field<int8_t> { static const bool value = false; };

template <> struct is_ring_field<double> { static const bool value = true; };

template <> struct is_ring_field<float> { static const bool value = true; };

// Trait of totally ordered set

template <typename T> struct is_totally_ordered {
  static const bool value = false;
};

template <> struct is_totally_ordered<int64_t> {
  static const bool value = true;
};

template <> struct is_totally_ordered<int32_t> { static const bool value = true; };

template <> struct is_totally_ordered<int16_t> { static const bool value = true; };

template <> struct is_totally_ordered<int8_t> {
  static const bool value = true;
};

template <> struct is_totally_ordered<double> {
  static const bool value = true;
};

template <> struct is_totally_ordered<float> {
  static const bool value = true;
};

// is floating arithmetic

template <typename T> struct is_float_arithmetic {
  static const bool value = false;
};

template <> struct is_float_arithmetic<float> {
  static const bool value = true;
};

template <> struct is_float_arithmetic<double> {
  static const bool value = true;
};

// Trait definition for underlying ring

template <typename T> struct underlying_totally_ordered_ring {};

template <> struct underlying_totally_ordered_ring<int64_t> {
  typedef int64_t real_type;
};

template <> struct underlying_totally_ordered_ring<int32_t> {
  typedef int32_t real_type;
};

template <> struct underlying_totally_ordered_ring<int16_t> {
  typedef int16_t real_type;
};

template <> struct underlying_totally_ordered_ring<int8_t> {
  typedef int8_t real_type;
};

// clang-format off
#endif  // SRC_NUMBER_TEMPLATETRAITS_H_
// clang-format on
