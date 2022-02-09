#ifndef DEFINE_TEMPLATE_TRAITS_H
#define DEFINE_TEMPLATE_TRAITS_H

// Trait for specific types.
// the std::is_same<T,mpz_class> is not adequate because it requires the type
// mpz_class to be known in the scope.

template <typename T>
struct is_mpq_class {
  static const bool value = false;
};

template <typename T>
struct is_mpz_class {
  static const bool value = false;
};



// is_euclidean_domain

template <typename T>
struct is_euclidean_domain {
  static const bool value = false;
};

template <>
struct is_euclidean_domain<short> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<long> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<int> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<long long> {
  static const bool value = true;
};

// overlying_field block

template<typename T>
struct overlying_field {
};

template<>
struct overlying_field<double> {
  typedef double field_type;
};

template<>
struct overlying_field<float> {
  typedef float field_type;
};


// underlying_ring

template<typename T>
struct underlying_ring {
};

template<>
struct underlying_ring<int> {
  typedef int ring_type;
};

template<>
struct underlying_ring<long> {
  typedef long ring_type;
};




// Trait definition for subset of integers

template <typename T>
struct is_implementation_of_Z {
};

template<>
struct is_implementation_of_Z<double> {
  static const bool value = false;
};

template<>
struct is_implementation_of_Z<float> {
  static const bool value = false;
};

template<>
struct is_implementation_of_Z<int> {
  static const bool value = true;
};

template<>
struct is_implementation_of_Z<long> {
  static const bool value = true;
};

// Trait definition for exactness

template <typename T>
struct is_exact_arithmetic {
};

template<>
struct is_exact_arithmetic<double> {
  static const bool value = false;
};

template<>
struct is_exact_arithmetic<float> {
  static const bool value = false;
};

// Trait definition for fields

template <typename T>
struct is_ring_field {
};

template <>
struct is_ring_field<short> {
  static const bool value = false;
};

template <>
struct is_ring_field<long> {
  static const bool value = false;
};

template <>
struct is_ring_field<int> {
  static const bool value = false;
};

template <>
struct is_ring_field<long long> {
  static const bool value = false;
};

template <>
struct is_ring_field<double> {
  static const bool value = true;
};

template <>
struct is_ring_field<float> {
  static const bool value = true;
};

// Trait of totally ordered set

template <typename T>
struct is_totally_ordered {
  static const bool value = false;
};

template <>
struct is_totally_ordered<short> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<long> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<int> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<long long> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<double> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<float> {
  static const bool value = true;
};

// is floating arithmetic

template <typename T>
struct is_float_arithmetic {
  static const bool value = false;
};

template<>
struct is_float_arithmetic<float> {
  static const bool value = true;
};

template<>
struct is_float_arithmetic<double> {
  static const bool value = true;
};

// Trait definition for underlying ring

template<typename T>
struct underlying_totally_ordered_ring {
};

template<>
struct underlying_totally_ordered_ring<int> {
  typedef int real_type;
};

template<>
struct underlying_totally_ordered_ring<long> {
  typedef long real_type;
};


#endif
