#ifndef TYPE_CONVERSION_INCLUDE
#define TYPE_CONVERSION_INCLUDE

//All the definitions of special fields are in other include.
//Nothing of this should depend on GMP or MPREAL or FLINT or whatever.
//
// All GMP are in NumberTheory.h
// All mpreal are in mpreal_related.h

#include <math.h>
#include <type_traits>

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

// is mpreal

template <typename T>
struct is_mpreal {
  static const bool value = false;
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

// type mappings for the rings and fields:
template<typename T>
struct overlying_field {
};

template<typename T>
struct underlying_ring {
};






// Trait definition for fields.

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

// is double type

template <typename T>
struct is_double_type {
  static const bool value = false;
};

template <>
struct is_double_type<double> {
  static const bool value = true;
};

// is double type

template <typename T>
struct is_int_type {
  static const bool value = false;
};

template <>
struct is_int_type<int> {
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

// Trait definition for is_mpq

template <typename T>
struct is_mpq_class {
  static const bool value = false;
};


// basic numeric code

int IntFloor(double const& x)
{
  return int(floor(x));
}

// Trait definition for is_mpz

template <typename T>
struct is_mpz_class {
  static const bool value = false;
};




template<typename T>
T T_abs(T const& eVal)
{
  if (eVal > 0)
    return eVal;
  T fVal= - eVal;
  return fVal;
}

void NearestInteger_double_int(double const& xI, int & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_z=int(xRnd_d);
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](int const& u) -> double {
    double diff = double(u) - xI;
    return T_abs(diff);
  };
  double err=GetErr(xRnd_z);
  //  std::cerr << "err=" << err << "\n";
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      int xTest = xRnd_z + shift;
      double TheErr=GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest << " TheErr=" << TheErr << "\n";
      if (TheErr < err) {
	IsOK=false;
	xRnd_z=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_z;
}


void TYPE_CONVERSION(double const& a1, double & a2)
{
  a2 = a1;
}

void TYPE_CONVERSION(double const& a1, int & a2)
{
  a2 = int(a1);
}

void TYPE_CONVERSION(double const& a1, long & a2)
{
  a2 = long(a1);
}

void TYPE_CONVERSION(int const& a1, double & a2)
{
  a2 = double(a1);
}

void TYPE_CONVERSION(long const& a1, double & a2)
{
  a2 = double(a1);
}




template<typename To>
void NearestInteger_double_To(double const& xI, To & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_i=int(xRnd_d);
  To xRnd_To=xRnd_i;
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](To const& u) -> double {
    double u_doubl;
    TYPE_CONVERSION(u, u_doubl);
    double diff = u_doubl - xI;
    return T_abs(diff);
  };
  double err=GetErr(xRnd_To);
  //  std::cerr << "err=" << err << "\n";
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      To xTest = xRnd_To + shift;
      double TheErr=GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest << " TheErr=" << TheErr << "\n";
      if (TheErr < err) {
	IsOK=false;
	xRnd_To=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_To;
}





template<typename To, typename Ti>
inline typename std::enable_if<(not is_mpreal<Ti>::value) && (not is_double_type<Ti>::value),To>::type UniversalNearestInteger(Ti const& a)
{
  To ret;
  NearestInteger(a, ret);
  //  std::cerr << "Temp_common a=" << a << " ret=" << ret << "\n";
  return ret;
}


template<typename To, typename Ti>
inline typename std::enable_if<is_double_type<Ti>::value,To>::type UniversalNearestInteger(Ti const& a)
{
  To ret;
  NearestInteger_double_To<To>(a, ret);
  //  std::cerr << "Temp_common(I) a=" << a << " ret=" << ret << "\n";
  return ret;
}


template<typename T1, typename T2>
T1 UniversalTypeConversion(T2 const& a)
{
  T1 ret;
  TYPE_CONVERSION(a, ret);
  return ret;
}




#endif
