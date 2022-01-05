#ifndef TYPE_CONVERSION_INCLUDE
#define TYPE_CONVERSION_INCLUDE

//All the definitions of special fields are in other include.
//Nothing of this should depend on GMP or MPREAL or FLINT or whatever.
//
// All GMP are in NumberTheory.h
// All mpreal are in mpreal_related.h

#include <math.h>
#include <type_traits>
#include "ExceptionEnding.h"
#include "TemplateTraits.h"

void NearestInteger_double_int(double const& xI, int & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_z=int(xRnd_d);
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](int const& u) -> double {
    double diff = double(u) - xI;
    if (diff < 0)
      return -diff;
    return diff;
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


inline void TYPE_CONVERSION(double const& a1, double & a2)
{
  a2 = a1;
}

inline void TYPE_CONVERSION(double const& a1, uint8_t & a2)
{
  a2 = uint8_t(a1);
}

inline void TYPE_CONVERSION(double const& a1, int & a2)
{
  a2 = int(a1);
}

inline void TYPE_CONVERSION(double const& a1, long & a2)
{
  a2 = long(a1);
}

/*
inline void TYPE_CONVERSION(double const& a1, int64_t & a2)
{
  a2 = int64_t(a1);
}
*/

inline void TYPE_CONVERSION(int const& a1, double & a2)
{
  a2 = double(a1);
}

inline void TYPE_CONVERSION(uint8_t const& a1, double & a2)
{
  a2 = double(a1);
}

/*
inline void TYPE_CONVERSION(int64_t const& a1, double & a2)
{
  a2 = double(a1);
}
*/

inline void TYPE_CONVERSION(long const& a1, double & a2)
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
    if (diff < 0)
      return -diff;
    return diff;
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
inline To UniversalFloorScalarInteger(Ti const& a)
{
  To ret;
  FloorInteger(a, ret);
  return ret;
}

template<typename To, typename Ti>
inline To UniversalCeilScalarInteger(Ti const& a)
{
  To ret;
  CeilInteger(a, ret);
  return ret;
}



template<typename To, typename Ti>
inline typename std::enable_if<(not is_mpreal<Ti>::value) && (not is_double_type<Ti>::value),To>::type UniversalNearestScalarInteger(Ti const& a)
{
  To ret;
  NearestInteger(a, ret);
  return ret;
}


template<typename To, typename Ti>
inline typename std::enable_if<is_double_type<Ti>::value,To>::type UniversalNearestScalarInteger(Ti const& a)
{
  To ret;
  NearestInteger_double_To<To>(a, ret);
  return ret;
}



template<typename T1, typename T2>
T1 UniversalScalarConversion(T2 const& a)
{
  T1 ret;
  try {
    TYPE_CONVERSION(a, ret);
  }
  catch (ConversionException & e) {
    std::cerr << "ConversionError e=" << e.val << "\n";
    throw TerminalException{1};
  }
  return ret;
}



template<typename T1, typename T2>
std::pair<bool,T1> UniversalScalarConversionCheck(T2 const& a)
{
  T1 ret;
  try {
    TYPE_CONVERSION(a, ret);
  }
  catch (ConversionException & e) {
    return {false,ret};
  }
  return {true,ret};
}







#endif
