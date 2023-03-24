// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_FP_H_
#define SRC_NUMBER_FP_H_

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <boost/math/special_functions/round.hpp>
#include "rational.h"
//#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
//#include "TypeConversion.h"

/*
  A finite prime field arithmetic class supposed.
  Typical use would be for example Fp<long, 2147389441> which should be pretty fast.
 */
template <typename Tint, Tint P> struct Fp {
private:
  Tint num;

public:
  // Constructors
  Fp() : num(0) {}
  Fp(Tint const &x) : num(x) { reduce(); }
  // Assignment operators
  Fp<Tint, P> operator=(Tint const &u) {
    // assignment operator from int
    num = u;
    reduce();
    return *this;
  }
  Tint &get_num() { return num; }
  const Tint &get_const_num() const { return num; }

private:
  // A few internal functions.
  static Tint comp_gcd(Tint const &m, Tint const &n) {
    Tint f = m;
    if (m < 0)
      f = -m;
    Tint g = n;
    if (n < 0)
      g = -n;
    Tint h;
    while (g != 0) {
      h = g;
      g = f % g;
      f = h;
    }
    return f;
  }
  // assume m is reduced
  static Tint mod_inverse(Tint const &m) {
    Tint t = 0; Tint newt = 1;
    Tint r = P; Tint newr = m; 
    Tint q, tmp;
    while (newr != 0) {
      q = r / newr;
      tmp = t;
      t = newt;
      newt = tmp - q * newt;
      tmp = r;
      r = newr;
      newr = tmp - q * newr;
    }
    if (r > 1)
      return 0;
    if (t < 0)
      t = t+P;
    return t;
  }
  void reduce() {
    num = num % P;
    if( num < 0 )
      num += P;
  }

public:
  Rational<Tint> rational_lift() {
    // we create the integer lattice (a,b) s.t. a = b*num mod P 
    // has basis (num, 1) and (P,0)
    // use lagrange reduction to find small (a,b)
    Tint a0 = P;
    Tint a1 = 0;
    Tint b0 = num;
    Tint b1 = 1;
    using Tfloat = long double;
    Tfloat A = P*P;
    Tfloat B = num*num+1;
    while(true) {
      Tfloat n = Tfloat(a0*b0 + a1*b1);
      Tint r = Tint(boost::math::lround(n/B));
      Tfloat T = A-2*Tfloat(r)*n + Tfloat(r)*Tfloat(r)*B;
      if( T+0.5 >= B ) {
        // we are done
        if( b1 < 0 )
          return Rational<Tint>(-b0, -b1);
        return Rational<Tint>(b0,b1);
      }
      Tint tmp0 = a0 - r*b0;
      Tint tmp1 = a1 - r*b1;
      a0 = b0;
      a1 = b1;
      b0 = tmp0;
      b1 = tmp1;
      A = B;
      B = T;
    }
  }
  
  void operator/=(Fp<Tint,P> const &x) {
    num = num / x;
  }
  void operator+=(Fp<Tint,P> const &x) {
    num = num+x;
  }
  void operator-=(Fp<Tint,P> const &x) {
    num = num-x;
  }
  friend Fp<Tint,P> operator+(Fp<Tint,P> const &x,
                                  Fp<Tint,P> const &y) {
    Fp<Tint, P> z;
    Tint zz = x.num + y.num;
    if (zz > P) zz -= P;
    z.num = zz;
    return z;
  }
  friend Fp<Tint,P> operator-(Fp<Tint,P> const &x,
                                  Fp<Tint,P> const &y) {
    Fp<Tint, P> z;
    Tint zz = x.num - y.num;
    if( zz < 0 ) zz += P;
    z.num = zz;
    return z;
  }
  friend Fp<Tint,P> operator/(Fp<Tint,P> const &x,
                                  Fp<Tint,P> const &y) {
    Tint den = mod_inverse(y.num);  
    Fp<Tint,P> z(x.num * den);
    return z;
  }
  void operator*=(Fp<Tint,P> const &x) {
    num = num * x;
  }
  friend Fp<Tint,P> operator*(Fp<Tint,P> const &x,
                                  Fp<Tint,P> const &y) {
    Fp<Tint, P> z(x.num * y.num);
    return z;
  }
  friend std::ostream &operator<<(std::ostream &os, Fp<Tint,P> const &v) {
    return os << v.num << " (mod " << P << ")";
  }
  friend bool operator==(Fp<Tint,P> const &x, Fp<Tint,P> const &y) {
    return (x.num == y.num);
  }
  friend bool operator!=(Fp<Tint,P> const &x, Fp<Tint,P> const &y) {
    return (x.num != y.num);
  }
  friend bool operator!=(Fp<Tint,P> const &x, Tint const &y) {
    return x.num != y;
  }
};

template <typename Tint, Tint P> struct is_euclidean_domain<Fp<Tint,P>> {
  static const bool value = true;
};

template <typename Tint, Tint P> struct is_exact_arithmetic<Fp<Tint,P>> {
  static const bool value = is_exact_arithmetic<Tint>::value;
};

template <typename Tint, Tint P> struct is_implementation_of_Z<Fp<Tint, P>> {
  static const bool value = false;
};

template <typename Tint, Tint P> struct is_implementation_of_Q<Fp<Tint,P>> {
  static const bool value = false;
};

template <typename Tint, Tint P> struct is_ring_field<Fp<Tint,P>> {
  static const bool value = true;
};

template <typename Tint, Tint P> struct is_totally_ordered<Fp<Tint,P>> {
  static const bool value = false;
};

template <typename Tint, Tint P> struct underlying_ring<Fp<Tint,P>> {
  typedef Tint ring_type;
};

// hashing function

namespace std {
template <typename T, T P> struct hash<Fp<T,P>> {
  std::size_t operator()(const Fp<T,P> &x) const {
    return std::hash<T>()(x.get_const_num());
  }
};
// clang-format off
}  // namespace std
// clang-format on

// The conversion tools (int)
/*
inline void TYPE_CONVERSION(stc<Rational<int>> const &a1, Rational<int> &a2) {
  a2 = a1.val;
}

template <typename T>
void TYPE_CONVERSION_Rational_T(stc<Rational<T>> const &a1, T &a2) {
  const T &den = a1.val.get_const_den();
  if (den != 1) {
    std::string str_err =
        "The denominator should be 1. It is den = " + std::to_string(den);
    throw ConversionException{str_err};
  }
  a2 = a1.val.get_const_num();
}

inline void TYPE_CONVERSION(stc<Rational<int>> const &a1, int &a2) {
  TYPE_CONVERSION_Rational_T<int>(a1, a2);
}

inline void TYPE_CONVERSION(stc<int> const &a1, Rational<int> &a2) {
  a2 = a1.val;
}

// The conversion tools (long)

inline void TYPE_CONVERSION(stc<Rational<long>> const &a1, Rational<long> &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<Rational<long>> const &a1, long &a2) {
  TYPE_CONVERSION_Rational_T<long>(a1, a2);
}

inline void TYPE_CONVERSION(stc<long> const &a1, Rational<long> &a2) {
  a2 = a1.val;
}

// Obtention of denominators

template <typename Tint>
inline Rational<Tint> GetDenominator(Rational<Tint> const &x) {
  return x.get_const_den();
}

template <typename Tint> inline Tint GetDenominator_z(Rational<Tint> const &x) {
  return x.get_const_den();
}

template <typename Tint>
void ScalingInteger_Kernel(stc<Rational<Tint>> const &x, Tint &x_ret) {
  x_ret = x.val.get_const_den();
}

// Floor / Ceil / Nearest operations

template <typename Tint>
Rational<Tint> FractionalPart(Rational<Tint> const &x) {
  Tint res = ResInt(x.get_const_num(), x.get_const_den());
  Rational<Tint> fr(res, x.get_const_den());
  return fr;
}

template <typename Tint>
inline void FloorInteger(Rational<Tint> const &xI, Rational<Tint> &xO) {
  Rational<Tint> fr = FractionalPart(xI);
  xO = xI - fr;
}

template <typename Tint>
inline void CeilInteger(Rational<Tint> const &xI, Rational<Tint> &xO) {
  Rational<Tint> fr = FractionalPart(xI);
  if (fr == 0) {
    xO = xI;
  } else {
    xO = 1 + xI - fr;
  }
}

template <typename Tint>
Rational<Tint> NearestInteger_rni(Rational<Tint> const &a) {
  Rational<Tint> fr = FractionalPart(a);
  Rational<Tint> eDiff1 = fr;
  Rational<Tint> eDiff2 = 1 - fr;
  Rational<Tint> RetVal = a - fr;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return 1 + RetVal;
  }
}

template <typename Tint>
inline void NearestInteger(Rational<Tint> const &xI, Rational<Tint> &xO) {
  xO = NearestInteger_rni(xI);
}

template <typename Tint>
inline void NearestInteger(Rational<Tint> const &xI, Tint &xO) {
  Rational<Tint> xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_num();
}

// Serialization stuff

namespace boost::serialization {

template <class Archive, typename T>
inline void serialize(Archive &ar, Rational<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("rational_num", val.get_num());
  ar &make_nvp("rational_den", val.get_den());
}
*/
// clang-format off
//}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_FP_H_
// clang-format on
