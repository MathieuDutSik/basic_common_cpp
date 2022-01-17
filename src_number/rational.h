#ifndef DEFINE_RATIONAL_H
#define DEFINE_RATIONAL_H

#include "TemplateTraits.h"
#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include "ExceptionEnding.h"


/*
  A rational class supposed to build integers.
  Typical use would be for example Rational<long> which should be pretty fast.
 */
template<typename Tint>
struct Rational {
private:
  Tint num;
  Tint den;
public:
  Rational() : num(0), den(1) {
  }
  Rational(Tint const& x) : num(x), den(1) {
  }
  Rational(Tint const& num, Tint const& den) : num(num), den(den) {
  }
  Rational(Rational<Tint> const& x) : num(x.num), den(x.den) {
  }
  Rational<Tint> operator=(Tint const& u) { // assignment operator from int
    num = u;
    den = 1;
    return *this;
  }
  Rational<Tint> operator=(Rational<Tint> const& x) { // assignment operator
    num = x.num;
    den = x.den;
    return *this;
  }
  const Tint& get_num() const {
    return num;
  }
  const Tint& get_den() const {
    return den;
  }
private: // A few internal functions.
  static Tint comp_gcd(Tint const& m, Tint const& n) {
    Tint f = m;
    if (m < 0) f = -m;
    Tint g = n;
    if (n < 0) g = -n;
    Tint h;
    while (g != 0) {
      h = g;
      g = f % g;
      f = h;
    }
    return f;
  }
  void gcd_reduction() {
    Tint gcd = comp_gcd(num,den);
    num = num / gcd;
    den = den / gcd;
  }
public:
  void operator/=(Rational<Tint> const&x) {
    if (x.num > 0) {
      num *= x.den;
      den *= x.num;
    } else {
      num *= -x.den;
      den *= -x.num;
    }
    gcd_reduction();
  }
  void operator/=(Tint const&x) {
    if (x > 0) {
      den *= x;
    } else {
      num = - num;
      den *= -x;
    }
    gcd_reduction();
  }
  void operator+=(Rational<Tint> const&x) {
    Tint gcd = comp_gcd(den, x.den);
    Tint new_den = den * x.den / gcd;
    num = num * (x.den / gcd) + x.num * (den/gcd);
    den = new_den;
    gcd_reduction(); // Yes, it is needed: example 1/2 + 1/2
  }
  void operator-=(Rational<Tint> const&x) {
    Tint gcd = comp_gcd(den, x.den);
    Tint new_den = den * x.den / gcd;
    num = num * (x.den / gcd) - x.num * (den/gcd);
    den = new_den;
    gcd_reduction(); // Yes, it is needed: example 1/2 - 1/2
  }
  friend Rational<Tint> operator+(Rational<Tint> const&x, Rational<Tint> const&y) {
    Rational<Tint> z;
    Tint gcd = Rational<Tint>::comp_gcd(x.den, x.den);
    z.den = y.den * x.den / gcd;
    z.num = y.num * (x.den / gcd) + x.num * (y.den/gcd);
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator+(Tint const&x, Rational<Tint> const&y) {
    Rational<Tint> z;
    z.num = x * y.den + y.num;
    z.den = y.den;
    return z;
  }
  friend Rational<Tint> operator+(Rational<Tint> const&x, Tint const&y) {
    Rational<Tint> z;
    z.num = x.num + y * x.den;
    z.den = x.den;
    return z;
  }
  friend Rational<Tint> operator-(Rational<Tint> const&x, Rational<Tint> const&y) {
    Rational<Tint> z;
    Tint gcd = Rational<Tint>::comp_gcd(x.den, x.den);
    z.den = x.den * y.den / gcd;
    z.num = x.num * (y.den / gcd) - y.num * (x.den/gcd);
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator-(Rational<Tint> const& x) {
    Rational<Tint> z;
    z.num = -x.num;
    z.den = x.den;
    return z;
  }
  friend Rational<Tint> operator/(Tint const&x, Rational<Tint> const&y) {
    Rational<Tint> z;
    if (y.num > 0) {
      z.num = x * y.den;
      z.den = y.num;
    } else {
      z.num = -x * y.den;
      z.den = -y.num;
    }
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator/(Rational<Tint> const&x, Rational<Tint> const&y) {
    Rational<Tint> z;
    if (y.num > 0) {
      z.num = x.num * y.den;
      z.den = x.den * y.num;
    } else {
      z.num = -x.num * y.den;
      z.den = -x.den * y.num;
    }
    z.gcd_reduction();
    return z;
  }
  void operator*=(Rational<Tint> const& x) {
    num = num * x.num;
    den = den * x.den;
    gcd_reduction();
  }
  friend Rational<Tint> operator*(Rational<Tint> const&x, Rational<Tint> const&y) {
    Rational<Tint> z;
    z.num = x.num * y.num;
    z.den = x.den * y.den;
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator*(Tint const&x, Rational<Tint> const&y) {
    Rational<Tint> z;
    z.num = x * y.num;
    z.den = y.den;
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator*(Rational<Tint> const&x, Tint const&y) {
    Rational<Tint> z;
    z.num = x.num * y;
    z.den = x.den;
    z.gcd_reduction();
    return z;
  }
  friend std::ostream& operator<<(std::ostream& os, Rational<Tint> const &v) {
    if (v.den == 1)
      return os << v.num;
    else
      return os << v.num << "/" << v.den;
  }
  friend std::istream& operator>>(std::istream &is, Rational<Tint> &v) {
    char c;
    std::string s;
    size_t miss_val = std::numeric_limits<size_t>::max();
    size_t pos_slash = miss_val;
    size_t pos = 0;
    while(true) {
      is >> c;
      if (c != ' ' && c != '\n') {
        s += c;
        pos++; // First character cannot be a slash
        break;
      }
    }
    while(true) {
      if (is.eof())
        break;
      is >> c;
      if (c == ' ' || c == '\n')
        break;
      s += c;
      if (c == '/')
        pos_slash = pos;
      pos++;
    }
    //
    if (pos_slash == miss_val) {
      std::istringstream isN(s);
      isN >> v.num;
      v.den = 1;
    } else {
      std::string sN = s.substr(0,pos_slash);
      std::string sD = s.substr(pos_slash+1, s.size() - 1 - pos_slash);
      std::istringstream isN(sN);
      isN >> v.num;
      std::istringstream isD(sD);
      isD >> v.den;
    }
    return is;
  }
  friend bool operator==(Rational<Tint> const&x, Rational<Tint> const&y) {
    return (x.num == y.num) && (x.den == y.den);
  }
  friend bool operator!=(Rational<Tint> const&x, Rational<Tint> const&y) {
    return (x.num != y.num) || (x.den != y.den);
  }
  friend bool operator!=(Rational<Tint> const&x, Tint const&y) {
    if (x.den > 1)
      return true;
    return x.num != y;
  }
  friend bool operator>=(Rational<Tint> const& x, Rational<Tint> const& y) {
    // x >= y is equivalent to x_n * y_d >= y_n * x_d
    return x.num * y.den >= y.num * x.den;
  }
  friend bool operator<=(Rational<Tint> const& x, Rational<Tint> const& y) {
    return x.num * y.den <= y.num * x.den;
  }
  friend bool operator<=(Rational<Tint> const& x, Tint const& y) {
    return x.num <= y * x.den;
  }
  friend bool operator>(Rational<Tint> const& x, Rational<Tint> const& y) {
    return x.num * y.den > y.num * x.den;
  }
  friend bool operator<(Rational<Tint> const& x, Rational<Tint> const& y) {
    return x.num * y.den < y.num * x.den;
  }
  friend bool operator<(Rational<Tint> const& x, Tint const& y) {
    return x.num < y * x.den;
  }
};

namespace std {
  template<typename Tint>
  std::string to_string(const Rational<Tint>& e_val)
  {
    std::stringstream s;
    s << e_val;
    std::string converted(s.str());
    return converted;
  };
}


template<typename Tint>
struct is_euclidean_domain<Rational<Tint>> {
  static const bool value = true;
};

template<typename Tint>
struct is_exact_arithmetic<Rational<Tint>> {
  static const bool value = is_exact_arithmetic<Tint>::value;
};

template<typename Tint>
struct is_implementation_of_Z<Rational<Tint>> {
  static const bool value = false;
};

template<typename Tint>
struct is_ring_field<Rational<Tint>> {
  static const bool value = true;
};

template<typename Tint>
struct is_totally_ordered<Rational<Tint>> {
  static const bool value = true;
};

template<typename Tint>
struct underlying_ring<Rational<Tint>> {
  typedef Tint ring_type;
};


template<typename Tint>
struct underlying_totally_ordered_ring<Rational<Tint>> {
  typedef Rational<Tint> real_type;
};

template<typename Tint>
inline Rational<Tint> T_NormGen(Rational<Tint> const& x)
{
  return T_abs(x);
}

template<typename Tint>
inline std::pair<Rational<Tint>,Rational<Tint>> ResQuoInt_kernel(Rational<Tint> const& a, Rational<Tint> const& b)
{
  // a = a_n / a_d
  // b = b_n / b_d
  // a = res + q * b  with 0 <= res < |b|
  // equivalent to
  // a_n / a_d = res + q * (b_n / b_d)
  // equivalent to
  // a_n * b_d = res * a_d * b_d + (q * a_d) * b_n
  using Tf= Rational<Tint>;
  Tint a_n = a.get_num();
  Tint b_n = b.get_num();
  Tint a_d = a.get_den();
  Tint b_d = b.get_den();
  Tint a1 = a_n * b_d;
  Tint b1 = a_d * b_n;
  Tint q = a1 / b1;
  Tf res = a - q * b;
  int sign;
  Tf b_abs;
  if (b < 0) {
    sign = -1;
    b_abs = -b;
  } else {
    sign = 1;
    b_abs = b;
  }
  while(true) {
    if (res < 0) {
      res += b_abs;
      q -= sign;
    } else {
      if (res >= b_abs) {
        res -= b_abs;
        q += sign;
      } else {
        if (res + q * b != a) {
          std::cerr << "Some error somewhere\n";
          std::cerr << "a=" << a << " b=" << b << "\n";
          std::cerr << "res=" << res << " q=" << q << "\n";
          Tf val = res + q * b;
          Tf val1 = q * b;
          Tf val2 = res + val1;
          std::cerr << "v=" << val << " val1=" << val1 << " val2=" << val2 << "\n";
          throw TerminalException{1};
        }
        return {res,q}; // Conversion of q from Tint to Rational<Tint> occurring here.
      }
    }
  }
}

template<typename Tint>
inline Rational<Tint> ResInt(Rational<Tint> const& a, Rational<Tint> const& b)
{
  return ResQuoInt_kernel(a, b).first;
}

template<typename Tint>
inline Rational<Tint> QuoInt(Rational<Tint> const& a, Rational<Tint> const& b)
{
  return ResQuoInt_kernel(a, b).second;
}


// OVerlying fields, the original motivation

template<>
struct overlying_field<int> {
  typedef Rational<int> field_type;
};

template<>
struct overlying_field<short> {
  typedef Rational<short> field_type;
};

template<>
struct overlying_field<long> {
  typedef Rational<long> field_type;
};

template<>
struct overlying_field<Rational<int>> {
  typedef Rational<int> field_type;
};

template<>
struct overlying_field<Rational<short>> {
  typedef Rational<short> field_type;
};

template<>
struct overlying_field<Rational<long>> {
  typedef Rational<long> field_type;
};

// The conversion tools (int)

inline void TYPE_CONVERSION(Rational<int> const& a1, Rational<int> & a2)
{
  a2 = a1;
}

inline void TYPE_CONVERSION(Rational<int> const& a1, int & a2)
{
  const int& den = a1.get_den();
  if (den != 1) {
    std::string str_err = "The denominator should be 1. It is den = " + std::to_string(den);
    throw ConversionException{str_err};
  }
  a2 = a1.get_num();
}

inline void TYPE_CONVERSION(int const& a1, Rational<int> & a2)
{
  a2 = a1;
}

// The conversion tools (long)

inline void TYPE_CONVERSION(Rational<long> const& a1, Rational<long> & a2)
{
  a2 = a1;
}

inline void TYPE_CONVERSION(Rational<long> const& a1, long & a2)
{
  const long& den = a1.get_den();
  if (den != 1) {
    std::string str_err = "The denominator should be 1. It is den = " + std::to_string(den);
    throw ConversionException{str_err};
  }
  a2 = a1.get_num();
}

inline void TYPE_CONVERSION(long const& a1, Rational<long> & a2)
{
  a2 = a1;
}

// Obtention of denominators

template<typename Tint>
inline Rational<Tint> GetDenominator(Rational<Tint> const& x)
{
  return x.get_den();
}

template<typename Tint>
inline Tint GetDenominator_z(Rational<Tint> const& x)
{
  return x.get_den();
}

// Fllor / Ceil / Nearest operations

template<typename Tint>
Rational<Tint> FractionalPart(Rational<Tint> const& x)
{
  Tint res = ResInt(x.num, x.den);
  Rational<Tint> fr(res, x.den);
  return fr;
}

template<typename Tint>
inline void FloorInteger(Rational<Tint> const& xI, Rational<Tint> & xO)
{
  Rational<Tint> fr = FractionalPart(xI);
  xO = xI - fr;
}

template<typename Tint>
inline void CeilInteger(Rational<Tint> const& xI, Rational<Tint> & xO)
{
  Rational<Tint> fr = FractionalPart(xI);
  if (fr == 0) {
    xO = xI;
  } else {
    xO = 1 + xI - fr;
  }
}

template<typename Tint>
Rational<Tint> NearestInteger_rni(Rational<Tint> const& a)
{
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

template<typename Tint>
inline void NearestInteger(Rational<Tint> const& xI, Rational<Tint> & xO)
{
  xO = NearestInteger_rni(xI);
}

template<typename Tint>
inline void NearestInteger(Rational<Tint> const& xI, Tint & xO)
{
  Rational<Tint> xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_num();
}




#endif
