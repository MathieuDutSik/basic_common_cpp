// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_QUADFIELD_H_
#define SRC_NUMBER_QUADFIELD_H_

#include "NumberTheory.h"
#include "Temp_common.h"

template <typename Tinp, int d> class QuadField {
public:
  using Tresidual = Tinp;
private:
  using T = Tinp;
  T a;
  T b;

public:
  T& get_a() {
    return a;
  }
  T& get_b() {
    return b;
  }
  const T& get_const_a() const{
    return a;
  }
  const T& get_const_b() const {
    return b;
  }

  // Constructor
  QuadField() : a(0), b(0) {}
  QuadField(T const &u) : a(u), b(0) {}
  QuadField(T const &_a, T const& _b) : a(_a), b(_b) {}
  QuadField(QuadField<T, d> const &x) : a(x.a), b(x.b) {}
  //  QuadField<T,d>& operator=(QuadField<T,d> const&); // assignment operator
  //  QuadField<T,d>& operator=(T const&); // assignment operator from T
  //  QuadField<T,d>& operator=(int const&); // assignment operator from T
  // assignment operator from int
  QuadField<T, d> operator=(T const &u) {
    a = u;
    b = 0;
    return *this;
  }
  // assignment operator
  QuadField<T, d> operator=(QuadField<T, d> const &x) {
    a = x.a;
    b = x.b;
    return *this;
  }
  //
  // Arithmetic operators below:
  void operator+=(QuadField<T, d> const &x) {
    a += x.a;
    b += x.b;
  }
  void operator-=(QuadField<T, d> const &x) {
    a -= x.a;
    b -= x.b;
  }
  void operator/=(QuadField<T, d> const &x) {
    T disc  = x.a * x.a - d * x.b * x.b;
    T a_new = (a * x.a - d * b * x.b) / disc;
    b       = (b * x.a -     a * x.b) / disc;
    a = a_new;
  }
  friend QuadField<T, d> operator+(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z.a = x.a + y.a;
    z.b = x.b + y.b;
    return z;
  }
  friend QuadField<T, d> operator-(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z.a = x.a - y.a;
    z.b = x.b - y.b;
    return z;
  }
  friend QuadField<T, d> operator-(QuadField<T, d> const &x, T const &y) {
    QuadField<T, d> z;
    z.a = x.a - y;
    z.b = x.b;
    return z;
  }
  friend QuadField<T, d> operator-(QuadField<T, d> const &x) {
    QuadField<T, d> z;
    z.a = -x.a;
    z.b = -x.b;
    return z;
  }
  friend QuadField<T, d> operator/(T const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    T disc = y.a * y.a - d * y.b * y.b;
    z.a = x * y.a / disc;
    z.b = -x * y.b / disc;
    return z;
  }
  friend QuadField<T, d> operator/(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    QuadField<T, d> z;
    T disc = y.a * y.a - d * y.b * y.b;
    z.a = (x.a * y.a - d * x.b * y.b) / disc;
    z.b = (x.b * y.a - x.a * y.b) / disc;
    return z;
  }
  double get_d() const {
    T hA, hB;
    double hA_d, hB_d, doubl_d;
    hA_d = a.get_d();
    hB_d = b.get_d();
    doubl_d = static_cast<double>(d);
    return hA_d + sqrt(doubl_d) * hB_d;
  }
  void operator*=(QuadField<T, d> const &x) {
    T hA = a * x.a + d * b * x.b;
    b    = a * x.b + b * x.a;
    a = hA;
  }
  friend QuadField<T, d> operator*(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z.a = x.a * y.a + d * x.b * y.b;
    z.b = x.a * y.b + x.b * y.a;
    return z;
  }
  friend QuadField<T, d> operator*(T const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z.a = x * y.a;
    z.b = x * y.b;
    return z;
  }
  friend std::ostream &operator<<(std::ostream &os, QuadField<T, d> const &v) {
    return os << "(" << v.a << "," << v.b << ")";
  }
  friend std::istream &operator>>(std::istream &is, QuadField<T, d> &v) {
    char c;
    std::string s;
    size_t miss_val = std::numeric_limits<size_t>::max();
    size_t pos;
    // First skipping the spaces
    while(true) {
      is.get(c);
      if (c != ' ' && c != '\n') {
        s += c;
        break;
      }
    }
    // Second reading the data till a space or end.
    size_t pos_comma = miss_val;
    while(true) {
      if (is.eof()) {
        break;
      }
      is.get(c);
      if (c == ' ' || c == '\n') {
        break;
      }
      s += c;
      if (c == ',')
        pos_comma = pos;
      pos++;
    }
    // Now parsing the data
    if (pos_comma == miss_val) {
      std::cerr << "Failed to find the comma\n";
      throw TerminalException{1};
    }
    std::string sA = s.substr(0, pos_comma);
    std::string sB = s.substr(pos_comma+1, s.size() - 2 - pos_comma);
    std::istringstream isA(sA);
    isA >> v.a;
    std::istringstream isB(sB);
    isB >> v.b;
    return is;
  }
  friend bool operator==(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    if (x.a != y.a)
      return false;
    if (x.b != y.b)
      return false;
    return true;
  }
  friend bool operator!=(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    if (x.a != y.a)
      return true;
    if (x.b != y.b)
      return true;
    return false;
  }
  friend bool operator!=(QuadField<T, d> const &x, T const &y) {
    if (x.a != y)
      return true;
    if (x.b != 0)
      return true;
    return false;
  }
  friend bool IsNonNegative(QuadField<T, d> const &x) {
    if (x.a == 0 && x.b == 0)
      return true;
    if (x.a >= 0 && x.b >= 0)
      return true;
    if (x.a <= 0 && x.b <= 0)
      return false;
    T disc = x.a * x.a - d * x.b * x.b;
    if (disc > 0) {
      if (x.a >= 0 && x.b <= 0)
        return true;
      if (x.a <= 0 && x.b >= 0)
        return false;
    } else {
      if (x.a >= 0 && x.b <= 0)
        return false;
      if (x.a <= 0 && x.b >= 0)
        return true;
    }
    std::cerr << "Major errors in the code\n";
    return false;
  }
  friend bool operator>=(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = x - y;
    return IsNonNegative(z);
  }
  friend bool operator<=(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = y - x;
    return IsNonNegative(z);
  }
  friend bool operator<=(QuadField<T, d> const &x, T const &y) {
    QuadField<T, d> z;
    z = y - x;
    return IsNonNegative(z);
  }
  friend bool operator>(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = x - y;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  friend bool operator<(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = y - x;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  friend bool operator<(QuadField<T, d> const &x, T const &y) {
    QuadField<T, d> z;
    z = y - x;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
};

template<typename T, int d>
struct overlying_field<QuadField<T,d>> { typedef QuadField<typename QuadField<T,d>::Tresidual,d> field_type; };

// Note that the underlying ring is not unique, there are many possibiliies actually
// but we can represent only one in our scheme.
template<typename T, int d>
struct underlying_ring<QuadField<T,d>> { typedef QuadField<typename QuadField<T,d>::Tresidual,d> ring_type; };


template <typename T, int d>
inline void TYPE_CONVERSION(stc<QuadField<T, d>> const &eQ, double &eD) {
  eD = eQ.val.get_d();
}

template <typename T, int d> struct is_totally_ordered<QuadField<T, d>> {
  static const bool value = true;
};

template <typename T, int d> struct is_ring_field<QuadField<T, d>> {
  static const bool value = is_ring_field<T>::value;
};

template <typename T, int d> struct is_exact_arithmetic<QuadField<T, d>> {
  static const bool value = true;
};

// Hashing function

template <typename T, int d>
struct is_implementation_of_Z<QuadField<T,d>> {
  static const bool value = false;
};

template <typename T, int d>
struct is_implementation_of_Q<QuadField<T,d>> {
  static const bool value = false;
};

// Hashing function

namespace std {
  template <typename T,int d> struct hash<QuadField<T,d>> {
    std::size_t operator()(const QuadField<T,d> &x) const {
      auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      	seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      size_t seed = std::hash<int>()(d);
      size_t e_hash1 = std::hash<T>()(x.get_const_a());
      size_t e_hash2 = std::hash<T>()(x.get_const_b());
      combine_hash(seed, e_hash1);
      combine_hash(seed, e_hash2);
      return seed;
    }
  };
}


// Local typing info

template<typename T> struct is_quad_field {};

template<typename T, int d> struct is_quad_field<QuadField<T,d>> { static const bool value = false; };

// Some functionality

template<typename T, int d>
bool IsInteger(QuadField<T,d> const& x) {
  if (x.get_b() != 0)
    return false;
  return IsInteger(x.get_const_a());
}


// The conversion tools (int)


template<typename T1, typename T2, int d>
inline void TYPE_CONVERSION(stc<QuadField<T1,d>> const &x1, QuadField<T2,d> &x2) {
  stc<T1> a1 { x1.val.get_const_a() };
  stc<T1> b1 { x1.val.get_const_b() };
  TYPE_CONVERSION(a1, x2.get_a());
  TYPE_CONVERSION(b1, x2.get_b());
}

template<typename T1, typename T2, int d>
inline void TYPE_CONVERSION(stc<QuadField<T1,d>> const &x1, double &x2) {
  stc<T1> a1 { x1.val.get_const_a() };
  stc<T1> b1 { x1.val.get_const_b() };
  double a2, b2;
  TYPE_CONVERSION(a1, a2);
  TYPE_CONVERSION(b1, b2);
  x2 = a2 + sqrt(d) * b2;
}

template<typename T1, typename T2, int d>
inline typename std::enable_if<not is_quad_field<T2>::value, void>::type TYPE_CONVERSION(stc<QuadField<T1,d>> const &x1, T2 &x2) {
  if (x1.val.b != 0) {
    std::string str = "Conversion error for quadratic field";
    throw ConversionException{str};
  }
  stc<T1> a1 { x1.val.get_const_a() };
  TYPE_CONVERSION(a1, x2);
}

// Serialization stuff

namespace boost::serialization {

template <class Archive, typename T, int d>
inline void serialize(Archive &ar, QuadField<T,d> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("quadfield_a", val.get_a());
  ar &make_nvp("quadfield_b", val.get_b());
}

}

// Turning into something rational

template<typename Tring, typename Tquad>
void ScalingInteger_Kernel(stc<Tquad> const& x, Tring& x_res) {
  using Tfield = typename Tquad::Tresidual;
  Tfield const& a = x.val.get_const_a();
  Tfield const& b = x.val.get_const_b();
  x_res = LCMpair(GetDenominator_z(a), GetDenominator_z(b));
}


// clang-format off
#endif  // SRC_NUMBER_QUADFIELD_H_
// clang-format on
