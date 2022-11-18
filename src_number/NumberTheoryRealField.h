// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_REALFIELD_H_
#define SRC_NUMBER_REALFIELD_H_

#include "NumberTheory.h"
#include "Temp_common.h"

// For general real field, 


template<typename Tfield>
struct class HelperClassRealField {
private:
  using T = Tfield;
  HelperClassRealField(std::string const& eFile) {
    std::ifstream is(eFile);
    is >> deg;
    for (size_t u=0; u<=deg; u++) {
      T val;
      is >> val;
      Pminimal.push_back(val);
    }
    size_t n_approx;
    is >> n_approx;
    for (size_t i=0; i<n_approx; i++) {
      T val, err;
      is >> val;
      is >> err;
      SequenceApproximant.push_back({val,err});
    }
    // Compute a few stuff for inverse and products.
    ...
  }
  std::vector<T> FindQuotient(std::vector<T> const& num, std::vector<T> const& den) const {
  }
  std::vector<T> FindQuotient(T const& num, std::vector<T> const& den) const {
  }
  std::vector<T> FindInverse(std::vector<T> const& x) const {
    
  }
  std::vector<T> FindProduct(std::vector<T> const& a, std::vector<T> const& b) const {
  }
  bool IsStrictlyPositive(std::vector<T> const& x) const {
    auto 
  }
  double evaluate_as_double(std::vector<T> const& x) const {
  }
private:
  size_t deg;
  std::vector<T> Pminimal;
  std::vector<std::pair<T,T>> SequenceApproximant;
  MyMatrix<T> Amat;
};


std::map<int,HelperClassRealField<mpq_class>> list_helper;


void insert_helper(int i_field, std::string const& eFile)
{
  list_helper[i_field] = HelperClassRealField<mpq_class>(eFile);
}


template <int i_field> class RealField {
public:
  using T = mpq_class;
  using Tresidual = T;
private:
  std::vector<T> a;

public:
  std::vector<T>& get_seq() {
    return a;
  }
  const std::vector<T>& get_const_seq() const {
    return a;
  }

  // Note: We are putting "int" as argument here because we want to do the comparison with the
  // stuff like x > 0 or x = 1.
  // For the type "rational<T>" we had to forbid that because this lead to erroneous conversion
  // of say int64_t to int with catastrophic loss of precision.
  // But for the QuadField<T> the loss of precision does not occur because T is typically mpq_class.
  // or some other type that does not convert to integers easily. And at the same time the natural
  // conversion of int to int64_t allows the comparison x > 0 and equality set x = 1 to work despite
  // the lack of a operator=(int const& u)

  // Constructor
  RealField() {
    size_t deg = list_helper.at(i_field).deg;
    a = std::vector<T>(deg+1,0);
  }
  RealField(int const &u) {
    size_t deg = list_helper.at(i_field).deg;
    a = std::vector<T>(deg+1,0);
    a[0] = u;
  }
  RealField(T const &u) {
    size_t deg = list_helper.at(i_field).deg;
    a = std::vector<T>(deg+1,0);
    a[0] = u;
  }
  RealField(std::vector<T> const &_a) : a(_a) {}
  RealField(Realield<i_field> const &x) : a(x.a) {}
  //  QuadField<T,d>& operator=(QuadField<T,d> const&); // assignment operator
  //  QuadField<T,d>& operator=(T const&); // assignment operator from T
  //  QuadField<T,d>& operator=(int const&); // assignment operator from T
  // assignment operator from int
  RealField<i_field> operator=(int const &u) {
    for (auto &x : a) {
      x = 0;
    }
    a[0] = u;
    return *this;
  }
  // assignment operator
  RealField<i_field> operator=(RealField<i_field> const &x) {
    a = x.a;
    return *this;
  }
  //
  // Arithmetic operators below:
  void operator+=(RealField<i_field> const &x) {
    size_t len = a.size();
    for (size_t u=0; u<len; u++)
      a[u] += x.a[u];
  }
  void operator-=(RealField<i_field> const &x) {
    size_t len = a.size();
    for (size_t u=0; u<len; u++)
      a[u] -= x.a[u];
  }
  void operator/=(RealField<i_field> const &x) {
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    a = hcrf.FindQuotient(a, x.a);
  }
  friend RealField<i_field> operator+(RealField<i_field> const &x,
                                      RealField<T, d> const &y) {
    size_t len = x.a.size();
    std::vector<T> V(len);
    for (size_t u=0; u<len; u++)
      V[u] = x.a[u] + y.a[u];
    return RealField<i_field>(V);
  }
  friend RealField<i_field> operator-(RealField<i_field> const &x,
                                      RealField<i_field> const &y) {
    size_t len = x.a.size();
    std::vector<T> V(len);
    for (size_t u=0; u<len; u++)
      V[u] = x.a[u] - y.a[u];
    return RealField<i_field>(V);
  }
  friend RealField<i_field> operator-(RealField<i_field> const &x, int const &y) {
    std::vector<T> V = x.a;
    V[0] -= y;
    return RealField<i_field>(V);
  }
  friend RealField<i_field> operator-(RealField<i_field> const &x) {
    size_t len = x.a.size();
    std::vector<T> V(len);
    for (size_t u=0; u<len; u++)
      V[u] = -x.a[u];
    return RealField<i_field>(V);
  }
  friend RealField<i_field> operator/(int const &x, RealField<i_field> const &y) {
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    std::vector<T> V = hcrf.FindQuotient(x, y.a);
    return RealField<i_field>(V);
  }
  friend RealField<i_field> operator/(RealField<i_field> const &x,
                                      RealField<i_field> const &y) {
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    std::vector<T> V = hcrf.FindQuotient(x.a, y.a);
    return RealField<i_field>(V);
  }
  double get_d() const {
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    return hcrf.evaluate_as_double(a);
  }
  void operator*=(RealField<i_field> const &x) {
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    a = hcrf.ComputeProduct(a, x.a);
  }
  friend RealField<i_field> operator*(RealField<i_field> const &x,
                                      RealField<i_field> const &y) {
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    std::vector<T> V = gcrf.ComputeProduct(x.a, y.a);
    return RealField<i_field>(V);
  }
  friend RealField<i_field> operator*(int const &x, RealField<i_field> const &y) {
    size_t len = y.a.size();
    std::vector<T> V(len);
    for (size_t u=0; u<len; u++)
      V[u] = x * y.a[u];
    return RealField<i_field>(V);
  }
  friend std::ostream &operator<<(std::ostream &os, QuadField<T, d> const &v) {
    os << "(";
    size_t len = v.a.size();
    for (size_t u=0; u<len; u++) {
      if (u>0)
        os << ",";
      os << v.a[u];
    }
    os << ")";
    return os;
  }
  friend std::istream &operator>>(std::istream &is, QuadField<T, d> &v) {
    char c;
    std::string s;
    size_t miss_val = std::numeric_limits<size_t>::max();
    size_t pos = 0;
    // First skipping the spaces
    std::vector<size_t> ListPosComma;
    while(true) {
      is.get(c);
      if (c != ' ' && c != '\n') {
        s += c;
        pos++;
        break;
      }
    }
    // Second reading the data till a space or end.
    while(true) {
      if (is.eof()) {
        break;
      }
      is.get(c);
      if (c == ' ' || c == '\n') {
        break;
      }
      if (c == ',')
        ListPosComma.push_back(pos);
      s += c;
      pos++;
    }
    // Now parsing the data
    size_t len = ListPosComma.size() + 1;
    std::vector<T> V(len);
    for (size_t u=0; u<len; u++) {
      size_t pos_first, pos_last;
      if (u == 0) {
        
      }

      
    }

    
    std::string sA = s.substr(1, pos_comma);
    std::string sB = s.substr(pos_comma+1, s.size() - 2 - pos_comma);
    std::istringstream isA(sA);
    isA >> v.a;
    std::istringstream isB(sB);
    isB >> v.b;
    //    std::cerr << "v.a=" << v.a << "  v.b=" << v.b << "\n";
    //    std::cerr << "-------------------\n";
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
  friend bool operator!=(QuadField<T, d> const &x, int const &y) {
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
  friend bool operator>=(QuadField<T, d> const &x, int const &y) {
    QuadField<T, d> z;
    z = x - y;
    return IsNonNegative(z);
  }
  friend bool operator<=(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = y - x;
    return IsNonNegative(z);
  }
  friend bool operator<=(QuadField<T, d> const &x, int const &y) {
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
  friend bool operator>(QuadField<T, d> const &x, int const &y) {
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
  friend bool operator<(QuadField<T, d> const &x, int const &y) {
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
