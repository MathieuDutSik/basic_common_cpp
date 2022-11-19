// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_REALFIELD_H_
#define SRC_NUMBER_REALFIELD_H_

#include "NumberTheory.h"
#include "Temp_common.h"
#include "MAT_Matrix.h"

// For general real field.
// We need to use more sophisticated and slower algorithms than quadratic fields.
// (A) linear algebra is needed.
// (B) analysis is needed for deciding signs.

template<typename Tfield>
struct class HelperClassRealField {
private:
  using T = Tfield;
  HelperClassRealField(std::vector<T> const& _Pminimal, double const& _val_double, std::vector<std::pair<T,T>> const& l_approx) : Pminimal(_Pminimal), val_double(_val_double) {
    deg = PMinimal.size() - 1;
    for (auto & e_approx : l_approx) {
      T val = e_approx.first;
      T err= e_approx.second;
      if (val < 0) {
        std::cerr << "We require that the value val is positive. val=" << val << "\n";
        std::cerr << "This is arbitrary of us, but this is our choice\n";
        throw TerminalException{1};
      }
      if (err < 0) {
        std::cerr << "The error err should be positive. err=" << err << "\n";
        throw TerminalException{1};
      }
      if (err > val) {
        std::cerr << "The error err should be smaller than val. val=" << val << " err=" << err << "\n";
        throw TerminalException{1};
      }
      T val_low = val - err;
      T val_upp = val + err;
      std::vector<T> l_pow_low;
      std::vector<T> l_pow_low;
      T pow_low = val_low;
      T pow_upp = val_upp;
      for (int i=1; i<deg; i++) {
        l_pow_low.push_back(pow_low);
        l_pow_upp.push_back(pow_upp);
        pow_low *= val_low;
        pow_upp *= val_upp;
      }
      SequenceApproximant.push_back({l_pow_low,l_pow_upp});
    }
  }
  HelperClassRealField(std::string const& eFile) {
    std::ifstream is(eFile);
    is >> deg;
    std::vector<T> _Pminimal;
    for (size_t u=0; u<=deg; u++) {
      T val;
      is >> val;
      _Pminimal.push_back(val);
    }
    //
    double _val_double;
    is >> _val_double;
    //
    std::vector<std::pair<T,T>> l_approx;
    size_t n_approx;
    is >> n_approx;
    for (size_t u=0; u<n_approx; u++) {
      T val, err;
      is >> val;
      is >> err;
      l_approx.push_back({val,err});
    }
    HelperClassRealField(_Pminimal, _val_double, l_approx);
  }
private:
  void SetMatrix(MyMatrix<T> & M, std::vector<T> const& den) const {
    for (int i_col=0; i_col<deg; i_col++) {
      M(0,i_col) = den[i_col];
    }
    for (int i_row=1; i_row<deg; i_row++) {
      M(i_row,0) = 0;
      for (int i_col=1; i_col<deg; i_col++)
        M(i_row,i_col) = M(i_row-1,i_col-1);
      T val = M(i_row-1,deg-1);
      for (int i_col=0; i_col<deg; i_col++)
        M(i_row,i_col) += val * ExprXdeg[i_col];
    }
  }
  std::vector<T> GetSolution(MyMatrix<T> const& M, std::vector<T> const& num) {
    MyVector<T> w(deg);
    for (int i=0; i<deg; i++)
      w(i) = num[i];
    std::optional<MyVector<T>> opt = SolutionMat(M, w);
    if (!opt) {
      std::cerr << "Failed to solve the linear system\n";
      throw TerminalException{1};
    }
    MyVector<T> const& eSol = *opt;
    std::vector<T> V(deg);
    for (int u=0; u<deg; u++)
      V[u] = eSol(u);
    return V;
  }
  std::vector<T> FindQuotient(std::vector<T> const& num, std::vector<T> const& den) const {
    MyMatrix<T> M(deg,deg);
    SetMatrix(M, den);
    return GetSolution(M, num);
  }
  std::vector<T> FindQuotient(T const& num, std::vector<T> const& den) const {
    std::vector<T> num_V(deg,0);
    num_V[0] = num;
    MyMatrix<T> M(deg,deg);
    SetMatrix(M, den);
    return GetSolution(M, num_V);
  }
  std::vector<T> FindInverse(std::vector<T> const& x) const {
    std::vector<T> num_V(deg,0);
    num_V[0] = 1;
    MyMatrix<T> M(deg,deg);
    SetMatrix(M, den);
    return GetSolution(M, num_V);
  }
  std::vector<T> ComputeProduct(std::vector<T> const& a, std::vector<T> const& b) const {
    MyMatrix<T> M(deg,deg);
    SetMatrix(M, b);
    std::vector<T> V(deg, 0);
    for (int i=0; i<deg; i++) {
      for (int j=0; j<deg; j++)
        V[j] += a[i] * M(i,j);
    }
    return V;
  }
  bool IsStrictlyPositive(std::vector<T> const& x) const {
    // x is assumed to be non-zero
    auto get_bounds=[&](std::pair<std::vector<T>, std::vector<T>> const& epair) -> std::pair<T,T> {
      T val_low = x[0];
      T val_upp = x[0];
      for (int i=1; i<deg; i++) {
        T val = x[i];
        if (val > 0) {
          val_low += val * epair.first[i-1];
          val_upp += val * epair.second[i-1];
        }
        if (val < 0) {
          val_low += val * epair.second[i-1];
          val_upp += val * epair.first[i-1];
        }
      }
      return {val_low,val_upp};
    };
    for (auto & epair : SequenceApproximant) {
      auto pair_bound = get_bounds(epair);
      T const& val_low = pair_bound.first;
      T const& val_upp = pair_bound.second;
      if (val_upp <= 0) {
        return false;
      }
      if (val_low >= 0) {
        return true;
      }
    }
    std::cerr << "Failed to find an approximant that allows to conclude, please produce better approximants\n";
    throw TerminalException{1};
  }
  double evaluate_as_double(std::vector<T> const& x) const {
    double ret_val = 0;
    double pow_double = 1.0;
    for (int i=0; i<deg; i++) {
      double coeff = UniversalScalarConversion<double,T>(x[i]);
      ret_val += coeff * pow_double;
      pow_double *= val_double;
    }
    return ret_val;
  }
private:
  size_t deg;
  std::vector<T> Pminimal;
  std::vector<T> ExprXdeg;
  double val_double;
  std::vector<std::pair<std::vector<T>,std::vector<T>>> SequenceApproximant;
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
  bool IsZero(std::vector<T> const& V) const {
    for (auto & val : V)
      if (val != 0)
        return false;
    return true;
  }
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
    std::vector<T> V = hcrf.ComputeProduct(x.a, y.a);
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
        pos_first = 1;
      } else {
        pos_first = ListPosComma[u-1]+1;
      }
      if (u == len-1) {
        pos_last = s.size() - 1;
      } else {
        pos_last = ListPosComma[u];
      }
      std::string sRed = s.substr(pos_first, pos_last - pos_frist);
      T val;
      std::istringstream isRed(sRed);
      isRed >> val;
      V.push_back(val);
    }
    return RealField<i_field>(V);
  }
  friend bool operator==(RealField<i_field> const &x, RealField<i_field> const &y) {
    size_t deg = x.a.size();
    for (size_t u=0; u<deg; u++)
      if (x.a[u] != y.a[u])
        return false;
    return true;
  }
  friend bool operator!=(RealField<i_field> const &x, RealField<i_field> const &y) {
    size_t deg = x.a.size();
    for (size_t u=0; u<deg; u++)
      if (x.a[u] != y.a[u])
        return true;
    return false;
  }
  friend bool operator!=(RealField<i_field> const &x, int const &y) {
    size_t deg = x.a.size();
    for (size_t u=1; u<deg; u++)
      if (x.a[u] != 0)
        return true;
    return x.a[0] != y;
  }
  friend bool IsNonNegative(RealField<i_field> const &x) {
    if (IsZero(x.a))
      return true;
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    return hcrf.IsStrictlyPositive(x.a);
  }
  friend bool operator>=(RealField<i_field> const &x, RealField<i_field> const &y) {
    size_t deg = x.a.size();
    std::vector<T> V(len);
    for (size_t u=0; u<deg; u++) {
      T val = x.a[u] = y.a[u];
      V[u] = val;
    }
    if (IsZero(V))
      return true;
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    return hcrf.IsStrictlyPositive(V);
  }
  friend bool operator>=(RealField<i_field> const &x, int const &y) {
    std::vector<T> V = x.a;
    V[0] -= y;
    return IsNonNegative(RealField<i_field>(V));
  }
  friend bool operator<=(RealField<i_field> const &x, RealField<i_field> const &y) {
    RealField<i_field> z = y - x;
    return IsNonNegative(z);
  }
  friend bool operator<=(RealField<i_field> const &x, int const &y) {
    RealField<i_field> z = y - x;
    return IsNonNegative(z);
  }
  friend bool operator>(RealField<i_field> const &x, RealField<i_field> const &y) {
    RealField<i_field> z = x - y;
    if (IsZero(z.a))
      return false;
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    return hcrf.IsStrictlyPositive(z.a);
  }
  friend bool operator>(RealField<i_field> const &x, int const &y) {
    RealField<T, d> z = x - y;
    if (IsZero(z.a))
      return false;
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    return hcrf.IsStrictlyPositive(z.a);
  }
  friend bool operator<(RealField<i_field> const &x, RealField<i_field> const &y) {
    RealField<i_field> z = y - x;
    if (IsZero(z.a))
      return false;
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    return hcrf.IsStrictlyPositive(z.a);
  }
  friend bool operator<(RealField<i_field> const &x, int const &y) {
    RealField<i_field> z = y - x;
    if (IsZero(z.a))
      return false;
    HelperClassRealField<T> const& hcrf = list_helper.at(i_field);
    return hcrf.IsStrictlyPositive(z.a);
  }
};

// For this construction we cannot hope to handle rings and fields nicely

template<int i_field>
struct overlying_field<RealField<i_field>> { typedef RealField<i_field> field_type; };

// Note that the underlying ring is not unique, there are many possibiliies actually
// but we can represent only one in our scheme.
template<int i_field>
struct underlying_ring<RealField<i_field>> { typedef RealField<i_field> ring_type; };


template <int i_field>
inline void TYPE_CONVERSION(stc<RealField<i_field>> const &eQ, double &eD) {
  eD = eQ.val.get_d();
}

template <int i_field> struct is_totally_ordered<RealField<i_field>> {
  static const bool value = true;
};

template <int i_field> struct is_ring_field<RealField<i_field>> {
  static const bool value = true;
};

template <int i_field> struct is_exact_arithmetic<RealField<i_field>> {
  static const bool value = true;
};

// Hashing function

template <int i_field>
struct is_implementation_of_Z<RealField<i_field>> {
  static const bool value = false;
};

template <int i_field>
struct is_implementation_of_Q<RealField<i_field>> {
  static const bool value = false;
};

// Hashing function

namespace std {
  template <int i_field> struct hash<RealField<i_field>> {
    std::size_t operator()(const RealField<i_field> &x) const {
      auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      	seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      size_t seed = 0x9e2479b9;
      std::vector<mpq_class> V = x.get_const_seq();
      for (auto & val : V) {
        size_t e_hash = std::hash<mpq_class>()(val);
        combine_hash(seed, e_hash);
      }
      return seed;
    }
  };
}


// Local typing info

template<typename T> struct is_real_algebraic_field {};

template<int i_field> struct is_real_algebraic_field<RealField<i_field>> { static const bool value = true; };

// Some functionality

template<int i_field>
bool IsInteger(RealField<i_field> const& x) {
  std::vector<mpq_class> V = x.get_const_seq();
  size_t len = V.size();
  for (size_t u=1; u<len; u++)
    if (V[u] != 0)
      return false;
  return IsInteger(V[0]);
}

// The conversion tools (int)

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

template <class Archive, int i_field>
inline void serialize(Archive &ar, RealField<i_field> &val,
                      [[maybe_unused]] const unsigned int version) {
  std::vector<T> & V = val.get_seq();
  for (auto & val : V)
    ar &make_nvp("quadfield_a", val);
}

}

// Turning into something rational

template<typename Tring, typename Tquad>
void ScalingInteger_Kernel(stc<Tquad> const& x, Tring& x_res) {
  using Tfield = typename Tquad::Tresidual;
  std::vector<Tring> V;
  for (auto & val : x.val.get_const_seq())
    V.push_back(GetDenominator_z(val));
  x_res = LCMlist(V);
}


// clang-format off
#endif  // SRC_NUMBER_REALFIELD_H_
// clang-format on
