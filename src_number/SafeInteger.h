#ifndef SRC_NUMBER_SAFEINTEGER_H_
#define SRC_NUMBER_SAFEINTEGER_H_

#include "BasicNumberTypes.h"

#define MAX_INT64_PROD 2147483647
#define MAX_INT64_SUM  4611686018427387903


void check_prod_int64(int64_t val) {
  if (val > MAX_INT64_PROD || val < -MAX_INT64_PROD) {
    std::cerr << "Safety check triggerred for product operation of int64_t\n";
    std::cerr << "val=" << val << " MAX_INT64_PROD=" << MAX_INT64_PROD << "\n";
    throw SafeIntException{val};
  }
}

void check_sum_int64(int64_t val) {
  if (val > MAX_INT64_SUM || val < -MAX_INT64_SUM) {
    std::cerr << "Safety check triggerred for sum operation of int64_t\n";
    std::cerr << "val=" << val << " MAX_INT64_SUM=" << MAX_INT64_SUM << "\n";
    throw SafeIntException{val};
  }
}

void check_reasonableness(std::string oper, int64_t val) {
  //  int64_t limit = 1000;
  //  int64_t limit = MAX_INT64_PROD;
  int64_t limit = MAX_INT64_SUM;
  if (T_abs(val) > limit) {
    std::cerr << "We have val=" << val << " int(val)=" << int(val) << "\n";
    std::cerr << "which does not seem reasonable at oper=" << oper << "\n";
    throw TerminalException{1};
  }
}


struct SafeInt64 {
private:
  int64_t val;
public:
  SafeInt64() : val(0) {
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "Default constructor for SafeInt64 val=" << val << "\n";
#endif
  }
  SafeInt64(int64_t const &x) : val(x) {
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "Int64_t constructor for SafeInt64 val=" << x << "\n";
#endif
  }
  SafeInt64(SafeInt64 const &x) : val(x.val) {
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "SafeInt64 constructor for SafeInt64 val=" << x << "\n";
#endif
  }
  double get_d() const {
    return static_cast<double>(val);
  }
  // Assignment operators
  SafeInt64 operator=(SafeInt64 const &u) {
    // assignment operator from int
    val = u.val;
    return *this;
  }
  int64_t & get_val() { return val; }
  const int64_t & get_const_val() const { return val; }
  void operator+=(SafeInt64 const &x) {
    check_sum_int64(val);
    check_sum_int64(x.val);
    val += x.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+", val);
#endif
  }
  void operator-=(SafeInt64 const &x) {
    check_sum_int64(val);
    check_sum_int64(x.val);
    val -= x.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator-=", val);
#endif
  }
  friend SafeInt64 operator+(SafeInt64 const &x, SafeInt64 const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x.val + y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+(SafeInt64,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator+(int64_t const &x, SafeInt64 const &y) {
    check_sum_int64(x);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x + y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+(int64_t,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator+(SafeInt64 const &x, int64_t const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y);
    SafeInt64 z;
    z.val = x.val + y;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+(SafeInt64,int64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator-(SafeInt64 const &x, SafeInt64 const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x.val - y.val;
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "x.val=" << x.val << " y.val=" << y.val << "\n";
    check_reasonableness("operator-(SafeInt64,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator-(SafeInt64 const &x) {
    SafeInt64 z;
    z.val = -x.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator-(SafeInt64)", z.val);
#endif
    return z;
  }
  SafeInt64 operator++() {
    check_sum_int64(val);
    val++;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator++()", val);
#endif
    return *this;
  }
  SafeInt64 operator++(int) {
    SafeInt64 tmp = *this;
    check_sum_int64(val);
    val++;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator++(int)", val);
#endif
    return tmp;
  }
  SafeInt64 operator--() {
    check_sum_int64(val);
    val--;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator--()", val);
#endif
    return *this;
  }
  SafeInt64 operator--(int) {
    SafeInt64 tmp = *this;
    check_sum_int64(val);
    val--;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator--(int)", val);
#endif
    return tmp;
  }
  void operator*=(SafeInt64 const &x) {
    check_sum_int64(x.val);
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "val=" << val << "\n";
#endif
    val *= x.val;
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "x.val=" << x.val << "\n";
    check_reasonableness("operator*=()", val);
#endif
  }
  friend SafeInt64 operator*(SafeInt64 const &x,
                             SafeInt64 const &y) {
    check_prod_int64(x.val);
    check_prod_int64(y.val);
    SafeInt64 z;
    z.val = x.val * y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator*(SafeInt64_t,SafeInt64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator/(SafeInt64 const &x,
                             SafeInt64 const &y) {
    SafeInt64 z;
    z.val = x.val / y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator/(SafeInt64_t,SafeInt64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator%(SafeInt64 const &x,
                             SafeInt64 const &y) {
    SafeInt64 z;
    z.val = x.val % y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator/(SafeInt64_t,SafeInt64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator*(int64_t const &x, SafeInt64 const &y) {
    check_prod_int64(x);
    check_prod_int64(y.val);
    SafeInt64 z;
    z.val = x * y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator*(int64_t,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator*(SafeInt64 const &x, int64_t const &y) {
    check_prod_int64(x.val);
    check_prod_int64(y);
    SafeInt64 z;
    z.val = x.val * y;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator*(SafeInt64,int64_t)", z.val);
#endif
    return z;
  }
  friend std::ostream &operator<<(std::ostream &os, SafeInt64 const &v) {
    return os << v.val;
  }
  friend std::istream &operator>>(std::istream &is, SafeInt64 &v) {
    is >> v.val;
    return is;
  }
  friend bool operator==(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val == y.val;
  }
  friend bool operator!=(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val != y.val;
  }
  friend bool operator!=(SafeInt64 const &x, int64_t const &y) {
    return x.val != y;
  }
  friend bool operator>=(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val >= y.val;
  }
  friend bool operator<=(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val <= y.val;
  }
  friend bool operator<=(SafeInt64 const &x, int64_t const &y) {
    return x.val <= y;
  }
  friend bool operator>(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val > y.val;
  }
  friend bool operator<(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val < y.val;
  }
  friend bool operator<(SafeInt64 const &x, int64_t const &y) {
    return x.val < y;
  }
};





inline SafeInt64 GetDenominator_z([[maybe_unused]] SafeInt64 const &x) { return SafeInt64(1); }

inline void QUO_INT(stc<SafeInt64> const &a, stc<SafeInt64> const &b, SafeInt64 & q) {
  const int64_t& a_int = a.val.get_const_val();
  const int64_t& b_int = b.val.get_const_val();
  int64_t quot = QuoInt_C_integer<int64_t>(a_int, b_int);
  q = SafeInt64(quot);
}

#include "QuoIntFcts.h"

void ResInt_Kernel(SafeInt64 const &a, SafeInt64 const &b, SafeInt64 & res) {
  const int64_t& a_int = a.get_const_val();
  const int64_t& b_int = b.get_const_val();
  res.get_val() = ResInt_C_integer<int64_t>(a_int, b_int);
}






#endif
