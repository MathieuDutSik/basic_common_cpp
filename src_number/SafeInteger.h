#ifndef SRC_NUMBER_SAFEINTEGER_H_
#define SRC_NUMBER_SAFEINTEGER_H_


#define MAX_INT64_PROD 2147483647
#define MAX_INT64_SUM 4611686018427387903


void check_prod_int64(int64_t val) {
  if (val > MAX_INT64_PROD || val < -MAX_INT64_PROD) {
    std::cerr << "Error in product operation\n";
    throw TerminalException{1};
  }
}

void check_sum_int64(int64_t val) {
  if (val > MAX_INT64_SUM || val < -MAX_INT64_SUM) {
    std::cerr << "Error in sum operation\n";
    throw TerminalException{1};
  }
}


struct SafeInt64 {
private:
  int64_t val;
public:
  SafeInt64() : val(0) {}
  SafeInt64(int64_t const &x) : val(x) {}
  SafeInt64(SafeInt64 const &x) : val(x.val) {}
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
  }
  void operator-=(SafeInt64 const &x) {
    check_sum_int64(val);
    check_sum_int64(x.val);
    val -= x.val;
  }
  friend SafeInt64 operator+(SafeInt64 const &x, SafeInt64 const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x.val + y.val;
    return z;
  }
  friend SafeInt64 operator+(int64_t const &x, SafeInt64 const &y) {
    check_sum_int64(x);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x + y.val;
    return z;
  }
  friend SafeInt64 operator+(SafeInt64 const &x, int64_t const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y);
    SafeInt64 z;
    z.val = x.val + y;
    return z;
  }
  friend SafeInt64 operator-(SafeInt64 const &x, SafeInt64 const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x.val + y.val;
    return z;
  }
  friend SafeInt64 operator-(SafeInt64 const &x) {
    SafeInt64 z;
    z.val = -x.val;
    return z;
  }
  SafeInt64 operator++() {
    check_sum_int64(val);
    val++;
    return *this;
  }
  SafeInt64 operator++(int) {
    SafeInt64 tmp = *this;
    check_sum_int64(val);
    val++;
    return tmp;
  }
  SafeInt64 operator--() {
    check_sum_int64(val);
    val--;
    return *this;
  }
  SafeInt64 operator--(int) {
    SafeInt64 tmp = *this;
    check_sum_int64(val);
    val--;
    return tmp;
  }
  void operator*=(SafeInt64 const &x) {
    check_sum_int64(x.val);
    val *= x.val;
  }
  friend SafeInt64 operator*(SafeInt64 const &x,
                             SafeInt64 const &y) {
    check_prod_int64(x.val);
    check_prod_int64(y.val);
    SafeInt64 z;
    z.val = x.val * y.val;
    return z;
  }
  friend SafeInt64 operator*(int64_t const &x, SafeInt64 const &y) {
    check_prod_int64(x);
    check_prod_int64(y.val);
    SafeInt64 z;
    z.val = x * y.val;
    return z;
  }
  friend SafeInt64 operator*(SafeInt64 const &x, int64_t const &y) {
    check_prod_int64(x.val);
    check_prod_int64(y);
    SafeInt64 z;
    z.val = x.val * y;
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

SafeInt64 QuoInt(SafeInt64 const &a, SafeInt64 const &b) {
  const int64_t& a_int = a.get_const_val();
  const int64_t& b_int = b.get_const_val();
  int64_t quot = QuoInt_C_integer<int64_t>(a_int, b_int);
  return SafeInt64(quot);
}

void ResInt_Kernel(SafeInt64 const &a, SafeInt64 const &b, SafeInt64 & res) {
  const int64_t& a_int = a.get_const_val();
  const int64_t& b_int = b.get_const_val();
  res.get_val() = ResInt_C_integer<int64_t>(a_int, b_int);
}






#endif
