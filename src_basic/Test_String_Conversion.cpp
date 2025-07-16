// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Unit tests for string conversion function optimizations (Item 1)
#include "Basic_string.h"
#include <iostream>
#include <cassert>
#include <string>
#include <cmath>
#include <stdexcept>

// Test helper function for exact string comparison
void test_assert_string(const std::string& actual, const std::string& expected, const std::string& test_name) {
  if (actual != expected) {
    std::cerr << "FAILED: " << test_name << std::endl;
    std::cerr << "Expected: \"" << expected << "\"" << std::endl;
    std::cerr << "Actual:   \"" << actual << "\"" << std::endl;
    exit(1);
  } else {
    std::cerr << "PASSED: " << test_name << std::endl;
  }
}

// Test helper function for integer comparison
void test_assert_int(int actual, int expected, const std::string& test_name) {
  if (actual != expected) {
    std::cerr << "FAILED: " << test_name << std::endl;
    std::cerr << "Expected: " << expected << std::endl;
    std::cerr << "Actual:   " << actual << std::endl;
    exit(1);
  } else {
    std::cerr << "PASSED: " << test_name << std::endl;
  }
}

// Test helper function for approximate double comparison
void test_assert_double_approx(double actual, double expected, double tolerance, const std::string& test_name) {
  if (std::abs(actual - expected) > tolerance) {
    std::cerr << "FAILED: " << test_name << std::endl;
    std::cerr << "Expected: " << expected << " (±" << tolerance << ")" << std::endl;
    std::cerr << "Actual:   " << actual << std::endl;
    exit(1);
  } else {
    std::cerr << "PASSED: " << test_name << std::endl;
  }
}

// Tests for IntToString function
void test_IntToString() {
  std::cerr << "Testing IntToString function..." << std::endl;

  // Test basic integers
  test_assert_string(IntToString(0), "0", "IntToString: zero");
  test_assert_string(IntToString(1), "1", "IntToString: positive single digit");
  test_assert_string(IntToString(-1), "-1", "IntToString: negative single digit");
  test_assert_string(IntToString(123), "123", "IntToString: positive multi-digit");
  test_assert_string(IntToString(-123), "-123", "IntToString: negative multi-digit");

  // Test edge cases
  test_assert_string(IntToString(2147483647), "2147483647", "IntToString: max int");
  test_assert_string(IntToString(-2147483648), "-2147483648", "IntToString: min int");

  // Test powers of 10
  test_assert_string(IntToString(10), "10", "IntToString: 10");
  test_assert_string(IntToString(100), "100", "IntToString: 100");
  test_assert_string(IntToString(1000), "1000", "IntToString: 1000");
  test_assert_string(IntToString(-1000), "-1000", "IntToString: -1000");

  std::cerr << "IntToString tests completed." << std::endl << std::endl;
}

// Tests for StringToInt function
void test_StringToInt() {
  std::cerr << "Testing StringToInt function..." << std::endl;

  // Test basic conversions
  test_assert_int(StringToInt("0"), 0, "StringToInt: zero");
  test_assert_int(StringToInt("1"), 1, "StringToInt: positive single digit");
  test_assert_int(StringToInt("-1"), -1, "StringToInt: negative single digit");
  test_assert_int(StringToInt("123"), 123, "StringToInt: positive multi-digit");
  test_assert_int(StringToInt("-123"), -123, "StringToInt: negative multi-digit");

  // Test edge cases
  test_assert_int(StringToInt("2147483647"), 2147483647, "StringToInt: max int");
  test_assert_int(StringToInt("-2147483648"), -2147483648, "StringToInt: min int");

  // Test powers of 10
  test_assert_int(StringToInt("10"), 10, "StringToInt: 10");
  test_assert_int(StringToInt("100"), 100, "StringToInt: 100");
  test_assert_int(StringToInt("1000"), 1000, "StringToInt: 1000");
  test_assert_int(StringToInt("-1000"), -1000, "StringToInt: -1000");

  // Test leading zeros
  test_assert_int(StringToInt("001"), 1, "StringToInt: leading zeros");
  test_assert_int(StringToInt("0123"), 123, "StringToInt: multiple leading zeros");

  std::cerr << "StringToInt tests completed." << std::endl << std::endl;
}

// Tests for DoubleTo4dot2f function
void test_DoubleTo4dot2f() {
  std::cerr << "Testing DoubleTo4dot2f function..." << std::endl;

  // Test basic doubles
  test_assert_string(DoubleTo4dot2f(0.0), "0.00", "DoubleTo4dot2f: zero");
  test_assert_string(DoubleTo4dot2f(1.0), "1.00", "DoubleTo4dot2f: one");
  test_assert_string(DoubleTo4dot2f(-1.0), "-1.00", "DoubleTo4dot2f: negative one");
  test_assert_string(DoubleTo4dot2f(1.23), "1.23", "DoubleTo4dot2f: exact 2 decimals");
  test_assert_string(DoubleTo4dot2f(1.234), "1.23", "DoubleTo4dot2f: rounded down");
  test_assert_string(DoubleTo4dot2f(1.236), "1.24", "DoubleTo4dot2f: rounded up");

  // Test edge cases
  test_assert_string(DoubleTo4dot2f(123.456), "123.46", "DoubleTo4dot2f: large number");
  test_assert_string(DoubleTo4dot2f(-123.456), "-123.46", "DoubleTo4dot2f: negative large");
  test_assert_string(DoubleTo4dot2f(0.001), "0.00", "DoubleTo4dot2f: very small");
  test_assert_string(DoubleTo4dot2f(0.009), "0.01", "DoubleTo4dot2f: small rounded up");

  std::cerr << "DoubleTo4dot2f tests completed." << std::endl << std::endl;
}

// Tests for DoubleTo4dot1f function
void test_DoubleTo4dot1f() {
  std::cerr << "Testing DoubleTo4dot1f function..." << std::endl;

  // Test basic doubles
  test_assert_string(DoubleTo4dot1f(0.0), "0.0", "DoubleTo4dot1f: zero");
  test_assert_string(DoubleTo4dot1f(1.0), "1.0", "DoubleTo4dot1f: one");
  test_assert_string(DoubleTo4dot1f(-1.0), "-1.0", "DoubleTo4dot1f: negative one");
  test_assert_string(DoubleTo4dot1f(1.2), "1.2", "DoubleTo4dot1f: exact 1 decimal");
  test_assert_string(DoubleTo4dot1f(1.23), "1.2", "DoubleTo4dot1f: rounded down");
  test_assert_string(DoubleTo4dot1f(1.26), "1.3", "DoubleTo4dot1f: rounded up");

  // Test edge cases
  test_assert_string(DoubleTo4dot1f(123.456), "123.5", "DoubleTo4dot1f: large number");
  test_assert_string(DoubleTo4dot1f(-123.456), "-123.5", "DoubleTo4dot1f: negative large");
  test_assert_string(DoubleTo4dot1f(0.01), "0.0", "DoubleTo4dot1f: very small");
  test_assert_string(DoubleTo4dot1f(0.09), "0.1", "DoubleTo4dot1f: small rounded up");

  std::cerr << "DoubleTo4dot1f tests completed." << std::endl << std::endl;
}

// Tests for DoubleToString function
void test_DoubleToString() {
  std::cerr << "Testing DoubleToString function..." << std::endl;

  // Test basic doubles
  test_assert_string(DoubleToString(0.0), "0.000000000", "DoubleToString: zero");
  test_assert_string(DoubleToString(1.0), "1.000000000", "DoubleToString: one");
  test_assert_string(DoubleToString(-1.0), "-1.000000000", "DoubleToString: negative one");
  test_assert_string(DoubleToString(1.123456789), "1.123456789", "DoubleToString: exact 9 decimals");
  test_assert_string(DoubleToString(1.1234567890123), "1.123456789", "DoubleToString: truncated");

  // Test edge cases
  test_assert_string(DoubleToString(123.456), "123.456000000", "DoubleToString: large number");
  test_assert_string(DoubleToString(-123.456), "-123.456000000", "DoubleToString: negative large");
  test_assert_string(DoubleToString(0.000000001), "0.000000001", "DoubleToString: very small");

  std::cerr << "DoubleToString tests completed." << std::endl << std::endl;
}

// Round-trip tests (IntToString -> StringToInt)
void test_roundtrip_int() {
  std::cerr << "Testing integer round-trip conversions..." << std::endl;

  int test_values[] = {0, 1, -1, 123, -123, 1000, -1000, 2147483647, -2147483648};

  for (int original : test_values) {
    std::string str_val = IntToString(original);
    int converted_back = StringToInt(str_val);
    test_assert_int(converted_back, original, "Round-trip: " + std::to_string(original));
  }

  std::cerr << "Integer round-trip tests completed." << std::endl << std::endl;
}

// Performance comparison test
void test_performance() {
  std::cerr << "Running performance verification..." << std::endl;

  // Test IntToString performance with various integers
  int test_ints[] = {0, 1, -1, 123, -123, 1000000, -1000000};
  for (int val : test_ints) {
    std::string result = IntToString(val);
    std::cerr << "IntToString(" << val << ") = \"" << result << "\"" << std::endl;
  }

  // Test double conversion performance
  double test_doubles[] = {0.0, 1.0, -1.0, 123.456, -123.456};
  for (double val : test_doubles) {
    std::string result1 = DoubleTo4dot2f(val);
    std::string result2 = DoubleTo4dot1f(val);
    std::string result3 = DoubleToString(val);
    std::cerr << "Double conversions for " << val << ": "
              << result1 << ", " << result2 << ", " << result3 << std::endl;
  }

  std::cerr << "Performance tests completed." << std::endl << std::endl;
}

int main() {
  std::cerr << "=== Unit Tests for String Conversion Optimizations (Item 1) ===" << std::endl << std::endl;

  try {
    test_IntToString();
    test_StringToInt();
    test_DoubleTo4dot2f();
    test_DoubleTo4dot1f();
    test_DoubleToString();
    test_roundtrip_int();
    test_performance();

    std::cerr << "=== ALL TESTS PASSED ===" << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Test failed with exception: " << e.what() << std::endl;
    return 1;
  }
}
