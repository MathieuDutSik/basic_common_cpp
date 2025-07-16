// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Unit tests for string optimization functions (Correction 5)
#include "Basic_string.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <string>

// Test helper function
void test_assert(bool condition, const std::string& test_name) {
  if (!condition) {
    std::cerr << "FAILED: " << test_name << std::endl;
    exit(1);
  } else {
    std::cerr << "PASSED: " << test_name << std::endl;
  }
}

// Test helper for vector comparison
void test_vector_equal(const std::vector<std::string>& actual,
                      const std::vector<std::string>& expected,
                      const std::string& test_name) {
  bool equal = (actual.size() == expected.size());
  if (equal) {
    for (size_t i = 0; i < actual.size(); i++) {
      if (actual[i] != expected[i]) {
        equal = false;
        break;
      }
    }
  }

  if (!equal) {
    std::cerr << "FAILED: " << test_name << "\n";
    std::cerr << "Expected: ";
    for (const auto& s : expected) std::cerr << "\"" << s << "\" ";
    std::cerr << "\n";
    std::cerr << "Actual:   ";
    for (const auto& s : actual) std::cerr << "\"" << s << "\" ";
    std::cerr << "\n";
    exit(1);
  } else {
    std::cerr << "PASSED: " << test_name << "\n";
  }
}

// Tests for IsFullyNumeric function
void test_IsFullyNumeric() {
  std::cerr << "Testing IsFullyNumeric function..." << std::endl;

  // Test valid numeric strings
  test_assert(IsFullyNumeric("123"), "IsFullyNumeric: simple integer");
  test_assert(IsFullyNumeric("123.45"), "IsFullyNumeric: decimal number");
  test_assert(IsFullyNumeric("0"), "IsFullyNumeric: zero");
  test_assert(IsFullyNumeric("0.0"), "IsFullyNumeric: zero decimal");
  test_assert(IsFullyNumeric("999.999"), "IsFullyNumeric: large decimal");
  test_assert(IsFullyNumeric(" 123 "), "IsFullyNumeric: with spaces");
  test_assert(IsFullyNumeric("1 2 3"), "IsFullyNumeric: spaced digits");
  test_assert(IsFullyNumeric("   "), "IsFullyNumeric: only spaces");
  test_assert(IsFullyNumeric(""), "IsFullyNumeric: empty string");
  test_assert(IsFullyNumeric("123."), "IsFullyNumeric: trailing dot");
  test_assert(IsFullyNumeric(".123"), "IsFullyNumeric: leading dot");

  // Test invalid numeric strings
  test_assert(!IsFullyNumeric("1.2.3"), "IsFullyNumeric: multiple dots");
  test_assert(!IsFullyNumeric("abc"), "IsFullyNumeric: letters only");
  test_assert(!IsFullyNumeric("123abc"), "IsFullyNumeric: mixed alphanumeric");
  test_assert(!IsFullyNumeric("12.3.4.5"), "IsFullyNumeric: too many dots");
  test_assert(!IsFullyNumeric("12a3"), "IsFullyNumeric: letter in middle");
  test_assert(!IsFullyNumeric("12-3"), "IsFullyNumeric: minus sign");
  test_assert(!IsFullyNumeric("12+3"), "IsFullyNumeric: plus sign");
  test_assert(!IsFullyNumeric("12,3"), "IsFullyNumeric: comma");
  test_assert(!IsFullyNumeric("12$3"), "IsFullyNumeric: special character");
  test_assert(!IsFullyNumeric("12\n3"), "IsFullyNumeric: newline");
  test_assert(!IsFullyNumeric("12\t3"), "IsFullyNumeric: tab");

  // Edge cases
  test_assert(IsFullyNumeric("0000"), "IsFullyNumeric: leading zeros");
  test_assert(IsFullyNumeric("00.00"), "IsFullyNumeric: zeros with decimal");
  test_assert(!IsFullyNumeric("12e3"), "IsFullyNumeric: scientific notation");
  test_assert(!IsFullyNumeric("12E3"), "IsFullyNumeric: scientific notation uppercase");

  std::cerr << "IsFullyNumeric tests completed." << std::endl << std::endl;
}

// Tests for STRING_SplitCharNb function
void test_STRING_SplitCharNb() {
  std::cerr << "Testing STRING_SplitCharNb function..." << std::endl;

  // Test basic number/character splitting
  test_vector_equal(STRING_SplitCharNb("abc123def"),
                   {"abc", "123", "def"},
                   "STRING_SplitCharNb: basic mixed");

  test_vector_equal(STRING_SplitCharNb("123"),
                   {"123"},
                   "STRING_SplitCharNb: numbers only");

  test_vector_equal(STRING_SplitCharNb("abc"),
                   {"abc"},
                   "STRING_SplitCharNb: letters only");

  test_vector_equal(STRING_SplitCharNb(""),
                   {},
                   "STRING_SplitCharNb: empty string");

  test_vector_equal(STRING_SplitCharNb("a1b2c3"),
                   {"a", "1", "b", "2", "c", "3"},
                   "STRING_SplitCharNb: alternating");

  test_vector_equal(STRING_SplitCharNb("111aaa222bbb"),
                   {"111", "aaa", "222", "bbb"},
                   "STRING_SplitCharNb: groups");

  test_vector_equal(STRING_SplitCharNb("1"),
                   {"1"},
                   "STRING_SplitCharNb: single digit");

  test_vector_equal(STRING_SplitCharNb("a"),
                   {"a"},
                   "STRING_SplitCharNb: single letter");

  // Test with minus sign (should be treated as number)
  test_vector_equal(STRING_SplitCharNb("-123abc"),
                   {"-123", "abc"},
                   "STRING_SplitCharNb: negative number");

  test_vector_equal(STRING_SplitCharNb("abc-123"),
                   {"abc", "-123"},
                   "STRING_SplitCharNb: minus after letters");

  test_vector_equal(STRING_SplitCharNb("-"),
                   {"-"},
                   "STRING_SplitCharNb: minus only");

  test_vector_equal(STRING_SplitCharNb("--123"),
                   {"--123"},
                   "STRING_SplitCharNb: double minus");

  // Test consecutive numbers and letters
  test_vector_equal(STRING_SplitCharNb("123456789"),
                   {"123456789"},
                   "STRING_SplitCharNb: long number");

  test_vector_equal(STRING_SplitCharNb("abcdefghij"),
                   {"abcdefghij"},
                   "STRING_SplitCharNb: long letters");

  // Test edge cases
  test_vector_equal(STRING_SplitCharNb("0a0"),
                   {"0", "a", "0"},
                   "STRING_SplitCharNb: zeros separated");

  test_vector_equal(STRING_SplitCharNb("9z9"),
                   {"9", "z", "9"},
                   "STRING_SplitCharNb: nines separated");

  // Test with special characters (should be treated as letters)
  test_vector_equal(STRING_SplitCharNb("123.456"),
                   {"123", ".", "456"},
                   "STRING_SplitCharNb: decimal point");

  test_vector_equal(STRING_SplitCharNb("123 456"),
                   {"123", " ", "456"},
                   "STRING_SplitCharNb: space separator");

  test_vector_equal(STRING_SplitCharNb("123+456"),
                   {"123", "+", "456"},
                   "STRING_SplitCharNb: plus separator");

  std::cerr << "STRING_SplitCharNb tests completed." << std::endl << std::endl;
}

// Performance comparison test (informal)
void test_performance() {
  std::cerr << "Running basic performance verification..." << std::endl;

  // Create a large string to test performance
  std::string large_str;
  for (int i = 0; i < 1000; i++) {
    large_str += "123.456.789 abcdefghijk ";
  }

  // Test IsFullyNumeric on large string
  bool result1 = IsFullyNumeric(large_str);
  std::cerr << "IsFullyNumeric on large string: " << (result1 ? "true" : "false") << std::endl;

  // Test STRING_SplitCharNb on large string
  std::vector<std::string> result2 = STRING_SplitCharNb(large_str);
  std::cerr << "STRING_SplitCharNb on large string: " << result2.size() << " parts" << std::endl;

  std::cerr << "Performance tests completed." << std::endl << std::endl;
}

int main() {
  std::cerr << "=== Unit Tests for String Optimization (Correction 5) ===" << std::endl << std::endl;

  try {
    test_IsFullyNumeric();
    test_STRING_SplitCharNb();
    test_performance();

    std::cerr << "=== ALL TESTS PASSED ===" << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Test failed with exception: " << e.what() << std::endl;
    return 1;
  }
}
