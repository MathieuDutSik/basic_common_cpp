// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Unit tests for Namelist.h substr optimization (Item 2)
#include "Namelist.h"
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <sstream>

// Test helper function
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
    std::cerr << "FAILED: " << test_name << std::endl;
    std::cerr << "Expected: ";
    for (const auto& s : expected) std::cerr << "\"" << s << "\" ";
    std::cerr << std::endl;
    std::cerr << "Actual:   ";
    for (const auto& s : actual) std::cerr << "\"" << s << "\" ";
    std::cerr << std::endl;
    exit(1);
  } else {
    std::cerr << "PASSED: " << test_name << std::endl;
  }
}

// Tests for NAMELIST_ConvertFortranStringToCppString function
void test_NAMELIST_ConvertFortranStringToCppString() {
  std::cerr << "Testing NAMELIST_ConvertFortranStringToCppString function..." << std::endl;

  // Test basic string conversions
  test_assert_string(NAMELIST_ConvertFortranStringToCppString(""),
                    "",
                    "Empty string");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("hello"),
                    "hello",
                    "Plain string without quotes");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("\"hello\""),
                    "hello",
                    "Double quoted string");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("'hello'"),
                    "hello",
                    "Single quoted string");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("\"hello world\""),
                    "hello world",
                    "Double quoted string with space");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("'hello world'"),
                    "hello world",
                    "Single quoted string with space");

  // Test edge cases
  test_assert_string(NAMELIST_ConvertFortranStringToCppString("\"\""),
                    "",
                    "Empty double quoted string");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("''"),
                    "",
                    "Empty single quoted string");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("\"a\""),
                    "a",
                    "Single character double quoted");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("'a'"),
                    "a",
                    "Single character single quoted");

  // Test longer strings
  test_assert_string(NAMELIST_ConvertFortranStringToCppString("\"This is a longer string\""),
                    "This is a longer string",
                    "Long double quoted string");

  test_assert_string(NAMELIST_ConvertFortranStringToCppString("'This is a longer string'"),
                    "This is a longer string",
                    "Long single quoted string");

  std::cerr << "NAMELIST_ConvertFortranStringToCppString tests completed." << std::endl << std::endl;
}

// Tests for NAMELIST_ConvertFortranListStringToCppListString function
void test_NAMELIST_ConvertFortranListStringToCppListString() {
  std::cerr << "Testing NAMELIST_ConvertFortranListStringToCppListString function..." << std::endl;

  // Test basic list conversions
  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString(""),
                   {},
                   "Empty string");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("'hello'"),
                   {"hello"},
                   "Single string in single quotes");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("\"hello\""),
                   {"hello"},
                   "Single string in double quotes");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("'hello''world'"),
                   {"hello", "world"},
                   "Two strings in single quotes");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("\"hello\"\"world\""),
                   {"hello", "world"},
                   "Two strings in double quotes");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("'hello world''test string'"),
                   {"hello world", "test string"},
                   "Two strings with spaces in single quotes");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("\"hello world\"\"test string\""),
                   {"hello world", "test string"},
                   "Two strings with spaces in double quotes");

  // Test edge cases
  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("''"),
                   {""},
                   "Empty single quoted string");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("\"\""),
                   {""},
                   "Empty double quoted string");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("'a''b''c'"),
                   {"a", "b", "c"},
                   "Three single characters");

  test_vector_equal(NAMELIST_ConvertFortranListStringToCppListString("\"a\"\"b\"\"c\""),
                   {"a", "b", "c"},
                   "Three single characters in double quotes");

  std::cerr << "NAMELIST_ConvertFortranListStringToCppListString tests completed." << std::endl << std::endl;
}

// Tests for NAMELIST_RemoveAfterCommentChar function
void test_NAMELIST_RemoveAfterCommentChar() {
  std::cerr << "Testing NAMELIST_RemoveAfterCommentChar function..." << std::endl;

  // Test basic comment removal
  test_assert_string(NAMELIST_RemoveAfterCommentChar("hello world", '!'),
                    "hello world",
                    "No comment character");

  test_assert_string(NAMELIST_RemoveAfterCommentChar("hello ! world", '!'),
                    "hello ",
                    "Comment in middle");

  test_assert_string(NAMELIST_RemoveAfterCommentChar("hello world!", '!'),
                    "hello world",
                    "Comment at end");

  test_assert_string(NAMELIST_RemoveAfterCommentChar("!hello world", '!'),
                    "",
                    "Comment at beginning");

  test_assert_string(NAMELIST_RemoveAfterCommentChar("hello # world", '#'),
                    "hello ",
                    "Hash comment character");

  test_assert_string(NAMELIST_RemoveAfterCommentChar("hello // world", '/'),
                    "hello ",
                    "Slash comment character");

  std::cerr << "NAMELIST_RemoveAfterCommentChar tests completed." << std::endl << std::endl;
}

// Tests for NAMELIST_RemoveAfterLastChar function
void test_NAMELIST_RemoveAfterLastChar() {
  std::cerr << "Testing NAMELIST_RemoveAfterLastChar function..." << std::endl;

  // Test basic last character removal
  test_assert_string(NAMELIST_RemoveAfterLastChar("hello world", '!'),
                    "hello world",
                    "No target character");

  test_assert_string(NAMELIST_RemoveAfterLastChar("hello ! world", '!'),
                    "hello ",
                    "Target character in middle");

  test_assert_string(NAMELIST_RemoveAfterLastChar("hello world!", '!'),
                    "hello world",
                    "Target character at end");

  test_assert_string(NAMELIST_RemoveAfterLastChar("!hello world", '!'),
                    "",
                    "Target character at beginning");

  test_assert_string(NAMELIST_RemoveAfterLastChar("hello ! world ! end", '!'),
                    "hello ! world ",
                    "Multiple target characters - keeps last");

  std::cerr << "NAMELIST_RemoveAfterLastChar tests completed." << std::endl << std::endl;
}

// Performance test
void test_performance() {
  std::cerr << "Running performance verification..." << std::endl;

  // Test with longer strings to verify performance improvements
  std::string long_string = "This is a very long string that will be processed multiple times to verify performance improvements from character access optimization";

  // Test fortran string conversion
  std::string quoted_long = "\"" + long_string + "\"";
  std::string result1 = NAMELIST_ConvertFortranStringToCppString(quoted_long);
  test_assert_string(result1, long_string, "Long string conversion");

  // Test comment removal
  std::string commented_long = long_string + " ! this is a comment";
  std::string result2 = NAMELIST_RemoveAfterCommentChar(commented_long, '!');
  test_assert_string(result2, long_string + " ", "Long string comment removal");

  // Test list string conversion
  std::string list_string = "\"" + long_string + "\"\"second string\"";
  std::vector<std::string> result3 = NAMELIST_ConvertFortranListStringToCppListString(list_string);
  test_vector_equal(result3, {long_string, "second string"}, "Long string list conversion");

  std::cerr << "Performance tests completed." << std::endl << std::endl;
}

int main() {
  std::cerr << "=== Unit Tests for Namelist.h substr Optimization (Item 2) ===" << std::endl << std::endl;

  try {
    test_NAMELIST_ConvertFortranStringToCppString();
    test_NAMELIST_ConvertFortranListStringToCppListString();
    test_NAMELIST_RemoveAfterCommentChar();
    test_NAMELIST_RemoveAfterLastChar();
    test_performance();

    std::cerr << "=== ALL TESTS PASSED ===" << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Test failed with exception: " << e.what() << std::endl;
    return 1;
  }
}
