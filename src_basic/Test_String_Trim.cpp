// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Unit tests for STRING_RemoveSpacesBeginningEnd optimization (Item 4)
#include "Basic_string.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <string>

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

// Tests for STRING_RemoveSpacesBeginningEnd function
void test_STRING_RemoveSpacesBeginningEnd() {
  std::cerr << "Testing STRING_RemoveSpacesBeginningEnd function..." << std::endl;

  // Test basic trimming
  test_assert_string(STRING_RemoveSpacesBeginningEnd("hello"),
                    "hello",
                    "No spaces to remove");

  test_assert_string(STRING_RemoveSpacesBeginningEnd(" hello"),
                    "hello",
                    "Leading space only");

  test_assert_string(STRING_RemoveSpacesBeginningEnd("hello "),
                    "hello",
                    "Trailing space only");

  test_assert_string(STRING_RemoveSpacesBeginningEnd(" hello "),
                    "hello",
                    "Both leading and trailing spaces");

  // Test multiple spaces
  test_assert_string(STRING_RemoveSpacesBeginningEnd("  hello  "),
                    "hello",
                    "Multiple leading and trailing spaces");

  test_assert_string(STRING_RemoveSpacesBeginningEnd("   hello   world   "),
                    "hello   world",
                    "Multiple spaces, preserve internal spaces");

  // Test edge cases
  test_assert_string(STRING_RemoveSpacesBeginningEnd(""),
                    "",
                    "Empty string");

  test_assert_string(STRING_RemoveSpacesBeginningEnd(" "),
                    "",
                    "Single space");

  test_assert_string(STRING_RemoveSpacesBeginningEnd("   "),
                    "",
                    "Multiple spaces only");

  test_assert_string(STRING_RemoveSpacesBeginningEnd("a"),
                    "a",
                    "Single character");

  test_assert_string(STRING_RemoveSpacesBeginningEnd(" a "),
                    "a",
                    "Single character with spaces");

  // Test with internal spaces preserved
  test_assert_string(STRING_RemoveSpacesBeginningEnd("hello world"),
                    "hello world",
                    "Internal space preserved");

  test_assert_string(STRING_RemoveSpacesBeginningEnd(" hello world "),
                    "hello world",
                    "Internal space preserved with trimming");

  test_assert_string(STRING_RemoveSpacesBeginningEnd("  a  b  c  "),
                    "a  b  c",
                    "Multiple internal spaces preserved");

  // Test with longer strings
  test_assert_string(STRING_RemoveSpacesBeginningEnd("   This is a longer string with multiple words   "),
                    "This is a longer string with multiple words",
                    "Long string with multiple words");

  // Test with special characters
  test_assert_string(STRING_RemoveSpacesBeginningEnd(" !@#$%^&*() "),
                    "!@#$%^&*()",
                    "Special characters");

  test_assert_string(STRING_RemoveSpacesBeginningEnd("  123 456 789  "),
                    "123 456 789",
                    "Numbers with spaces");

  // Test extreme cases
  std::string long_spaces(100, ' ');
  test_assert_string(STRING_RemoveSpacesBeginningEnd(long_spaces),
                    "",
                    "Long string of spaces only");

  std::string long_string = long_spaces + "content" + long_spaces;
  test_assert_string(STRING_RemoveSpacesBeginningEnd(long_string),
                    "content",
                    "Long spaces around content");

  std::cerr << "STRING_RemoveSpacesBeginningEnd tests completed." << std::endl << std::endl;
}

// Performance comparison test
void test_performance() {
  std::cerr << "Running performance verification..." << std::endl;

  // Create strings with various space patterns
  std::string test_strings[] = {
    "   short   ",
    "                    medium length string                    ",
    std::string(1000, ' ') + "long string with many spaces" + std::string(1000, ' '),
    "no spaces",
    "     ",
    ""
  };

  for (const auto& test_str : test_strings) {
    std::string result = STRING_RemoveSpacesBeginningEnd(test_str);
    std::cerr << "Input length: " << test_str.length()
              << ", Output length: " << result.length() << std::endl;
  }

  std::cerr << "Performance tests completed." << std::endl << std::endl;
}

int main() {
  std::cerr << "=== Unit Tests for STRING_RemoveSpacesBeginningEnd Optimization (Item 4) ===" << std::endl << std::endl;

  try {
    test_STRING_RemoveSpacesBeginningEnd();
    test_performance();

    std::cerr << "=== ALL TESTS PASSED ===" << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Test failed with exception: " << e.what() << std::endl;
    return 1;
  }
}
