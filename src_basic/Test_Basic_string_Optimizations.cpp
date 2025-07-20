// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Unit tests for Basic_string.h optimizations (substr, reserve, char access)
#include "Basic_string.h"
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <chrono>

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

// Test helper for boolean comparison
void test_assert_bool(bool actual, bool expected, const std::string& test_name) {
  if (actual != expected) {
    std::cerr << "FAILED: " << test_name << std::endl;
    std::cerr << "Expected: " << (expected ? "true" : "false") << std::endl;
    std::cerr << "Actual:   " << (actual ? "true" : "false") << std::endl;
    exit(1);
  } else {
    std::cerr << "PASSED: " << test_name << std::endl;
  }
}

// Test helper for integer comparison
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

// Tests for STRING_IsStringReduceToSpace function
void test_STRING_IsStringReduceToSpace() {
  std::cerr << "Testing STRING_IsStringReduceToSpace function..." << std::endl;

  // Test basic cases
  test_assert_bool(STRING_IsStringReduceToSpace(""), true, "Empty string");
  test_assert_bool(STRING_IsStringReduceToSpace(" "), true, "Single space");
  test_assert_bool(STRING_IsStringReduceToSpace("   "), true, "Multiple spaces");
  test_assert_bool(STRING_IsStringReduceToSpace("a"), false, "Single character");
  test_assert_bool(STRING_IsStringReduceToSpace("hello"), false, "Word");
  test_assert_bool(STRING_IsStringReduceToSpace(" a "), false, "Character with spaces");
  test_assert_bool(STRING_IsStringReduceToSpace("hello world"), false, "Multiple words");
  test_assert_bool(STRING_IsStringReduceToSpace("                    "), true, "Many spaces");

  std::cerr << "STRING_IsStringReduceToSpace tests completed." << std::endl << std::endl;
}

// Tests for STRING_GetCharPositionInString function
void test_STRING_GetCharPositionInString() {
  std::cerr << "Testing STRING_GetCharPositionInString function..." << std::endl;

  // Test basic cases
  test_assert_int(STRING_GetCharPositionInString("hello", 'h'), 0, "First character");
  test_assert_int(STRING_GetCharPositionInString("hello", 'o'), 4, "Last character");
  test_assert_int(STRING_GetCharPositionInString("hello", 'l'), 2, "First occurrence of repeated char");
  test_assert_int(STRING_GetCharPositionInString("hello", 'x'), -1, "Character not found");
  test_assert_int(STRING_GetCharPositionInString("", 'a'), -1, "Empty string");
  test_assert_int(STRING_GetCharPositionInString("a", 'a'), 0, "Single character match");
  test_assert_int(STRING_GetCharPositionInString("a", 'b'), -1, "Single character no match");
  test_assert_int(STRING_GetCharPositionInString("hello world", ' '), 5, "Space character");
  test_assert_int(STRING_GetCharPositionInString("hello world", 'w'), 6, "Character after space");

  std::cerr << "STRING_GetCharPositionInString tests completed." << std::endl << std::endl;
}

// Tests for StringSubstitution function
void test_StringSubstitution() {
  std::cerr << "Testing StringSubstitution function..." << std::endl;

  // Test basic substitution
  std::vector<std::pair<std::string, std::string>> subst1 = {{"name", "John"}};
  test_assert_string(StringSubstitution("Hello $name!", subst1), "Hello John!", "Basic substitution");

  // Test multiple substitutions
  std::vector<std::pair<std::string, std::string>> subst2 = {{"name", "John"}, {"age", "30"}};
  test_assert_string(StringSubstitution("$name is $age years old", subst2), "John is 30 years old", "Multiple substitutions");

  // Test no substitution
  test_assert_string(StringSubstitution("Hello world", subst1), "Hello world", "No substitution needed");

  // Test dollar without match
  test_assert_string(StringSubstitution("Price is $50", subst1), "Price is $50", "Dollar without match");

  // Test empty string
  test_assert_string(StringSubstitution("", subst1), "", "Empty string");

  std::cerr << "StringSubstitution tests completed." << std::endl << std::endl;
}

// Tests for StringVectorStringGAP function
void test_StringVectorStringGAP() {
  std::cerr << "Testing StringVectorStringGAP function..." << std::endl;

  // Test basic cases
  test_assert_string(StringVectorStringGAP({}), "[]", "Empty vector");
  test_assert_string(StringVectorStringGAP({"a"}), "[a]", "Single element");
  test_assert_string(StringVectorStringGAP({"a", "b"}), "[a,b]", "Two elements");
  test_assert_string(StringVectorStringGAP({"hello", "world", "test"}), "[hello,world,test]", "Multiple elements");
  test_assert_string(StringVectorStringGAP({"", "a", ""}), "[,a,]", "Empty strings");

  std::cerr << "StringVectorStringGAP tests completed." << std::endl << std::endl;
}

// Tests for StringNumber function
void test_StringNumber() {
  std::cerr << "Testing StringNumber function..." << std::endl;

  // Test basic cases
  test_assert_string(StringNumber(0, 1), "0", "Zero with 1 digit");
  test_assert_string(StringNumber(1, 1), "1", "One with 1 digit");
  test_assert_string(StringNumber(1, 3), "001", "One with 3 digits");
  test_assert_string(StringNumber(123, 3), "123", "Three digits exact");
  test_assert_string(StringNumber(45, 4), "0045", "Two digits with 4 padding");
  test_assert_string(StringNumber(999, 3), "999", "Max 3 digits");

  std::cerr << "StringNumber tests completed." << std::endl << std::endl;
}

// Tests for DropSpace function
void test_DropSpace() {
  std::cerr << "Testing DropSpace function..." << std::endl;

  // Test basic cases
  test_assert_string(DropSpace("hello world"), "helloworld", "Basic space removal");
  test_assert_string(DropSpace("  hello  world  "), "helloworld", "Leading/trailing spaces");
  test_assert_string(DropSpace("hello"), "hello", "No spaces");
  test_assert_string(DropSpace(""), "", "Empty string");
  test_assert_string(DropSpace("   "), "", "Only spaces");
  test_assert_string(DropSpace("a b c d"), "abcd", "Multiple single spaces");
  test_assert_string(DropSpace("a    b"), "ab", "Multiple consecutive spaces");

  std::cerr << "DropSpace tests completed." << std::endl << std::endl;
}

// Tests for STRING_Replace function
void test_STRING_Replace() {
  std::cerr << "Testing STRING_Replace function..." << std::endl;

  // Test basic replacement
  test_assert_string(STRING_Replace("hello world", " ", "_"), "hello_world", "Basic replacement");
  test_assert_string(STRING_Replace("hello world hello", "hello", "hi"), "hi world hi", "Multiple replacements");
  test_assert_string(STRING_Replace("hello", "world", "hi"), "hello", "No matches");
  test_assert_string(STRING_Replace("", "a", "b"), "", "Empty string");
  test_assert_string(STRING_Replace("aaa", "aa", "b"), "ba", "Overlapping matches");
  test_assert_string(STRING_Replace("hello", "hello", "world"), "world", "Full replacement");

  std::cerr << "STRING_Replace tests completed." << std::endl << std::endl;
}

// Tests for STRING_Split function
void test_STRING_Split() {
  std::cerr << "Testing STRING_Split function..." << std::endl;

  // Test basic splitting
  test_vector_equal(STRING_Split("a,b,c", ","), {"a", "b", "c"}, "Basic comma split");
  test_vector_equal(STRING_Split("hello world", " "), {"hello", "world"}, "Space split");
  test_vector_equal(STRING_Split("hello", ","), {"hello"}, "No separator");
  test_vector_equal(STRING_Split("", ","), {}, "Empty string");
  test_vector_equal(STRING_Split("a,,b", ","), {"a", "", "b"}, "Empty elements");
  test_vector_equal(STRING_Split("a::b::c", "::"), {"a", "b", "c"}, "Multi-character separator");

  std::cerr << "STRING_Split tests completed." << std::endl << std::endl;
}

// Performance test
void test_performance() {
  std::cerr << "Running performance verification..." << std::endl;

  // Test with longer strings to verify performance improvements
  std::string long_string = "This is a very long string that will be processed multiple times to verify performance improvements from character access optimization and string reserve calls";
  
  // Test STRING_IsStringReduceToSpace performance
  std::string space_string(1000, ' ');
  bool result1 = STRING_IsStringReduceToSpace(space_string);
  test_assert_bool(result1, true, "Long space string");

  // Test STRING_GetCharPositionInString performance
  std::string search_string = long_string + long_string + long_string;
  int pos = STRING_GetCharPositionInString(search_string, 'z');
  test_assert_int(pos, -1, "Long string character search");

  // Test StringSubstitution performance
  std::vector<std::pair<std::string, std::string>> subst = {{"long", "SHORT"}, {"string", "STR"}};
  std::string result2 = StringSubstitution(long_string, subst);
  std::cerr << "StringSubstitution result length: " << result2.length() << std::endl;

  // Test DropSpace performance
  std::string spaced_string = "a b c d e f g h i j k l m n o p q r s t u v w x y z";
  std::string result3 = DropSpace(spaced_string);
  test_assert_string(result3, "abcdefghijklmnopqrstuvwxyz", "Long spaced string");

  // Test STRING_Split performance
  std::vector<std::string> result4 = STRING_Split(spaced_string, " ");
  test_assert_int(static_cast<int>(result4.size()), 26, "Long string split");

  std::cerr << "Performance tests completed." << std::endl << std::endl;
}

int main() {
  std::cerr << "=== Unit Tests for Basic_string.h Optimizations ===\n" << std::endl;

  try {
    test_STRING_IsStringReduceToSpace();
    test_STRING_GetCharPositionInString();
    test_StringSubstitution();
    test_StringVectorStringGAP();
    test_StringNumber();
    test_DropSpace();
    test_STRING_Replace();
    test_STRING_Split();
    test_performance();

    std::cerr << "=== ALL TESTS PASSED ===\n" << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Test failed with exception: " << e.what() << std::endl;
    return 1;
  }
}
