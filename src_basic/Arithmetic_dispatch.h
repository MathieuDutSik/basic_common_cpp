// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_ARITHMETIC_DISPATCH_H_
#define SRC_BASIC_ARITHMETIC_DISPATCH_H_

// clang-format off
#include "ExceptionsFunc.h"
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  The table of the arithmetics a program accepts, each with the code that
  runs it. An optional arithmetic is added inside the #ifdef that compiles
  it, so the list printed in the usage text and in the error message is
  always the list of the arithmetics actually compiled in:

    ArithmeticDispatch dispatch;
    dispatch.add("mpq_class", [&]() { process<mpq_class>(...); });
  #ifdef ENABLE_FLINT_SUPPORT
    dispatch.add("fmpq_class", [&]() { process<fmpq_class>(...); });
  #endif
    ...
    std::cerr << "arith: " << dispatch.choices() << "\n";   // usage text
    ...
    dispatch.run(arith);   // throws, listing the choices, if not found
*/
struct ArithmeticDispatch {
  std::vector<std::pair<std::string, std::function<void()>>> entries;
  void add(std::string const &name, std::function<void()> const &f) {
    entries.emplace_back(name, f);
  }
  // The comma separated list of the accepted names.
  std::string choices() const {
    std::string str;
    for (auto &entry : entries) {
      if (!str.empty())
        str += ", ";
      str += entry.first;
    }
    return str;
  }
  void run(std::string const &name) const {
    for (auto &entry : entries) {
      if (entry.first == name) {
        entry.second();
        return;
      }
    }
    std::cerr << "Failed to find a matching entry for arith=" << name << "\n";
    std::cerr << "Allowed values: " << choices() << "\n";
    throw TerminalException{1};
  }
};

// clang-format off
#endif  // SRC_BASIC_ARITHMETIC_DISPATCH_H_
// clang-format on
