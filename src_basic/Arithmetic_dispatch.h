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

  An entry can carry a description, printed by print_descriptions, and a
  parametrized entry is matched by its prefix, the rest of the name being
  passed to its function (e.g. "RealAlgebraic=" followed by a file name).
  The registration of the common families of arithmetics is in
  src_number/Arithmetic_dispatch_types.h.
*/
struct ArithmeticDispatch {
  struct Entry {
    std::string name;  // the name, or the prefix of a parametrized entry
    std::string param; // the name of the parameter, empty if none
    std::string description;
    std::function<void(std::string const &)> f;
  };
  std::vector<Entry> entries;
  void add(std::string const &name, std::function<void()> const &f) {
    add(name, "", f);
  }
  void add(std::string const &name, std::string const &description,
           std::function<void()> const &f) {
    entries.push_back(
        {name, "", description, [f](std::string const &) -> void { f(); }});
  }
  // A parametrized entry, accepted as prefix + value.
  void add_prefix(std::string const &prefix, std::string const &param,
                  std::string const &description,
                  std::function<void(std::string const &)> const &f) {
    entries.push_back({prefix, param, description, f});
  }
  // The comma separated list of the accepted names.
  std::string choices() const {
    std::string str;
    for (auto &entry : entries) {
      if (!str.empty())
        str += ", ";
      str += entry.name + entry.param;
    }
    return str;
  }
  // One line per accepted name with its description.
  void print_descriptions(std::ostream &os) const {
    for (auto &entry : entries) {
      os << entry.name << entry.param;
      if (!entry.description.empty())
        os << " : " << entry.description;
      os << "\n";
    }
  }
  void run(std::string const &name) const {
    for (auto &entry : entries) {
      if (entry.param.empty()) {
        if (entry.name == name) {
          entry.f("");
          return;
        }
      } else {
        if (name.size() > entry.name.size() &&
            name.compare(0, entry.name.size(), entry.name) == 0) {
          entry.f(name.substr(entry.name.size()));
          return;
        }
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
