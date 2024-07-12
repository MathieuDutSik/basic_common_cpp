// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_INPUTOUTPUT_H_
#define SRC_NUMBER_INPUTOUTPUT_H_

// clang-format off
#include "Basic_functions.h"
#include "ExceptionsFunc.h"
#include <vector>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>
// clang-format on

template <typename T>
void WriteVectorFromRealAlgebraicString(std::ostream &os,
                                        std::vector<T> const &V) {
  bool DoSomething = false;
  bool is_first = true;
  for (size_t u = 0; u < V.size(); u++) {
    T const &val = V[u];
    if (val != 0) {
      DoSomething = true;
      if (u == 0) {
        os << val;
      } else {
        if (val > 0) {
          if (!is_first)
            os << "+";
          if (val != 1) {
            os << val << "*";
          }
        } else {
          if (val == -1) {
            os << "-";
          } else {
            os << val << "*";
          }
        }
        if (u == 1) {
          os << "x";
        }
        if (u > 1) {
          os << "x^" << u;
        }
      }
      is_first = false;
    }
  }
  if (!DoSomething) {
    os << "0";
  }
}

template <typename T>
std::vector<T> ReadVectorFromRealAlgebraicString(std::istream &is,
                                                 size_t const &deg) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  char c;
  std::string s;
  // First skipping the spaces
  while (true) {
    is.get(c);
    if (c != ' ' && c != '\n') {
      s += c;
      break;
    }
  }
  // Second reading the data till a space or end.
  while (true) {
    if (is.eof()) {
      break;
    }
    is.get(c);
    // If number of characters read is 0, then we have reached eof.
    if (is.gcount() == 0) {
      break;
    }
    if (c == ' ' || c == '\n') {
      break;
    }
    s += c;
  }
  // Now parsing the data
  std::vector<T> V(deg, 0);
  std::vector<size_t> W;
  for (size_t u = 1; u < s.size(); u++) {
    std::string echar = s.substr(u, 1);
    if (echar == "+" || echar == "-") {
      W.push_back(u);
    }
  }
  auto eval_expo = [](std::string const &s_expo) -> size_t {
    if (s_expo == "x")
      return 1;
    std::string s1 = s_expo.substr(0, 2);
    if (s1 != "x^") {
      std::cerr << "We should have an expression of the form x^n with n the "
                   "exponent\n";
      std::cerr << "s_expo=" << s_expo << " s1=" << s1 << "\n";
      throw TerminalException{1};
    }
    std::string s2 = s_expo.substr(2, s_expo.size() - 2);
    return ParseScalar<size_t>(s2);
  };
  auto get_position = [&](std::string const &s,
                          std::string const &echar) -> size_t {
    size_t len = s.size();
    for (size_t u = 0; u < len; u++) {
      std::string c = s.substr(u, 1);
      if (c == echar) {
        return u;
      }
    }
    return miss_val;
  };
  auto eval_scalar = [](std::string const &sc) -> T {
    if (sc == "") {
      return T(1);
    }
    if (sc == "-") {
      return T(-1);
    }
    if (sc == "+") {
      return T(1);
    }
    if (sc.substr(0, 1) == "+") {
      size_t s_len = sc.size();
      return ParseScalar<T>(sc.substr(1, s_len - 1));
    }
    return ParseScalar<T>(sc);
  };
  auto eval = [&](std::string const &sb) -> std::pair<T, size_t> {
    size_t lenb = sb.size();
    size_t pos_x = get_position(sb, "x");
    if (pos_x == miss_val) {
      return {eval_scalar(sb), 0};
    }
    size_t pos_mult = get_position(sb, "*");
    size_t pos_div = get_position(sb, "/");
    std::optional<T> div;
    if (pos_div != miss_val && pos_div > pos_x) {
      div = ParseScalar<T>(sb.substr(pos_div + 1, lenb - 1 - pos_div));
    }
    std::string sbred = sb.substr(0, pos_div);
    auto get_coeff = [&]() -> T {
      if (pos_mult == miss_val) {
        return eval_scalar(sbred.substr(0, pos_x));
      } else {
        if (pos_mult + 1 != pos_x) {
          std::cerr << "We expect the string to have a *x component\n";
          throw TerminalException{1};
        }
        return eval_scalar(sbred.substr(0, pos_mult));
      }
    };
    T coeff = get_coeff();
    if (div) {
      coeff /= *div;
    }
    size_t lenbred = sbred.size();
    size_t expo = eval_expo(sbred.substr(pos_x, lenbred - pos_x));
    return {coeff, expo};
  };
  for (size_t w = 0; w <= W.size(); w++) {
    size_t pos_first = 0, pos_last = s.size();
    if (w > 0) {
      pos_first = W[w - 1];
    }
    if (w < W.size()) {
      pos_last = W[w];
    }
    std::string sb = s.substr(pos_first, pos_last - pos_first);
    std::pair<T, size_t> ep = eval(sb);
    if (ep.second >= deg) {
      std::cerr << "sb=" << sb << "\n";
      std::cerr << "We found a term of the form x^{" << ep.second
                << "} while deg=" << deg << "\n";
      throw TerminalException{1};
    }
    V[ep.second] = ep.first;
  }
  return V;
}

// clang-format off
#endif  // SRC_NUMBER_INPUTOUTPUT_H_
// clang-format on
