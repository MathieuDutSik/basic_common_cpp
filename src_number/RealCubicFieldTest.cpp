// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "approximation.h"
// clang-format on

// The quantity is 2*cos(2*pi/7)
// The minimal polynomial is X^3 + X^2 - 2X - 1

int main() {
  try {
    using T_rat = mpq_class;
    std::string eFile = "../Examples/CubicField/CubicFieldDisc_49";
    HelperClassRealField<T_rat> hcrf(eFile);
    int const idx_discriminant_49 = 1;
    insert_helper_real_algebraic_field(idx_discriminant_49, hcrf);
    using T = RealField<idx_discriminant_49>;
    std::cerr << "STEP 1: Initial setup done\n";
    //
    T x1 = 42;
    std::cerr << "x1=" << x1 << "\n";
    std::vector<T_rat> V2{1, 2, 4};
    T x2(V2);
    std::vector<T_rat> V3{1, 3, 1};
    T x3(V3);
    std::cerr << "x2=" << x2 << "\n";
    T x2inv = -1 / x2;
    std::cerr << "x2inv=" << x2inv << "\n";
    T minus_x2 = -x2;
    std::cerr << "minus_x2=" << minus_x2 << "\n";
    T x2_3 = x2 * x3;
    std::cerr << "x2_3=" << x2_3 << "\n";
    std::cerr << "STEP 2: Product and quotient checked\n";
    //
    std::ostringstream convert;
    convert << x2inv;
    std::string strI = convert.str();
    std::cerr << "strI=" << strI << "X\n";
    T x2inv_cp;
    std::istringstream(strI) >> x2inv_cp;
    if (x2inv != x2inv_cp) {
      std::cerr << "Error in the conversion T -> string -> T\n";
      std::cerr << "INFO: x2inv=" << x2inv << "\n";
      std::cerr << "INFO: strI=" << strI << "\n";
      std::cerr << "INFO: x2inv_cp=" << x2inv_cp << "\n";
      throw TerminalException{1};
    }
    std::cerr << "STEP 3: conversion to string and back ckecked\n";
    //
    T thr = 1 / T(100000);
    std::cerr << "thr=" << thr << "\n";
    T approx = find_approximation_dichotomy(x2inv, thr);
    double x2inv_d = UniversalScalarConversion<double, T>(x2inv);
    double approx_d = UniversalScalarConversion<double, T>(approx);
    std::cerr << "x2inv_d=" << x2inv_d << " approx_d=" << approx_d << "\n";
    std::cerr << "STEP 4: comparison operations\n";
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
}
