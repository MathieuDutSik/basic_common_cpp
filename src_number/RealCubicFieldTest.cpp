// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryRealField.h"
#include "NumberTheory.h"


// The quantity is 2*cos(2*pi/7)
// The minimal polynomial is X^3 + X^2 - 2X - 1

int main(int argc, char *argv[]) {
  try {
    using T_rat = mpq_class;
    std::string eFile = "../Examples/CubicField/CubicFieldDisc_49";
    HelperClassRealField<T_rat> hcrf(eFile);
    //    std::cerr << "hcrf.deg=" << hcrf.deg << " |ExprXdeg|=" << hcrf.ExprXdeg.size() << "\n";
    int const idx_discriminant_49 = 1;
    insert_helper_real_algebraic_field(idx_discriminant_49, hcrf);
    //    print_all_helpers(579);
    //
    using T = RealField<idx_discriminant_49>;
    //    print_all_helpers(1453);
    T x1 = 42;
    std::cerr << "x1=" << x1 << "\n";
    std::vector<T_rat> V2{1, 2, 4};
    T x2(V2);
    std::cerr << "x2=" << x2 << "\n";
    T x2inv = 1 / x2;
    //    print_all_helpers(1924);
    std::cerr << "x2inv=" << x2inv << "\n";
    //
    std::ostringstream convert;
    convert << x2inv;
    std::string strI = convert.str();
    T x2inv_cp;
    std::istringstream(strI) >> x2inv_cp;
    if (x2inv != x2inv_cp) {
      std::cerr << "Error in the conversion T -> string -> T\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
}
