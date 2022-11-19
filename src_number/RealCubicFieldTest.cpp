// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryRealField.h"
#include "NumberTheory.h"


// The quantity is 2*cos(2*pi/7)
// The minimal polynomial is X^3 + X^2 - 2X - 1

int main(int argc, char *argv[]) {
  try {
    using T_rat = mpq_class;
    std::vector<T_rat> Pminimal{-1, 2, 1, 1};
    double val_double = 1.2469796037174672;
    std::vector<std::pair<T_rat,T_rat>> l_approx = {
      { 722/579 , 1/262144 },
      { 715371/573683 , 1/549755813888 },
      { 22572427/18101681 , 1/4503599627370496 },
      { 6957966595/5579855977 , 1/36893488147419103232 },
      { 137870231707/110563341450 , 1/302231454903657293676544 },
      { 44236337402433/35474788256806 , 1/9903520314283042199192993792 },
      { 78206434740138735/62716691200875701 , 1/5192296858534827628530496329220096 },
      { 1995859401875560391/1600554969724885645 , 1/10633823966279326983230456482242756608 },
      { 383121711408695944835/307239757784764293999 , 1/348449143727040986586495598010130648530944 }};
    HelperClassRealField<mpq_class> hcrf(Pminimal, val_double, l_approx);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
