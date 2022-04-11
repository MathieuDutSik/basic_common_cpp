#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "TypeConversion.h"

int main(int argc, char *argv[]) {
  //  using T=int;
  //  using T=long;
  //  using T=long;
  using T = boost::multiprecision::cpp_int;
  //  using T=mpq_class;
  using Tint = mpz_class;
  try {
    size_t n_error = 0;
    auto TestCons = [&](int a, int b) -> void {
      if (b != 0) {
        T a_T = UniversalScalarConversion<T, int>(a);
        T b_T = UniversalScalarConversion<T, int>(b);
        Tint a_cont = UniversalScalarConversion<Tint, int>(a);
        Tint b_cont = UniversalScalarConversion<Tint, int>(b);
        //
        T res_T = ResInt(a_T, b_T);
        int res1 = UniversalScalarConversion<int, T>(res_T);
        //
        Tint res_cont = ResInt(a_cont, b_cont);
        int res2 = UniversalScalarConversion<int, Tint>(res_cont);
        if (res1 != res2) {
          std::cerr << "Inconsistency a=" << a << " b=" << b << " res1=" << res1
                    << " res2=" << res2 << "\n";
          n_error++;
        }
        //
        T quo_T = QuoInt(a_T, b_T);
        int quo1 = UniversalScalarConversion<int, T>(quo_T);
        //
        Tint quo_cont = QuoInt(a_cont, b_cont);
        int quo2 = UniversalScalarConversion<int, Tint>(quo_cont);
        if (quo1 != quo2) {
          std::cerr << "Inconsistency a=" << a << " b=" << b << " quo1=" << quo1
                    << " quo2=" << quo2 << "\n";
          n_error++;
        }
      }
    };
    int nb = 100;
    int siz = 10000;
    for (int i = 0; i < nb; i++) {
      std::cerr << "i=" << i << "/" << nb << "\n";
      int a = rand() % (2 * siz + 1) - siz;
      int b = rand() % (2 * siz + 1) - siz;
      TestCons(a, b);
    }
    for (int a = -10; a < 10; a++)
      for (int b = -10; b < 10; b++)
        TestCons(a, b);
    std::cerr << "n_error=" << n_error << "\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
