#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
int main(int argc, char *argv[])
{
  //  using T=int;
  using T=long;
  try {
    int nb = 100;
    int siz = 10000;
    for (int i=0; i<nb; i++) {
      std::cerr << "i=" << i << "/" << nb << "\n";
      T a = rand() % (2*siz + 1) - siz;
      T b = rand() % (2*siz + 1) - siz;
      if (b != 0) {
        T res1 = ResInt_Generic<T>(a, b);
        T res2 = ResInt(a,b);
        if (res1 != res2) {
          std::cerr << "Inconsistency a=" << a << " b=" << b << " res1=" << res1 << " res2=" << res2 << "\n";
        }
        //
        T quo1 = QuoInt_Generic<T>(a, b);
        T quo2 = QuoInt(a, b);
        if (quo1 != quo2) {
          std::cerr << "Inconsistency a=" << a << " b=" << b << " quo1=" << quo1 << " quo2=" << quo2 << "\n";
        }
      }
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
