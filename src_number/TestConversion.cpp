#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "TypeConversion.h"
#include "TypeConversionFinal.h"


int main(int argc, char *argv[])
{
  using T1=boost::multiprecision::cpp_int;
  using T2=mpz_class;

  T1 val1;
  T2 val2;
  /*
  std::string str1 = "43";
  std::istringstream is(str1);
  is >> val1;
  std::cerr << "val1=" << val1;
  */

  /*
  val1 = 43;
  TYPE_CONVERSION_STRING(val1, val2);
  std::cerr << "val1=" << val1 << " val2=" << val2 << "\n";
  */

  val1 = 43;
  val2 = UniversalScalarConversion<T2,T1>(val1);
  std::cerr << "val1=" << val1 << " val2=" << val2 << "\n";
}
