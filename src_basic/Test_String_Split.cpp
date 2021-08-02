#include "Basic_string.h"

void test(const std::vector<std::string>& LStr1, const std::vector<std::string>& LStr2)
{
  if (LStr1.size() != LStr2.size()) {
    std::cerr << "Different lengths.\n";
    std::cerr << "|LStr1=" << LStr1.size() << " |LStr2|=" << LStr2.size() << "\n";
    exit(1);
  }
  for (size_t i=0; i<LStr1.size(); i++)
    if (LStr1[i] != LStr2[i]) {
      std::cerr << "Different strings at i=" << i << " str1=" << " str1=" << LStr1[i] << " str2=" << LStr2[i] << "\n";
      exit(1);
    }
}



int main()
{
  std::string eStrA, eStrB;
  //
  eStrA="A B C D", eStrB=" ";
  test(STRING_Split(eStrA, eStrB), {"A", "B", "C", "D"});
  //
  eStrA="Bonjour Mathieu et Maja", eStrB=" ";
  test(STRING_Split(eStrA, eStrB), {"Bonjour", "Mathieu", "et", "Maja"});
  //
  eStrA="A  B", eStrB=" ";
  test(STRING_Split(eStrA, eStrB), {"A", "B"});
  //
  eStrA="ABC", eStrB=" ";
  test(STRING_Split(eStrA, eStrB), {"ABC"});
  //
  eStrA=" ABC ", eStrB="ABC";
  test(STRING_Split(eStrA, eStrB), {" ", " "});
  //
  eStrA=" ABC ", eStrB=" ";
  test(STRING_Split(eStrA, eStrB), {"ABC"});
}
