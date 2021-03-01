#include "COMB_Combinatorics.h"

void PrintVector(std::vector<int> const& V)
{
  std::cerr << "[";
  for (size_t i=0; i<V.size(); i++) {
    if (i>0)
      std::cerr << ",";
    std::cerr << " " << V[i];
  }
  std::cerr << " ]";
}


int main()
{
  BlockCppIterator BLK(5,2);
  // First algorithm
  std::cerr << "The first type of algorithm\n";
  BlockCppIterator::const_iterator iter = BLK.cbegin();
  while(iter != BLK.cend()) {
    PrintVector(*iter);
    std::cerr << "\n";
    iter++;
  }

  // Second algorithm
  std::cerr << "The second type of algorithm\n";
  for (auto const& eVect : BLK) {
    PrintVector(eVect);
    std::cerr << "\n";
  }
}

