#include "COMB_Combinatorics_elem.h"

template <typename Tidx> void PrintVector(std::vector<Tidx> const &V) {
  std::cerr << "[";
  for (size_t i = 0; i < V.size(); i++) {
    if (i > 0)
      std::cerr << ",";
    std::cerr << " " << V[i];
  }
  std::cerr << " ]";
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Program is used as\n";
    std::cerr << "Test_SetCppIterator [dim] [siz]\n";
    throw TerminalException{1};
  }
  int dim, size;
  sscanf(argv[1], "%d", &dim);
  sscanf(argv[2], "%d", &size);

  SetCppIterator SCI(dim, size);
  // First algorithm
  std::cerr << "The first type of algorithm\n";
  SetCppIterator::const_iterator iter = SCI.cbegin();
  while (iter != SCI.cend()) {
    PrintVector(*iter);
    std::cerr << "\n";
    iter++;
  }

  // Second algorithm
  std::cerr << "The second type of algorithm\n";
  for (auto const &eVect : SCI) {
    PrintVector(eVect);
    std::cerr << "\n";
  }
}
