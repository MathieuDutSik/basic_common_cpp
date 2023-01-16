// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_file.h"

int main() {
  std::vector<std::string> l_prog = {"ppl_lcdd", "cat", "normaliz"};
  for (auto &e_prog : l_prog) {
    bool test = IsProgramInPath(e_prog);
    std::cerr << "e_prog=" << e_prog << " test=" << test << "\n";
  }
}
