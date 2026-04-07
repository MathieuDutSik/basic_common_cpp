// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_file.h"

int main() {
  try {
    if (!IsProgramInPath("ls")) {
      std::cerr << "Expected ls to be available in PATH\n";
      throw TerminalException{1};
    }
    if (IsProgramInPath("lsx")) {
      std::cerr << "Expected lsx to be unavailable in PATH\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of Test_program_in_path\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_program_in_path\n";
    exit(e.eVal);
  }
}
