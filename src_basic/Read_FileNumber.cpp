// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "basic_datafile.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "Read_FileNumber [FileNb]\n";
      std::cerr << "\n";
      std::cerr << "FileNb : the file written by the FileNumber type\n";
      throw TerminalException{1};
    }
    std::string FileNb = argv[1];
    size_t val = FileNumber_Read(FileNb);
    std::cout << "val=" << val << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
}
