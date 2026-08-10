// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "basic_datafile.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Read_FileBool [FileBool]\n";
      std::cerr << "or\n";
      std::cerr << "Read_FileBool [FileBool] [n_ent]\n";
      std::cerr << "\n";
      std::cerr << "FileBool : the file written by the FileBool type\n";
      std::cerr << "n_ent    : the number of bits to read. If missing, all\n";
      std::cerr << "           the bits of the file are read, that is 8\n";
      std::cerr << "           times the number of bytes of the file\n";
      throw TerminalException{1};
    }
    std::string FileBl = argv[1];
    if (!IsExistingFile(FileBl)) {
      std::cerr << "The file FileBl=" << FileBl << " is missing\n";
      throw TerminalException{1};
    }
    size_t n_byte = static_cast<size_t>(std::filesystem::file_size(FileBl));
    size_t n_ent_max = 8 * n_byte;
    size_t n_ent = n_ent_max;
    if (argc == 3) {
      n_ent = ParseScalar<size_t>(argv[2]);
      if (n_ent > n_ent_max) {
        std::cerr << "The file FileBl=" << FileBl << " contains only "
                  << n_ent_max << " bits but n_ent=" << n_ent << "\n";
        throw TerminalException{1};
      }
    }
    std::vector<uint8_t> l_status = FileBool_FullRead(FileBl, n_ent, std::cerr);
    size_t n_false = 0;
    size_t n_true = 0;
    for (auto &val : l_status) {
      if (val == 0) {
        n_false += 1;
      } else {
        n_true += 1;
      }
    }
    std::cout << "n_ent=" << n_ent << " n_false=" << n_false
              << " n_true=" << n_true << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
}
