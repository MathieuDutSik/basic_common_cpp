// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "basic_datafile.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "Read_FileFace [FileFace] [siz]\n";
      std::cerr << "or\n";
      std::cerr << "Read_FileFace [FileFace] [siz] [n_face]\n";
      std::cerr << "\n";
      std::cerr << "FileFace : the file written by the FileFace type\n";
      std::cerr << "siz      : the number of bits of a single face\n";
      std::cerr << "n_face   : the number of faces to read. If missing, all\n";
      std::cerr << "           the faces contained in the file are read\n";
      throw TerminalException{1};
    }
    std::string FileFc = argv[1];
    if (!IsExistingFile(FileFc)) {
      std::cerr << "The file FileFc=" << FileFc << " is missing\n";
      throw TerminalException{1};
    }
    size_t siz = ParseScalar<size_t>(argv[2]);
    if (siz == 0) {
      std::cerr << "The size siz of a face must be positive\n";
      throw TerminalException{1};
    }
    size_t n_bit = 8 * static_cast<size_t>(std::filesystem::file_size(FileFc));
    size_t n_face_max = n_bit / siz;
    size_t n_face = n_face_max;
    if (argc == 4) {
      n_face = ParseScalar<size_t>(argv[3]);
      if (n_face > n_face_max) {
        std::cerr << "The file FileFc=" << FileFc << " contains only "
                  << n_face_max << " faces of size siz=" << siz
                  << " but n_face=" << n_face << "\n";
        throw TerminalException{1};
      }
    }
    FileFace ff(FileFc, siz, n_face);
    std::vector<size_t> l_len(n_face);
    for (size_t i_face = 0; i_face < n_face; i_face++) {
      Face f = ff.getface(i_face);
      l_len[i_face] = f.count();
    }
    CollectedResult<size_t> rec = Collected(l_len);
    std::cout << "n_face=" << n_face << " siz=" << siz << "\n";
    for (size_t u = 0; u < rec.LVal.size(); u++) {
      std::cout << "len=" << rec.LVal[u] << " mult=" << rec.LMult[u] << "\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
}
