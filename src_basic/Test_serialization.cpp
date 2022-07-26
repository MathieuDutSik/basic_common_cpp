// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Temp_common.h"

int main() {
  std::vector<uint8_t> eV;
  for (uint8_t u=0; u<10; u++)
    eV.push_back(u);
  //
  std::string filename = "/tmp/Archive_std_vector_uint8";
  // save data to archive
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive
    oa << eV;
    // archive and stream closed when destructors are called
  }
  // load data from archive
  std::vector<uint8_t> fV;
  {
    // create and open an archive for input
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    // read class state from archive
    ia >> fV;
    // archive and stream closed when destructors are called
  }
  //
  if (eV.size() != fV.size()) {
    std::cerr << "Error are different\n";
  }
  for (size_t u=0; u<eV.size(); u++)
    if (eV[u] != fV[u])
      std::cerr << "eV[u] != fV[u] and so different\n";

}
