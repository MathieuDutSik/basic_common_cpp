// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Boost_bitset_kernel.h"
#include "MAT_Matrix.h"

int main() {
  size_t n = 10;
  size_t n_face = 20;
  vectface vf1(n);
  for (size_t i_face = 0; i_face < n_face; i_face++) {
    Face f(n);
    for (size_t i = 0; i < n; i++)
      f[i] = random() % 2;
    vf1.push_back(f);
  }
  std::string filename = "/tmp/VF_filename.boost_archive";

  //
  // Writing the data
  //

  // save data to archive
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive
    oa << vf1;
    // archive and stream closed when destructors are called
  }

  auto get_vf = [&]() -> vectface {
    vectface f;
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> f;
    return f;
  };
  vectface vf2 = get_vf();
  if (vf1 != vf2) {
    std::cerr << "ERROR: vf1 and vf2 not matching\n";
  } else {
    std::cerr << "Normal termination of the program\n";
  }
}
