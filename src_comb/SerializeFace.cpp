// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Boost_bitset_kernel.h"
#include "MAT_Matrix.h"

int main() {

  size_t n = 20;
  Face f1(n);
  for (size_t i = 0; i < n; i++)
    f1[i] = random() % 2;
  std::string filename = "/tmp/Face_filename.boost_archive";

  //
  // Writing the data
  //

  // save data to archive
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive
    oa << f1;
    // archive and stream closed when destructors are called
  }

  auto get_vf = [&]() -> Face {
    Face f;
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> f;
    return f;
  };
  Face f2 = get_vf();
  bool is_equal = true;
  if (f1.size() != f2.size()) {
    std::cerr << "f1 and f2 have DIFFERENT length\n";
    is_equal = false;
  }
  for (size_t u = 0; u < f1.size(); u++)
    if (f1[u] != f2[u]) {
      std::cerr << "f1 and f2 differ at u=" << u << "\n";
      is_equal = false;
    }
  if (is_equal) {
    std::cerr << "f1 and f2 are equal\n";
  } else {
    std::cerr << "f1 and f2 are NOT equal\n";
  }
}
