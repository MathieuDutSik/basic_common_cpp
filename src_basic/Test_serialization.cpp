// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Temp_common.h"

struct DataTest {
  std::vector<uint8_t> V;
  std::unordered_set<int> set;
  std::unordered_map<int, int> map;
};

namespace boost::serialization {
  template <class Archive>
  inline void serialize(Archive &ar, DataTest &eRec,
                        [[maybe_unused]] const unsigned int version) {
    ar &make_nvp("V", eRec.V);
    ar &make_nvp("set", eRec.set);
    ar &make_nvp("map", eRec.map);
  }
}

bool operator==(DataTest const &x, DataTest const &y) {
  if (x.V != y.V) {
    return false;
  }
  if (x.set != y.set) {
    return false;
  }
  if (x.map != y.map) {
    return false;
  }
  return true;
}

bool operator!=(DataTest const &x, DataTest const &y) {
  return !(x == y);
}

DataTest GenerateEntry() {
  std::vector<uint8_t> eV;
  for (uint8_t u = 0; u < 10; u++)
    eV.push_back(u);
  std::unordered_set<int> set;
  std::unordered_map<int, int> map;
  for (int u = 0; u < 10; u++) {
    int val1 = 2*u + 1;
    int val2 = 3*u + 4;
    set.insert(val1);
    map.insert(std::make_pair(val1, val2));
  }
  return {eV, set, map};
}


int main() {
  DataTest test = GenerateEntry();
  //
  std::string filename = "/tmp/Archive_std_vector_uint8";
  // save data to archive
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive
    oa << test;
    // archive and stream closed when destructors are called
  }
  // load data from archive
  DataTest test_read;
  {
    // create and open an archive for input
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    // read class state from archive
    ia >> test_read;
    // archive and stream closed when destructors are called
  }
  //
  if (test != test_read) {
    std::cerr << "Error are different\n";
  }
}
