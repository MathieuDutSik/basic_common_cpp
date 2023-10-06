// clang-format off
#include "sparse_map.h"
#include "robin_map.h"
#include "hopscotch_map.h"
#include <unordered_map>
#include "Temp_common.h"
// clang-format on


template<typename T_cont>
void test_container() {
  T_cont map;
  for (int u=0; u<10; u++) {
    map[u] = u;
  }
  for (int u=0; u<10; u++) {
    map[u] += u;
  }
  for (auto& kv: map) {
    int k = kv.first;
    int v = kv.second;
    if (2*k != v) {
      std::cerr << "Inconsistency, k=" << k << " v=" << v << "\n";
      throw TerminalException{1};
    }
  }
}



int main() {
  test_container<std::unordered_map<int,int>>();
  test_container<tsl::sparse_map<int,int>>();
  test_container<tsl::robin_map<int,int>>();
  test_container<tsl::hopscotch_map<int,int>>();
}
