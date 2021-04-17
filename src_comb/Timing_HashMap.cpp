#include "Boost_bitset.h"
#include "Boost_bitset.h"


/*
  The memory usage is shown via
  /usr/bin/time -v ./Timing_HashMap   |& grep resident

  The resulting memory usage is:
  ---unordered_map : 1035472
  ---sparse-map    : 773676
  ---robin-map     : 2624436
  ---hopscotch-map : 1634788
 */


#define UNORDERED_MAP
//#define TSL_SPARSE_MAP
//#define TSL_ROBIN_MAP
//#define TSL_HOPSCOTCH_MAP

#ifdef UNORDERED_MAP
# include <unordered_map>
# define MAP std::unordered_map
#endif

#ifdef TSL_SPARSE_MAP
# include "sparse_map.h"
# define MAP tsl::sparse_map
#endif

#ifdef TSL_ROBIN_MAP
# include "robin_map.h"
# define MAP tsl::robin_map
#endif

#ifdef TSL_HOPSCOTCH_MAP
# include "hopscotch_map.h"
# define MAP tsl::hopscotch_map
#endif




int main()
{
  int N = 10 * 1000 * 1000;
  int m=100;
  MAP<Face, int> eMap;
  //
  for (int i=0; i<N; i++) {
    Face f(m);
    for (int j=0; j<m; j++) {
      int rnd = rand() % 2;
      f[j] = rnd;
    }
    eMap[f] = i;
  }

}
