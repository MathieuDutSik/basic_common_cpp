#ifndef INCLUDE_COMBINATORICS_ELEM
#define INCLUDE_COMBINATORICS_ELEM

#include "Temp_common.h"

std::vector<int> BinomialStdvect_First(int const& k)
{
  std::vector<int> eVect(k);
  for (int i=0; i<k; i++)
    eVect[i]=i;
  return eVect;
}


bool BinomialStdvect_Increment(int const&n, int const&k, std::vector<int> & Tvect)
{
  Tvect[0]++;
  int xy2=1;
  while ((xy2 < k) && (Tvect[xy2-1] >= Tvect[xy2])) {
    Tvect[xy2]++;
    xy2++;
  }
  if (xy2 != 1) {
    for (int xy1=0; xy1<xy2-1; xy1++)
      Tvect[xy1]=xy1;
  }
  if (Tvect[k-1] == n)
    return false;
  return true;
}






#endif


