#ifndef FACTORIZATIONS_INCLUDE
#define FACTORIZATIONS_INCLUDE

#include <vector>
#include <map>



template<typename T>
std::pair<bool, T> rho_pollard_factorize(T const& number)
{
  T count;
  T x_fixed = 2, x = 2, size = 2, factor, diff;
  do {
    count = size;
    do {
      x = ResInt( x * x + 1, number);
      diff = x - x_fixed;
      if (diff < 0)
        diff = -diff;
      factor = gcd(diff, number);
    } while (--count && (factor == 1));
    size *= 2;
    x_fixed = x;
  } while (factor == 1);
  if (factor == number) {
    return {false, -1};
  } else {
    return {true, factor};
  }
}


template<typename T>
std::vector<T> successive_division_factorize(T const& N)
{
  T pos = 2;
  while(true) {
    T res = ResInt(N, pos);
    if (res == 0) {
      T quot = QuoInt(N, pos);
      if (quot > 1) {
        std::vector<T> eVect = successive_division_factorize(quot);
        eVect.push_back(pos);
        return eVect;
      }
      return {pos};
    }
    pos++;
    if (pos * pos > N)
      break;
  }
  return {N};
}





template<typename T>
bool successive_division_isprime(T const& N)
{
  T pos = 2;
  while(true) {
    T res = ResInt(N, pos);
    if (res == 0)
      return false;
    pos++;
    if (pos * pos > N)
      break;
  }
  return true;
}


template<typename T>
bool IsPrime(const T& N)
{
  std::pair<bool, T> epair = rho_pollard_factorize(N);
  if (epair.first) {
    return false;
  } else {
    return successive_division_isprime(N);
  }
}


template<typename T>
std::vector<T> FactorsInt(T const& N)
{
  std::pair<bool, T> epair = rho_pollard_factorize(N);
  if (epair.first) {
    T fact1 = epair.second;
    T fact2 = QuoInt(N, fact1);
    std::vector<T> ListPrime = FactorsInt(fact1);
    std::vector<T> V2 = FactorsInt(fact2);
    ListPrime.insert(ListPrime.end(), V2.begin(), V2.end());
    return ListPrime;
  } else {
    return successive_division_factorize(N);
  }
}


template<typename T>
std::vector<T> GetAllFactors(std::map<T,int> const& eMap)
{
  std::vector<T> LVal = {1};
  for (auto & kv : eMap) {
    std::vector<T> NewVal;
    T ePow = 1;
    T mult = kv.first;
    for (int i=0; i<=kv.second; i++) {
      for (auto & eVal : LVal)
        NewVal.push_back(ePow * eVal);
      ePow *= mult;
    }
    LVal = NewVal;
  }
  return LVal;
}


template<typename T>
std::vector<T> GetAllFactors(T const& N)
{
  std::vector<T> LFact = FactorsInt(N);
  std::map<T,int> eMap;
  for (auto & eVal : LFact)
    eMap[eVal]++;
  return GetAllFactors(eMap);
}




#endif
