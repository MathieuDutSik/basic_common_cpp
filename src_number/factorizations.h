// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_FACTORIZATIONS_H_
#define SRC_NUMBER_FACTORIZATIONS_H_

// clang-format off
#include "TemplateTraits.h"
#include <map>
#include <utility>
#include <vector>
// clang-format on

template <typename T>
std::pair<bool, T> rho_pollard_factorize(T const &number) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
#ifdef DEBUG_FACTORIZATION
  if (number == 0) {
    std::cerr << "FACT: Trying to factorize a number equal to zero is probably not what you had in mind\n";
    throw TerminalException{1};
  }
#endif
  T count;
  T x_fixed = 2, x = 2, size = 2, factor, diff;
  do {
    count = size;
    do {
      ResInt_Kernel(x * x + 1, number, x);
      diff = x - x_fixed;
      if (diff < 0)
        diff = -diff;
      factor = GcdPair(diff, number);
    } while (--count > 0 && factor == 1);
    size *= 2;
    x_fixed = x;
  } while (factor == 1);
  if (factor == number) {
    return {false, -1};
  } else {
    return {true, factor};
  }
}

template <typename T> std::vector<T> successive_division_factorize(T const &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  T pos = 2;
  while (true) {
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

template <typename T> bool successive_division_isprime(T const &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  if (N == 2) {
    return true;
  }
  T pos = 2;
  while (true) {
    T res = ResInt(N, pos);
    if (res == 0)
      return false;
    pos++;
    if (pos * pos > N)
      break;
  }
  return true;
}

template <typename T> bool IsPrime(const T &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  std::pair<bool, T> epair = rho_pollard_factorize(N);
  if (epair.first) {
    return false;
  } else {
    return successive_division_isprime(N);
  }
}

template <typename T> std::vector<T> Kernel_FactorsInt(T const &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  if (N == 1)
    return {};
  std::pair<bool, T> epair = rho_pollard_factorize(N);
  if (epair.first) {
    T fact1 = epair.second;
    T fact2 = QuoInt(N, fact1);
    std::vector<T> ListPrime = Kernel_FactorsInt(fact1);
    std::vector<T> V2 = Kernel_FactorsInt(fact2);
    ListPrime.insert(ListPrime.end(), V2.begin(), V2.end());
    return ListPrime;
  } else {
    return successive_division_factorize(N);
  }
}

template <typename T>
inline typename std::enable_if<is_implementation_of_Z<T>::value,
                               std::vector<T>>::type
FactorsInt(T const &N) {
  return Kernel_FactorsInt(N);
}

template <typename T>
inline typename std::enable_if<!is_implementation_of_Z<T>::value,
                               std::vector<T>>::type
FactorsInt(T const &N) {
  using Tint = typename underlying_ring<T>::ring_type;
  Tint N_int = UniversalScalarConversion<Tint, T>(N);
  std::vector<Tint> LFact_int = Kernel_FactorsInt(N_int);
  std::vector<T> LFact;
  for (auto &val_i : LFact_int) {
    T val = UniversalScalarConversion<T, Tint>(val_i);
    LFact.push_back(val);
  }
  return LFact;
}

template <typename T> std::map<T, size_t> FactorsIntMap(T const &N) {
  std::vector<T> vect = FactorsInt(N);
  std::map<T, size_t> map;
  for (auto &eV : vect) {
    map[eV] += 1;
  }
  return map;
}

// 389258285913703137976972453410674957488110405 2265346989680430163430848057126071 10374024096391
// Factorizing is expensive, so if you know some number that could get you some factors, then use them.
template<typename T> std::map<T, size_t> FactorsIntMap_help(T const& N, std::vector<T> const& list_help) {
  std::vector<T> ListN{N};
#ifdef DEBUG_FACTORIZATION
  std::cerr << "FACT: N=" << N << "\n";
#endif
  // Using the helper factors to break N into some factors,
  for (auto & help : list_help) {
#ifdef DEBUG_FACTORIZATION
    std::cerr << "help=" << help << "\n";
#endif
    std::vector<T> NewListN;
    for (auto & eN : ListN) {
      T gcd = GcdPair(eN, help);
#ifdef DEBUG_FACTORIZATION
      std::cerr << "  eN=" << eN << " gcd=" << gcd << "\n";
#endif
      if (gcd != 1 && gcd != eN) {
        NewListN.push_back(gcd);
        NewListN.push_back(eN / gcd);
      } else {
        NewListN.push_back(eN);
      }
    }
    ListN = NewListN;
  }
#ifdef DEBUG_FACTORIZATION
  std::cerr << "FACT: ListN=";
  for (auto & eN: ListN) {
    std::cerr << " " << eN;
  }
  std::cerr << "\n";
#endif
  struct result {
    size_t u;
    size_t v;
    T gcd;
  };
  std::map<T, size_t> MapN;
  for (auto & val : ListN) {
    MapN[val] += 1;
  }
#ifdef DEBUG_FACTORIZATION
  std::cerr << "FACT: 1 : MapN=";
  for (auto & kv: MapN) {
    std::cerr << " (" << kv.first << "," << kv.second << ")";
  }
  std::cerr << "\n";
#endif
  auto get_list=[](std::map<T, size_t> const& m) -> std::vector<std::pair<T,size_t>> {
    std::vector<std::pair<T, size_t>> v;
    for (auto & kv: m) {
      v.push_back({kv.first, kv.second});
    }
    return v;
  };
  auto get_map=[](std::vector<std::pair<T, size_t>> const& v) -> std::map<T,size_t> {
    std::map<T, size_t> m;
    for (auto & pair: v) {
      m[pair.first] += pair.second;
    }
    return m;
  };
  // Detect some factors
  auto get_pair=[](std::vector<std::pair<T,size_t>> const& V) -> std::optional<result> {
    size_t len = V.size();
    for (size_t u=0; u<len; u++) {
      for (size_t v=u+1; v<len; v++) {
        T gcd = GcdPair(V[u].first, V[v].first);
        if (gcd != 1) {
          result res{u, v, gcd};
          return res;
        }
      }
    }
    return {};
  };
  while(true) {
    std::vector<std::pair<T,size_t>> NewListN = get_list(MapN);
    std::optional<result> opt = get_pair(NewListN);
    if (opt) {
      result res = *opt;
      std::vector<std::pair<T,size_t>> NewListN2;
      size_t u = res.u;
      size_t v = res.v;
      T gcd = res.gcd;
      for (size_t i=0; i<NewListN.size(); i++) {
        if (i != u && i != v) {
          NewListN2.push_back(NewListN[i]);
        } else {
          T val1 = NewListN[i].first / gcd;
          std::vector<T> V{gcd, val1};
          for (auto & val : V) {
            if (val != 1) {
              NewListN2.push_back({val, NewListN[i].second});
            }
          }
        }
      }
      MapN = get_map(NewListN2);
#ifdef DEBUG_FACTORIZATION
      std::cerr << "FACT: 2 : MapN=";
      for (auto & kv: MapN) {
        std::cerr << " (" << kv.first << "," << kv.second << ")";
      }
      std::cerr << "\n";
#endif
    } else {
      break;
    }
  }
  std::map<T, size_t> map_ret;
  for (auto & kv1 : MapN) {
    T const& eN = kv1.first;
    size_t const& mult = kv1.second;
#ifdef DEBUG_FACTORIZATION
    std::cerr << "FACT: eN=" << eN << " mult=" << mult << "\n";
#endif
    std::map<T, size_t> map = FactorsIntMap(eN);
    for (auto & kv2 : map) {
      map_ret[kv2.first] += kv2.second * mult;
    }
  }
  return map_ret;
}



template <typename T>
std::vector<T> GetAllFactors(std::map<T, int> const &eMap) {
  std::vector<T> LVal = {1};
  for (auto &kv : eMap) {
    std::vector<T> NewVal;
    T ePow = 1;
    T mult = kv.first;
    for (int i = 0; i <= kv.second; i++) {
      for (auto &eVal : LVal)
        NewVal.push_back(ePow * eVal);
      ePow *= mult;
    }
    LVal = NewVal;
  }
  return LVal;
}

template <typename T> std::vector<T> GetAllFactors(T const &N) {
  std::vector<T> LFact = FactorsInt(N);
  std::map<T, int> eMap;
  for (auto &eVal : LFact)
    eMap[eVal]++;
  return GetAllFactors(eMap);
}

// clang-format off
#endif  // SRC_NUMBER_FACTORIZATIONS_H_
// clang-format on
