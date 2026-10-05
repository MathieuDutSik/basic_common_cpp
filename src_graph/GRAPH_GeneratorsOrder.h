// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GRAPH_GRAPH_GENERATORSORDER_H_
#define SRC_GRAPH_GRAPH_GENERATORSORDER_H_

// clang-format off
#include <cstddef>
#include <map>
#include <utility>
#include <vector>
// clang-format on

// The generators of the automorphism group of a graph, restricted to its
// first n_last vertices, with the order of the group when the graph program
// gives it exactly (has_order): the product of k^m over the pairs (k, m) of
// order_factors. It is the order of the group the generators generate when
// the restriction to the first vertices is faithful, and a multiple of it
// otherwise.
template <typename TidxG> struct GraphGeneratorsOrder {
  std::vector<std::vector<TidxG>> ListGen;
  bool has_order;
  std::vector<std::pair<size_t, size_t>> order_factors;
};

// The order factors from a list of integer factors (with repetitions).
template <typename Tfact>
std::vector<std::pair<size_t, size_t>>
GetGraphOrderFactors(std::vector<Tfact> const &factors) {
  std::map<size_t, size_t> map;
  for (auto &k : factors)
    if (k > 1)
      map[static_cast<size_t>(k)]++;
  return std::vector<std::pair<size_t, size_t>>(map.begin(), map.end());
}

// The order of the group in the integer type Tint, which has to hold it.
template <typename Tint, typename TidxG>
Tint GetGraphGroupOrder(GraphGeneratorsOrder<TidxG> const &x) {
  Tint order(1);
  for (auto &[k, m] : x.order_factors) {
    Tint k_t(k);
    for (size_t i = 0; i < m; i++)
      order *= k_t;
  }
  return order;
}

// clang-format off
#endif  // SRC_GRAPH_GRAPH_GENERATORSORDER_H_
// clang-format on
