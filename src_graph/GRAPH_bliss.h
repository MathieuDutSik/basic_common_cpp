// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GRAPH_GRAPH_BLISS_H_
#define SRC_GRAPH_GRAPH_BLISS_H_

#include "ExceptionsFunc.h"
#include "defs.hh"
#include "graph.hh"
#include "partition.hh"
#include "timer.hh"
#include "utils.hh"
#include <iostream>
#include <string>
#include <utility>
#include <vector>

bliss::Graph *ReadGraphFromFile(FILE *f, unsigned int &nof_vertices) {
  unsigned int nof_edges;
  int ret;
  bliss::Graph *g = 0;
  ret = fscanf(f, "%u %u\n", &nof_vertices, &nof_edges);
  if (ret != 1) {
    std::cerr << "fscanf error while reading graph 1\n";
    throw TerminalException{1};
  }
  g = new bliss::Graph(nof_vertices);
  for (int i = 0; i < static_cast<int>(nof_vertices); i++) {
    unsigned int color;
    ret = fscanf(f, "%u\n", &color);
    if (ret != 1) {
      std::cerr << "fscanf error while reading graph 2\n";
      throw TerminalException{1};
    }
    g->change_color(i, color);
  }
  for (int iEdge = 0; iEdge < static_cast<int>(nof_edges); iEdge++) {
    int a, b;
    ret = fscanf(f, "%u %u\n", &a, &b);
    if (ret != 1) {
      std::cerr << "fscanf error while reading graph 3\n";
      throw TerminalException{1};
    }
    g->add_edge(a - 1, b - 1);
  }
  return g;
}

template <typename Tgr> bliss::Graph GetBlissGraphFromGraph(Tgr const &eGR) {
  using T_bliss = unsigned int;
  T_bliss nbVert = T_bliss(eGR.GetNbVert());
  bliss::Graph g(nbVert);
  for (T_bliss iVert = 0; iVert < nbVert; iVert++) {
    T_bliss eColor = 0;
    if (eGR.GetHasVertexColor())
      eColor = T_bliss(eGR.GetColor(iVert));
    g.change_color(iVert, eColor);
  }
  for (T_bliss iVert = 0; iVert < nbVert - 1; iVert++)
    for (T_bliss jVert = iVert + 1; jVert < nbVert; jVert++)
      if (eGR.IsAdjacent(iVert, jVert))
        g.add_edge(iVert, jVert);
  return g;
}

template <typename Tgr>
bool IsIsomorphicGraph(Tgr const &eGR1, Tgr const &eGR2) {
  if (eGR1.GetNbVert() != eGR2.GetNbVert())
    return false;
  size_t nof_vertices = eGR1.GetNbVert();
  //
  bliss::Graph g1 = GetBlissGraphFromGraph(eGR1);
  bliss::Graph g2 = GetBlissGraphFromGraph(eGR2);
  bliss::Stats stats;
  //
  const unsigned int *cl1;
  cl1 = g1.canonical_form(stats);
  const unsigned int *cl2;
  cl2 = g2.canonical_form(stats);
  std::vector<size_t> clR2(nof_vertices);
  for (size_t i = 0; i < nof_vertices; i++)
    clR2[cl2[i]] = i;
  std::vector<size_t> TheEquivExp(nof_vertices);
  for (size_t iVert = 0; iVert < nof_vertices; iVert++)
    TheEquivExp[iVert] = clR2[cl1[iVert]];
  for (size_t iVert = 0; iVert < nof_vertices; iVert++) {
    size_t jVert = TheEquivExp[iVert];
    if (eGR1.GetColor(iVert) != eGR2.GetColor(jVert))
      return false;
  }
  for (size_t iVert1 = 0; iVert1 < nof_vertices; iVert1++) {
    size_t iVert2 = TheEquivExp[iVert1];
    for (size_t jVert1 = 0; jVert1 < nof_vertices; jVert1++) {
      size_t jVert2 = TheEquivExp[jVert1];
      if (eGR1.IsAdjacent(iVert1, jVert1) != eGR2.IsAdjacent(iVert2, jVert2))
        return false;
    }
  }
  return true;
}

template <typename Tgr> std::string GetCanonicalForm_string(Tgr const &eGR) {
  size_t nof_vertices = eGR.GetNbVert();
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  //
  const unsigned int *cl;
  cl = g.canonical_form(stats);
  std::vector<size_t> clR(nof_vertices);
  for (size_t i = 0; i < nof_vertices; i++)
    clR[cl[i]] = i;
  //
  std::string strRet;
  for (size_t iVert = 0; iVert < nof_vertices; iVert++) {
    size_t iVertCan = clR[iVert];
    if (eGR.GetHasVertexColor()) {
      size_t eColor = eGR.GetColor(iVertCan);
      strRet += " " + std::to_string(eColor);
    }
    //
    for (size_t jVert = 0; jVert < nof_vertices; jVert++) {
      size_t jVertCan = clR[jVert];
      bool eVal_b = eGR.IsAdjacent(iVertCan, jVertCan);
      strRet += " " + std::to_string(eVal_b);
    }
  }
  return strRet;
}

template <typename Tgr, typename TidxC>
std::vector<TidxC> BLISS_GetCanonicalOrdering(Tgr const &eGR) {
  size_t nof_vertices = eGR.GetNbVert();
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  //
  const unsigned int *cl;
  cl = g.canonical_form(stats);
  std::vector<TidxC> vectD(nof_vertices);
  for (size_t i = 0; i < nof_vertices; i++)
    vectD[i] = cl[i];
  return vectD;
}

struct RecParam {
  size_t n_last;
  std::vector<std::vector<unsigned int>> LGen;
};

template <typename TidxG>
std::vector<std::vector<TidxG>>
get_generators(bliss::Graph &g, bliss::Stats &stats, size_t n_last) {
  std::vector<std::vector<TidxG>> ListGen;
  std::vector<TidxG> eGen(n_last);
  std::function<void(unsigned int n, const unsigned int *aut)> report =
      [&]([[maybe_unused]] unsigned int n, const unsigned int *aut) -> void {
    for (size_t i = 0; i < n_last; i++) {
      eGen[i] = aut[i];
    }
    ListGen.push_back(eGen);
  };
  g.find_automorphisms(stats, report);
  return ListGen;
}

template <typename Tgr, typename TidxG>
std::vector<std::vector<TidxG>> BLISS_GetListGenerators(Tgr const &eGR,
                                                        size_t const &n_last) {
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  return get_generators<TidxG>(g, stats, n_last);
}

template <typename Tgr, typename TidxC, typename TidxG>
std::pair<std::vector<TidxC>, std::vector<std::vector<TidxG>>>
BLISS_GetCanonicalOrdering_ListGenerators(Tgr const &eGR,
                                          size_t const &n_last) {
  size_t nof_vertices = eGR.GetNbVert();
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  //
  const unsigned int *cl;
  cl = g.canonical_form(stats);
  std::vector<TidxC> vectD(nof_vertices);
  for (size_t i = 0; i < nof_vertices; i++)
    vectD[i] = cl[i];
  //
  std::vector<std::vector<TidxG>> ListGen =
      get_generators<TidxG>(g, stats, n_last);
  //
  return {std::move(vectD), std::move(ListGen)};
}

// clang-format off
#endif  // SRC_GRAPH_GRAPH_BLISS_H_
// clang-format on
