#ifndef INCLUDE_GRAPH_BLISS
#define INCLUDE_GRAPH_BLISS

#include <string>
#include <iostream>
#include "ExceptionEnding.h"
#include "defs.hh"
#include "graph.hh"
#include "partition.hh"
#include "timer.hh"
#include "utils.hh"


bliss::Graph* ReadGraphFromFile(FILE *f, unsigned int &nof_vertices)
{
  unsigned int nof_edges;
  int ret;
  bliss::Graph *g =0;
  ret=fscanf(f, "%u %u\n", &nof_vertices, &nof_edges);
  if (ret != 1) {
    std::cerr << "fscanf error while reading graph 1\n";
    throw TerminalException{1};
  }
  g = new bliss::Graph(nof_vertices);
  for (int i=0; i<int(nof_vertices); i++) {
    unsigned int color;
    ret=fscanf(f, "%u\n", &color);
    if (ret != 1) {
      std::cerr << "fscanf error while reading graph 2\n";
      throw TerminalException{1};
    }
    g->change_color(i, color);
  }
  for (int iEdge=0; iEdge<int(nof_edges); iEdge++) {
    int a, b;
    ret=fscanf(f, "%u %u\n", &a, &b);
    if (ret != 1) {
      std::cerr << "fscanf error while reading graph 3\n";
      throw TerminalException{1};
    }
    g->add_edge(a-1, b-1);
  }
  return g;
}


// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
static inline void report_aut_void(void* param, const unsigned int n, const unsigned int* aut)
{

}
#pragma GCC diagnostic pop


template<typename Tgr>
bliss::Graph GetBlissGraphFromGraph(Tgr const& eGR)
{
  int nbVert=eGR.GetNbVert();
  bliss::Graph g(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int eColor = 0;
    if (eGR.GetHasVertexColor())
      eColor = eGR.GetColor(iVert);
    g.change_color(iVert, eColor);
  }
  for (int iVert=0; iVert<nbVert-1; iVert++)
    for (int jVert=iVert+1; jVert<nbVert; jVert++)
      if (eGR.IsAdjacent(iVert,jVert))
	g.add_edge(iVert,jVert);
  return g;
}


template<typename Tgr>
bool IsIsomorphicGraph(Tgr const& eGR1, Tgr const& eGR2)
{
  if (eGR1.GetNbVert() != eGR2.GetNbVert())
    return false;
  int nof_vertices = eGR1.GetNbVert();
  //
  bliss::Graph g1 = GetBlissGraphFromGraph(eGR1);
  bliss::Graph g2 = GetBlissGraphFromGraph(eGR2);
  bliss::Stats stats;
  //
  const unsigned int* cl1;
  cl1=g1.canonical_form(stats, &report_aut_void, stderr);
  const unsigned int* cl2;
  cl2=g2.canonical_form(stats, &report_aut_void, stderr);
  std::vector<int> clR2(nof_vertices);
  for (int i=0; i<nof_vertices; i++)
    clR2[cl2[i]]=i;
  std::vector<int> TheEquivExp(nof_vertices);
  for (int iVert=0; iVert<nof_vertices; iVert++)
    TheEquivExp[iVert]=clR2[cl1[iVert]];
  for (int iVert=0; iVert<nof_vertices; iVert++) {
    int jVert=TheEquivExp[iVert];
    if (eGR1.GetColor(iVert) != eGR2.GetColor(jVert))
      return false;
  }
  for (int iVert1=0; iVert1<nof_vertices; iVert1++) {
    int iVert2=TheEquivExp[iVert1];
    for (int jVert1=0; jVert1<nof_vertices; jVert1++) {
      int jVert2=TheEquivExp[jVert1];
      if (eGR1.IsAdjacent(iVert1,jVert1) != eGR2.IsAdjacent(iVert2,jVert2) )
        return false;
    }
  }
  return true;
}

template<typename Tgr>
std::string GetCanonicalForm_string(Tgr const& eGR)
{
  int nof_vertices = eGR.GetNbVert();
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  //
  const unsigned int* cl;
  cl=g.canonical_form(stats, &report_aut_void, stderr);
  std::vector<int> clR(nof_vertices);
  for (int i=0; i<nof_vertices; i++)
    clR[cl[i]]=i;
  //
  std::string strRet;
  for (int iVert=0; iVert<nof_vertices; iVert++) {
    int iVertCan = clR[iVert];
    if (eGR.GetHasVertexColor()) {
      int eColor = eGR.GetColor(iVertCan);
      strRet += " " + std::to_string(eColor);
    }
    //
    for (int jVert=0; jVert<nof_vertices; jVert++) {
      int jVertCan = clR[jVert];
      bool eVal_b = eGR.IsAdjacent(iVertCan, jVertCan);
      strRet += " " + std::to_string(eVal_b);
    }
  }
  return strRet;
}


template<typename Tgr, typename Tidx>
std::vector<unsigned int> BLISS_GetCanonicalOrdering(Tgr const& eGR)
{
  int nof_vertices = eGR.GetNbVert();
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  //
  const unsigned int* cl;
  cl=g.canonical_form(stats, &report_aut_void, stderr);
  std::vector<Tidx> vectD(nof_vertices);
  for (int i=0; i<nof_vertices; i++)
    vectD[i] = cl[i];
  return vectD;
}

struct RecParam {
  int n_last;
  std::vector<std::vector<unsigned int>> LGen;
};


// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
static inline void report_aut_vectvectint(void* param, const unsigned int n, const unsigned int* aut)
{
  RecParam* rec_param = (RecParam*)param;
  int n_last = rec_param->n_last;
  std::vector<unsigned int> eVect(n_last);
  for (int i=0; i<n_last; i++)
    eVect[i] = aut[i];
  rec_param->LGen.push_back(eVect);
}
#pragma GCC diagnostic pop



template<typename Tgr, typename Tidx>
std::vector<std::vector<Tidx>> BLISS_GetListGenerators(Tgr const& eGR, int n_last)
{
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  RecParam rec_param;
  rec_param.n_last = n_last;
  RecParam* rec_param_ptr = &rec_param;
  g.find_automorphisms(stats, &report_aut_vectvectint, (void *)rec_param_ptr);
  std::vector<std::vector<Tidx>> ListGen;
  for (auto & eList : rec_param.LGen) {
    std::vector<Tidx> V(n_last);
    for (int i=0; i<n_last; i++)
      V[i] = eList[i];
    ListGen.push_back(V);
  }
  return ListGen;
}


template<typename Tgr, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>>  BLISS_GetCanonicalOrdering_ListGenerators(Tgr const& eGR, int n_last)
{
  int nof_vertices = eGR.GetNbVert();
  bliss::Graph g = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  //
  const unsigned int* cl;
  cl=g.canonical_form(stats, &report_aut_void, stderr);
  std::vector<Tidx> vectD(nof_vertices);
  for (int i=0; i<nof_vertices; i++)
    vectD[i] = cl[i];
  //
  RecParam rec_param;
  rec_param.n_last = n_last;
  RecParam* rec_param_ptr = &rec_param;
  g.find_automorphisms(stats, &report_aut_vectvectint, (void *)rec_param_ptr);
  std::vector<std::vector<Tidx>> ListGen;
  for (auto & eList : rec_param.LGen) {
    std::vector<Tidx> V(n_last);
    for (int i=0; i<n_last; i++)
      V[i] = eList[i];
    ListGen.push_back(V);
  }
  //
  return {std::move(vectD), std::move(ListGen)};
}




#endif
