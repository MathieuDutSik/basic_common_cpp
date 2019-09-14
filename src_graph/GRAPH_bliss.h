#ifndef INCLUDE_GRAPH_BLISS
#define INCLUDE_GRAPH_BLISS

#include <string>
#include "defs.hh"
#include "graph.hh"
#include "partition.hh"
#include "timer.hh"
#include "utils.hh"




static inline void report_aut_void(void* param, const unsigned int n, const unsigned int* aut)
{

}


template<typename Tgr>
bliss::Graph GetBlissGraphFromGraph(Tgr const& eGR)
{
  int nbVert=eGR.GetNbVert();
  bliss::Graph g(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int eColor=eGR.GetColor(iVert);
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
  bliss::Graph g2 = GetBlissGraphFromGraph(eGR);
  bliss::Stats stats;
  //
  const unsigned int* cl;
  cl=g1.canonical_form(stats, &report_aut_void, stderr);
  std::vector<int> clR(nof_vertices);
  for (int i=0; i<nof_vertices; i++)
    clR[cl[i]]=i;
  //
  std::string strRet;
  for (int iVert=0; iVert<nof_vertices; iVert++) {
    int iVertCan = clR[iVert];
    if (eGR.HasColor()) {
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




#endif
