#ifndef INCLUDE_HYPERGRAPH_BASIC
#define INCLUDE_HYPERGRAPH_BASIC


struct HypergraphBitset {
  int GetNbVert() const
  {
    return nbVert;
  }
  int GetNbHyperedge() const
  {
    return nbHyperedge;
  }
  std::vector<int> GetHyperedge(int const& iHyper) const
  {
    std::vector<int> LEnt;
    for (int iVert=0; iVert<nbVert; iVert++) {
      int val = ListEntry(iHyper * nbVert + iVert);
      if (val == 1)
        LEnt.push_back(iVert);
    }
    return LEnt;
  }
  std::vector<int> GetContainingHyperedge(int const& iVert) const
  {
    std::vector<int> LEnt;
    for (int iHyper=0; iHyper<nbHyperedge; iHyper++) {
      int val = ListEntry(iHyper * nbVert + iVert);
      if (val == 1)
        LEnt.push_back(iHyper);
    }
    return LEnt;
  }
private:
  int nbVert;
  int nbHyperedge;
  Face ListEntry;
}



#endif
