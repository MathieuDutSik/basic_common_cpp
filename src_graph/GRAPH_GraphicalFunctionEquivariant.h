// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GRAPH_GRAPH_GRAPHICALFUNCTIONEQUIVARIANT_H_
#define SRC_GRAPH_GRAPH_GRAPHICALFUNCTIONEQUIVARIANT_H_

#include "GRAPH_GraphicalFunctions.h"
#include "GRP_GroupFct.h"
#include <vector>

template <typename Tgr>
int ComputeDiameterEquivariant(Tgr const &GR, TheGroupFormat const &TheGRP) {
  int nbVert = GR.GetNbVert();
  Face eList(nbVert);
  for (int iVert = 0; iVert < nbVert; iVert++)
    eList[iVert] = 1;
  std::vector<Face> vvO = DecomposeOrbitPoint(TheGRP, eList);
  std::function<int(int const &)> LocalDiameter = [&](int const &nPt) -> int {
    int localDiameter = 0;
    Face eStatus(nbVert);
    for (int iVert = 0; iVert < nbVert; iVert++)
      eStatus[iVert] = 0;
    int nbDone = 0;
    std::vector<int> ListPoint{nPt};
    while (true) {
      std::vector<int> NewListPoint;
      for (auto &zPt : ListPoint) {
        nbDone++;
        for (const size_t &hPt : GR.Adjacency(zPt)) {
          if (eStatus[hPt] == 0) {
            eStatus[hPt] = 1;
            NewListPoint.push_back(hPt);
          }
        }
      }
      if (nbDone == nbVert)
        return localDiameter;
      if (NewListPoint.size() == 0)
        return -1;
      ListPoint = NewListPoint;
      localDiameter++;
    }
  };
  int TheDiameter = 0;
  for (auto &eOrb : vvO) {
    int ePt = eOrb[0];
    int localDiameter = LocalDiameter(ePt);
    if (localDiameter == -1) {
      return -1;
    }
    if (localDiameter > TheDiameter)
      TheDiameter = localDiameter;
  }
  return TheDiameter;
}

#endif // SRC_GRAPH_GRAPH_GRAPHICALFUNCTIONEQUIVARIANT_H_
