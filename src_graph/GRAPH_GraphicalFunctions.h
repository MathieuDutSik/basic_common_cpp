// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GRAPH_GRAPH_GRAPHICALFUNCTIONS_H_
#define SRC_GRAPH_GRAPH_GRAPHICALFUNCTIONS_H_

#include "COMB_Combinatorics.h"
#include "COMB_Vectors.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "MAT_Matrix.h"
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

struct GraphFunctional {
public:
  GraphFunctional() = delete;
  GraphFunctional(
      size_t const &inpNbVert,
      std::function<bool(size_t const &, size_t const &)> const &eFct)
      : nbVert(inpNbVert), f(eFct) {
    HasVertexColor = false;
  }
  ~GraphFunctional() {}
  GraphFunctional(GraphFunctional const &eG) {
    nbVert = eG.GetNbVert();
    f = eG.GetFCT();
    HasVertexColor = eG.GetHasVertexColor();
    fColor = eG.GetFColor();
  }
  GraphFunctional operator=(GraphFunctional const &eG) {
    nbVert = eG.GetNbVert();
    f = eG.GetFCT();
    HasVertexColor = eG.GetHasVertexColor();
    fColor = eG.GetFColor();
    return *this;
  }
  // lighter stuff
  size_t GetNbVert() const { return nbVert; }
  std::function<bool(size_t const &, size_t const &)> GetFCT() const {
    return f;
  }
  bool GetHasVertexColor() const { return HasVertexColor; }
  std::function<size_t(size_t const &)> GetFColor() const { return fColor; }
  //
  void SetFColor(std::function<size_t(size_t const &)> const &inpFColor) {
    HasVertexColor = true;
    fColor = inpFColor;
  }
  std::vector<size_t> Adjacency(int const &iVert) const {
    std::vector<size_t> retList;
    for (size_t jVert = 0; jVert < nbVert; jVert++)
      if (f(iVert, jVert))
        retList.push_back(jVert);
    return retList;
  }
  bool IsAdjacent(size_t const &iVert, size_t const &jVert) const {
    return f(iVert, jVert);
  }
  size_t GetColor(size_t const &iVert) const { return fColor(iVert); }

private:
  size_t nbVert;
  std::function<bool(size_t const &, size_t const &)> f;
  bool HasVertexColor;
  std::function<size_t(size_t const &)> fColor;
};

template <typename T> struct is_functional_graph_class {
  static const bool value = false;
};

template <> struct is_functional_graph_class<GraphFunctional> {
  static const bool value = true;
};

template <typename Tgr>
void GRAPH_PrintOutput(std::ostream &os, Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  os << "nbVert=" << nbVert << "\n";
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    std::vector<size_t> LVert = GR.Adjacency(iVert);
    os << "iVert=" << iVert << " |Adj|=" << LVert.size() << " Adj =";
    for (auto &eVert : LVert)
      os << " " << eVert;
    os << "\n";
  }
}

template <typename Tgr>
void GRAPH_PrintOutputGAP_vertex_colored(std::string const &eFile,
                                         Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  std::vector<std::vector<size_t>> ListEdges;
  std::ofstream os(eFile);
  os << "local ListAdjacency, ThePartition;\n";
  os << "ListAdjacency:=[";
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    if (iVert > 0)
      os << ",\n";
    std::vector<size_t> LVert = GR.Adjacency(iVert);
    WriteStdVectorGAP(os, LVert);
  }
  os << "];\n";
  std::vector<size_t> ListColor = GR.GetListVertexColor();
  size_t nbBlock = VectorMax(ListColor) + 1;
  std::vector<std::vector<size_t>> ListBlock(nbBlock);
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    size_t eColor = ListColor[iVert];
    ListBlock[eColor].push_back(iVert + 1);
  }
  os << "ThePartition:=[";
  for (size_t iBlock = 0; iBlock < nbBlock; iBlock++) {
    if (iBlock > 0)
      os << ",";
    WriteStdVectorGAP(os, ListBlock[iBlock]);
  }
  os << "];\n";
  os << "return [ListAdjacency, ThePartition];\n";
}

template <typename Tgr>
void GRAPH_PrintOutputGAP(std::ostream &os, Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  std::vector<std::vector<size_t>> ListEdges;
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    std::vector<size_t> LVert = GR.Adjacency(iVert);
    for (auto &eVert : LVert)
      if (eVert > iVert)
        ListEdges.push_back({iVert + 1, eVert + 1});
  }
  size_t nbEdge = ListEdges.size();
  os << "local ListEdges, eEdge, GRA;\n";
  os << "ListEdges:=[";
  for (size_t iEdge = 0; iEdge < nbEdge; iEdge++) {
    if (iEdge > 0)
      os << ",";
    os << "[" << ListEdges[iEdge][0] << "," << ListEdges[iEdge][1] << "]";
  }
  os << "];\n";
  os << "GRA:=NullGraph(Group(()), " << nbVert << ");\n";
  os << "for eEdge in ListEdges\n";
  os << "do\n";
  os << "  AddEdgeOrbit(GRA, eEdge);\n";
  os << "  AddEdgeOrbit(GRA, Reversed(eEdge));\n";
  os << "od;\n";
  os << "return GRA;\n";
}

template <typename Tgr> void GRAPH_Write(std::ostream &os, const Tgr &eGR) {
  size_t n = eGR.GetNbVert();
  bool HasVertexColor = eGR.GetHasVertexColor();
  os << "n=" << n << " HasColor=" << HasVertexColor << "\n";
  if (HasVertexColor) {
    os << "ListColor =";
    for (size_t i = 0; i < n; i++) {
      os << " " << eGR.GetColor(i);
    }
    os << "\n";
  }
  for (size_t i = 0; i < n; i++) {
    os << i << " :";
    for (auto &eAdj : eGR.Adjacency(i)) {
      os << " " << eAdj;
    }
    os << "\n";
  }
}

template <typename Tgr> Tgr GRAPH_Read(std::istream &is) {
  if (!is.good()) {
    std::cerr << "GRAPH_Read operation failed because stream is not valied\n";
    throw TerminalException{1};
  }
  size_t nbVert;
  is >> nbVert;
  Tgr eGR(nbVert);
  int SetHasColor_i;
  is >> SetHasColor_i;
  bool SetHasColor = SetHasColor_i;
  eGR.SetHasColor(SetHasColor);
  if (SetHasColor) {
    for (size_t iVert = 0; iVert < nbVert; iVert++) {
      size_t eColor;
      is >> eColor;
      eGR.SetColor(iVert, eColor);
    }
  }
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    size_t eDeg;
    is >> eDeg;
    for (size_t i = 0; i < eDeg; i++) {
      size_t eAdj;
      is >> eAdj;
      if (eAdj >= nbVert) {
        std::cerr << "The vertex iVert=" << iVert << " of degree " << eDeg
                  << "\n";
        std::cerr << "At index i=" << i << " the adjacent vertex is " << eAdj
                  << "\n";
        std::cerr << "But this is not adequate as the total number of vertices "
                     "is nbVert="
                  << nbVert << "\n";
        throw TerminalException{1};
      }
      eGR.AddAdjacent(iVert, eAdj);
    }
  }
  return eGR;
}

template <typename Tgr> Tgr GRAPH_ReadFile(std::string const&file) {
  std::ifstream GRAfs(file);
  return GRAPH_Read<Tgr>(GRAfs);
}

template <typename Tgr>
MyMatrix<size_t> ShortestPathDistanceMatrix(Tgr const &GR) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  size_t nbVert = GR.GetNbVert();
  MyMatrix<size_t> DistMat(nbVert, nbVert);
  for (size_t iVert = 0; iVert < nbVert; iVert++)
    for (size_t jVert = 0; jVert < nbVert; jVert++)
      DistMat(iVert, jVert) = miss_val;
  for (size_t iVert = 0; iVert < nbVert; iVert++)
    DistMat(iVert, iVert) = 0;
  while (true) {
    bool IsFinished = true;
    for (size_t iVert = 0; iVert < nbVert; iVert++)
      for (size_t jVert = 0; jVert < nbVert; jVert++)
        if (DistMat(iVert, jVert) != miss_val) {
          std::vector<size_t> LLAdj = GR.Adjacency(jVert);
          for (size_t &eAdj : LLAdj) {
            size_t CandDist = DistMat(iVert, jVert) + 1;
            if (DistMat(iVert, eAdj) > CandDist) {
              DistMat(iVert, eAdj) = CandDist;
              IsFinished = false;
            }
          }
        }
    if (IsFinished)
      break;
  }
  return DistMat;
}

template <typename Tgr>
std::vector<std::vector<size_t>>
GRAPH_FindAllShortestPath(Tgr const &GR, size_t const &x, size_t const &y) {
  MyMatrix<size_t> DistMat = ShortestPathDistanceMatrix(GR);
  size_t eDist = DistMat(x, y);
  if (eDist == std::numeric_limits<size_t>::max())
    return {};
  std::vector<std::vector<size_t>> ListPath{{x}};
  for (size_t i = 1; i <= eDist; i++) {
    std::vector<std::vector<size_t>> NewListPath;
    for (auto &ePath : ListPath) {
      size_t CurrentVert = ePath[ePath.size() - 1];
      std::vector<size_t> LAdj = GR.Adjacency(CurrentVert);
      for (auto &u : LAdj) {
        if (DistMat(u, y) == eDist - i) {
          std::vector<size_t> fPath = ePath;
          fPath.push_back(u);
          NewListPath.push_back(fPath);
        }
      }
    }
    ListPath = NewListPath;
  }
  return ListPath;
}

template <typename Tgr> size_t Diameter(Tgr const &GR) {
  MyMatrix<size_t> DistMat = ShortestPathDistanceMatrix(GR);
  if (DistMat(0, 0) == std::numeric_limits<size_t>::max())
    return 1;
  size_t TheDiam = 0;
  size_t nbVert = GR.GetNbVert();
  for (size_t iVert = 0; iVert < nbVert; iVert++)
    for (size_t jVert = 0; jVert < nbVert; jVert++) {
      size_t eDist = DistMat(iVert, jVert);
      if (eDist > TheDiam)
        TheDiam = eDist;
    }
  return TheDiam;
}

template <typename Tgr> bool IsSimpleGraph(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  for (size_t iVert = 0; iVert < nbVert; iVert++)
    if (GR.IsAdjacent(iVert, iVert))
      return false;
  return true;
}

template <typename Tgr> bool IsSymmetricGraph(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  MyMatrix<int8_t> M = ZeroMatrix<int8_t>(nbVert, nbVert);
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    std::vector<size_t> LAdj = GR.Adjacency(iVert);
    for (auto &eAdj : LAdj) {
      M(iVert, eAdj) = 1;
    }
  }
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    for (size_t jVert = iVert + 1; jVert < nbVert; jVert++) {
      if (M(iVert, jVert) != M(jVert, iVert)) {
        return false;
      }
    }
  }
  return true;
}

template <typename Tgr>
std::vector<std::vector<size_t>> GetEdgeSet(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  /*
  int nbEdge=0;
  for (int iVert=0; iVert<nbVert-1; iVert++)
    for (int jVert=iVert+1; jVert<nbVert; jVert++)
      if (GR.IsAdjacent(iVert,jVert) == 1)
      nbEdge++;
  MyMatrix<int> Edges(nbEdge,2);*/
  std::vector<std::vector<size_t>> Edges;
  for (size_t iVert = 0; iVert < nbVert - 1; iVert++)
    for (size_t jVert = iVert + 1; jVert < nbVert; jVert++)
      if (GR.IsAdjacent(iVert, jVert))
        Edges.push_back({iVert, jVert});
  return Edges;
}

template <typename Tgr>
std::vector<std::vector<size_t>> GRAPH_FindAllCycles(Tgr const &GR) {
  std::vector<std::vector<size_t>> Edges = GetEdgeSet(GR);
  size_t nbEdge = Edges.size();
  size_t nbVert = GR.GetNbVert();
  std::unordered_set<std::vector<size_t>> ListCycle;
  auto FuncInsertCycle = [&](std::vector<size_t> const &eL) -> void {
    std::vector<size_t> eLcan = MinimumDihedralOrbit(eL);
    ListCycle.insert(eLcan);
  };
  for (size_t iEdge = 0; iEdge < nbEdge; iEdge++) {
    GraphBitset GRred(nbVert);
    for (size_t jEdge = 0; jEdge < nbEdge; jEdge++)
      if (iEdge != jEdge) {
        size_t eVert = Edges[jEdge][0];
        size_t fVert = Edges[jEdge][1];
        GRred.AddAdjacent(eVert, fVert);
        GRred.AddAdjacent(fVert, eVert);
      }
    size_t x = Edges[iEdge][0];
    size_t y = Edges[iEdge][1];
    std::vector<std::vector<size_t>> ListPath =
        GRAPH_FindAllShortestPath(GRred, x, y);
    for (auto &eCycle : ListPath)
      FuncInsertCycle(eCycle);
  }
  std::vector<std::vector<size_t>> ListCycleR;
  for (auto &eCycle : ListCycle)
    ListCycleR.push_back(eCycle);
  return ListCycleR;
}

template <typename Tgr>
bool IsClique(Tgr const &GR, std::vector<size_t> const &eList) {
  size_t len = eList.size();
  for (size_t i = 0; i < len - 1; i++)
    for (size_t j = i + 1; j < len; j++)
      if (!GR.IsAdjacent(eList[i], eList[j]))
        return false;
  return true;
}

template <typename Tgr> std::vector<size_t> StartingCell(Tgr const &GR) {
  std::vector<size_t> A = GR.Adjacency(0);
  size_t eC = A[0];
  std::vector<size_t> Adj = IntersectionVect(GR.Adjacency(eC), A);
  size_t r = Adj.size();
  std::vector<size_t> TT;
  if (r == 0) {
    return {0, eC};
  } else {
    if (r == 1) {
      size_t Elt = Adj[0];
      size_t h = IntersectionVect(GR.Adjacency(Elt), GR.Adjacency(0)).size();
      size_t k = IntersectionVect(GR.Adjacency(Elt), GR.Adjacency(eC)).size();
      if (h == 1 && k == 1) {
        return {0, eC, Elt};
      } else {
        if (h > 1) {
          TT = {Elt, 0};
        } else {
          TT = {Elt, eC};
        }
      }
    } else {
      TT = {0, eC};
    }
  }
  Adj = IntersectionVect(GR.Adjacency(TT[0]), GR.Adjacency(TT[1]));
  r = Adj.size();
  std::vector<size_t> Lodd;
  size_t nbVert = GR.GetNbVert();
  for (size_t &i : Adj) {
    std::vector<size_t> VertSet(nbVert);
    for (size_t iVert = 0; iVert < nbVert; iVert++)
      VertSet[iVert] = iVert;
    std::vector<size_t> hSet{i, TT[0], TT[1]};
    std::vector<size_t> SE = DifferenceVect(VertSet, hSet);
    bool test = true;
    for (size_t &j : SE) {
      size_t oddness = IntersectionVect<size_t>(GR.Adjacency(j), hSet).size();
      if (oddness == 1 || oddness == 3)
        test = false;
    }
    if (!test)
      Lodd.push_back(i);
  }
  size_t s = Lodd.size();
  if (r == 2 && s == 0) {
    return {Adj[0], TT[0], TT[1]};
  } else {
    if (s == r || s == r - 1) {
      Lodd.push_back(TT[0]);
      Lodd.push_back(TT[1]);
      if (!IsClique(GR, Lodd))
        return {};
      return Lodd;
    } else {
      return {};
    }
  }
}

template <typename Tgr>
std::vector<std::vector<size_t>> SpanningTree(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  std::vector<int> ListStatus(nbVert, 0);
  ListStatus[0] = 1;
  std::vector<size_t> ListActiveVert{0};
  std::vector<std::vector<size_t>> TheSpann;
  while (true) {
    std::vector<size_t> NewListActiveVert;
    for (auto &eVert : ListActiveVert) {
      std::vector<size_t> LLadj = GR.Adjacency(eVert);
      for (auto &fVert : LLadj) {
        if (ListStatus[fVert] == 0) {
          ListStatus[fVert] = 1;
          TheSpann.push_back(VectorAsSet<size_t>({eVert, fVert}));
          NewListActiveVert.push_back(fVert);
        }
      }
    }
    if (NewListActiveVert.size() == 0)
      break;
    ListActiveVert = NewListActiveVert;
  }
  return TheSpann;
}

template <typename Tgr>
std::vector<std::vector<int>> InverseLineGraphConnected(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  if (nbVert == 1) {
    return {{0, 1}};
  }
  std::vector<size_t> eCell = StartingCell(GR);
  if (eCell.size() == 0) {
    std::cerr << "Failed to find starting cell\n";
    throw TerminalException{1};
  }
  std::cerr << "eCell=";
  WriteVectorInt_GAP(std::cerr, eCell);
  std::cerr << "\n";
  std::vector<size_t> eCellS = VectorAsSet(eCell);
  if (eCell[0] == std::numeric_limits<size_t>::max()) {
    return {{-1}};
  }
  std::vector<std::vector<size_t>> P = {eCell};
  std::vector<std::vector<size_t>> TotalEdge(nbVert);
  for (size_t iVert = 0; iVert < nbVert; iVert++)
    TotalEdge[iVert] = GR.Adjacency(iVert);
  for (auto &eVert : eCell)
    TotalEdge[eVert] = DifferenceVect(TotalEdge[eVert], eCell);
  size_t PreTot = 0;
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    size_t eSize = TotalEdge[iVert].size();
    PreTot += eSize;
  }
  size_t Tot = PreTot / 2;
  std::cerr << "Tot=" << Tot << "\n";
  auto GetAdding = [&]() -> size_t {
    for (size_t &iVert : eCellS)
      if (TotalEdge[iVert].size() > 0)
        return iVert;
    std::cerr << "Error in GetAdding\n";
    throw TerminalException{1};
  };
  while (true) {
    std::cerr << "Tot=" << Tot << "\n";
    if (Tot == 0)
      break;
    size_t Adding = GetAdding();
    std::cerr << "Adding=" << Adding << "\n";
    eCell = UnionVect(TotalEdge[Adding], {Adding});
    for (auto &eVert : eCell)
      for (auto &fVert : eCell)
        if (eVert != fVert)
          if (PositionVect(TotalEdge[eVert], fVert) == -1)
            return {{-1}};
    P.push_back(eCell);
    eCellS = VectorAsSet(UnionVect(eCellS, eCell));
    for (size_t &iV : eCell)
      TotalEdge[iV] = DifferenceVect(TotalEdge[iV], eCell);
    size_t Csiz = eCell.size();
    size_t nbEdge = Csiz * (Csiz - 1) / 2;
    std::cerr << "Before Csiz=" << Csiz << " nbEdge=" << nbEdge
              << " Tot=" << Tot << "\n";
    Tot -= nbEdge;
    std::cerr << "After Tot=" << Tot << "\n";
  }
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    int nb = 0;
    for (auto &eP : P) {
      int pos = PositionVect(eP, iVert);
      if (pos != -1)
        nb++;
    }
    if (nb == 1)
      P.push_back({iVert});
  }
  std::vector<std::vector<int>> Label(nbVert);
  size_t Plen = P.size();
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    for (size_t iP = 0; iP < Plen; iP++) {
      int pos = PositionVect(P[iP], iVert);
      if (pos != -1)
        Label[iVert].push_back(static_cast<int>(iP));
    }
  }
  return Label;
}

template <typename T>
void Print_VectorVector(std::ostream &os,
                        std::vector<std::vector<T>> const &TheList) {
  size_t siz = TheList.size();
  os << "|TheList|=" << siz << "\n";
  for (size_t i = 0; i < siz; i++) {
    os << "i=" << i << " V=";
    for (auto &eVal : TheList[i])
      os << " " << eVal;
    os << "\n";
  }
}

template <typename Tgr>
std::vector<std::vector<int>> InverseLineGraph(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  std::vector<std::vector<size_t>> ListConn = ConnectedComponents_set(GR);
  size_t nbConn = ListConn.size();
  auto Hsize = [&](std::vector<std::vector<int>> const &ListListSet) -> int {
    int eMax = 0;
    for (auto &eListSet : ListListSet)
      for (auto &eVal : eListSet)
        if (eVal > eMax)
          eMax = eVal;
    return eMax + 1;
  };
  std::vector<std::vector<int>> ListLabel(nbVert);
  int TotShift = 0;
  for (size_t iConn = 0; iConn < nbConn; iConn++) {
    std::vector<size_t> eConn = ListConn[iConn];
    size_t sizConn = eConn.size();
    GraphBitset GRind = InducedSubgraph<GraphBitset, Tgr>(GR, eConn);
    std::vector<std::vector<int>> PartLabel = InverseLineGraphConnected(GRind);
    int eFirstFirst = PartLabel[0][0];
    if (eFirstFirst == -1)
      return {{-1}};
    for (size_t iVertConn = 0; iVertConn < sizConn; iVertConn++) {
      size_t eVert = eConn[iVertConn];
      std::vector<int> eLabel;
      for (int &eVal : PartLabel[iVertConn])
        eLabel.push_back(eVal + TotShift);
      ListLabel[eVert] = eLabel;
    }
    int eShift = Hsize(PartLabel);
    TotShift += eShift;
  }
  return ListLabel;
}

MyMatrix<int> CreateEmbedding(int const &StartPoint,
                              std::vector<std::vector<size_t>> const &ListEdge,
                              std::vector<std::vector<int>> const &ListLabel) {
  int nbEdge = static_cast<int>(ListEdge.size());
  //
  size_t nbVert = 0;
  for (auto &eEdge : ListEdge)
    for (auto &eVal : eEdge)
      if (eVal > nbVert)
        nbVert = eVal;
  nbVert++;
  //
  int nbLabel = 0;
  for (auto &eLabel : ListLabel)
    for (auto &eVal : eLabel)
      if (eVal > nbLabel)
        nbLabel = eVal;
  nbLabel++;
  //
  MyMatrix<int> Embedding(nbVert, nbLabel);
  for (int i = 0; i < nbLabel; i++)
    Embedding(StartPoint, i) = 0;
  std::vector<int> ListStatus(nbVert, 0);
  ListStatus[StartPoint] = 1;
  while (true) {
    bool IsFinished = true;
    for (int iEdge = 0; iEdge < nbEdge; iEdge++) {
      std::vector<size_t> eEdge = ListEdge[iEdge];
      int eStat1 = ListStatus[eEdge[0]];
      int eStat2 = ListStatus[eEdge[1]];
      int RelPos = -1;
      int eVert1 = -1, eVert2 = -1;
      if (eStat1 == 0 && eStat2 == 1) {
        RelPos = 0;
        eVert1 = static_cast<int>(eEdge[1]);
        eVert2 = static_cast<int>(eEdge[0]);
      }
      if (eStat2 == 0 && eStat1 == 1) {
        RelPos = 0;
        eVert1 = static_cast<int>(eEdge[0]);
        eVert2 = static_cast<int>(eEdge[1]);
      }
      if (RelPos != -1) {
        IsFinished = false;
        for (int i = 0; i < nbLabel; i++)
          Embedding(eVert2, i) = Embedding(eVert1, i);
        for (auto &eVal : ListLabel[iEdge]) {
          Embedding(eVert2, eVal) = 1 - Embedding(eVert1, eVal);
        }
        ListStatus[eVert2] = 1;
      }
    }
    if (IsFinished)
      break;
  }
  return Embedding;
}

template <typename Tgr>
std::vector<std::vector<size_t>> EnumerationClique(Tgr const &GR,
                                                   size_t const &kSiz) {
  size_t nbVert = GR.GetNbVert();
  std::vector<std::vector<size_t>> eList;
  if (kSiz == 1) {
    for (size_t iVert = 0; iVert < nbVert; iVert++)
      eList.push_back({iVert});
  } else {
    for (auto &eLowList : EnumerationClique(GR, kSiz - 1)) {
      size_t eFirst = eLowList[kSiz - 2] + 1;
      for (size_t i = eFirst; i < nbVert; i++) {
        bool test = true;
        for (auto &eVal : eLowList)
          if (!GR.IsAdjacent(eVal, i))
            test = false;
        if (test) {
          std::vector<size_t> NewList = UnionVect(eLowList, {i});
          eList.push_back(NewList);
        }
      }
    }
  }
  return eList;
}

int L1_distance(std::vector<int> const &V1, std::vector<int> const &V2) {
  size_t siz = V1.size();
  assert(siz == V2.size());
  int dist = 0;
  for (size_t i = 0; i < siz; i++)
    if (V1[i] != V2[i])
      dist++;
  return dist;
}

template <typename Tgr>
std::vector<MyMatrix<int>> GRAPH_S_Embedding(Tgr const &GR, size_t const &s_sz,
                                             size_t const &MaxIter,
                                             size_t &iter) {
  int s = static_cast<int>(s_sz);
  iter = 0;
  if (!IsSimpleGraph(GR)) {
    std::cerr << "The graph should be simple\n";
    return {};
  }
  MyMatrix<int> M =
      UniversalMatrixConversion<int, size_t>(ShortestPathDistanceMatrix(GR));
  std::cerr << "M=\n";
  WriteMatrix(std::cerr, M);
  if (M(0, 0) == -1) {
    std::cerr << "The graph should be connected\n";
    return {};
  }
  std::vector<int> SetZero{0};
  std::vector<int> SetOne{1};
  std::vector<int> SetTwo{2};
  std::cerr << "****Start treating " << s << "-embedding\n";
  size_t n = GR.GetNbVert();
  std::vector<std::vector<size_t>> ListEdge = GetEdgeSet(GR);
  size_t nbEdge = ListEdge.size();
  std::cerr << "nbEdge=" << nbEdge << "\n";
  size_t DiamGraph = Diameter(GR);
  std::cerr << "DiamGraph=" << DiamGraph << "\n";
  auto EdgeDist = [&](size_t const &iEdge, size_t const &jEdge) -> int {
    size_t a, b, c, d, swp;
    a = ListEdge[iEdge][0];
    b = ListEdge[iEdge][1];
    c = ListEdge[jEdge][0];
    d = ListEdge[jEdge][1];
    if (M(a, c) > s || M(a, d) > s || M(b, c) > s || M(b, d) > s)
      return -400;
    std::vector<int> ListDist;
    for (size_t u = 0; u < n; u++)
      if (T_abs(M(a, u) - M(b, u)) == 1 && T_abs(M(c, u) - M(d, u)) == 1 &&
          M(a, u) <= s && M(b, u) <= s && M(c, u) <= s && M(d, u) <= s) {
        if (M(u, b) == M(u, a) - 1) {
          swp = a;
          a = b;
          b = swp;
        }
        if (M(u, d) == M(u, c) - 1) {
          swp = c;
          c = d;
          d = swp;
        }
        int DeltaE = 1 + M(u, a) - M(u, b);
        int DeltaF = 1 + M(u, c) - M(u, d);
        if (DeltaE != 0 || DeltaF != 0) {
          std::cerr << "Not correct u distances\n";
          throw TerminalException{1};
        }
        int eDist = -M(b, d) + M(a, d) + M(b, c) - M(a, c);
        ListDist.push_back(eDist);
      }
    std::vector<int> ListDistRed = VectorAsSet(ListDist);
    if (ListDistRed.size() == 0) {
      return -400;
    }
    if (ListDistRed.size() > 1) {
      std::cerr << "We find |ListDist|=" << ListDistRed.size() << "ListDist=";
      WriteVectorInt_GAP(std::cerr, ListDistRed);
      std::cerr << "\n";
    }
    for (auto &eVal : ListDistRed)
      if (eVal >= 0)
        return eVal;
    return -400;
  };
  std::vector<std::vector<size_t>> ConSet(nbEdge);
  for (size_t iEdge = 0; iEdge < nbEdge; iEdge++)
    ConSet[iEdge] = {iEdge};
  for (size_t iEdge = 0; iEdge < nbEdge; iEdge++) {
    std::cerr << "iEdge=" << iEdge << " e=" << ListEdge[iEdge][0] << ","
              << ListEdge[iEdge][1] << "\n";
  }
  MyMatrix<std::vector<int>> MCE(nbEdge, nbEdge);
  int NbPrev = 0;
  int NbUnsolved = 0;
  std::vector<std::vector<size_t>> ListPairNegative;
  for (size_t iEdge = 0; iEdge < nbEdge - 1; iEdge++)
    for (size_t jEdge = iEdge + 1; jEdge < nbEdge; jEdge++) {
      int val = EdgeDist(iEdge, jEdge);
      std::vector<int> LVal;
      std::cerr << "iEdge=" << iEdge << " jEdge=" << jEdge << " val=" << val
                << "\n";
      if (val == -400) {
        LVal = {0, 1, 2};
        NbUnsolved++;
      } else {
        LVal = {val};
        if (val > 2 || val < -2) {
          std::cerr << "Logical error for the distances\n";
          throw TerminalException{1};
        }
        if (val < 0) {
          ListPairNegative.push_back({iEdge, jEdge});
        }
      }
      MCE(iEdge, jEdge) = LVal;
      MCE(jEdge, iEdge) = LVal;
      NbPrev += 2;
    }
  std::cerr << "NbUnsolved=" << NbUnsolved << "\n";
  std::cerr << "NbPrev    =" << NbPrev << "\n";
  std::cerr << "DiamGraph=" << DiamGraph << "\n";
  size_t nbNegative = ListPairNegative.size();
  if (nbNegative > 0) {
    std::cerr
        << "Found some pairs of edges with intersection of negative size\n";
    std::cerr << "|ListPairNegative|=" << nbNegative << "\n";
    for (size_t iNeg = 0; iNeg < nbNegative; iNeg++) {
      size_t iEdge = ListPairNegative[iNeg][0];
      size_t jEdge = ListPairNegative[iNeg][1];
      int eDist = EdgeDist(iEdge, jEdge);
      size_t a = ListEdge[iEdge][0];
      size_t b = ListEdge[iEdge][1];
      size_t c = ListEdge[jEdge][0];
      size_t d = ListEdge[jEdge][1];
      std::cerr << "i=" << iNeg << " e={" << a << "," << b << "} f={" << c
                << "," << d << "} eDist=" << eDist << "\n";
    }
    throw TerminalException{1};
  }
  if (NbPrev == NbUnsolved) {
    std::cerr << "Nothing computable. Error actually\n";
    return {};
  }
  auto GetGraphIdentityEdge =
      [&ListEdge](std::vector<std::vector<size_t>> const &LEdge,
                  MyMatrix<std::vector<int>> const &MCEwork) -> GraphBitset {
    size_t nbEdgeLoc = LEdge.size();
    GraphBitset Aspec(nbEdgeLoc);
    for (size_t iEdge = 0; iEdge < nbEdgeLoc - 1; iEdge++)
      for (size_t jEdge = iEdge + 1; jEdge < nbEdgeLoc; jEdge++) {
        int iPos = PositionVect(ListEdge, LEdge[iEdge]);
        int jPos = PositionVect(ListEdge, LEdge[jEdge]);
        if (MCEwork(iPos, jPos) == std::vector<int>({2})) {
          Aspec.AddAdjacent(iEdge, jEdge);
          Aspec.AddAdjacent(jEdge, iEdge);
        }
      }
    return Aspec;
  };
  auto IsCoherentEdgeValues =
      [&ListEdge](std::vector<std::vector<size_t>> const &Cspec,
                  std::vector<std::vector<size_t>> const &LEdge,
                  MyMatrix<std::vector<int>> const &MCEwork) -> bool {
    size_t nbComp = Cspec.size();
    for (size_t iComp = 0; iComp < nbComp - 1; iComp++)
      for (size_t jComp = iComp + 1; jComp < nbComp; jComp++) {
        std::vector<size_t> eComp = Cspec[iComp];
        std::vector<size_t> fComp = Cspec[jComp];
        size_t siz1 = eComp.size();
        size_t siz2 = fComp.size();
        bool IsFirst = true;
        std::vector<int> CommonVal;
        for (size_t iElt = 0; iElt < siz1; iElt++)
          for (size_t jElt = 0; jElt < siz2; jElt++) {
            std::vector<size_t> eEdge = LEdge[eComp[iElt]];
            std::vector<size_t> fEdge = LEdge[fComp[jElt]];
            int iPos = PositionVect(ListEdge, eEdge);
            int jPos = PositionVect(ListEdge, fEdge);
            std::vector<int> eVal = MCEwork(iPos, jPos);
            if (IsFirst) {
              CommonVal = eVal;
              IsFirst = false;
            } else {
              if (CommonVal != eVal)
                return false;
            }
          }
      }
    return true;
  };
  auto GraphIntersectionOne =
      [&ListEdge](std::vector<std::vector<size_t>> const &Cspec,
                  std::vector<std::vector<size_t>> const &LEdge,
                  MyMatrix<std::vector<int>> const &MCEwork) -> GraphBitset {
    size_t nbComp = Cspec.size();
    GraphBitset eGR(nbComp);
    for (size_t iComp = 0; iComp < nbComp - 1; iComp++)
      for (size_t jComp = iComp + 1; jComp < nbComp; jComp++) {
        std::vector<size_t> eComp = Cspec[iComp];
        std::vector<size_t> fComp = Cspec[jComp];
        std::vector<size_t> eEdge = LEdge[eComp[0]];
        std::vector<size_t> fEdge = LEdge[fComp[0]];
        int iPos = PositionVect(ListEdge, eEdge);
        int jPos = PositionVect(ListEdge, fEdge);
        std::vector<int> eVal = MCEwork(iPos, jPos);
        if (eVal == std::vector<int>({1})) {
          eGR.AddAdjacent(iComp, jComp);
          eGR.AddAdjacent(jComp, iComp);
        }
      }
    return eGR;
  };
  auto ExtendListLabel = [](std::vector<std::vector<size_t>> const &Cspec,
                            std::vector<std::vector<int>> const &PartListLabel)
      -> std::vector<std::vector<int>> {
    size_t nbVert = 0;
    for (auto &eList : Cspec)
      for (auto &eVal : eList)
        if (eVal > nbVert)
          nbVert = eVal;
    nbVert++;
    //
    std::vector<std::vector<int>> ListLabel(nbVert);
    size_t nbComp = Cspec.size();
    for (size_t iComp = 0; iComp < nbComp; iComp++) {
      std::vector<int> eLabel = PartListLabel[iComp];
      for (auto &iVert : Cspec[iComp])
        ListLabel[iVert] = eLabel;
    }
    return ListLabel;
  };
  auto CheckEmbedding = [&](MyMatrix<int> const &eEmbed) -> bool {
    size_t TheDim = eEmbed.cols();
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++) {
        std::vector<int> V1(TheDim);
        std::vector<int> V2(TheDim);
        for (size_t iC = 0; iC < TheDim; iC++) {
          V1[iC] = eEmbed(i, iC);
          V2[iC] = eEmbed(j, iC);
        }
        int dist1 = M(i, j);
        int dist2 = L1_distance(V1, V2);
        if (dist1 <= s) {
          if (2 * dist1 != dist2) {
            std::cerr << "Fail embedding at i=" << i << " j=" << j << "\n";
            return false;
          }
        }
      }
    return true;
  };
  if (DiamGraph == s_sz) {
    std::cerr << "Direct embed, step 0\n";
    std::vector<std::vector<size_t>> ET = SpanningTree(GR);
    std::cerr << "Direct embed, step 1\n";
    GraphBitset Aspec = GetGraphIdentityEdge(ET, MCE);
    std::cerr << "Direct embed, step 2\n";
    std::vector<std::vector<size_t>> Cspec = ConnectedComponents_set(Aspec);
    std::cerr << "Direct embed, step 3\n";
    for (std::vector<size_t> &eList : Cspec)
      if (!IsClique(Aspec, eList)) {
        std::cerr << "Not a clique. So no embedding\n";
        return {};
      }
    std::cerr << "Direct embed, step 4\n";
    if (!IsCoherentEdgeValues(Cspec, ET, MCE)) {
      std::cerr << "Non-coherency of values\n";
      return {};
    }
    std::cerr << "Direct embed, step 5\n";
    GraphBitset Dspec = GraphIntersectionOne(Cspec, ET, MCE);
    std::cerr << "Direct embed, step 6\n";
    std::vector<std::vector<int>> PartListLabel = InverseLineGraph(Dspec);
    std::cerr << "Direct embed, step 7\n";
    if (PartListLabel[0][0] == -1) {
      std::cerr << "Not a line graph\n";
      return {};
    }
    std::cerr << "Direct embed, step 8\n";
    std::vector<std::vector<int>> ListLabel =
        ExtendListLabel(Cspec, PartListLabel);
    std::cerr << "Direct embed, step 9\n";
    MyMatrix<int> Embedding = CreateEmbedding(0, ET, ListLabel);
    bool res = CheckEmbedding(Embedding);
    if (!res) {
      std::cerr << "Maybe some bug in the code\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<int>> ListEmbedding{Embedding};
    std::cerr << "|ListEmbedding|=" << ListEmbedding.size() << "\n";
    return ListEmbedding;
  }
  auto RefinementBlock1 =
      [&ListEdge, &GetGraphIdentityEdge](MyMatrix<std::vector<int>> &MCEwork,
                                         int &nbOper) -> bool {
    GraphBitset Aspec = GetGraphIdentityEdge(ListEdge, MCEwork);
    std::vector<std::vector<size_t>> Cspec = ConnectedComponents_set(Aspec);
    for (auto &eConn : Cspec) {
      size_t lenConn = eConn.size();
      for (size_t i = 0; i < lenConn - 1; i++)
        for (size_t j = i + 1; j < lenConn; j++) {
          size_t iEdge = eConn[i];
          size_t jEdge = eConn[j];
          std::vector<int> eVal = MCEwork(iEdge, jEdge);
          int pos = PositionVect(eVal, 2);
          if (pos == -1) {
            return false;
          }
          if (eVal != std::vector<int>({2})) {
            MCEwork(iEdge, jEdge) = {2};
            MCEwork(jEdge, iEdge) = {2};
            nbOper++;
          }
        }
    }
    return true;
  };
  auto RefinementBlock2 =
      [&ListEdge, &GetGraphIdentityEdge](MyMatrix<std::vector<int>> &MCEwork,
                                         int &nbOper) -> bool {
    GraphBitset Aspec = GetGraphIdentityEdge(ListEdge, MCEwork);
    std::vector<std::vector<size_t>> Cspec = ConnectedComponents_set(Aspec);
    size_t nbConn = Cspec.size();
    for (size_t iConn = 0; iConn < nbConn - 1; iConn++)
      for (size_t jConn = iConn + 1; jConn < nbConn; jConn++) {
        std::vector<int> ListPoss{0, 1, 2};
        std::vector<size_t> eConn = Cspec[iConn];
        std::vector<size_t> fConn = Cspec[jConn];
        size_t iSize = eConn.size();
        size_t jSize = fConn.size();
        for (size_t i = 0; i < iSize; i++)
          for (size_t j = 0; j < jSize; j++) {
            size_t iEdge = eConn[i];
            size_t jEdge = fConn[j];
            std::vector<int> eVal = MCEwork(iEdge, jEdge);
            ListPoss = IntersectionVect(ListPoss, eVal);
          }
        if (ListPoss.size() == 0) {
          return false;
        }
        for (size_t i = 0; i < iSize; i++)
          for (size_t j = 0; j < jSize; j++) {
            size_t iEdge = eConn[i];
            size_t jEdge = fConn[j];
            std::vector<int> eVal = MCEwork(iEdge, jEdge);
            if (eVal != ListPoss) {
              MCEwork(iEdge, jEdge) = ListPoss;
              MCEwork(jEdge, iEdge) = ListPoss;
              nbOper++;
            }
          }
      }
    return true;
  };
  // Iterative refinement triangles
  // If M(a,b) = M(a,c) = M(b,c) = 1
  // TYPE 1:
  //   We have a=UV, b=UW, c=VW
  //   If we have a fourth edge with M(a,d)=1 then
  //   d=Ux with x not equal to V.
  //   ---If x is not V or W then M(b,d)=1 and M(b,c)=0 and there can be more
  //   than one such
  //   ---If x=W then M(b,d)=2 and M(c,d)=1
  //   So possible configurations: [1,1,0] and [1,1,2] or [0,0,0]
  // TYPE 2:
  //   We have a=UR, b=US, c=UT
  //   If we have a fourth edge with M(a,d)=1 then
  //   If d=Ux with x not R,S,T then pattern is [1,1,1]
  //   If d=Rx with x not U,S,T then pattern is [1,0,0]
  //   If d=RS or combination   then pattern is [1,1,0]
  //   Other possible patterns: [2,1,1] and [0,0,0]
  // In summary:
  //    TYPE 1: [1,1,0],                   [2,1,1], [0,0,0]
  //    TYPE 2: [1,1,1], [1,0,0], [1,1,0], [2,1,1], [0,0,0]
  //                    Pattern [1,1,0] can occur at most 3 times.
  auto RefinementTriangles = [&ListEdge, &GetGraphIdentityEdge,
                              &GraphIntersectionOne, &SetOne, &SetZero,
                              &SetTwo](MyMatrix<std::vector<int>> &MCEwork,
                                       int &nbOper) -> bool {
    GraphBitset Aspec = GetGraphIdentityEdge(ListEdge, MCEwork);
    std::vector<std::vector<size_t>> Cspec = ConnectedComponents_set(Aspec);
    GraphBitset Dspec = GraphIntersectionOne(Cspec, ListEdge, MCEwork);
    size_t nbVertD = Dspec.GetNbVert();
    std::vector<size_t> VertSet(nbVertD);
    for (size_t iVert = 0; iVert < nbVertD; iVert++)
      VertSet[iVert] = iVert;
    std::vector<std::vector<size_t>> ListTriple = EnumerationClique(Dspec, 3);
    for (auto &eTriple : ListTriple) {
      bool Type1feasible = true;
      bool Type2feasible = true;
      int nb110 = 0;
      for (auto &d : DifferenceVect(VertSet, eTriple)) {
        size_t dEdge = Cspec[d][0];
        int nb1 = 0;
        int nb0 = 0;
        for (size_t i = 0; i < 3; i++) {
          size_t jEdge = Cspec[eTriple[i]][0];
          if (MCEwork(dEdge, jEdge) == SetOne)
            nb1++;
          if (MCEwork(dEdge, jEdge) == SetZero)
            nb0++;
        }
        if (nb1 == 2 && nb0 == 1)
          nb110++;
      }
      if (nb110 > 3)
        Type2feasible = false;
      for (int i = 0; i < 3; i++) {
        size_t iNext = NextIdx(3, i);
        size_t iPrev = PrevIdx(3, i);
        size_t iEdge = Cspec[eTriple[i]][0];
        size_t iEdgeNext = Cspec[eTriple[iNext]][0];
        size_t iEdgePrev = Cspec[eTriple[iPrev]][0];
        for (auto &d : DifferenceVect(VertSet, eTriple)) {
          size_t dEdge = Cspec[d][0];
          //
          int pos1 = PositionVect(MCEwork(dEdge, iEdge), 1);
          int pos2 = PositionVect(MCEwork(dEdge, iEdge), 2);
          if (MCEwork(dEdge, iEdgeNext) == SetOne &&
              MCEwork(dEdge, iEdgePrev) == SetOne) {
            if (pos1 != -1 && !Type2feasible) {
              std::vector<int> NewVal =
                  DifferenceVect(MCEwork(dEdge, iEdge), SetOne);
              MCEwork(dEdge, iEdge) = NewVal;
              MCEwork(iEdge, dEdge) = NewVal;
              nbOper++;
            }
            // value is necessarily 1 or 0 and this forbids Type1
            if (pos2 == -1) {
              Type1feasible = false;
            }
          }
          for (int j = 0; j < 2; j++) {
            size_t iEdge1, iEdge2;
            if (j == 0) {
              iEdge1 = iEdgePrev;
              iEdge2 = iEdgeNext;
            } else {
              iEdge1 = iEdgeNext;
              iEdge2 = iEdgePrev;
            }
            if (MCEwork(dEdge, iEdge1) == SetOne &&
                MCEwork(dEdge, iEdge2) == SetZero) {
              if (pos2 != -1) {
                std::vector<int> NewVal =
                    DifferenceVect(MCEwork(dEdge, iEdge), SetTwo);
                MCEwork(dEdge, iEdge) = NewVal;
                MCEwork(iEdge, dEdge) = NewVal;
                nbOper++;
              }
              if (pos1 == -1 && !Type2feasible) {
                std::cerr << "Finding false 1\n";
                return false;
              }
              if (MCEwork(dEdge, iEdge) != SetOne && !Type2feasible) {
                MCEwork(dEdge, iEdge) = SetOne;
                MCEwork(iEdge, dEdge) = SetOne;
                nbOper++;
              }
            }
          }
          if (MCEwork(dEdge, iEdgePrev) == SetZero &&
              MCEwork(dEdge, iEdgeNext) == SetZero) {
            if (pos2 != -1) {
              std::vector<int> NewVal =
                  DifferenceVect(MCEwork(dEdge, iEdge), SetTwo);
              MCEwork(dEdge, iEdge) = NewVal;
              MCEwork(iEdge, dEdge) = NewVal;
              nbOper++;
            }
            if (!Type2feasible && pos1 != -1) {
              std::vector<int> NewVal =
                  DifferenceVect(MCEwork(dEdge, iEdge), SetOne);
              MCEwork(dEdge, iEdge) = NewVal;
              MCEwork(iEdge, dEdge) = NewVal;
              nbOper++;
            }
            if (MCEwork(dEdge, iEdge) == SetOne) {
              Type1feasible = false;
            }
          }
        }
      }
      if (!Type1feasible && !Type2feasible)
        return false;
    }
    return true;
  };
  auto IsAprioriFeasible =
      [&nbEdge](MyMatrix<std::vector<int>> const &eMat) -> bool {
    for (size_t iEdge = 0; iEdge < nbEdge - 1; iEdge++)
      for (size_t jEdge = iEdge + 1; jEdge < nbEdge; jEdge++) {
        size_t siz = eMat(iEdge, jEdge).size();
        if (siz == 0)
          return false;
      }
    return true;
  };
  auto IterativeRefinementMCE =
      [&](MyMatrix<std::vector<int>> &MCEwork) -> bool {
    bool test;
    while (true) {
      int nbOper = 0;
      test = RefinementBlock1(MCEwork, nbOper);
      std::cerr << "RefinementBlock1 test=" << test << " nbOper=" << nbOper
                << "\n";
      if (!test)
        return false;
      test = RefinementBlock2(MCEwork, nbOper);
      std::cerr << "RefinementBlock2 test=" << test << " nbOper=" << nbOper
                << "\n";
      if (!test)
        return false;
      test = RefinementTriangles(MCEwork, nbOper);
      std::cerr << "RefinementTriangles test=" << test << " nbOper=" << nbOper
                << "\n";
      if (!test)
        return false;
      test = IsAprioriFeasible(MCEwork);
      std::cerr << "IsAprioriFeasible test=" << test << "\n";
      if (!test)
        return false;
      if (nbOper == 0)
        break;
    }
    return true;
  };
  auto GetFirstNonZero =
      [&nbEdge](
          MyMatrix<std::vector<int>> const &MCEwork) -> std::vector<size_t> {
    for (size_t iEdge = 0; iEdge < nbEdge - 1; iEdge++)
      for (size_t jEdge = iEdge + 1; jEdge < nbEdge; jEdge++) {
        size_t siz = MCEwork(iEdge, jEdge).size();
        if (siz > 1)
          return {iEdge, jEdge};
      }
    std::cerr << "Fail to find an index\n";
    throw TerminalException{1};
  };
  auto IsFinalState =
      [&nbEdge](MyMatrix<std::vector<int>> const &MCEwork) -> bool {
    for (size_t iEdge = 0; iEdge < nbEdge - 1; iEdge++)
      for (size_t jEdge = iEdge + 1; jEdge < nbEdge; jEdge++) {
        size_t siz = MCEwork(iEdge, jEdge).size();
        if (siz > 1)
          return false;
      }
    return true;
  };
  auto GetEmbeddingMatrix =
      [&ListEdge, &GetGraphIdentityEdge, &GraphIntersectionOne,
       &ExtendListLabel,
       &CheckEmbedding](MyMatrix<std::vector<int>> const &MCEwork,
                        MyMatrix<int> &FinalEmbed) -> bool {
    std::cerr << "GetEmbeddingMatrix, step 1\n";
    std::cerr << "GetEmbeddingMatrix, step 2\n";
    GraphBitset Aspec = GetGraphIdentityEdge(ListEdge, MCEwork);
    std::cerr << "GetEmbeddingMatrix, step 3\n";
    std::vector<std::vector<size_t>> Cspec = ConnectedComponents_set(Aspec);
    std::cerr << "GetEmbeddingMatrix, step 4 |Cspec|=" << Cspec.size() << "\n";
    GraphBitset Dspec = GraphIntersectionOne(Cspec, ListEdge, MCEwork);
    std::vector<std::vector<int>> PartListLabel = InverseLineGraph(Dspec);
    std::cerr << "GetEmbeddingMatrix, step 6\n";
    if (PartListLabel[0][0] == -1) {
      std::cerr << "Embedding failure\n";
      return false;
    }
    std::vector<std::vector<int>> ListLabel =
        ExtendListLabel(Cspec, PartListLabel);
    std::cerr << "GetEmbeddingMatrix, step 7\n";
    FinalEmbed = CreateEmbedding(0, ListEdge, ListLabel);
    bool res = CheckEmbedding(FinalEmbed);
    return res;
  };
  struct OneLevel {
    std::vector<size_t> PairEdge;
    std::vector<int> ListPoss;
    int idx;
    MyMatrix<std::vector<int>> MCE;
  };
  std::vector<OneLevel> TheTree;
  auto GoUpNextInTree = [&]() -> bool {
    int idx;
    while (true) {
      if (TheTree.size() == 0)
        return false;
      size_t len = TheTree.size();
      std::cerr << "TheTree len=" << len << "\n";
      int lenPoss = static_cast<int>(TheTree[len - 1].ListPoss.size());
      std::cerr << "lenPoss=" << lenPoss << "\n";
      if (len == 1) {
        std::cerr << "ListPoss=";
        for (auto &eVal : TheTree[len - 1].ListPoss)
          std::cerr << " " << eVal;
        std::cerr << "\n";
      }
      idx = TheTree[len - 1].idx;
      std::cerr << "idx=" << idx << "\n";
      if (idx == lenPoss - 1) {
        TheTree.pop_back();
      } else {
        idx++;
        size_t iEdge = TheTree[len - 1].PairEdge[0];
        size_t jEdge = TheTree[len - 1].PairEdge[1];
        TheTree[len - 1].idx = idx;
        MyMatrix<std::vector<int>> MCEnew;
        if (len == 1) {
          MCEnew = MCE;
        } else {
          MCEnew = TheTree[len - 2].MCE;
        }
        if (len == 1) {
          int eVal = TheTree[len - 1].ListPoss[idx];
          std::cerr << "1: Now operating with ListPoss=" << eVal << "\n";
        }
        int eVal = TheTree[len - 1].ListPoss[idx];
        MCEnew(iEdge, jEdge) = {eVal};
        bool test = IterativeRefinementMCE(MCEnew);
        TheTree[len - 1].MCE = MCEnew;
        if (test)
          return true;
      }
    }
  };
  auto GetLatestMCE = [&TheTree, &MCE]() -> MyMatrix<std::vector<int>> {
    size_t len = TheTree.size();
    if (len == 0)
      return MCE;
    return TheTree[len - 1].MCE;
  };
  auto NextInTree = [&]() -> bool {
    MyMatrix<std::vector<int>> MCEstart = GetLatestMCE();
    bool test = IsFinalState(MCEstart);
    if (test)
      return GoUpNextInTree();
    std::vector<size_t> ePairEdge = GetFirstNonZero(MCEstart);
    size_t iEdge = ePairEdge[0];
    size_t jEdge = ePairEdge[1];
    std::cerr << "Assign iEdge=" << iEdge << " jEdge=" << jEdge << "\n";
    std::vector<int> ListPoss = MCEstart(iEdge, jEdge);
    int idx = 0;
    MCEstart(iEdge, jEdge) = {ListPoss[idx]};
    MCEstart(jEdge, iEdge) = {ListPoss[idx]};
    TheTree.push_back({ePairEdge, ListPoss, idx, MCEstart});
    size_t len = TheTree.size();
    bool test2 = IterativeRefinementMCE(TheTree[len - 1].MCE);
    if (!test2)
      return GoUpNextInTree();
    return true;
  };
  bool test = IterativeRefinementMCE(MCE);
  std::cerr << "IterativeRefinementMCE test=" << test << "\n";
  if (!test)
    return {};
  std::vector<MyMatrix<int>> ListEmbedding;
  size_t nbChoice = nbEdge * (nbEdge - 1) / 2;
  while (true) {
    MyMatrix<std::vector<int>> MCEstart = GetLatestMCE();
    bool testF = IsFinalState(MCEstart);
    std::cerr << "IsFinalState test=" << testF << "\n";
    if (testF) {
      MyMatrix<int> FinalEmbed;
      bool test2 = GetEmbeddingMatrix(MCEstart, FinalEmbed);
      std::cerr << "test2=" << test2 << "\n";
      if (test2)
        ListEmbedding.push_back(FinalEmbed);
    }
    std::cerr << "Passing by NextInTree nbEmbedding=" << ListEmbedding.size()
              << "\n";
    std::cerr << "|TheTree|=" << TheTree.size() << " nbChoice=" << nbChoice
              << " iter=" << iter << "\n";
    if (iter > MaxIter)
      break;
    bool test2 = NextInTree();
    if (!test2)
      break;
    iter++;
  }
  std::cerr << "|ListEmbedding|=" << ListEmbedding.size() << "\n";
  return ListEmbedding;
}

// clang-format off
#endif  // SRC_GRAPH_GRAPH_GRAPHICALFUNCTIONS_H_
// clang-format on
