// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GRAPH_GRAPH_GRAPHICALBASIC_H_
#define SRC_GRAPH_GRAPH_GRAPHICALBASIC_H_

// clang-format off
#include "MAT_Matrix.h"
#include "Temp_common.h"
#include "Boost_bitset.h"
#include "COMB_Combinatorics_elem.h"
#include <limits>
#include <string>
#include <utility>
#include <vector>
// clang-format on

struct GraphListAdj {
public:
  GraphListAdj() = delete;
  GraphListAdj(size_t const &inpNbVert) : nbVert(inpNbVert) {
    HasVertexColor = false;
    for (size_t iVert = 0; iVert < nbVert; iVert++)
      ListListAdj.push_back({});
  }
  GraphListAdj(std::vector<std::pair<size_t, size_t>> const &ListEdge,
               size_t const &inpNbVert) {
    HasVertexColor = false;
    nbVert = inpNbVert;
    for (size_t iVert = 0; iVert < nbVert; iVert++)
      ListListAdj.push_back({});
    //
    for (auto &eEdge : ListEdge) {
      size_t eVert1 = eEdge.first;
      size_t eVert2 = eEdge.second;
      ListListAdj[eVert1].push_back(eVert2);
      ListListAdj[eVert2].push_back(eVert1);
    }
  }
  GraphListAdj(MyMatrix<size_t> const &ListEdge, size_t const &inpNbVert) {
    HasVertexColor = false;
    nbVert = inpNbVert;
    for (size_t iVert = 0; iVert < nbVert; iVert++)
      ListListAdj.push_back({});
    //
    size_t nbEdge = ListEdge.rows();
    for (size_t iEdge = 0; iEdge < nbEdge; iEdge++) {
      size_t eVert1 = ListEdge(iEdge, 0);
      size_t eVert2 = ListEdge(iEdge, 1);
      ListListAdj[eVert1].push_back(eVert2);
      ListListAdj[eVert2].push_back(eVert1);
    }
  }
  ~GraphListAdj() {}
  GraphListAdj(GraphListAdj const &eG) {
    nbVert = eG.GetNbVert();
    ListListAdj = eG.GetListListAdj();
    HasVertexColor = eG.GetHasVertexColor();
    ListVertexColor = eG.GetListVertexColor();
  }
  GraphListAdj operator=(GraphListAdj const &eG) {
    nbVert = eG.GetNbVert();
    ListListAdj = eG.GetListListAdj();
    HasVertexColor = eG.GetHasVertexColor();
    ListVertexColor = eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  size_t GetNbAdjacent() const {
    size_t nb_adj = 0;
    for (auto &LAdj : ListListAdj)
      nb_adj += LAdj.size();
    return nb_adj;
  }
  size_t GetNbVert() const { return nbVert; }
  std::vector<std::vector<size_t>> GetListListAdj() const {
    return ListListAdj;
  }
  bool GetHasVertexColor() const { return HasVertexColor; }
  std::vector<size_t> GetListVertexColor() const { return ListVertexColor; }
  //
  void SetHasColor(bool const &TheVal) {
    if (TheVal == HasVertexColor)
      return;
    HasVertexColor = TheVal;
    if (TheVal)
      ListVertexColor = std::vector<size_t>(nbVert);
    if (!TheVal)
      ListVertexColor.clear();
  }
  void SetColor(size_t const &iVert, size_t const &eColor) {
    ListVertexColor[iVert] = eColor;
  }
  std::vector<size_t> Adjacency(size_t const &iVert) const {
    return ListListAdj[iVert];
  }
  void AddAdjacent(size_t const &iVert, size_t const &jVert) {
    ListListAdj[iVert].push_back(jVert);
  }
  void RemoveAdjacent(size_t const &iVert, size_t const &jVert) {
    // Not sure to find the right simple code for removing entry in std::vector
    std::vector<size_t> NewList;
    for (auto &eVert : ListListAdj[iVert])
      if (size_t(eVert) != jVert)
        NewList.push_back(eVert);
    ListListAdj[iVert] = NewList;
  }
  bool IsAdjacent(size_t const &iVert, size_t const &jVert) const {
    for (auto &eVert : ListListAdj[iVert])
      if (eVert == jVert)
        return true;
    return false;
  }
  size_t GetColor(size_t const &iVert) const {
    if (!HasVertexColor) {
      std::cerr << "Call to GetColor while HasVertexColor=false\n";
      throw TerminalException{1};
    }
    return ListVertexColor[iVert];
  }

private:
  size_t nbVert;
  std::vector<std::vector<size_t>> ListListAdj;
  bool HasVertexColor;
  std::vector<size_t> ListVertexColor;
};

struct GraphSparseImmutable {
public:
  GraphSparseImmutable() = delete;
  GraphSparseImmutable(size_t const &_nbVert,
                       std::vector<size_t> const &_ListStart,
                       std::vector<size_t> const &_ListListAdj)
      : nbVert(_nbVert), ListStart(_ListStart), ListListAdj(_ListListAdj) {
    HasVertexColor = false;
  }
  ~GraphSparseImmutable() {}
  GraphSparseImmutable(MyMatrix<size_t> const &ListEdge, size_t const &_nbVert)
      : nbVert(_nbVert) {
    HasVertexColor = false;
    size_t nbEdge = ListEdge.rows();
    std::vector<size_t> ListDeg(nbVert, 0);
    for (size_t iEdge = 0; iEdge < nbEdge; iEdge++)
      for (size_t i = 0; i < 2; i++) {
        size_t eVert = ListEdge(iEdge, i);
        ListDeg[eVert]++;
      }
    ListStart.resize(nbVert + 1, 0);
    size_t nbAdj = 0;
    for (size_t iVert = 0; iVert < nbVert; iVert++) {
      size_t eDeg = ListDeg[iVert];
      ListStart[iVert + 1] = ListStart[iVert] + eDeg;
      nbAdj += eDeg;
    }
    std::vector<size_t> ListPos(nbVert, 0);
    ListListAdj.resize(nbAdj, 0);
    for (size_t iEdge = 0; iEdge < nbEdge; iEdge++)
      for (size_t i = 0; i < 2; i++) {
        size_t j = 1 - i;
        size_t eVert = ListEdge(iEdge, i);
        size_t fVert = ListEdge(iEdge, j);
        size_t eStart = ListStart[eVert];
        size_t pos = ListPos[eVert];
        ListListAdj[eStart + pos] = fVert;
        ListPos[eVert]++;
      }
  }
  GraphSparseImmutable(GraphSparseImmutable const &eG) {
    nbVert = eG.GetNbVert();
    ListStart = eG.GetListStart();
    ListListAdj = eG.GetListListAdj();
    HasVertexColor = eG.GetHasVertexColor();
    ListVertexColor = eG.GetListVertexColor();
  }
  GraphSparseImmutable operator=(GraphSparseImmutable const &eG) {
    nbVert = eG.GetNbVert();
    ListStart = eG.GetListStart();
    ListListAdj = eG.GetListListAdj();
    HasVertexColor = eG.GetHasVertexColor();
    ListVertexColor = eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  size_t GetNbVert() const { return nbVert; }
  std::vector<size_t> GetListListAdj() const { return ListListAdj; }
  std::vector<size_t> GetListStart() const { return ListStart; }
  bool GetHasVertexColor() const { return HasVertexColor; }
  std::vector<size_t> GetListVertexColor() const { return ListVertexColor; }
  //
  void SetHasColor(bool const &TheVal) {
    if (TheVal == HasVertexColor) {
      return;
    }
    HasVertexColor = TheVal;
    if (TheVal)
      ListVertexColor = std::vector<size_t>(nbVert);
    if (!TheVal)
      ListVertexColor.clear();
  }
  void SetColor(size_t const &iVert, size_t const &eColor) {
    ListVertexColor[iVert] = eColor;
  }
  std::vector<size_t> Adjacency(size_t const &iVert) const {
    size_t eStart = ListStart[iVert];
    size_t eEnd = ListStart[iVert + 1];
    std::vector<size_t> TheRet;
    for (size_t i = eStart; i < eEnd; i++)
      TheRet.push_back(ListListAdj[i]);
    return TheRet;
  }
  bool IsAdjacent(size_t const &iVert, size_t const &jVert) const {
    size_t eStart = ListStart[iVert];
    size_t eEnd = ListStart[iVert + 1];
    for (size_t i = eStart; i < eEnd; i++)
      if (ListListAdj[i] == jVert)
        return true;
    return false;
  }
  size_t GetColor(size_t const &iVert) const {
    if (!HasVertexColor) {
      std::cerr << "Call to GetColor while HasVertexColor=false\n";
      throw TerminalException{1};
    }
    return ListVertexColor[iVert];
  }
  std::pair<std::vector<size_t>, std::vector<size_t>>
  Get_ListStart_ListListAdj() const {
    return {ListStart, ListListAdj};
  }
  size_t GetIndex(size_t const &eVert, size_t const &fVert) const {
    size_t eStart = ListStart[eVert];
    size_t eEnd = ListStart[eVert + 1];
    for (size_t i = eStart; i < eEnd; i++)
      if (ListListAdj[i] == fVert)
        return i;
    return std::numeric_limits<size_t>::max();
  }

private:
  size_t nbVert;
  std::vector<size_t> ListStart;
  std::vector<size_t> ListListAdj;
  bool HasVertexColor;
  std::vector<size_t> ListVertexColor;
};

template <> struct is_graphsparseimmutable_class<GraphSparseImmutable> {
  static const bool value = true;
};

template <typename Tgr>
std::vector<size_t> ConnectedComponents_vector(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::vector<size_t> ListStatus(nbVert, miss_val);
  size_t iStatus = 0;
  auto Assignment = [&](size_t const &iVert) -> void {
    std::vector<size_t> TheSet{iVert};
    ListStatus[iVert] = iStatus;
    while (true) {
      std::vector<size_t> NewSet;
      for (size_t &eVal : TheSet) {
        for (size_t &fVal : GR.Adjacency(eVal)) {
          if (ListStatus[fVal] == miss_val) {
            ListStatus[fVal] = iStatus;
            NewSet.push_back(fVal);
          }
        }
      }
      if (NewSet.size() == 0)
        break;
      TheSet = NewSet;
    }
    iStatus++;
  };
  for (size_t iVert = 0; iVert < nbVert; iVert++)
    if (ListStatus[iVert] == miss_val)
      Assignment(iVert);
  return ListStatus;
}

template<typename Tgr>
bool IsConnectedGraphMinusSubset(Tgr const &GR, std::vector<size_t> const& V) {
  size_t nbVert = GR.GetNbVert();
  size_t miss_val = std::numeric_limits<size_t>::max();
  Face f_out(nbVert);
  for (auto & eVert: V) {
    f_out[eVert] = 1;
  }
  auto get_out_vert=[&]() -> size_t {
    for (size_t u=0; u<nbVert; u++) {
      if (f_out[u] == 0) {
        return u;
      }
    }
    return miss_val;
  };
  size_t out_vert = get_out_vert();
  if (out_vert == miss_val) {
    return true;
  }
  Face fatt(nbVert);
  std::vector<size_t> Lvert(nbVert);
  Lvert[0] = out_vert;
  size_t vert_start = 0;
  size_t vert_end = 1;
  fatt[out_vert] = 1;
  while(true) {
    size_t pos_vert = vert_end;
    for (size_t u=vert_start; u<vert_end; u++) {
      for (auto & eAdj : GR.Adjacency(Lvert[u])) {
        if (f_out[eAdj] == 0 && fatt[eAdj] == 0) {
          Lvert[pos_vert] = eAdj;
          fatt[eAdj] = 1;
          pos_vert += 1;
        }
      }
    }
    vert_start = vert_end;
    vert_end = pos_vert;
    if (vert_start == vert_end) {
      break;
    }
  }
  size_t total_len = nbVert - V.size();
  if (vert_start == total_len) {
    return true;
  }
  return false;
}

template<typename Tgr>
bool IsKConnectedGraph(Tgr const &GR, size_t const& k) {
  size_t nbVert = GR.GetNbVert();
  SetCppIterator set(nbVert, k-1);
  for (auto & V : set) {
    std::vector<size_t> V2;
    for (auto & val : V) {
      V2.push_back(val);
    }
    bool test = IsConnectedGraphMinusSubset(GR, V2);
    if (!test) {
      return false;
    }
  }
  return true;
}


template <typename Tgr>
std::vector<std::vector<size_t>> ConnectedComponents_set(Tgr const &GR) {
  size_t nbVert = GR.GetNbVert();
  std::vector<size_t> ListStatus = ConnectedComponents_vector(GR);
  size_t nbConn = VectorMax(ListStatus) + 1;
  std::vector<std::vector<size_t>> ListConn(nbConn);
  for (size_t iVert = 0; iVert < nbVert; iVert++) {
    size_t iStatus = ListStatus[iVert];
    ListConn[iStatus].push_back(iVert);
  }
  return ListConn;
}

template <typename Tret, typename Tgr>
inline typename std::enable_if<is_graphsparseimmutable_class<Tret>::value,
                               Tret>::type
InducedSubgraph(Tgr const &GR, std::vector<int> const &eList) {
  size_t nbVert = GR.GetNbVert();
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::vector<size_t> ListStat(nbVert, miss_val);
  size_t nbVertRed = eList.size();
  for (size_t iVertRed = 0; iVertRed < nbVertRed; iVertRed++) {
    size_t iVert = eList[iVertRed];
    ListStat[iVert] = iVertRed;
  }
  std::vector<MyVector<size_t>> PreListEdge;
  for (size_t iVertRed = 0; iVertRed < nbVertRed; iVertRed++) {
    size_t iVert = eList[iVertRed];
    std::vector<size_t> LLadj = GR.Adjacency(iVert);
    for (size_t &jVert : LLadj) {
      size_t jVertRed = ListStat[jVert];
      if (jVertRed != miss_val && jVertRed > iVertRed) {
        MyVector<size_t> eVect(2);
        eVect(0) = iVertRed;
        eVect(1) = jVertRed;
        PreListEdge.push_back(eVect);
      }
    }
  }
  MyMatrix<size_t> ListEdge = MatrixFromVectorFamily(PreListEdge);
  return GraphSparseImmutable(ListEdge, nbVertRed);
}

template <typename Tgr>
std::string GetCanonicalForm_string(Tgr const &eGR,
                                    std::vector<unsigned int> const &cl) {
  size_t nof_vertices = eGR.GetNbVert();
  std::vector<size_t> clR(nof_vertices);
  for (size_t i = 0; i < nof_vertices; i++)
    clR[cl[i]] = i;

  std::string strRet;
  for (size_t iVert = 0; iVert < nof_vertices; iVert++) {
    size_t iVertCan = clR[iVert];
    strRet += std::to_string(iVert) + std::string(" : ");
    if (eGR.GetHasVertexColor()) {
      size_t eColor = eGR.GetColor(iVertCan);
      strRet += " " + std::to_string(eColor);
    }
    strRet += " : ";
    //
    for (size_t jVert = 0; jVert < nof_vertices; jVert++) {
      size_t jVertCan = clR[jVert];
      bool eVal_b = eGR.IsAdjacent(iVertCan, jVertCan);
      strRet += " " + std::to_string(eVal_b);
    }
  }
  return strRet;
}

// clang-format off
#endif  // SRC_GRAPH_GRAPH_GRAPHICALBASIC_H_
// clang-format on
