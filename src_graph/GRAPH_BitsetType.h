#ifndef SRC_GRAPH_GRAPH_BITSETTYPE_H_
#define SRC_GRAPH_GRAPH_BITSETTYPE_H_

#include "Boost_bitset.h"

struct GraphBitset {
public:
  GraphBitset() = delete;
  GraphBitset(size_t const &inpNbVert) : nbVert(inpNbVert) {
    size_t eProd = nbVert * nbVert;
#ifdef DEBUG
    std::cerr << "nbVert=" << nbVert << " eProd=" << eProd << "\n";
#endif
    LLAdj = boost::dynamic_bitset<>(eProd);
    HasVertexColor = false;
  }
  ~GraphBitset() {
    ~LLAdj;
    if (HasVertexColor)
      ListVertexColor.clear();
  }
  GraphBitset(GraphBitset const &eG) {
    nbVert = eG.GetNbVert();
    LLAdj = eG.GetBitset();
    HasVertexColor = eG.GetHasVertexColor();
    if (HasVertexColor)
      ListVertexColor = eG.GetListVertexColor();
  }
  GraphBitset operator=(GraphBitset const &eG) {
    nbVert = eG.GetNbVert();
    LLAdj = eG.GetBitset();
    HasVertexColor = eG.GetHasVertexColor();
    if (HasVertexColor)
      ListVertexColor = eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  size_t GetNbAdjacent() const { return LLAdj.count(); }
  size_t GetNbVert() const { return nbVert; }
  boost::dynamic_bitset<> GetBitset() const { return LLAdj; }
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
    std::vector<size_t> retList;
    for (size_t jVert = 0; jVert < nbVert; jVert++) {
      size_t idxMat = iVert + nbVert * jVert;
      if (LLAdj[idxMat])
        retList.push_back(jVert);
    }
    return retList;
  }
  void AddAdjacent(size_t const &iVert, size_t const &jVert) {
    size_t idxMat = iVert + nbVert * jVert;
    LLAdj.set(idxMat, 1);
  }
  void RemoveAdjacent(size_t const &iVert, size_t const &jVert) {
    size_t idxMat = iVert + nbVert * jVert;
    LLAdj.set(idxMat, 1);
  }
  bool IsAdjacent(size_t const &iVert, size_t const &jVert) const {
    size_t idxMat = iVert + nbVert * jVert;
    if (LLAdj[idxMat])
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
  boost::dynamic_bitset<> LLAdj;
  bool HasVertexColor = false;
  std::vector<size_t> ListVertexColor;
};

template <> struct is_graphbitset_class<GraphBitset> {
  static const bool value = true;
};

template <typename Tret, typename Tgr>
inline typename std::enable_if<is_graphbitset_class<Tret>::value, Tret>::type
InducedSubgraph(Tgr const &GR, std::vector<size_t> const &eList) {
  size_t nbVert = GR.GetNbVert();
  std::vector<size_t> ListStat(nbVert, std::numeric_limits<size_t>::max());
  size_t nbVertRed = eList.size();
  GraphBitset GRred(nbVertRed);
  for (size_t iVertRed = 0; iVertRed < nbVertRed; iVertRed++) {
    size_t iVert = eList[iVertRed];
    ListStat[iVert] = iVertRed;
  }
  for (size_t iVertRed = 0; iVertRed < nbVertRed; iVertRed++) {
    size_t iVert = eList[iVertRed];
    std::vector<size_t> LLadj = GR.Adjacency(iVert);
    for (size_t &jVert : LLadj) {
      size_t jVertRed = ListStat[jVert];
      if (jVertRed != std::numeric_limits<size_t>::max()) {
        GRred.AddAdjacent(iVertRed, jVertRed);
        GRred.AddAdjacent(jVertRed, iVertRed);
      }
    }
  }
  return GRred;
}

#endif // SRC_GRAPH_GRAPH_BITSETTYPE_H_
