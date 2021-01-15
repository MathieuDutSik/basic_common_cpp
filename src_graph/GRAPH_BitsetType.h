#ifndef INCLUDE_GRAPH_BITSET
#define INCLUDE_GRAPH_BITSET


#include "Boost_bitset.h"


struct GraphBitset {
public:
  GraphBitset() = delete;
  GraphBitset(size_t const& inpNbVert) : nbVert(inpNbVert)
  {
    size_t eProd = nbVert * nbVert;
#ifdef DEBUG
    std::cerr << "nbVert=" << nbVert << " eProd=" << eProd << "\n";
#endif
    LLAdj=boost::dynamic_bitset<>(eProd);
    HasVertexColor=false;
  }
  ~GraphBitset()
  {
    ~LLAdj;
    if (HasVertexColor)
      ListVertexColor.clear();
  }
  GraphBitset(GraphBitset const& eG)
  {
    nbVert=eG.GetNbVert();
    LLAdj=eG.GetBitset();
    HasVertexColor=eG.GetHasVertexColor();
    if (HasVertexColor)
      ListVertexColor=eG.GetListVertexColor();
  }
  GraphBitset operator=(GraphBitset const& eG)
  {
    nbVert=eG.GetNbVert();
    LLAdj=eG.GetBitset();
    HasVertexColor=eG.GetHasVertexColor();
    if (HasVertexColor)
      ListVertexColor=eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  size_t GetNbAdjacent() const
  {
    return LLAdj.count();
  }
  size_t GetNbVert() const
  {
    return nbVert;
  }
  boost::dynamic_bitset<> GetBitset() const
  {
    return LLAdj;
  }
  bool GetHasVertexColor() const
  {
    return HasVertexColor;
  }
  std::vector<int> GetListVertexColor() const
  {
    return ListVertexColor;
  }
  //
  void SetHasColor(bool const& TheVal)
  {
    if (TheVal == HasVertexColor) {
      return;
    }
    HasVertexColor=TheVal;
    if (TheVal)
      ListVertexColor = std::vector<int>(nbVert);
    if (!TheVal)
      ListVertexColor.clear();
  }
  void SetColor(size_t const& iVert, int const& eColor)
  {
    ListVertexColor[iVert]=eColor;
  }
  std::vector<int> Adjacency(size_t const& iVert) const
  {
    std::vector<int> retList;
    for (size_t jVert=0; jVert<nbVert; jVert++) {
      size_t idxMat=iVert + nbVert*jVert;
      if (LLAdj[idxMat])
	retList.push_back(jVert);
    }
    return retList;
  }
  void AddAdjacent(size_t const& iVert, size_t const& jVert)
  {
    size_t idxMat=iVert + nbVert*jVert;
    LLAdj.set(idxMat,1);
  }
  void RemoveAdjacent(size_t const& iVert, size_t const& jVert)
  {
    size_t idxMat=iVert + nbVert*jVert;
    LLAdj.set(idxMat,1);
  }
  bool IsAdjacent(size_t const& iVert, size_t const& jVert) const
  {
    size_t idxMat=iVert + nbVert*jVert;
    if (LLAdj[idxMat])
      return true;
    return false;
  }
  int GetColor(size_t const& iVert) const
  {
    if (!HasVertexColor) {
      std::cerr << "Call to GetColor while HasVertexColor=false\n";
      throw TerminalException{1};
    }
    return ListVertexColor[iVert];
  }
private:
  size_t nbVert;
  boost::dynamic_bitset<> LLAdj;
  bool HasVertexColor=false;
  std::vector<int> ListVertexColor;
};

template<>
struct is_graphbitset_class<GraphBitset> {
  static const bool value = true;
};


template<typename Tret, typename Tgr>
inline typename std::enable_if<is_graphbitset_class<Tret>::value,Tret>::type InducedSubgraph(Tgr const& GR, std::vector<int> const& eList)
{
  size_t nbVert=GR.GetNbVert();
  std::vector<int> ListStat(nbVert,-1);
  size_t nbVertRed=eList.size();
  GraphBitset GRred(nbVertRed);
  for (size_t iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    size_t iVert=eList[iVertRed];
    ListStat[iVert]=iVertRed;
  }
  for (size_t iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    size_t iVert=eList[iVertRed];
    std::vector<int> LLadj=GR.Adjacency(iVert);
    for (int & jVert : LLadj) {
      int jVertRed=ListStat[jVert];
      if (jVertRed != -1) {
	GRred.AddAdjacent(iVertRed, jVertRed);
	GRred.AddAdjacent(jVertRed, iVertRed);
      }
    }
  }
  return GRred;
}



#endif
