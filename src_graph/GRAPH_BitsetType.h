#ifndef INCLUDE_GRAPH_BITSET
#define INCLUDE_GRAPH_BITSET


#include "Boost_bitset.h"


struct GraphBitset {
public:
  GraphBitset() = delete;
  GraphBitset(int const& inpNbVert) : nbVert(inpNbVert)
  {
    LLAdj=boost::dynamic_bitset<>(nbVert*nbVert);
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
  int GetNbVert() const
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
      ListVertexColor=std::vector<int>(nbVert);
    if (!TheVal)
      ListVertexColor.clear();
  }
  void SetColor(int const& iVert, int const& eColor)
  {
    ListVertexColor[iVert]=eColor;
  }
  std::vector<int> Adjacency(int const& iVert) const
  {
    std::vector<int> retList;
    for (int jVert=0; jVert<nbVert; jVert++) {
      int idxMat=iVert + nbVert*jVert;
      if (LLAdj[idxMat])
	retList.push_back(jVert);
    }
    return retList;
  }
  void AddAdjacent(int const& iVert, int const& jVert)
  {
    int idxMat=iVert + nbVert*jVert;
    LLAdj.set(idxMat,1);
  }
  void RemoveAdjacent(int const& iVert, int const& jVert)
  {
    int idxMat=iVert + nbVert*jVert;
    LLAdj.set(idxMat,1);
  }
  bool IsAdjacent(int const& iVert, int const& jVert) const
  {
    int idxMat=iVert + nbVert*jVert;
    if (LLAdj[idxMat])
      return true;
    return false;
  }
  int GetColor(int const& iVert) const
  {
    if (!HasVertexColor) {
      std::cerr << "Call to GetColor while HasVertexColor=false\n";
      throw TerminalException{1};
    }
    return ListVertexColor[iVert];
  }
private:
  int nbVert;
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
  int nbVert=GR.GetNbVert();
  std::vector<int> ListStat(nbVert,-1);
  int nbVertRed=eList.size();
  GraphBitset GRred(nbVertRed);
  for (int iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    int iVert=eList[iVertRed];
    ListStat[iVert]=iVertRed;
  }
  for (int iVertRed=0; iVertRed<nbVertRed; iVertRed++) {
    int iVert=eList[iVertRed];
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
