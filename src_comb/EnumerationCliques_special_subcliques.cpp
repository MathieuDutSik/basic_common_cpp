#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

// This is a variant (INCOMPLETE CODE) of the clique enumeration
//
// We try to insert the cliques in a more ordered way. The idea
// is that if we insert one element, then we can just as well
// insert some other.
// In the case of 2-dimensional Keller-like clique packing, this
// would give you the packing done line by line.
// QUESTION: Is this compatible with the ordered enumeration?
// ANSWER: 1) If we insert set S = {v1, ...., vM}
//            then the minimum starting point at the next step
//            ought to be v1+1, not vM+1.
//         2) Anything less to consider? It would seem no.
//
// Structure of set of points to add.
// When considering the expansion process, we have a number of possibilities V
// F(x) = {y in V s.t. x adj y and for all z in V we have z adj x <=> z adj y}
// So, we can insert the set F(x) right away in the enumeration process. and not
// just x.
//
// Another optimization is that if a point x in V satisfies {y adj x for all y in V}
// then we can insert it right away.
//
// The algorithm needs to keep track of the insertions.
// We need a data structure for storing the partitions of point sets.
// or do we need to store the partition? We only need to keep track
// of which points have been considered in the process.
// So, we just need to have a vector of point status.


struct GraphType {
  int nbPoint;
  std::vector<std::vector<int>> LLAdj;
};

GraphType ReadGraphFile(std::string const& eFile)
{
  std::ifstream is(eFile);
  int nbPoint;
  is >> nbPoint;
  std::vector<std::vector<int>> LLAdj;
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    std::vector<int> LAdj;
    for (int jPoint=0; jPoint<nbPoint; jPoint++) {
      int eVal;
      is >> eVal;
      LAdj.push_back(eVal);
    }
    LLAdj.push_back(LAdj);
  }
  return {nbPoint, LLAdj};
}

struct GroupType {
  int nbElt;
  int nbPoint;
  std::vector<std::vector<int>> ARR;
};



GroupType ReadGroupFile(std::string const& eFile)
{
  std::ifstream is(eFile);
  int nbPoint, nbElt;
  is >> nbPoint;
  is >> nbElt;
  std::vector<std::vector<int>> ARR;
  for (int iElt=0; iElt<nbElt; iElt++) {
    std::vector<int> LAdj;
    for (int iPoint=0; iPoint<nbPoint; iPoint++) {
      int eVal;
      is >> eVal;
      LAdj.push_back(eVal);
    }
    ARR.push_back(LAdj);
  }
  return {nbElt, nbPoint, ARR};
}

struct VectRec {
  int len;
  std::vector<int> V;
};

void PrintVector(std::string const& eStr, VectRec const& eVect)
{
  std::cerr << eStr;
  for (int idx=0; idx<eVect.len; idx++)
    std::cerr << " " << eVect.V[idx];
  std::cerr << "\n";
}


bool IsMinimal(GroupType const& eGroup, VectRec const& eVect)
{
  std::vector<int> eVectOrd(eVect.len);
  for (int i=0; i<eVect.len; i++)
    eVectOrd[i] = eVect.V[i];
  std::sort(eVectOrd.begin(), eVectOrd.end());
  std::vector<int> ImageVect(eVect.len);
  auto IsCounterexample=[&](std::vector<int> const& uVect) -> bool {
    for (int i=0; i<eVect.len; i++) {
      if (uVect[i] < eVectOrd[i])
        return true;
      if (uVect[i] > eVectOrd[i]) // for the next item to be counterexample, it has to be equal
        return false;
    }
    return false;
  };
  for (int iElt=0; iElt<eGroup.nbElt; iElt++) {
    for (int i=0; i<eVect.len; i++) {
      int iImg=eGroup.ARR[iElt][eVect.V[i]];
      ImageVect[i]=iImg;
    }
    std::sort(ImageVect.begin(), ImageVect.end());
    if (IsCounterexample(ImageVect))
      return false;
  }
  return true;
}



// Now eVect at level i is no longer of length i. Can be larger.
struct OneLevel {
  VectRec eVect;
  int nbPossibilityTotal;
  //
  int nbCompletelyAdjacent;
  std::vector<int> ListCompletelyAdjacent;
  int nbComplex;
  std::vector<int> ListComplex;
  int CurrPos;
  std::vector<int> PointStatus;
};

struct FullChain {
  int CurrLevel;
  std::vector<int> ListPossibility;
  std::vector<OneLevel> ListLevel;
};


void SetListPoss(GraphType const& eGraph, FullChain & eChain, int const& iLevel)
{
  std::cerr << "SetListPoss, iLevel=" << iLevel << "\n";
  std::cerr << "  len=" << eChain.ListLevel[iLevel].eVect.len << "\n";
  std::cerr << "  V=";
  for (int idx=0; idx<eChain.ListLevel[iLevel].eVect.len; idx++)
    std::cerr << eChain.ListLevel[iLevel].eVect.V[idx] << " ";
  std::cerr << "\n";
  auto IsCorrect=[&](int const& iPoint) -> bool {
    for (int idx=0; idx<eChain.ListLevel[iLevel].eVect.len; idx++) {
      int jPoint = eChain.ListLevel[iLevel].eVect.V[idx];
      if (eGraph.LLAdj[iPoint][jPoint] == 0)
        return false;
    }
    return true;
  };
  int CurrPos = eChain.ListLevel[iLevel-1].CurrPos;
  int nbComplexPrev = eChain.ListLevel[iLevel-1].nbComplex;
  std::cerr << "  CurrPos=" << CurrPos << " nbComplexPrev=" << nbComplexPrev << "\n";
  // Using the previous set of list possibility.
  int nbPossibility=0;
  for (int idx=CurrPos + 1; idx<nbComplexPrev; idx++) {
    int iPoint = eChain.ListLevel[iLevel-1].ListComplex[idx];
    std::cerr << "  idx=" << idx << " iPoint=" << iPoint << "\n";
    if (eChain.ListLevel[iLevel-1].PointStatus[idx] == 0) {
      if (IsCorrect(iPoint)) {
        eChain.ListPossibility[nbPossibility] = iPoint;
        std::cerr << "      nbPossibility=" << nbPossibility << " iPoint=" << iPoint << "\n";
        nbPossibility++;
      }
    }
  }
  std::cerr << "  nbPossibility=" << nbPossibility << "\n";
  int nbCompletelyAdjacent=0;
  int nbComplex=0;
  for (int iPoss=0; iPoss<nbPossibility; iPoss++) {
    auto IsCompletelyAdjacent=[&](int const& iPoss, int const& ePoint) -> bool {
      for (int jPoss=0; jPoss<nbPossibility; jPoss++)
        if (iPoss != jPoss) {
          int fPoint = eChain.ListPossibility[jPoss];
          if (eGraph.LLAdj[ePoint][fPoint] == 0)
            return false;
        }
      return true;
    };
    int eVert = eChain.ListPossibility[iPoss];
    if (IsCompletelyAdjacent(iPoss, eVert)) {
      eChain.ListLevel[iLevel].ListCompletelyAdjacent[nbCompletelyAdjacent] = eVert;
      nbCompletelyAdjacent++;
    }
    else {
      eChain.ListLevel[iLevel].ListComplex[nbComplex] = eVert;
      eChain.ListLevel[iLevel].PointStatus[nbComplex] = 0;
      nbComplex++;
    }
  }
  int nbComplement=0;
  if (nbPossibility == 0) { // no need to compute if there is already one regular extension
    int iPointStart = eChain.ListLevel[iLevel-1].ListComplex[CurrPos];
    auto IsNotMaximal=[&]() -> int {
      for (int iPoint=0; iPoint<iPointStart; iPoint++)
        if (IsCorrect(iPoint))
          return 1;
      return 0;
    };
    nbComplement = IsNotMaximal();
  }
  std::cerr << "  nbComplex=" << nbComplex << " nbCompletelyAdjacent=" << nbCompletelyAdjacent << "\n";
  eChain.ListLevel[iLevel].nbPossibilityTotal = nbPossibility + nbComplement;
  eChain.ListLevel[iLevel].nbCompletelyAdjacent = nbCompletelyAdjacent;
  eChain.ListLevel[iLevel].nbComplex = nbComplex;
  eChain.ListLevel[iLevel].CurrPos = 0;
}


// At any step, the ListPoss must be initialized at the end
// of the process.
FullChain GetTotalFullLevel(int const& nbPoint)
{
  std::vector<OneLevel> ListLevel(nbPoint+1);
  for (int iPoint=0; iPoint<=nbPoint; iPoint++) {
    int len=0;
    std::vector<int> eVect(nbPoint, -400);
    VectRec V{len, eVect};
    int nbPossibilityTotal = -400;
    //
    int nbCompletelyAdjacent = -400;
    std::vector<int> ListCompletelyAdjacent(nbPoint, -400);
    int nbComplex = -400;
    std::vector<int> ListComplex(nbPoint, -400);
    int CurrPos=-400;
    std::vector<int> PointStatus(nbPoint);
    OneLevel eLevel{V, nbPossibilityTotal, nbCompletelyAdjacent, ListCompletelyAdjacent, nbComplex, ListComplex, CurrPos, PointStatus};
    ListLevel[iPoint] = eLevel;
  }
  for (int iPoint=0; iPoint<nbPoint; iPoint++)
    ListLevel[0].ListComplex[iPoint] = iPoint;
  ListLevel[0].nbComplex = nbPoint;
  ListLevel[0].nbCompletelyAdjacent = 0; // otherwise your input is really bad
  ListLevel[0].CurrPos = 0;
  std::vector<int> ListPossibility(nbPoint);
  return {0, ListPossibility, ListLevel};
}

// The vector eVect starts at startPos which already contains previous eVect and the completely
// adjacent points.
// Idea is that points in the same adajcency pattern can be considered.
void AssignationVect(GraphType const& eGraph, FullChain & eChain, int const& iLevel, int const& startPos, int const& CurrPos)
{
  int nbComplex=eChain.ListLevel[iLevel].nbComplex;
  int eVert = eChain.ListLevel[iLevel].ListComplex[CurrPos];
  auto CheckPoint=[&](int const& posIns) -> bool {
    int fVert=eChain.ListLevel[iLevel].ListComplex[posIns];
    if (eGraph.LLAdj[eVert][fVert] == 0)
      return false;
    for (int pos=0; pos<nbComplex; pos++) {
      if (pos != CurrPos && pos != posIns) {
        int gVert=eChain.ListLevel[iLevel].ListComplex[pos];
        if (eGraph.LLAdj[eVert][gVert] != eGraph.LLAdj[fVert][gVert])
          return false;
      }
    }
    return true;
  };
  int posIns=startPos;
  auto FuncInsert=[&](int const& idx_v) -> void {
    eChain.ListLevel[iLevel].PointStatus[idx_v] = 1;
    eChain.ListLevel[iLevel+1].eVect.V[posIns] = eChain.ListLevel[iLevel].ListComplex[idx_v];
    std::cerr << "Assigning posIns=" << posIns << " val=" << eChain.ListLevel[iLevel].ListComplex[idx_v] << "\n";
    posIns++;
  };
  FuncInsert(CurrPos);
  for (int idx=CurrPos+1; idx<eChain.ListLevel[iLevel].nbComplex; idx++)
    if (CheckPoint(idx))
      FuncInsert(idx);
  std::cerr << "AssignationVect, iLevel=" << iLevel << " posIns=" << posIns << "\n";
  eChain.ListLevel[iLevel+1].eVect.len = posIns;
}

bool GoUpNextInTree(GroupType const& eGroup, GraphType const& eGraph, FullChain & eChain)
{
  int iLevel=eChain.CurrLevel;
  std::cerr << "GoUpNextInTree, iLevel=" << iLevel << "\n";
  if (iLevel > 0) {
    int startPos = eChain.ListLevel[iLevel-1].eVect.len + eChain.ListLevel[iLevel-1].nbCompletelyAdjacent;
    int CurrPos = eChain.ListLevel[iLevel-1].CurrPos;
    while(true) {
      if (CurrPos == eChain.ListLevel[iLevel-1].nbComplex)
        break;
      if (eChain.ListLevel[iLevel-1].PointStatus[CurrPos] == 0) {
        std::cerr << "  GoUpNextInTree iLevel=" << iLevel << " CurrPos=" << CurrPos << "\n";
        AssignationVect(eGraph, eChain, iLevel-1, startPos, CurrPos);
        std::cerr << "  GoUpNextInTree iLevel=" << iLevel << " val=" << eChain.ListLevel[iLevel-1].PointStatus[CurrPos] << "\n";
        bool test = IsMinimal(eGroup, eChain.ListLevel[iLevel].eVect);
        if (test) {
          eChain.ListLevel[iLevel-1].CurrPos = CurrPos;
          SetListPoss(eGraph, eChain, iLevel);
          return true;
        }
      }
      CurrPos++;
    }
  }
  if (iLevel == 0)
    return false;
  eChain.CurrLevel = iLevel-1;
  return GoUpNextInTree(eGroup, eGraph, eChain);
}



bool NextInTree(GroupType const& eGroup, GraphType const& eGraph, FullChain & eChain)
{
  int iLevel=eChain.CurrLevel;
  int lenPrev = eChain.ListLevel[iLevel].eVect.len;
  for (int idx=0; idx<lenPrev; idx++)
    eChain.ListLevel[iLevel+1].eVect.V[idx] = eChain.ListLevel[iLevel].eVect.V[idx];
  int nbCompletelyAdjacent = eChain.ListLevel[iLevel].nbCompletelyAdjacent;
  for (int idx=0; idx<nbCompletelyAdjacent; idx++)
    eChain.ListLevel[iLevel+1].eVect.V[idx + lenPrev] = eChain.ListLevel[iLevel].ListCompletelyAdjacent[idx];
  int startPos = lenPrev + nbCompletelyAdjacent;
  std::cerr << "NextInTree, iLevel=" << iLevel << " lenPrev=" << lenPrev << " nbCompletelyAdjacent=" << nbCompletelyAdjacent << "\n";
  int FirstPos=0;
  std::cerr << "nbComplex=" << eChain.ListLevel[iLevel].nbComplex << "\n";
  while(true) {
    if (FirstPos == eChain.ListLevel[iLevel].nbComplex)
      break;
    if (eChain.ListLevel[iLevel].PointStatus[FirstPos] == 0) {
      AssignationVect(eGraph, eChain, iLevel, startPos, FirstPos);
      PrintVector("NextInTree eVect=", eChain.ListLevel[iLevel+1].eVect);
      bool test = IsMinimal(eGroup, eChain.ListLevel[iLevel+1].eVect);
      std::cerr << "test=" << test << "\n";
      if (test) {
        eChain.ListLevel[iLevel].CurrPos = FirstPos;
        SetListPoss(eGraph, eChain, iLevel+1);
        eChain.CurrLevel = iLevel + 1;
        return true;
      }
    }
    FirstPos++;
  }
  return GoUpNextInTree(eGroup, eGraph, eChain);
}



void DoEnumeration(GroupType const& eGroup, GraphType const& eGraph, std::string const& MaximalFile)
{
  bool IsFirst=true;
  std::ofstream os(MaximalFile);
  FullChain eChain = GetTotalFullLevel(eGraph.nbPoint);
  int nbIter = 0;
  os << "return [\n";
  while(true) {
    PrintVector("Level Main Loop, eVect=", eChain.ListLevel[eChain.CurrLevel].eVect);
    std::cerr << "nbPossibilityTotal=" << eChain.ListLevel[eChain.CurrLevel].nbPossibilityTotal << "\n";
    if (eChain.ListLevel[eChain.CurrLevel].nbPossibilityTotal == 0) {
      if (!IsFirst)
        os << ",\n";
      IsFirst=false;
      std::vector<int> eVect(eChain.ListLevel[eChain.CurrLevel].eVect.len);
      for (int idx=0; idx<eChain.ListLevel[eChain.CurrLevel].eVect.len; idx++)
        eVect[idx] = eChain.ListLevel[eChain.CurrLevel].eVect.V[idx] + 1;
      std::sort(eVect.begin(), eVect.end());
      os << "[";
      for (int idx=0; idx<eChain.ListLevel[eChain.CurrLevel].eVect.len; idx++) {
        if (idx > 0)
          os << ",";
        os << eVect[idx];
      }
      os << "]";
    }
    bool test = NextInTree(eGroup, eGraph, eChain);
    nbIter++;
    //    std::cerr << "test=" << test << " CurrLevel=" << eChain.CurrLevel << "\n";
    if (!test)
      break;
  }
  os << "];\n";
  std::cerr << "nbIter=" << nbIter << "\n";
}


int main(int argc, char* argv[])
{
  if (argc != 4) {
    std::cerr << "Program is used as\n";
    std::cerr << "EnumerationCliques [GraphFile] [GroupFile] [MaximalFile]\n";
    exit(1);
  }
  std::string GraphFile = argv[1];
  std::string GroupFile = argv[2];
  std::string MaximalFile = argv[3];

  GraphType eGraph = ReadGraphFile(GraphFile);
  GroupType eGroup = ReadGroupFile(GroupFile);
  //
  DoEnumeration(eGroup, eGraph, MaximalFile);
}
