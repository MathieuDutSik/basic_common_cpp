#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

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


struct TypeVectClique {
  int lenClique;
  std::vector<int> eVect;
};

struct TypeVectBRW {
  int len;
  std::vector<int> eVect;
};

struct TypeListGroupElement {
  int nbAdmissibleElt;
  std::vector<int> ListGroupPossibility;
};

struct OneLevel {
  TypeVectClique eVectClique;
  TypeVectBRW eVectBRW;
  TypeListGroupElement RecAdmissibleGroupElement;
};

struct FullChain {
  int CurrLevel;
  std::vector<OneLevel> ListLevel;
};


bool IsMinimal(GroupType const& eGroup, ArrMinimal const& eArr, TypeVectBRW const& eVectBRW, TypeListGroupElement & ReturnListGroupElement, TypeListGroupElement const& InputListGroupElement)
{
  int lenBRW=eVectBRW.len;
  int nbPoint=eGroup.nbPoint;
  // StatusPermutation is about how the permutation stands for later.
  // status=0 if it is not useful for further computation
  // status=1 if it may be useful later on.
  // status=2 if it proves non-minimality.
  auto StatusPermutation=[&](int const& iElt) -> int {
    for (int idx=0; idx<lenBRW; idx++) {
      eArr.VectCompar[idx]=2;
    }
    for (int idx=0; idx<lenBRW; idx++) {
      int iImg=eGroup.ARR[iElt][idx];
      if (iImg < lenBRW)
        eArr.VectCompar[iImg] = eVectBRW.eVect[idx];
    }
    // red is less than blue
    // we search for the lexicographically minimal element
    for (int idx=0; idx<lenBRW; idx++) {
      if (eVectBRW[idx] == 1) { // node is red
        if (eArr.VectCompar[iImg] == 2) // white, permutation may be useful
          return 1;
        if (eArr.VectCompar[iImg] == 0) // blue, permutation cannot help
          return 0;
        // last case is neutral, going forward
      }
      if (eVectBRW[idx] == 0) { // node is blue
        if (eArr.VectCompar[iImg] == 2) // white, permutation may be useful
          return 1;
        if (eArr.VectCompar[iImg] == 1) // red, permutation allows to conclude non-minimality
          return 2;
      }
    }
    return 1;
  };
  int iPoss=0;
  for (int idx=0; idx<InputListGroupElement.nbAdmissibleElt; idx++) {
    int iElt = InputListGroupElement.ListGroupPossibility[idx];
    int status = StatusPermutation(iElt);
    if (status == 2)
      return false;
    if (status == 1) {
      ReturnListGroupElement.ListGroupPossibility[iPoss]=iElt;
      iPoss++;
    }
  }
  ReturnListGroupElement.nbAdmissibleElt = iPoss;
  return true;
}








// At any step, the ListPoss must be initialized at the end
// of the process.
FullChain GetTotalFullLevel(int const& nbPoint, int const& nbElt)
{
  std::vector<OneLevel> ListLevel(nbPoint+1);
  std::vector<int> ListGroupPossibility(nbElt);
  for (int iElt=0; iElt<nbElt; iElt++)
    ListGroupPossibility[iElt]=iElt;
  TypeListGroupElement RecAdmissibleGroupElement{nbElt, ListGroupPossibility};
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    int sizClique=0;
    std::vector<int> eVect1(iPoint0, 0);
    TypeVectClique eVectClique{sizClique, eVect1};
    std::vector<int> eVect2(iPoint, 0);
    TypeVectBRW eVectBRW{iPoint,eVect2};
    OneLevel eLevel{eVectClique, eVectBRW, RecAdmissibleGroupElement};
    ListLevel[iPoint] = eLevel;
  }
  return {0, ListLevel};
}

bool IsAdmissibleVertex(GraphType const& eGraph, TypeVectClique const& eVectClique, int const& eVert)
{
  for (int idx=0; idx<eVectClique.sizClique; idx++) {
    int iVert = eVectClique.eVect[idx];
    if (eGraph.ARR[eVert][iVert] == 0)
      return false;
  }
  return true;
}
 
bool GoUpNextInTree(GroupType const& eGroup, GraphType const& eGraph, ArrMinimal & eArr, FullChain & eChain)
{
  int iLevel=eChain.CurrLevel;
  if (iLevel > 0) {
    if (eChain.ListLevel[iLevel].eVectBRW.eVect[iLevel-1] == 0) {
      if (IsAdmissibleVertex(eGraph, iLevel-1)) {
        eChain.ListLevel[iLevel].eVectBRW.eVect[iLevel-1]=1;
        if (IsMinimal(eGroup, eArr, eChain.ListLevel[iLevel].eVectBRW, eChain.ListLevel[iLevel].RecAdmissibleGroupElement, eChain.ListLevel[iLevel-1].RecAdmissibleGroupElement)) {
          int sizClique=eChain.ListLevel[iLevel].eVectClique.sizClique;
          eChain.ListLevel[iLevel].eVectClique.sizClique++;
          eChain.ListLevel[iLevel].eVectClique[sizClique] = iLevel-1;
          return true;
        }
      }
  }
  if (iLevel == 0)
    return false;
  eChain.CurrLevel = iLevel-1;
  //  std::cerr << "Before calling GoUpNextInTree 2\n";
  return GoUpNextInTree(eGroup, eGraph, eChain);
}


bool NextInTree(GroupType const& eGroup, GraphType const& eGraph, FullChain & eChain)
{
  int iLevel=eChain.CurrLevel;
  for (int idx=0; idx<iLevel; idx++) {
    eChain.ListLevel[iLevel+1].eVect[idx] = eChain.ListLevel[iLevel].eVect[idx];
  }
  for (int iPoss=0; iPoss<eChain.ListLevel[iLevel].nbPossibility; iPoss++) {
    eChain.ListLevel[iLevel+1].eVect[iLevel] = eChain.ListLevel[iLevel].ListPoss[iPoss];
    bool test = IsMinimal(eGroup, eChain.ListLevel[iLevel+1].eVect);
    if (test) {
      eChain.ListLevel[iLevel].CurrPos = iPoss;
      SetListPoss(eGraph, eChain, iLevel+1);
      eChain.CurrLevel = iLevel + 1;
      return true;
    }
  }
  return GoUpNextInTree(eGroup, eGraph, eChain);
}

void PrintLastLevel(FullChain const& eChain)
{
  int iLevel=eChain.CurrLevel;
  std::cerr << "iLevel=" << iLevel << " eVect=[";
  for (int idx=0; idx<iLevel; idx++) {
    if (idx>0)
      std::cerr << ",";
    std::cerr << eChain.ListLevel[iLevel].eVect[idx];
  }
}


void DoEnumeration(GroupType const& eGroup, GraphType const& eGraph, std::string const& MaximalFile)
{
  bool IsFirst=true;
  std::ofstream os(MaximalFile);
  FullChain eChain = GetTotalFullLevel(eGraph.nbPoint);
  os << "return [\n";
  while(true) {
    if (eChain.ListLevel[eChain.CurrLevel].nbPossibilityTotal == 0) {
      if (!IsFirst)
        os << ",\n";
      IsFirst=false;
      os << "[";
      for (int iPoint=0; iPoint<eChain.CurrLevel; iPoint++) {
        if (iPoint>0)
          os << ",";
        int eVal = eChain.ListLevel[eChain.CurrLevel].eVect[iPoint] + 1;
        os << eVal;
      }
      os << "]";
    }
    bool test = NextInTree(eGroup, eGraph, eChain);
    std::cerr << "test=" << test << " CurrLevel=" << eChain.CurrLevel << "\n";
    if (!test)
      break;
  }
  os << "];\n";
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
