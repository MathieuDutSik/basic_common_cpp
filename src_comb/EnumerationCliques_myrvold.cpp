// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

// This program uses the method of the manuscript of Wendy Myrvold and
// Patrick Fowler, "Fast Enumeration of all independent sets of a graph"
// It is more complicated but the group is considered in a more intelligent
// way.
//
// Idea of the iteration is that we keep the elements in a Blue, Red, White
// way. So, we go step by step making binary choice about the Blue, Red vertex
// nature.
// The group elements are kept if they can potentially be helpful in deciding
// the minimality of the generated cliques.
//
// It is a priori faster but it turned out to be slower than the ordered
// enumeration.

struct GraphType {
  int nbPoint;
  std::vector<std::vector<int>> LLAdj;
};

GraphType ReadGraphFile(std::string const &eFile) {
  std::ifstream is(eFile);
  int nbPoint;
  is >> nbPoint;
  std::vector<std::vector<int>> LLAdj;
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    std::vector<int> LAdj;
    for (int jPoint = 0; jPoint < nbPoint; jPoint++) {
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

GroupType ReadGroupFile(std::string const &eFile) {
  std::ifstream is(eFile);
  int nbPoint, nbElt;
  is >> nbPoint;
  is >> nbElt;
  std::vector<std::vector<int>> ARR;
  for (int iElt = 0; iElt < nbElt; iElt++) {
    std::vector<int> LAdj;
    for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
      int eVal;
      is >> eVal;
      LAdj.push_back(eVal);
    }
    ARR.push_back(LAdj);
  }
  return {nbElt, nbPoint, ARR};
}

struct TypeVectClique {
  int sizClique;
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

struct ArrMinimal {
  std::vector<int> VectCompar;
};

ArrMinimal GetRecordArrMinimal(int const &nbPoint) {
  std::vector<int> eVect(nbPoint, -1);
  return {eVect};
}

bool IsMinimal(GroupType const &eGroup, ArrMinimal &eArr,
               TypeVectBRW const &eVectBRW,
               TypeListGroupElement &ReturnListGroupElement,
               TypeListGroupElement const &InputListGroupElement) {
  //  std::cerr << "IsMinimal beginning\n";
  int lenBRW = eVectBRW.len;
  // StatusPermutation is about how the permutation stands for later.
  // status=0 if it is not useful for further computation
  // status=1 if it may be useful later on.
  // status=2 if it proves non-minimality.
  auto StatusPermutation = [&](int const &iElt) -> int {
    for (int idx = 0; idx < lenBRW; idx++) {
      eArr.VectCompar[idx] = 2;
    }
    for (int idx = 0; idx < lenBRW; idx++) {
      int iImg = eGroup.ARR[iElt][idx];
      if (iImg < lenBRW)
        eArr.VectCompar[iImg] = eVectBRW.eVect[idx];
    }
    // red is less than blue
    // we search for the lexicographically minimal element
    for (int idx = 0; idx < lenBRW; idx++) {
      if (eVectBRW.eVect[idx] == 1) {  // node is red
        if (eArr.VectCompar[idx] == 2) // white, permutation may be useful
          return 1;
        if (eArr.VectCompar[idx] == 0) // blue, permutation cannot help
          return 0;
        // last case is neutral, going forward
      }
      if (eVectBRW.eVect[idx] == 0) {  // node is blue
        if (eArr.VectCompar[idx] == 2) // white, permutation may be useful
          return 1;
        if (eArr.VectCompar[idx] ==
            1) // red, permutation allows to conclude non-minimality
          return 2;
      }
    }
    return 1;
  };
  int iPoss = 0;
  for (int idx = 0; idx < InputListGroupElement.nbAdmissibleElt; idx++) {
    int iElt = InputListGroupElement.ListGroupPossibility[idx];
    int status = StatusPermutation(iElt);
    if (status == 2) {
      //      std::cerr << "IsMinimal, Returning with false\n";
      return false;
    }
    if (status == 1) {
      ReturnListGroupElement.ListGroupPossibility[iPoss] = iElt;
      iPoss++;
    }
  }
  //  std::cerr << "IsMinimal, Exiting iPoss=" << iPoss << "\n";
  ReturnListGroupElement.nbAdmissibleElt = iPoss;
  return true;
}

// At any step, the ListPoss must be initialized at the end
// of the process.
FullChain GetTotalFullLevel(int const &nbPoint, int const &nbElt) {
  std::vector<OneLevel> ListLevel(nbPoint + 1);
  std::vector<int> ListGroupPossibility(nbElt);
  for (int iElt = 0; iElt < nbElt; iElt++)
    ListGroupPossibility[iElt] = iElt;
  TypeListGroupElement RecAdmissibleGroupElement{nbElt, ListGroupPossibility};
  for (int iPoint = 0; iPoint <= nbPoint; iPoint++) {
    int sizClique = 0;
    std::vector<int> eVect1(iPoint, 0);
    TypeVectClique eVectClique{sizClique, eVect1};
    std::vector<int> eVect2(iPoint, 0);
    TypeVectBRW eVectBRW{iPoint, eVect2};
    OneLevel eLevel{eVectClique, eVectBRW, RecAdmissibleGroupElement};
    ListLevel[iPoint] = eLevel;
  }
  return {0, ListLevel};
}

bool IsAdmissibleVertex(GraphType const &eGraph,
                        TypeVectClique const &eVectClique, int const &eVert) {
  for (int idx = 0; idx < eVectClique.sizClique; idx++) {
    int iVert = eVectClique.eVect[idx];
    if (eGraph.LLAdj[eVert][iVert] == 0)
      return false;
  }
  return true;
}

bool GoUpNextInTree(GroupType const &eGroup, GraphType const &eGraph,
                    ArrMinimal &eArr, FullChain &eChain) {
  int iLevel = eChain.CurrLevel;
  //  std::cerr << "Start of GoUpNextInTree iLevel=" << iLevel << "\n";
  if (iLevel > 0) {
    if (eChain.ListLevel[iLevel].eVectBRW.eVect[iLevel - 1] == 0) {
      if (IsAdmissibleVertex(eGraph, eChain.ListLevel[iLevel].eVectClique,
                             iLevel - 1)) {
        eChain.ListLevel[iLevel].eVectBRW.eVect[iLevel - 1] = 1;
        if (IsMinimal(eGroup, eArr, eChain.ListLevel[iLevel].eVectBRW,
                      eChain.ListLevel[iLevel].RecAdmissibleGroupElement,
                      eChain.ListLevel[iLevel - 1].RecAdmissibleGroupElement)) {
          int sizClique = eChain.ListLevel[iLevel].eVectClique.sizClique;
          //          std::cerr << "GoUpNextInTree sizClique=" << sizClique <<
          //          "\n";
          eChain.ListLevel[iLevel].eVectClique.sizClique++;
          //          std::cerr << "GoUpNextInTree ASSIGN
          //          eChain.ListLevel[iLevel].eVectClique.sizClique=" <<
          //          eChain.ListLevel[iLevel].eVectClique.sizClique << "\n";
          eChain.ListLevel[iLevel].eVectClique.eVect[sizClique] = iLevel - 1;
          return true;
        }
      }
    }
  }
  if (iLevel == 0)
    return false;
  eChain.CurrLevel = iLevel - 1;
  //  std::cerr << "Before calling GoUpNextInTree 2\n";
  return GoUpNextInTree(eGroup, eGraph, eArr, eChain);
}

bool NextInTree(GroupType const &eGroup, GraphType const &eGraph,
                ArrMinimal &eArr, FullChain &eChain) {
  int iLevel = eChain.CurrLevel;
  //  std::cerr << "NextInTree iLevel=" << iLevel << " nbPoint=" <<
  //  eGraph.nbPoint << "\n";
  int nbPoint = eGraph.nbPoint;
  if (iLevel == nbPoint)
    return GoUpNextInTree(eGroup, eGraph, eArr, eChain);
  int sizClique = eChain.ListLevel[iLevel].eVectClique.sizClique;
  eChain.ListLevel[iLevel + 1].eVectClique.sizClique = sizClique;
  //  std::cerr << "NextInTree After ASSIGN 1
  //  eChain.ListLevel[iLevel+1].eVectClique.sizClique=" <<
  //  eChain.ListLevel[iLevel+1].eVectClique.sizClique << " sizClique=" <<
  //  sizClique << "\n";
  for (int idx = 0; idx < sizClique; idx++) {
    eChain.ListLevel[iLevel + 1].eVectClique.eVect[idx] =
        eChain.ListLevel[iLevel].eVectClique.eVect[idx];
  }
  for (int idx = 0; idx < iLevel; idx++) {
    //    std::cerr << "BEFORE idx=" << idx << "\n";
    eChain.ListLevel[iLevel + 1].eVectBRW.eVect[idx] =
        eChain.ListLevel[iLevel].eVectBRW.eVect[idx];
    //    std::cerr << " AFTER idx=" << idx << "\n";
  }
  for (int iColor = 0; iColor < 2; iColor++) {
    if (iColor == 0 ||
        IsAdmissibleVertex(eGraph, eChain.ListLevel[iLevel].eVectClique,
                           iLevel)) {
      eChain.ListLevel[iLevel + 1].eVectBRW.eVect[iLevel] = iColor;
      if (IsMinimal(eGroup, eArr, eChain.ListLevel[iLevel + 1].eVectBRW,
                    eChain.ListLevel[iLevel + 1].RecAdmissibleGroupElement,
                    eChain.ListLevel[iLevel].RecAdmissibleGroupElement)) {
        eChain.CurrLevel = iLevel + 1;
        //        std::cerr << "iColor=" << iColor << " sizClique=" << sizClique
        //        << "\n";
        if (iColor == 1) {
          eChain.ListLevel[iLevel + 1].eVectClique.sizClique++;
          eChain.ListLevel[iLevel + 1].eVectClique.eVect[sizClique] = iLevel;
        }
        //        std::cerr << "NextInTree After ASSIGN 2
        //        eChain.ListLevel[iLevel+1].eVectClique.sizClique=" <<
        //        eChain.ListLevel[iLevel+1].eVectClique.sizClique << "\n";
        return true;
      }
    }
  }
  return GoUpNextInTree(eGroup, eGraph, eArr, eChain);
}

void PrintLastLevel(FullChain const &eChain) {
  int iLevel = eChain.CurrLevel;
  std::cerr << "iLevel=" << iLevel << " eVect=[";
  for (int idx = 0; idx < iLevel; idx++) {
    if (idx > 0)
      std::cerr << ",";
    std::cerr << eChain.ListLevel[iLevel].eVectClique.eVect[idx];
  }
}

bool IsMaximalClique(GraphType const &eGraph, FullChain const &eChain) {
  int nbPoint = eGraph.nbPoint;
  int sizClique = eChain.ListLevel[nbPoint].eVectClique.sizClique;
  auto IsVertexOK = [&](int const &eVert) -> bool {
    for (int idx = 0; idx < sizClique; idx++) {
      int fVert = eChain.ListLevel[nbPoint].eVectClique.eVect[idx];
      if (eGraph.LLAdj[fVert][eVert] == 0)
        return false;
    }
    return true;
  };
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    if (eChain.ListLevel[nbPoint].eVectBRW.eVect[iPoint] ==
        0) { // a blue vertex so candidate for extension
      if (IsVertexOK(iPoint))
        return false;
    }
  }
  return true;
}

void DoEnumeration(GroupType const &eGroup, GraphType const &eGraph,
                   std::string const &MaximalFile) {
  bool IsFirst = true;
  std::ofstream os(MaximalFile);
  int nbPoint = eGraph.nbPoint;
  ArrMinimal eArr = GetRecordArrMinimal(nbPoint);
  FullChain eChain = GetTotalFullLevel(nbPoint, eGroup.nbElt);
  int nbIter = 0;
  int nbIterFinal = 0;
  os << "return [\n";
  while (true) {
    if (eChain.CurrLevel == nbPoint) {
      nbIterFinal++;
      if (IsMaximalClique(eGraph, eChain)) {
        std::cerr << "|STAB|="
                  << eChain.ListLevel[nbPoint]
                         .RecAdmissibleGroupElement.nbAdmissibleElt
                  << "\n";
        if (!IsFirst)
          os << ",\n";
        IsFirst = false;
        int sizClique = eChain.ListLevel[nbPoint].eVectClique.sizClique;
        os << "[";
        for (int idx = 0; idx < sizClique; idx++) {
          int eVert = eChain.ListLevel[nbPoint].eVectClique.eVect[idx] + 1;
          if (idx > 0)
            os << ",";
          os << eVert;
        }
        os << "]";
      }
    }
    nbIter++;
    bool test = NextInTree(eGroup, eGraph, eArr, eChain);
    if (!test)
      break;
  }
  os << "];\n";
  std::cerr << "nbIter=" << nbIter << " nbIterFinal=" << nbIterFinal << "\n";
}

int main(int argc, char *argv[]) {
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
