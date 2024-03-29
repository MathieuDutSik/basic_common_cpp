// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

struct GraphType {
  int nbPoint;
  std::vector<int> LLAdj;
};

GraphType ReadGraphFile(std::string const &eFile) {
  std::ifstream is(eFile);
  int nbPoint;
  is >> nbPoint;
  std::vector<int> LLAdj;
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    for (int jPoint = 0; jPoint < nbPoint; jPoint++) {
      int eVal;
      is >> eVal;
      LLAdj.push_back(eVal);
    }
  }
  return {nbPoint, LLAdj};
}

struct GroupType {
  int nbElt;
  int nbPoint;
  std::vector<int> ARR;
};

GroupType ReadGroupFile(std::string const &eFile) {
  std::ifstream is(eFile);
  int nbPoint, nbElt;
  is >> nbPoint;
  is >> nbElt;
  std::vector<int> ARR;
  for (int iElt = 0; iElt < nbElt; iElt++) {
    for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
      int eVal;
      is >> eVal;
      ARR.push_back(eVal);
    }
  }
  return {nbElt, nbPoint, ARR};
}

bool IsMinimal(GroupType const &eGroup, std::vector<int> const &eVect) {
  int len = eVect.size();
  int nbPoint = eGroup.nbPoint;
  std::vector<int> ImageVect(len);
  auto IsCounterexample = [&](std::vector<int> const &uVect) -> bool {
    for (int i = 0; i < len; i++) {
      if (uVect[i] < eVect[i])
        return true;
      // for the next item to be counterexample, it has to be equal
      if (uVect[i] > eVect[i])
        return false;
    }
    return false;
  };
  for (int iElt = 0; iElt < eGroup.nbElt; iElt++) {
    for (int i = 0; i < len; i++) {
      int iImg = eGroup.ARR[eVect[i] + iElt * nbPoint];
      ImageVect[i] = iImg;
    }
    std::sort(ImageVect.begin(), ImageVect.end());
    if (IsCounterexample(ImageVect))
      return false;
  }
  return true;
}

struct OneLevel {
  std::vector<int> eVect;
  int nbPossibility;
  int nbPossibilityTotal;
  int CurrPos;
  std::vector<int> ListPoss;
};

struct FullChain {
  int CurrLevel;
  std::vector<OneLevel> ListLevel;
};

void SetListPoss(GraphType const &eGraph, FullChain &eChain,
                 int const &iLevel) {
  int nbPoint = eGraph.nbPoint;
  /*
  std::cerr << "iLevel=" << iLevel << " eVect=[";
  for (int idx=0; idx<iLevel; idx++) {
    if (idx>0)
      std::cerr << ",";
    int jPoint = eChain.ListLevel[iLevel].eVect[idx];
    std::cerr << jPoint;
  }
  std::cerr << "];\n"; */
  auto IsCorrect = [&](int const &iPoint) -> bool {
    for (int idx = 0; idx < iLevel; idx++) {
      int jPoint = eChain.ListLevel[iLevel].eVect[idx];
      //      std::cerr << "IsCorrect idx=" << idx << " jPoint=" << jPoint << "
      //      val=" << eGraph.LLAdj[iPoint][jPoint] << "\n";
      if (eGraph.LLAdj[jPoint + iPoint * nbPoint] == 0)
        return false;
    }
    return true;
  };
  int nbPossibility = 0;
  int nbComplement = 0;
  int iPointStart = 0;
  if (iLevel > 0) {
    iPointStart = eChain.ListLevel[iLevel].eVect[iLevel - 1] + 1;
  }
  for (int iPoint = iPointStart; iPoint < nbPoint; iPoint++) {
    if (IsCorrect(iPoint)) {
      //      std::cerr << "Inserting iPoint=" << iPoint << "\n";
      eChain.ListLevel[iLevel].ListPoss[nbPossibility] = iPoint;
      nbPossibility++;
    }
  }
  if (nbPossibility == 0) {
    // no need to compute if there is already one reular extension
    for (int iPoint = 0; iPoint < iPointStart; iPoint++)
      if (IsCorrect(iPoint))
        nbComplement++;
  }
  eChain.ListLevel[iLevel].nbPossibility = nbPossibility;
  eChain.ListLevel[iLevel].nbPossibilityTotal = nbPossibility + nbComplement;
  eChain.ListLevel[iLevel].CurrPos = 0;
}

// At any step, the ListPoss must be initialized at the end
// of the process.
FullChain GetTotalFullLevel(int const &nbPoint) {
  std::vector<OneLevel> ListLevel(nbPoint + 1);
  for (int iPoint = 0; iPoint <= nbPoint; iPoint++) {
    std::vector<int> eVect(iPoint, -1);
    int nbPoss = -1;
    int CurrPos = -1;
    std::vector<int> ListPoss(nbPoint);
    OneLevel eLevel{eVect, nbPoss, nbPoss, CurrPos, ListPoss};
    ListLevel[iPoint] = eLevel;
  }
  for (int iPoint = 0; iPoint < nbPoint; iPoint++)
    ListLevel[0].ListPoss[iPoint] = iPoint;
  ListLevel[0].nbPossibility = nbPoint;
  ListLevel[0].CurrPos = 0;
  return {0, ListLevel};
}

bool GoUpNextInTree(GroupType const &eGroup, GraphType const &eGraph,
                    FullChain &eChain) {
  int iLevel = eChain.CurrLevel;
  if (iLevel > 0) {
    int CurrPos = eChain.ListLevel[iLevel - 1].CurrPos;
    int nbPossibility = eChain.ListLevel[iLevel - 1].nbPossibility;
    /*
    std::cerr << "CP: iLevel-1=" << (iLevel-1) << " CurrPos=" << CurrPos << "
    nbPossibility=" << nbPossibility << "\n"; std::cerr << "CP: eVect=["; for
    (int idx=0; idx<iLevel-1; idx++) { if (idx > 0) std::cerr << ","; std::cerr
    << eChain.ListLevel[iLevel-1].eVect[idx];
    }
    std::cerr << "]\n"; */
    for (int iPoss = CurrPos + 1; iPoss < nbPossibility; iPoss++) {
      eChain.ListLevel[iLevel].eVect[iLevel - 1] =
          eChain.ListLevel[iLevel - 1].ListPoss[iPoss];
      bool test = IsMinimal(eGroup, eChain.ListLevel[iLevel].eVect);
      //      std::cerr << "CurrPos=" << CurrPos << " iPoss=" << iPoss << "
      //      test=" << test << " ePoint=" <<
      //      eChain.ListLevel[iLevel-1].ListPoss[iPoss] << "\n";
      if (test) {
        //        std::cerr << "Assigning iLevel-1=" << (iLevel-1) << "
        //        CurrPos=" << iPoss << "\n";
        eChain.ListLevel[iLevel - 1].CurrPos = iPoss;
        SetListPoss(eGraph, eChain, iLevel);
        return true;
      }
    }
  }
  if (iLevel == 0)
    return false;
  eChain.CurrLevel = iLevel - 1;
  //  std::cerr << "Before calling GoUpNextInTree 2\n";
  return GoUpNextInTree(eGroup, eGraph, eChain);
}

bool NextInTree(GroupType const &eGroup, GraphType const &eGraph,
                FullChain &eChain) {
  int iLevel = eChain.CurrLevel;
  for (int idx = 0; idx < iLevel; idx++)
    eChain.ListLevel[iLevel + 1].eVect[idx] =
        eChain.ListLevel[iLevel].eVect[idx];
  for (int iPoss = 0; iPoss < eChain.ListLevel[iLevel].nbPossibility; iPoss++) {
    eChain.ListLevel[iLevel + 1].eVect[iLevel] =
        eChain.ListLevel[iLevel].ListPoss[iPoss];
    bool test = IsMinimal(eGroup, eChain.ListLevel[iLevel + 1].eVect);
    if (test) {
      eChain.ListLevel[iLevel].CurrPos = iPoss;
      SetListPoss(eGraph, eChain, iLevel + 1);
      eChain.CurrLevel = iLevel + 1;
      return true;
    }
  }
  //  std::cerr << "Before calling GoUpNextInTree 1\n";
  return GoUpNextInTree(eGroup, eGraph, eChain);
}

void PrintLastLevel(FullChain const &eChain) {
  int iLevel = eChain.CurrLevel;
  std::cerr << "iLevel=" << iLevel << " eVect=[";
  for (int idx = 0; idx < iLevel; idx++) {
    if (idx > 0)
      std::cerr << ",";
    std::cerr << eChain.ListLevel[iLevel].eVect[idx];
  }
  std::cerr << "] ListPoss=[";
  for (int iPoss = 0; iPoss < eChain.ListLevel[iLevel].nbPossibility; iPoss++) {
    if (iPoss > 0)
      std::cerr << ",";
    std::cerr << eChain.ListLevel[iLevel].ListPoss[iPoss];
  }
  std::cerr << "]\n";
}

void DoEnumeration(GroupType const &eGroup, GraphType const &eGraph,
                   std::string const &MaximalFile) {
  bool IsFirst = true;
  std::ofstream os(MaximalFile);
  FullChain eChain = GetTotalFullLevel(eGraph.nbPoint);
  int nbIter = 0;
  os << "return [\n";
  while (true) {
    //    PrintLastLevel(eChain);
    if (eChain.ListLevel[eChain.CurrLevel].nbPossibilityTotal == 0) {
      if (!IsFirst)
        os << ",\n";
      IsFirst = false;
      os << "[";
      for (int iPoint = 0; iPoint < eChain.CurrLevel; iPoint++) {
        if (iPoint > 0)
          os << ",";
        int eVal = eChain.ListLevel[eChain.CurrLevel].eVect[iPoint] + 1;
        os << eVal;
      }
      os << "]";
    }
    bool test = NextInTree(eGroup, eGraph, eChain);
    nbIter++;
    //    std::cerr << "test=" << test << " CurrLevel=" << eChain.CurrLevel <<
    //    "\n";
    if (!test)
      break;
  }
  os << "];\n";
  std::cerr << "nbIter=" << nbIter << "\n";
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
