#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

struct GraphType {
  int nbPoint;
  std::vector<std::vector<int>> LLAdj;
};

struct GroupType {
  int nbElt;
  int nbPoint;
  std::vector<std::vector<int>> ARR;
};

struct ArrMinimal {
  std::vector<int> VectCompar;
};

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


bool IsMinimal(GroupType const& eGroup, ArrMinimal & eArr, std::vector<int> const& eVect)
{
  int len=eVect.size();
  auto IsCounterexample=[&]() -> bool {
    for (int i=0; i<len; i++) {
      if (eArr.VectCompar[i] < eVect[i])
        return true;
      if (eArr.VectCompar[i] > eVect[i]) // for the next item to be counterexample, it has to be equal
        return false;
    }
    return false;
  };
  for (int iElt=0; iElt<eGroup.nbElt; iElt++) {
    for (int i=0; i<len; i++) {
      int iImg=eGroup.ARR[iElt][eVect[i]];
      eArr.VectCompar[i]=iImg;
    }
    auto iter = eArr.VectCompar.begin();
    std::sort(iter, iter+len);
    if (IsCounterexample())
      return false;
  }
  return true;
}

ArrMinimal GetRecordArrMinimal(int const& nbPoint)
{
  std::vector<int> V(nbPoint,0);
  return {V};
}



void SetListPoss(GraphType const& eGraph, FullChain & eChain, int const& iLevel)
{
  int nbPoint=eGraph.nbPoint;
  /*
  std::cerr << "iLevel=" << iLevel << " eVect=[";
  for (int idx=0; idx<iLevel; idx++) {
    if (idx>0)
      std::cerr << ",";
    int jPoint = eChain.ListLevel[iLevel].eVect[idx];
    std::cerr << jPoint;
  }
  std::cerr << "];\n"; */
  auto IsCorrect=[&](int const& iPoint) -> bool {
    for (int idx=0; idx<iLevel; idx++) {
      int jPoint = eChain.ListLevel[iLevel].eVect[idx];
      //      std::cerr << "IsCorrect idx=" << idx << " jPoint=" << jPoint << " val=" << eGraph.LLAdj[iPoint][jPoint] << "\n";
      if (eGraph.LLAdj[iPoint][jPoint] == 0)
        return false;
    }
    return true;
  };
  int nbPossibility=0;
  int nbComplement=0;
  int iPointStart=0;
  if (iLevel > 0) {
    iPointStart = eChain.ListLevel[iLevel].eVect[iLevel-1] + 1;
  }
  for (int iPoint=iPointStart; iPoint<nbPoint; iPoint++) {
    if (IsCorrect(iPoint)) {
      //      std::cerr << "Inserting iPoint=" << iPoint << "\n";
      eChain.ListLevel[iLevel].ListPoss[nbPossibility] = iPoint;
      nbPossibility++;
    }
  }
  if (nbPossibility == 0) { // no need to compute if there is already one reular extension
    auto IsNotMaximal=[&]() -> int {
      for (int iPoint=0; iPoint<iPointStart; iPoint++)
        if (IsCorrect(iPoint))
          return 1;
      return 0;
    };
    nbComplement = IsNotMaximal();
  }
  eChain.ListLevel[iLevel].nbPossibility = nbPossibility;
  eChain.ListLevel[iLevel].nbPossibilityTotal = nbPossibility + nbComplement;
  eChain.ListLevel[iLevel].CurrPos = 0;
}


// At any step, the ListPoss must be initialized at the end
// of the process.
FullChain GetTotalFullLevel(int const& nbPoint)
{
  std::vector<OneLevel> ListLevel(nbPoint+1);
  for (int iPoint=0; iPoint<=nbPoint; iPoint++) {
    std::vector<int> eVect(iPoint, -1);
    int nbPoss=-1;
    int CurrPos=-1;
    std::vector<int> ListPoss(nbPoint);
    OneLevel eLevel{eVect, nbPoss, nbPoss, CurrPos, ListPoss};
    ListLevel[iPoint] = eLevel;
  }
  for (int iPoint=0; iPoint<nbPoint; iPoint++)
    ListLevel[0].ListPoss[iPoint] = iPoint;
  ListLevel[0].nbPossibility = nbPoint;
  ListLevel[0].CurrPos = 0;
  return {0, ListLevel};
}


bool GoUpNextInTree(GroupType const& eGroup, GraphType const& eGraph, ArrMinimal & eArr, FullChain & eChain)
{
  int iLevel=eChain.CurrLevel;
  if (iLevel > 0) {
    int CurrPos=eChain.ListLevel[iLevel-1].CurrPos;
    int nbPossibility=eChain.ListLevel[iLevel-1].nbPossibility;
    /*
    std::cerr << "CP: iLevel-1=" << (iLevel-1) << " CurrPos=" << CurrPos << " nbPossibility=" << nbPossibility << "\n";
    std::cerr << "CP: eVect=[";
    for (int idx=0; idx<iLevel-1; idx++) {
      if (idx > 0)
        std::cerr << ",";
      std::cerr << eChain.ListLevel[iLevel-1].eVect[idx];
    }
    std::cerr << "]\n"; */
    for (int iPoss=CurrPos+1; iPoss<nbPossibility; iPoss++) {
      eChain.ListLevel[iLevel].eVect[iLevel-1] = eChain.ListLevel[iLevel-1].ListPoss[iPoss];
      bool test = IsMinimal(eGroup, eArr, eChain.ListLevel[iLevel].eVect);
      //      std::cerr << "CurrPos=" << CurrPos << " iPoss=" << iPoss << " test=" << test << " ePoint=" << eChain.ListLevel[iLevel-1].ListPoss[iPoss] << "\n";
      if (test) {
        //        std::cerr << "Assigning iLevel-1=" << (iLevel-1) << " CurrPos=" << iPoss << "\n";
        eChain.ListLevel[iLevel-1].CurrPos = iPoss;
        SetListPoss(eGraph, eChain, iLevel);
        return true;
      }
    }
  }
  if (iLevel == 0)
    return false;
  eChain.CurrLevel = iLevel-1;
  //  std::cerr << "Before calling GoUpNextInTree 2\n";
  return GoUpNextInTree(eGroup, eGraph, eArr, eChain);
}


bool NextInTree(GroupType const& eGroup, GraphType const& eGraph, ArrMinimal & eArr, FullChain & eChain)
{
  int iLevel=eChain.CurrLevel;
  for (int idx=0; idx<iLevel; idx++)
    eChain.ListLevel[iLevel+1].eVect[idx] = eChain.ListLevel[iLevel].eVect[idx];
  for (int iPoss=0; iPoss<eChain.ListLevel[iLevel].nbPossibility; iPoss++) {
    eChain.ListLevel[iLevel+1].eVect[iLevel] = eChain.ListLevel[iLevel].ListPoss[iPoss];
    bool test = IsMinimal(eGroup, eArr, eChain.ListLevel[iLevel+1].eVect);
    if (test) {
      eChain.ListLevel[iLevel].CurrPos = iPoss;
      SetListPoss(eGraph, eChain, iLevel+1);
      eChain.CurrLevel = iLevel + 1;
      return true;
    }
  }
  //  std::cerr << "Before calling GoUpNextInTree 1\n";
  return GoUpNextInTree(eGroup, eGraph, eArr, eChain);
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
  std::cerr << "] ListPoss=[";
  for (int iPoss=0; iPoss<eChain.ListLevel[iLevel].nbPossibility; iPoss++) {
    if (iPoss > 0)
      std::cerr << ",";
    std::cerr << eChain.ListLevel[iLevel].ListPoss[iPoss];
  }
  std::cerr << "]\n";
}


void DoEnumeration(GroupType const& eGroup, GraphType const& eGraph, std::string const& MaximalFile)
{
  bool IsFirst=true;
  std::ofstream os(MaximalFile);
  int nbPoint = eGraph.nbPoint;
  FullChain eChain = GetTotalFullLevel(nbPoint);
  ArrMinimal eArr = GetRecordArrMinimal(nbPoint);
  int nbIter = 0;
  os << "return [\n";
  while(true) {
    //    PrintLastLevel(eChain);
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
    bool test = NextInTree(eGroup, eGraph, eArr, eChain);
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
