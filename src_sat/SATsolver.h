#ifndef SRC_SAT_SATSOLVER_H_
#define SRC_SAT_SATSOLVER_H_

#include "Boost_bitset.h"
#include "Temp_common.h"
#include "minisat/core/Solver.h"

struct SATformulation {
  int nbVar;
  std::vector<std::vector<int>> ListExpr;
};


bool IsFeasibleExpression(std::vector<int> eList, Face const &eFace) {
  bool IsMatching = false;
  for (auto &x : eList) {
    int pos;
    bool val;
    if (x > 0) {
      pos = x - 1;
      val = true;
    } else {
      pos = -1 - x;
      val = false;
    }
    if (eFace[pos] == 0 && !val)
      IsMatching = true;
    if (eFace[pos] == 1 && val)
      IsMatching = true;
  }
  return IsMatching;
}

void PrintFeasibilityInfo(SATformulation const &eExpr, Face const &eFace) {
  int nbCond = eExpr.ListExpr.size();
  bool IsFeasible = true;
  for (int iCond = 0; iCond < nbCond; iCond++) {
    std::vector<int> eList = eExpr.ListExpr[iCond];
    bool IsMatching = IsFeasibleExpression(eList, eFace);
    if (!IsMatching) {
      std::cerr << "Condition iCond=" << iCond << " violated. eList=";
      WriteStdVector(std::cerr, eList);
      IsFeasible = false;
    }
  }
  std::cerr << "IsFeasible=" << IsFeasible << "\n";
}

bool IsFeasible(SATformulation const &eExpr, Face const &eFace) {
  for (auto &eList : eExpr.ListExpr) {
    bool IsMatching = IsFeasibleExpression(eList, eFace);
    if (!IsMatching) {
      return false;
    }
  }
  return true;
}

void SATenumeration(SATformulation const &eExpr, int const &MAX_ITER,
                    std::function<int(Face const &)> const &f) {
  Minisat::Solver S;
  int verb = 2;
  S.verbosity = verb;
  //  Solver *solver;
  //  solver = &S;
  int nbVar = eExpr.nbVar;
  std::cerr << "nbVar=" << nbVar << "\n";
  std::vector<Minisat::Var> ListVar;
  for (int iVar = 0; iVar < nbVar; iVar++)
    ListVar.push_back(S.newVar());
  auto GetLit = [&](int const &x) -> Minisat::Lit {
    int pos;
    bool val;
    if (x > 0) {
      pos = x - 1;
      val = true;
    } else {
      pos = -1 - x;
      val = false;
    }
    //    std::cerr << "x=" << x << " pos=" << pos << "\n";
    return Minisat::mkLit(ListVar[pos], val);
  };
  std::cerr << "nbClause=" << eExpr.ListExpr.size() << "\n";
  int TotalSiz = 0;
  for (auto &eList : eExpr.ListExpr) {
    int siz = eList.size();
    TotalSiz += siz;
    //    std::cerr << "siz=" << siz << "\n";
    if (siz == 2) {
      Minisat::Lit lit1 = GetLit(eList[0]);
      Minisat::Lit lit2 = GetLit(eList[1]);
      //      std::cerr << "Before addClause 1\n";
      S.addClause(lit1, lit2);
      //      std::cerr << "After addClause 1\n";
    } else {
      Minisat::vec<Minisat::Lit> ListLit(siz);
      for (int i = 0; i < siz; i++)
        ListLit[i] = GetLit(eList[i]);
      //      std::cerr << "Before addClause 2\n";
      S.addClause_(ListLit);
      //      std::cerr << "After addClause 2\n";
    }
  }
  std::cerr << "TotalSiz=" << TotalSiz << "\n";
  int iter = 0;
  while (S.solve()) {
    Face eFace(nbVar);
    for (int i = 0; i < nbVar; i++)
      eFace[i] = 0;
    for (int i = 0; i < nbVar; i++) {
      using namespace Minisat;
      if (S.modelValue(i) ==
          l_False) // Curiously, this is how we should do it, it seems.
        eFace[i] = 1;
    }
    int terminate = f(eFace);
    if (terminate == 1)
      break;
    Minisat::vec<Minisat::Lit> blocking_clause(nbVar);
    for (int i = 0; i < nbVar; i++) {
      if (eFace[i] == 1)
        blocking_clause[i] = Minisat::mkLit(ListVar[i], false);
      else
        blocking_clause[i] = Minisat::mkLit(ListVar[i], true);
    }
    S.addClause(blocking_clause);
    iter++;
    if (MAX_ITER == -1)
      if (iter > MAX_ITER)
        break;
  }
  std::cerr << "iter=" << iter << "\n";
}

#endif
