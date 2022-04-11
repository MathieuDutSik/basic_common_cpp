#include "SATsolver.h"
#include "Temp_common.h"

template <typename Tgr>
std::pair<bool, std::vector<int>> GetColoringOrFail(Tgr eGR,
                                                    int const &nbColor) {
  Minisat::Solver S;
  int nbVert = eGR.GetNbVert();
  std::vector<Minisat::Var> ListVar;
  int nbVar = nbColor * nbVert;
  for (int iVar = 0; iVar < nbVar; iVar++)
    ListVar.push_back(S.newVar());
  //
  // At least one color should be selected for each vertex
  //
  for (int iVert = 0; iVert < nbVert; iVert++) {
    Minisat::vec<Minisat::Lit> one_color_select(nbColor);
    for (int iColor = 0; iColor < nbColor; iColor++) {
      int pos = nbColor * iVert + iColor;
      Minisat::Lit cond = Minisat::mkLit(ListVar[pos], true);
      one_color_select[iColor] = cond;
    }
    S.addClause(one_color_select);
  }
  //
  // For each vertex at most one color should be selected
  //
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (int iC1 = 0; iC1 < nbColor; iC1++)
      for (int iC2 = iC1 + 1; iC2 < nbColor; iC2++) {
        int pos1 = nbColor * iVert + iC1;
        int pos2 = nbColor * iVert + iC2;
        Minisat::Lit cond1 = Minisat::mkLit(ListVar[pos1], false);
        Minisat::Lit cond2 = Minisat::mkLit(ListVar[pos2], false);
        S.addClause(cond1, cond2);
      }
  }
  //
  // For each edge we should not have the same color on both
  //
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (auto &jVert : eGR.Adjacency(iVert)) {
      for (int iColor = 0; iColor < nbColor; iColor++) {
        int pos1 = nbColor * iVert + iColor;
        int pos2 = nbColor * jVert + iColor;
        Minisat::Lit cond1 = Minisat::mkLit(ListVar[pos1], false);
        Minisat::Lit cond2 = Minisat::mkLit(ListVar[pos2], false);
        S.addClause(cond1, cond2);
      }
    }
  }
  //
  // Solving the system
  //
  if (!S.solve())
    return {false, {}};
  //
  // Extracting the solution
  //
  std::vector<int> V(nbVert, -1);
  for (int iVert = 0; iVert < nbVert; iVert++)
    for (int iColor = 0; iColor < nbColor; iColor++) {
      int pos = nbColor * iVert + iColor;
      using namespace Minisat;
      if (S.modelValue(pos) ==
          l_False) { // Curiously, this is how we should do it, it seems.
        if (V[iVert] != -1) {
          std::cerr << "iVert=" << iVert << " has already been assigned\n";
          throw TerminalException{1};
        }
        V[iVert] = iColor;
      }
    }
  for (int iVert = 0; iVert < nbVert; iVert++) {
    if (V[iVert] == -1) {
      std::cerr << "iVert=" << iVert << " has not been assigned\n";
      throw TerminalException{1};
    }
  }
  for (int iVert = 0; iVert < nbVert; iVert++) {
    for (auto &jVert : eGR.Adjacency(iVert)) {
      if (V[iVert] == V[jVert]) {
        std::cerr << "Adjacent vertices iVert=" << iVert << " jVert=" << jVert
                  << " have same color\n";
        throw TerminalException{1};
      }
    }
  }
  return {true, std::move(V)};
}
