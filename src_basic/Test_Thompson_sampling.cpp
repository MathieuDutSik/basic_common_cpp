// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Heuristic_fct.h"

int main() {
  //  using T = mpq_class;
  using T = T_uint64_t;
  FullNamelist eFull = NAMELIST_GetStandard_RecursiveDualDescription();
  ThompsonSamplingHeuristic<T> TSH(std::cerr, eFull);

  std::map<std::string,T> TheCand;
  std::string choice = TSH.GetEvaluation(TheCand);
  double result = 2.3;
  TSH.SubmitResult(result);
}
