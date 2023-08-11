// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYMULTITYPE_H_
#define SRC_NUMBER_NUMBERTHEORYMULTITYPE_H_

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
// clang-format on




template<typename F, typename... Targs>
void process_by_numeric_type(std::stirng const& arith, F f, Targs... args)
{
  if (arith == "rational") {
    using T = mpq_class;
    return f<T>(args...);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return f<T>(args...);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return f<T>(args...);
  }
  std::optional<std::string> opt_realalgebraic =
      get_postfix(arith, "RealAlgebraic=");
  if (opt_realalgebraic) {
    std::string const &FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    using T_rat = mpq_class;
    HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return f<T>(args...);
  }
  std::cerr << "Failed to find a matching arithmetic for f\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYMULTITYPE_H_
// clang-format on
