// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORY_H_
#define SRC_NUMBER_NUMBERTHEORY_H_

// NumberTheoryGmp.h comes first so that the gmp types are declared by the time
// NumberTheoryCommon.h reaches TypeConversionFinal.h, whose mixed gmp/boost
// conversions are guarded on SRC_NUMBER_NUMBERTHEORYGMP_H_. In the other order
// the include guard of TypeConversionFinal.h fires before gmp is known and
// that block is skipped for good.
// clang-format off
#include "NumberTheoryGmp.h"
#include "NumberTheoryCommon.h"
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORY_H_
// clang-format on
