// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORY_H_
#define SRC_NUMBER_NUMBERTHEORY_H_

// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
// The pooled GMP allocator. Inert unless the build defines GMP_POOL, in which
// case it installs itself once at startup (see gmp_pool_allocator.h). Included
// here so every program using the number types picks it up when enabled.
#include "gmp_pool_allocator.h"
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORY_H_
// clang-format on
