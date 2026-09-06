#!/bin/bash
# Tests of the real algebraic arithmetic of src_number/NumberTheoryRealField.h:
# the field RealField and the order RealRing = Z[x] spanned by the powers of
# the generator, which is what underlying_ring<RealField> returns.
#
#   * Test_RealCubicField:        the field itself.
#   * Test_RealRing:              the ring itself -- arithmetic, ordering, the
#                                 exact division and its failure mode, the
#                                 rejection of a non-monic minimal polynomial,
#                                 the conversions and the string round trip.
#   * Test_RealRingConsistency:   the field and the ring code paths return the
#                                 same answers -- determinant, product, rank,
#                                 adjugate, the canonicalizations, the subset
#                                 solver, the scaling into the ring, the
#                                 ordering and the division. These are the
#                                 dispatches that see, in the ring, a type that
#                                 is neither a field nor a euclidean domain.
#   * Bench_real_ring:            the field against the ring on the same
#                                 matrices, printing the timing comparison and
#                                 checking that both compute the same
#                                 determinant. The timings are reported, never
#                                 asserted: they vary with the runner.
#
# All of it runs on the cubic field of discriminant 49 described by
# CubicFieldDisc_49 in this directory (the generator is 2*cos(2*pi/7), of
# minimal polynomial X^3 + X^2 - 2X - 1, which is monic).
#
# Usage:  ./run_test.sh [n_trials] [bench_size]      (defaults: 50 8)
#
# Honours the same environment variables as src_matrix/Makefile
# (CXX, GMP_INCDIR, BOOST_INCDIR, EIGEN_PATH, GMP_CXX_LINK); Homebrew
# defaults are used when they are unset.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/../.."

: "${CXX:=g++}"
: "${GMP_INCDIR:=/opt/homebrew/include}"
: "${BOOST_INCDIR:=/opt/homebrew/include}"
: "${EIGEN_PATH:=/opt/homebrew/include/eigen3}"
: "${GMP_LIBDIR:=/opt/homebrew/lib}"
: "${GMP_CXX_LINK:=-L${GMP_LIBDIR} -lgmp -lgmpxx}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

CXXFLAGS="-std=c++20 -Wall -Wextra -O3 -I${ROOT}/src_basic -I${ROOT}/src_number -I${ROOT}/src_matrix -I${ROOT}/src_comb -I${GMP_INCDIR} -I${BOOST_INCDIR} -I${EIGEN_PATH}"
LDFLAGS="-lm ${GMP_CXX_LINK} -pthread"

N_TRIALS="${1:-50}"
BENCH_SIZE="${2:-8}"

# The programs locate CubicFieldDisc_49 by walking up from the current
# directory, so they are run from the repository root.
cd "$ROOT"

for prog in src_number/Test_RealCubicField src_number/Test_RealRing \
            src_matrix/Test_RealRingConsistency src_number/Bench_real_ring; do
  name="$(basename "$prog")"
  echo "Building $name ..."
  "$CXX" $CXXFLAGS "$ROOT/$prog.cpp" -o "$WORK/$name" $LDFLAGS
done

echo
echo "===== Test_RealCubicField ====="
"$WORK/Test_RealCubicField"

echo
echo "===== Test_RealRing ====="
"$WORK/Test_RealRing"

echo
echo "===== Test_RealRingConsistency ($N_TRIALS trials per section) ====="
"$WORK/Test_RealRingConsistency" "$N_TRIALS"

echo
echo "===== Bench_real_ring (RealField against RealRing) ====="
"$WORK/Bench_real_ring" "$BENCH_SIZE" 20 20

echo
echo "All the real algebraic field and ring tests passed"
