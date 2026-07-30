#!/bin/bash
# Correctness tests of the fraction-free adjugate machinery of
# src_matrix/MAT_MatrixInverse.h and MAT_MatrixInt.h:
#   * AdjugateDeterminant: A adj(A) = adj(A) A = det(A) I, determinant
#     matching DeterminantMat, row-swap path, singular rejection;
#   * InverseFractionFreeLU on unimodular matrices;
#   * SelectIndependentRowsRing against the field selection;
#   * the ring CanonicalizeOrderedMatrix against the computation over
#     the overlying field.
#
# This builds and drives src_matrix/Test_AdjugateDeterminant over the
# integer rings, the rational fields and several matrix sizes. Any
# failure aborts the script with a non-zero exit status.
#
# Usage:  ./run_test.sh [n1 n2 ...]        (default sizes: 4 6 8)
#
# Honours the same environment variables as src_matrix/Makefile
# (CXX, GMP_INCDIR, BOOST_INCDIR, EIGEN_PATH, GMP_CXX_LINK); Homebrew
# defaults are used when they are unset.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/../.."
SRC="$ROOT/src_matrix"

: "${CXX:=g++}"
: "${GMP_INCDIR:=/opt/homebrew/include}"
: "${BOOST_INCDIR:=/opt/homebrew/include}"
: "${EIGEN_PATH:=/opt/homebrew/include/eigen3}"
: "${GMP_LIBDIR:=/opt/homebrew/lib}"
: "${GMP_CXX_LINK:=-L${GMP_LIBDIR} -lgmp -lgmpxx}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

CXXFLAGS="-std=c++20 -Wall -Wextra -O3 -DSANITY_CHECK -I${SRC} -I${ROOT}/src_basic -I${ROOT}/src_number -I${ROOT}/src_comb -I${GMP_INCDIR} -I${BOOST_INCDIR} -I${EIGEN_PATH}"
LDFLAGS="-lm ${GMP_CXX_LINK} -pthread"

echo "Building Test_AdjugateDeterminant ..."
"$CXX" $CXXFLAGS "$SRC/Test_AdjugateDeterminant.cpp" -o "$WORK/Test_AdjugateDeterminant" $LDFLAGS

SIZES="${*:-4 6 8}"
for n in $SIZES; do
  for arith in mpz_class safe_integer boost_cpp_int mpq_class safe_rational; do
    echo "----- n = $n arith = $arith -----"
    "$WORK/Test_AdjugateDeterminant" "$arith" "$n"
  done
done
echo "All the adjugate determinant tests passed"
