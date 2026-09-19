#!/bin/bash
# Tests of the modular matrix operations of src_matrix/MAT_MatrixMod.h,
# and of RecSolutionMatMod in particular: the membership of a vector in a
# space given modulo p.
#
# The program Test_MatrixMod is self checking. For each modulus in
# 2, 3, 5, 7, 11, 101 it draws spaces of every shape up to 7 x 6, with
# spanning families that are free and families that are not, and decides
# membership twice: once through the structure, once through the rank
# computed by SelectRowColMatMod, which shares the elimination but not
# the reformulation by equations. It also checks that a combination of
# the rows is accepted, that is_containing_m is the conjunction of
# has_solution_v over the rows, that a shift of the entries by large
# multiples of the modulus changes nothing, and the two extreme spaces
# (the full one and the zero one).
#
# The three arithmetics are run because they fail differently: on the
# unbounded ones a mistake is a wrong answer, whereas SafeInt64 also
# catches the intermediate values leaving the range of a 64 bit integer.
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

CXXFLAGS="-std=c++20 -Wall -Wextra -O3 -I${SRC} -I${ROOT}/src_basic -I${ROOT}/src_number -I${ROOT}/src_comb -I${GMP_INCDIR} -I${BOOST_INCDIR} -I${EIGEN_PATH}"
LDFLAGS="-lm ${GMP_CXX_LINK} -pthread"

echo "Building Test_MatrixMod ..."
"$CXX" $CXXFLAGS "$SRC/Test_MatrixMod.cpp" \
  -o "$WORK/Test_MatrixMod" $LDFLAGS

# The number of spaces tried per modulus. The default of the program is
# 42; MATRIXMOD_QUICK=1 cuts it down for a faster run.
if [ "${MATRIXMOD_QUICK:-0}" = "1" ]; then
  n_iter=8
else
  n_iter=42
fi

for arith in safe_integer mpz_class boost_cpp_int; do
  echo "----- arith = $arith -----"
  # The program reports its own failures and exits non-zero; set -e then
  # stops the script. The output is kept so that a failure is readable.
  if ! "$WORK/Test_MatrixMod" "$arith" "$n_iter" >"$WORK/out.txt" 2>&1; then
    echo "FAILURE with $arith"
    cat "$WORK/out.txt"
    exit 1
  fi
  grep "TESTMOD: random vectors" "$WORK/out.txt"
  grep "TESTMOD: all the" "$WORK/out.txt"
done
echo "All the modular matrix tests passed"
