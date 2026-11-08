#!/bin/bash
# Correctness of the Eigen matrix products (src_matrix/Test_EigenFlintProduct):
# every form of matrix-matrix and matrix-vector product against a triple
# loop. With ENABLE_FLINT_SUPPORT=1 this covers the products of fmpz_class /
# fmpq_class, which go to fmpz_mat_mul / fmpq_mat_mul
# (src_matrix/EigenFlintProduct.h); without it, the generic kernel for the
# gmp and machine types.
#
# Any failure aborts the script with a non-zero exit status.
#
# Honours the same environment variables as src_matrix/Makefile
# (CXX, GMP_INCDIR, BOOST_INCDIR, EIGEN_PATH, GMP_CXX_LINK, and with
# ENABLE_FLINT_SUPPORT=1, FLINT_INCDIR and FLINT_LINK); Homebrew defaults
# are used when they are unset.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/../.."

: "${CXX:=g++}"
: "${GMP_INCDIR:=/opt/homebrew/include}"
: "${BOOST_INCDIR:=/opt/homebrew/include}"
: "${EIGEN_PATH:=/opt/homebrew/include/eigen3}"
: "${GMP_LIBDIR:=/opt/homebrew/lib}"
: "${GMP_CXX_LINK:=-L${GMP_LIBDIR} -lgmp -lgmpxx}"
: "${BOOST_LINK:=-L${GMP_LIBDIR} -lboost_serialization}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

CXXFLAGS="-std=c++20 -Wall -Wextra -O3 -I${ROOT}/src_basic -I${ROOT}/src_number -I${ROOT}/src_matrix -I${ROOT}/src_comb -I${GMP_INCDIR} -I${BOOST_INCDIR} -I${EIGEN_PATH}"
LDFLAGS="-lm ${GMP_CXX_LINK} ${BOOST_LINK} -pthread"
if [ "${ENABLE_FLINT_SUPPORT:-0}" = "1" ]; then
  : "${FLINT_INCDIR:=${GMP_INCDIR}}"
  : "${FLINT_LINK:=-lflint}"
  CXXFLAGS="${CXXFLAGS} -DENABLE_FLINT_SUPPORT -I${FLINT_INCDIR}"
  LDFLAGS="${LDFLAGS} ${FLINT_LINK}"
  echo "With the flint types"
fi

echo "Building Test_EigenFlintProduct ..."
"$CXX" $CXXFLAGS "$ROOT/src_matrix/Test_EigenFlintProduct.cpp" -o "$WORK/Test_EigenFlintProduct" $LDFLAGS

echo
echo "===== Test_EigenFlintProduct ====="
"$WORK/Test_EigenFlintProduct"

echo
echo "All the matrix product tests passed"
