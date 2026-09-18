#!/bin/bash
# Smith normal form invariant factors on the boundary matrices of a chain
# complex.
#
# The matrices Matrix1_8.sparse .. Matrix7_8.sparse are the boundary maps
# of one cell complex, in the sparse format of ReadSparseMatrix. Their
# dimensions chain up (870 - 5253 - 11568 - 12119 - 6285 - 1477 - 126),
# they are between 0.06% and 1.8% dense, and their entries are almost all
# +-1. This is the regime the dense algorithms are worst at and the one
# the sparse path exists for: held densely, Matrix3_8 is 140 million cells
# for 81720 entries that carry a value.
#
# The test checks two things:
#   * the invariant factors are the expected ones, on every matrix and
#     for each arithmetic;
#   * the arithmetics agree with each other.
#
# The expected values are recorded below. They were cross-checked against
# the dense-input program (SmithNormalFormInvariant) on the matrices where
# that is feasible.
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

echo "Building SmithNormalFormInvariant_sparse ..."
"$CXX" $CXXFLAGS "$SRC/SmithNormalFormInvariant_sparse.cpp" \
  -o "$WORK/SmithNormalFormInvariant_sparse" $LDFLAGS

# The invariant factors of each matrix, as the multiset that the program
# prints: [value,multiplicity] in increasing order of the value.
expected_1="MultInv = [0,1] [1,862] [3,2] [6,3] [12,2]"
expected_2="MultInv = [0,869] [1,4384]"
expected_3="MultInv = [0,4384] [1,7179] [2,5]"
expected_4="MultInv = [0,1360] [1,4912] [2,10] [6,3]"
expected_5="MultInv = [0,122] [1,1354] [2,1]"
expected_6="MultInv = [0,4] [1,122]"
expected_7="MultInv ="

# Matrix3_8 and Matrix4_8 are the expensive ones (tens of seconds); they
# are covered by the default run since they are the cases that motivate
# the sparse path, but SNF_QUICK=1 restricts the run to the cheap ones.
if [ "${SNF_QUICK:-0}" = "1" ]; then
  matrices="1 2 5 6 7"
else
  matrices="1 2 3 4 5 6 7"
fi

for i in $matrices; do
  eval "expected=\$expected_$i"
  for arith in mpz_class boost_cpp_int; do
    echo "----- Matrix${i}_8 arith = $arith -----"
    obtained=$("$WORK/SmithNormalFormInvariant_sparse" "$arith" \
      "$HERE/Matrix${i}_8.sparse" 2>&1 | grep "^MultInv" || true)
    # grep strips nothing else; compare after collapsing the whitespace.
    obtained_n=$(echo "$obtained" | tr -s ' ' | sed 's/ *$//')
    expected_n=$(echo "$expected" | tr -s ' ' | sed 's/ *$//')
    if [ "$obtained_n" != "$expected_n" ]; then
      echo "MISMATCH on Matrix${i}_8 with $arith"
      echo "  expected: $expected_n"
      echo "  obtained: $obtained_n"
      exit 1
    fi
    echo "  $obtained_n"
  done
done
echo "All the Smith normal form computations gave the expected invariant factors"
