#!/bin/bash
# Correctness tests of the key functionality of src_matrix:
#   * Test_HermiteSmith: the row and column Hermite normal forms
#     (H = U M, unimodularity, echelon shape, canonicity) and the Smith
#     normal form (ROW M COL diagonal, divisibility chain, invariance of
#     the invariant factors under unimodular multiplication);
#   * Test_SolutionNullspaceInt: SolutionMat / SolutionIntMat (recovery
#     of planted solutions, rejection of unsolvable systems),
#     NullspaceIntMat / NullspaceIntTrMat (kernel, nullity, saturation),
#     GetZbasis (double inclusion) and IntersectionLattice;
#   * Test_RankMat: the fraction-free Bareiss rank over the rings
#     against the Gaussian elimination over the overlying field, with
#     planted low rank factorizations and rank invariances;
#   * Test_ScalingCanonic: RemoveFractionVector / RemoveFractionMatrix
#     (integrality, content one, positive multiplier),
#     NonUniqueScaleToIntegerVector, CanonicalizeVector and
#     ScalarCanonicalizationMatrix (invariance under positive rescaling).
#
# Any failure aborts the script with a non-zero exit status. The
# safe_integer runs use a smaller dimension: the intermediate entries of
# the normal form computations exceed the 64 bit range beyond it.
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

for prog in Test_HermiteSmith Test_SolutionNullspaceInt Test_ScalingCanonic Test_RankMat; do
  echo "Building $prog ..."
  "$CXX" $CXXFLAGS "$SRC/$prog.cpp" -o "$WORK/$prog" $LDFLAGS
done

for prog in Test_HermiteSmith Test_SolutionNullspaceInt; do
  for n in 5 7; do
    for arith in mpz_class boost_cpp_int; do
      echo "----- $prog n = $n arith = $arith -----"
      "$WORK/$prog" "$arith" "$n"
    done
  done
  echo "----- $prog n = 3 arith = safe_integer -----"
  "$WORK/$prog" safe_integer 3
done

for n in 5 8; do
  for arith in mpz_class boost_cpp_int mpq_class; do
    echo "----- Test_RankMat n = $n arith = $arith -----"
    "$WORK/Test_RankMat" "$arith" "$n"
  done
done
for arith in safe_integer safe_rational; do
  echo "----- Test_RankMat n = 5 arith = $arith -----"
  "$WORK/Test_RankMat" "$arith" 5
done

for n in 6 10; do
  for arith in mpq_class safe_rational boost_cpp_rational boost_mpq_rational; do
    echo "----- Test_ScalingCanonic n = $n arith = $arith -----"
    "$WORK/Test_ScalingCanonic" "$arith" "$n"
  done
done
echo "All the matrix functionality tests passed"
