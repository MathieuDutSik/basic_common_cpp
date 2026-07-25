#!/bin/bash
# Iterate the determinant / inverse benchmarks over the matrix dimension n.
#
# This does NOT add any new program: it builds and drives the existing
# src_matrix/DeterminantMat and src_matrix/Inverse, which already print the
# per-method timings. The script only generates input matrices and loops over n.
#
# For every n it generates:
#   * a UNIMODULAR (det = 1) integer matrix (A = L*U of unit-triangular integer
#     factors), so the inverse is integral and the 'integer' (mpz_class) path is
#     exercised -- the common det=1 use case;
#   * a general diagonally-dominant matrix, run with 'rational' (mpq_class).
#
# Inverse runs for every n (all methods are O(n^3)). DeterminantMat also runs an
# O(n!) permutation method, so it is only invoked for n <= 8.
#
# Usage:  ./run_benchmark.sh [n1 n2 ...]        (default sizes: 5 8 10 20 40 80)
#
# Honours the same environment variables as src_matrix/Makefile
# (CXX, GMP_INCDIR, BOOST_INCDIR, EIGEN_PATH, GMP_CXX_LINK); Homebrew defaults
# are used when they are unset.
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

CXXFLAGS="-std=c++20 -O3 -I${SRC} -I${ROOT}/src_basic -I${ROOT}/src_number -I${ROOT}/src_comb -I${GMP_INCDIR} -I${BOOST_INCDIR} -I${EIGEN_PATH}"
LDFLAGS="-lm ${GMP_CXX_LINK} -pthread"

echo "Building DeterminantMat and Inverse ..."
"$CXX" $CXXFLAGS "$SRC/DeterminantMat.cpp" -o "$WORK/DeterminantMat" $LDFLAGS
"$CXX" $CXXFLAGS "$SRC/Inverse.cpp"        -o "$WORK/Inverse"        $LDFLAGS

# gen_matrix <unimodular|general> <n>  ->  matrix file on stdout (ReadMatrix format)
gen_matrix() {
  python3 - "$1" "$2" <<'PY'
import sys
kind, n = sys.argv[1], int(sys.argv[2])
st = (0x2545F4914F6CDD1D ^ (n * 2654435761)) & 0xFFFFFFFFFFFFFFFF
def r(rng):
    global st
    st ^= (st << 13) & 0xFFFFFFFFFFFFFFFF
    st ^= st >> 7
    st ^= (st << 17) & 0xFFFFFFFFFFFFFFFF
    return st % (2 * rng + 1) - rng
if kind == "unimodular":
    # A = L * U with unit-triangular integer factors  ->  det(A) = 1
    L = [[1 if i == j else (r(3) if j < i else 0) for j in range(n)] for i in range(n)]
    U = [[1 if i == j else (r(3) if j > i else 0) for j in range(n)] for i in range(n)]
    A = [[sum(L[i][k] * U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
else:
    # diagonally dominant integer matrix  ->  invertible over Q
    A = [[r(5) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        A[i][i] = sum(abs(A[i][j]) for j in range(n)) + 1
out = ["%d %d" % (n, n)]
out += [" ".join(str(x) for x in row) for row in A]
print("\n".join(out))
PY
}

SIZES="${*:-5 8 10 20 40 80}"
for n in $SIZES; do
  gen_matrix unimodular "$n" > "$WORK/unimod_$n.txt"
  gen_matrix general    "$n" > "$WORK/general_$n.txt"

  echo
  echo "############################## n = $n ##############################"
  if [ "$n" -le 8 ]; then
    echo "----- DeterminantMat integer (unimodular, det=1) -----"
    "$WORK/DeterminantMat" integer  "$WORK/unimod_$n.txt"  2>&1 | grep -v "termination"
    echo "----- DeterminantMat rational (general) -----"
    "$WORK/DeterminantMat" rational "$WORK/general_$n.txt" 2>&1 | grep -v "termination"
  else
    echo "(DeterminantMat skipped for n>8: it runs an O(n!) permutation method)"
  fi
  echo "----- Inverse integer (unimodular, det=1) -----"
  "$WORK/Inverse" integer  "$WORK/unimod_$n.txt"  2>&1 | grep -v "termination"
  echo "----- Inverse rational (general) -----"
  "$WORK/Inverse" rational "$WORK/general_$n.txt" 2>&1 | grep -v "termination"
done
