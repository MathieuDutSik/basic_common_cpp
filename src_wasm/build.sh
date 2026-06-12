#!/bin/bash
set -euo pipefail

# Build the Wasm smoke test(s) in this directory and run them under node.
# Requires emscripten (emcc) on PATH. On macOS:  brew install emscripten

cd "$(dirname "$0")"

BOOST_INC="${BOOST_INC:-/opt/homebrew/opt/boost/include}"
# Any Eigen version works now that src_matrix/EigenBoostNumTraits.h
# overrides NumTraits<cpp_int>::Literal = int. Default to system Eigen.
EIGEN_INC="${EIGEN_INC:-/opt/homebrew/include/eigen3}"

if [ ! -d "$BOOST_INC" ]; then
  echo "Boost include dir not found at $BOOST_INC (set BOOST_INC env var)" >&2
  exit 1
fi
if [ ! -d "$EIGEN_INC" ]; then
  echo "Eigen include dir not found at $EIGEN_INC (set EIGEN_INC env var)" >&2
  exit 1
fi

ROOT="$(cd .. && pwd)"
INCLUDES=(
  -I"$ROOT/src_basic"
  -I"$ROOT/src_number"
  -I"$ROOT/src_matrix"
  -I"$ROOT/src_comb"
  -I"$ROOT/sparse-map/include/tsl"
  -I"$ROOT/robin-map/include/tsl"
  -I"$ROOT/hopscotch-map/include/tsl"
  # External: -idirafter (not -I) so the host include tree never shadows the
  # wasm sysroot. With BOOST_INC=/usr/include (Linux CI), plain -I would
  # inject the host glibc tree ahead of wasi-libc and libc++'s
  # `#include_next <stdlib.h>` would land on /usr/include/stdlib.h and fail
  # on glibc-only bits/libc-header-start.h. -idirafter appends these strictly
  # after the sysroot search path.
  -idirafter "$BOOST_INC"
  -idirafter "$EIGEN_INC"
)

CXXFLAGS=(
  -std=c++20
  -O2
  # Keep the generated JS glue readable (multi-line, comments preserved)
  # instead of being collapsed to a single minified line. Wasm-side
  # optimization is unaffected.
  --minify
  0
  -Wall
  -Wextra
  -Wno-deprecated-declarations
  # libc++ removed std::result_of in C++20; Eigen 3.3.9 still uses it.
  -D_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS
  -DINCLUDE_NUMBER_THEORY_BOOST_CPP_INT
)

mkdir -p build
shopt -s nullglob
sources=(Test_wasm_*.cpp)
shopt -u nullglob

if [ "${#sources[@]}" -eq 0 ]; then
  echo "No Test_wasm_*.cpp sources found." >&2
  exit 1
fi

for src in "${sources[@]}"; do
  name="${src%.cpp}"
  out="build/${name}.js"
  echo "==> Building $src -> $out"
  emcc "${CXXFLAGS[@]}" "${INCLUDES[@]}" "$src" -o "$out"
  echo "==> Running $out under node"
  node "$out"
done

echo "All Wasm tests passed."
