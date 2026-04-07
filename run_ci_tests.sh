#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

run_step() {
  local step_name="$1"
  shift
  echo
  echo "== $step_name =="
  "$@"
}

cd "$ROOT_DIR"

#run_step "Compile Sweep" ./compile.sh

run_step "Build Basic Tests" make -C src_basic \
  Test_timing \
  Test_PresenceProgram \
  Test_namelist \
  Test_String_Split \
  Test_serialization \
  Test_EmpiricalDistribution \
  Test_hash_table \
  Test_bool \
  Test_face \
  Test_String_Optimization \
  Test_String_Trim \
  Test_String_Conversion \
  Test_Namelist_Substr \
  Test_Basic_string_Optimizations \
  Test_filesystem \
  Test_program_in_path

#run_step "Test_Thompson_sampling" ./src_basic/Test_Thompson_sampling
run_step "Test_timing" ./src_basic/Test_timing
run_step "Test_PresenceProgram" ./src_basic/Test_PresenceProgram
run_step "Test_namelist" ./src_basic/Test_namelist
run_step "Test_String_Split" ./src_basic/Test_String_Split
run_step "Test_serialization" ./src_basic/Test_serialization
run_step "Test_EmpiricalDistribution" ./src_basic/Test_EmpiricalDistribution
run_step "Test_hash_table" ./src_basic/Test_hash_table
run_step "Test_bool" ./src_basic/Test_bool
run_step "Test_face" ./src_basic/Test_face
run_step "Test_String_Optimization" ./src_basic/Test_String_Optimization
run_step "Test_String_Trim" ./src_basic/Test_String_Trim
run_step "Test_String_Conversion" ./src_basic/Test_String_Conversion
run_step "Test_Namelist_Substr" ./src_basic/Test_Namelist_Substr
run_step "Test_Basic_string_Optimizations" ./src_basic/Test_Basic_string_Optimizations
run_step "Test_filesystem" ./src_basic/Test_filesystem
run_step "Test_program_in_path" ./src_basic/Test_program_in_path

run_step "Build Number Theory Tests" make -C src_number \
  Test_SquareRoot \
  Test_Factorize \
  Test_QuoInt \
  Test_TypeBoostGmp \
  Test_GetBit \
  Test_QuadraticResidue \
  Test_PrimeGenerator \
  Test_SequenceApproximant \
  Test_UnorderedMapMpzq \
  Test_ComputePairGcdDot \
  Test_RealCubicField \
  Test_QuadField \
  Test_Rational \
  Test_PracticalInf \
  Test_PrintStaticInfo

run_step "Test_SquareRoot" ./src_number/Test_SquareRoot
run_step "Test_Factorize" ./src_number/Test_Factorize check
run_step "Test_QuoInt" ./src_number/Test_QuoInt
run_step "Test_TypeBoostGmp" ./src_number/Test_TypeBoostGmp
run_step "Test_GetBit" ./src_number/Test_GetBit
run_step "Test_QuadraticResidue" ./src_number/Test_QuadraticResidue
run_step "Test_PrimeGenerator" ./src_number/Test_PrimeGenerator mpz 4000
run_step "Test_SequenceApproximant" ./src_number/Test_SequenceApproximant print
run_step "Test_UnorderedMapMpzq" ./src_number/Test_UnorderedMapMpzq
run_step "Test_ComputePairGcdDot int64_t" ./src_number/Test_ComputePairGcdDot int64_t
run_step "Test_ComputePairGcdDot mpz_class" ./src_number/Test_ComputePairGcdDot mpz_class
run_step "Test_ComputePairGcdDot SafeInt64" ./src_number/Test_ComputePairGcdDot SafeInt64
run_step "Test_ComputePairGcdDot boost_cpp_int" ./src_number/Test_ComputePairGcdDot boost_cpp_int
run_step "Test_RealCubicField" ./src_number/Test_RealCubicField
run_step "Test_QuadField" ./src_number/Test_QuadField
run_step "Test_Rational" ./src_number/Test_Rational
run_step "Test_PracticalInf" ./src_number/Test_PracticalInf
run_step "Test_PrintStaticInfo" ./src_number/Test_PrintStaticInfo

run_step "Build Matrix Tests" make -C src_matrix \
  Test_MatrixInverse \
  Test_PerformanceHNF \
  Test_HilbertMatrix \
  Test_NullspaceComputation \
  Test_SubspaceCompletion \
  Test_FindIsotropicMod

run_step "Test_MatrixInverse" ./src_matrix/Test_MatrixInverse mpz_class 10
run_step "Test_PerformanceHNF" ./src_matrix/Test_PerformanceHNF 10 10
run_step "Test_HilbertMatrix" ./src_matrix/Test_HilbertMatrix mpq_class 10
run_step "Test_NullspaceComputation" ./src_matrix/Test_NullspaceComputation
run_step "Test_SubspaceCompletion" ./src_matrix/Test_SubspaceCompletion 15 10

echo
echo "CI test sequence completed successfully."
