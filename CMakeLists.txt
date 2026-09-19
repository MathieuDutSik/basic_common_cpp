cmake_minimum_required(VERSION 3.30)

project(basic_common_cpp LANGUAGES CXX)

include(ExternalProject)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# On Windows (LLP64), Eigen::Index defaults to std::ptrdiff_t == `long long`,
# for which gmpxx provides no constructor, so Eigen-internal conversions such as
# RealScalar(Index) are ambiguous for mpz_class / mpq_class. Force Eigen::Index
# to `long`, which gmpxx accepts. On LP64 (Linux/macOS) std::ptrdiff_t already
# is `long`, so this is a no-op there.
if(WIN32)
  add_compile_definitions(EIGEN_DEFAULT_DENSE_INDEX_TYPE=long)
endif()

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS serialization)

pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
find_library(GMP_LIBRARY NAMES gmp REQUIRED)

# Optional flint arithmetic (fmpz_class / fmpq_class of
# src_number/NumberTheoryFlint.h and the matrix backend of
# src_matrix/MAT_MatrixFlint.h), the counterpart of the
# `make ENABLE_FLINT_SUPPORT=1` of the Makefiles.
#
# ENABLE_FLINT_SUPPORT drives the search: AUTO, the default, uses flint
# when it is found and builds without it otherwise; OFF ignores it even
# when installed; anything else, including the `1` that the Makefiles
# take, requires it and fails the configuration when it is missing.
set(
  ENABLE_FLINT_SUPPORT "AUTO"
  CACHE STRING
  "Flint arithmetic: AUTO (use when found), ON (require), OFF (ignore)"
)
set_property(CACHE ENABLE_FLINT_SUPPORT PROPERTY STRINGS AUTO ON OFF)

string(TOUPPER "${ENABLE_FLINT_SUPPORT}" FLINT_REQUEST)
if(FLINT_REQUEST STREQUAL "AUTO")
  set(FLINT_WANTED TRUE)
  set(FLINT_REQUIRED FALSE)
elseif(FLINT_REQUEST MATCHES "^(OFF|0|NO|FALSE|N|IGNORE|NOTFOUND)$")
  set(FLINT_WANTED FALSE)
  set(FLINT_REQUIRED FALSE)
else()
  set(FLINT_WANTED TRUE)
  set(FLINT_REQUIRED TRUE)
endif()

set(FLINT_USABLE FALSE)
if(FLINT_WANTED)
  # flint 3 ships a pkg-config file; older installations and some
  # distributions do not, hence the search by hand as a fallback.
  pkg_check_modules(FLINT IMPORTED_TARGET flint)
  if(FLINT_FOUND)
    set(FLINT_PROBE_INCLUDES "${FLINT_INCLUDE_DIRS}")
    set(FLINT_PROBE_LIBRARIES "${FLINT_LDFLAGS}")
    set(FLINT_ORIGIN "pkg-config, version ${FLINT_VERSION}")
  else()
    find_path(FLINT_INCLUDE_DIR NAMES flint/fmpz.h)
    find_library(FLINT_LIBRARY NAMES flint)
    if(FLINT_INCLUDE_DIR AND FLINT_LIBRARY)
      set(FLINT_PROBE_INCLUDES "${FLINT_INCLUDE_DIR}")
      set(FLINT_PROBE_LIBRARIES "${FLINT_LIBRARY}" "${GMP_LIBRARY}")
      set(FLINT_ORIGIN "${FLINT_LIBRARY}")
    endif()
  endif()

  # Being installed is not the same as being usable: the entry accessors
  # and the normal form routines below are what the two headers call, and
  # they are the parts that moved between the flint versions. A flint too
  # old to compile them is treated as absent rather than allowed to break
  # the build later. The outcome is cached as FLINT_COMPILES, so a flint
  # installed after a configuration that did not find one is picked up
  # only once that entry is removed from the cache.
  if(DEFINED FLINT_PROBE_INCLUDES)
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES ${FLINT_PROBE_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${FLINT_PROBE_LIBRARIES})
    check_cxx_source_compiles(
      "
      #include <flint/fmpz.h>
      #include <flint/fmpq.h>
      #include <flint/fmpz_mat.h>
      #include <flint/fmpq_mat.h>
      int main() {
        fmpz_t a;
        fmpz_init(a);
        fmpz_set_si(a, 3);
        fmpz_mat_t m, h;
        fmpz_mat_init(m, 2, 2);
        fmpz_mat_init(h, 2, 2);
        fmpz_set(fmpz_mat_entry(m, 0, 0), a);
        fmpz_mat_hnf(h, m);
        fmpz_mat_snf(h, m);
        fmpz_mat_clear(h);
        fmpz_mat_clear(m);
        fmpq_mat_t q;
        fmpq_mat_init(q, 2, 2);
        fmpq_mat_rref(q, q);
        fmpq_mat_clear(q);
        fmpz_clear(a);
        return 0;
      }
      "
      FLINT_COMPILES
    )
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(FLINT_COMPILES)
      set(FLINT_USABLE TRUE)
    endif()
  endif()
endif()

add_library(basic_common_flint INTERFACE)
if(FLINT_USABLE)
  target_compile_definitions(basic_common_flint INTERFACE ENABLE_FLINT_SUPPORT)
  if(FLINT_FOUND)
    target_link_libraries(basic_common_flint INTERFACE PkgConfig::FLINT)
  else()
    target_include_directories(basic_common_flint INTERFACE "${FLINT_INCLUDE_DIR}")
    target_link_libraries(basic_common_flint INTERFACE "${FLINT_LIBRARY}")
  endif()
  message(STATUS "flint found (${FLINT_ORIGIN}); building with ENABLE_FLINT_SUPPORT")
elseif(FLINT_REQUIRED)
  if(DEFINED FLINT_PROBE_INCLUDES)
    message(FATAL_ERROR
      "ENABLE_FLINT_SUPPORT=${ENABLE_FLINT_SUPPORT} was asked for but the "
      "flint found at ${FLINT_ORIGIN} does not compile the routines used by "
      "src_number/NumberTheoryFlint.h and src_matrix/MAT_MatrixFlint.h. "
      "See CMakeFiles/CMakeError.log for the compiler output.")
  else()
    message(FATAL_ERROR
      "ENABLE_FLINT_SUPPORT=${ENABLE_FLINT_SUPPORT} was asked for but flint "
      "was not found. Install it, or set ENABLE_FLINT_SUPPORT=OFF to build "
      "without it.")
  endif()
elseif(FLINT_WANTED)
  message(STATUS "flint not usable; building without the flint arithmetic")
else()
  message(STATUS "flint disabled by ENABLE_FLINT_SUPPORT=${ENABLE_FLINT_SUPPORT}")
endif()

set(NAUTY_PREFIX "${CMAKE_BINARY_DIR}/external/nauty")
set(NAUTY_INSTALL_DIR "${NAUTY_PREFIX}/install")
set(NAUTY_INCLUDE_DIR "${NAUTY_INSTALL_DIR}/include")
set(NAUTY_LIBRARY "${NAUTY_INSTALL_DIR}/lib/libnauty.a")

ExternalProject_Add(
  NautyDirectPackage
  PREFIX "${CMAKE_BINARY_DIR}/external"
  GIT_REPOSITORY https://github.com/MathieuDutSik/nauty
  GIT_TAG 2.8.9
  UPDATE_DISCONNECTED TRUE
  CONFIGURE_COMMAND ./configure --prefix=${NAUTY_INSTALL_DIR}
  BUILD_IN_SOURCE 1
  BUILD_BYPRODUCTS "${NAUTY_LIBRARY}"
)

add_library(basic_common_nauty STATIC IMPORTED GLOBAL)
set_target_properties(
  basic_common_nauty
  PROPERTIES
    IMPORTED_LOCATION "${NAUTY_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NAUTY_INCLUDE_DIR}"
)
add_dependencies(basic_common_nauty NautyDirectPackage)

set(BLISS_PREFIX "${CMAKE_BINARY_DIR}/external/bliss")
set(BLISS_INSTALL_DIR "${BLISS_PREFIX}/install")
set(BLISS_INCLUDE_DIR "${BLISS_INSTALL_DIR}/include")
set(BLISS_LIBRARY "${BLISS_INSTALL_DIR}/lib/libbliss.a")

ExternalProject_Add(
  BlissDirectPackage
  PREFIX "${BLISS_PREFIX}"
  URL "${CMAKE_SOURCE_DIR}/ExternalLib/bliss-0.77.zip"
  CMAKE_ARGS
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E make_directory ${BLISS_INCLUDE_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/src ${BLISS_INCLUDE_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${BLISS_INSTALL_DIR}/lib
    COMMAND ${CMAKE_COMMAND} -E copy <BINARY_DIR>/libbliss_static.a ${BLISS_LIBRARY}
  BUILD_BYPRODUCTS "${BLISS_LIBRARY}"
)

add_library(basic_common_bliss STATIC IMPORTED GLOBAL)
set_target_properties(
  basic_common_bliss
  PROPERTIES
    IMPORTED_LOCATION "${BLISS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${BLISS_INCLUDE_DIR}"
)
add_dependencies(basic_common_bliss BlissDirectPackage)

add_library(basic_common INTERFACE)
target_include_directories(
  basic_common
  INTERFACE
    "${CMAKE_SOURCE_DIR}/src_basic"
    "${CMAKE_SOURCE_DIR}/src_number"
    "${CMAKE_SOURCE_DIR}/src_matrix"
    "${CMAKE_SOURCE_DIR}/src_comb"
    "${CMAKE_SOURCE_DIR}/src_graph"
    "${CMAKE_SOURCE_DIR}/sparse-map/include/tsl"
    "${CMAKE_SOURCE_DIR}/robin-map/include/tsl"
    "${CMAKE_SOURCE_DIR}/hopscotch-map/include/tsl"
)
target_link_libraries(
  basic_common
  INTERFACE
    Eigen3::Eigen
    Boost::serialization
    Threads::Threads
    PkgConfig::GMPXX
    "${GMP_LIBRARY}"
    basic_common_flint
)

add_library(basic_common_graph_support INTERFACE)
target_link_libraries(
  basic_common_graph_support
  INTERFACE
    basic_common
    basic_common_nauty
    basic_common_bliss
)

function(add_basic_common_binary target source support_target)
  add_executable(${target} ${source} ${ARGN})
  target_link_libraries(${target} PRIVATE ${support_target})
endfunction()

set(
  BASIC_PROGRAMS
  Convert_Heuristic_to_Thompson_sampling
  Test_Thompson_sampling
  Fuzzing_Heuristic_and_Thompson_sampling
  Test_timing
  Test_PresenceProgram
  Test_namelist
  Test_String_Split
  Test_serialization
  Test_EmpiricalDistribution
  Analysis_EmpiricalDistribution
  Test_hash_table
  Test_bool
  Test_face
  Test_String_Optimization
  Test_String_Trim
  Test_String_Conversion
  Test_Namelist_Substr
  Test_Basic_string_Optimizations
  Test_filesystem
  Test_program_in_path
  Test_external_program
  Test_Heuristic
  Test_parse_split
  Read_FileNumber
  Read_FileBool
  Read_FileFace
)
foreach(target IN LISTS BASIC_PROGRAMS)
  add_basic_common_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_basic/${target}.cpp"
    basic_common
  )
endforeach()

set(
  NUMBER_PROGRAMS
  Test_SquareRoot
  Test_UnorderedMapMpzq
  Factorize
  Factorize_help
  Test_Factorize
  Test_PracticalInf
  Test_QuoInt
  Test_Operations
  Test_PrintStaticInfo
  Test_RealCubicField
  Test_QuadField
  Test_Rational
  Test_Format
  Test_TypeBoostGmp
  Test_GetBit
  Test_QuadraticResidue
  Test_PrimeGenerator
  Test_SequenceApproximant
  ReadValueString
  ComputeQuadraticResidue
  Test_ComputePairGcdDot
  Test_Conversion
  Test_Conversion8_16_32_64
  Test_Padic
  Test_FactorizeHelp
  Test_FloorCeilScalarInteger
  SequenceApproximant
  TEST_jet_number
  Test_Fp
)
foreach(target IN LISTS NUMBER_PROGRAMS)
  add_basic_common_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_number/${target}.cpp"
    basic_common
  )
endforeach()

set(
  MATRIX_PROGRAMS
  AnalysisSingMat
  NullspaceComputation
  ZbasisComputation
  RowHermiteNormalForm
  SerializeMatrix
  Test_MatrixInverse
  Test_PerformanceHNF
  ColHermiteNormalForm
  NullspaceIntTrMat
  NullspaceIntMat
  SolutionIntMat
  IntersectionLattice
  Test_HilbertMatrix
  Test_ProdMat
  Matrix_As_Set
  Test_NullspaceComputation
  Test_NullspaceTrMat
  ScalarCanonicalization
  MatrixSubsetSolver
  DeterminantMat
  DeterminantMatMod
  DeterminantMatHadamard
  ComputePointer
  NonUniqueScaleToInteger
  NullspaceComputationExperiment
  SelectRowColMatModHadamard
  SmithNormalForm
  Test_SubspaceCompletion
  RemoveFractionMatrix
  SolutionMatRepetitive
  Test_FindIsotropicMod
  Test_MatrixMod
  Test_EquaSpace
  RankMat
  Matrix_As_SetPosScal
  NullspaceMatMod
  IsotropicMatMod
  RescaleInverse
  SmithNormalFormInvariant
  ComputeChromaticLowerBounds
  Test_ReadListMat
  ComputeTranslationClasses
  Test_eigen_all_types
  Test_AdjugateDeterminant
  Test_HermiteSmith
  Test_SolutionNullspaceInt
  Test_ScalingCanonic
  Test_RankMat
  Test_SubsetSolver
  Test_SaturationCompletion
)
foreach(target IN LISTS MATRIX_PROGRAMS)
  add_basic_common_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_matrix/${target}.cpp"
    basic_common
  )
endforeach()

set(
  COMB_PROGRAMS
  EnumerationCliques_myrvold
  EnumerationCliques_ordered_enum
  EnumerationCliques_special_subcliques
  Block_CppIterator
  Timing_HashMap
  Vector_Face_iterator
  SerializeVectface
  SerializeFace
  IterFace
  Test_SetCppIterator
)
foreach(target IN LISTS COMB_PROGRAMS)
  add_basic_common_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_comb/${target}.cpp"
    basic_common
  )
endforeach()

set(
  WASM_PROGRAMS
  Test_wasm_cpp_int
  Test_wasm_matrix
  Test_wasm_ordered_containers
)
foreach(target IN LISTS WASM_PROGRAMS)
  add_basic_common_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_wasm/${target}.cpp"
    basic_common
  )
endforeach()

set(
  GRAPH_PROGRAMS
  GRAPH_EnumerateDistanceSembedding
  GRAPH_EnumerateShortCycles
  GRAPH_TestCanonicComputation
  GRAPH_CanonicComputation
  GRAPH_ComputeAutGroup
  TestPadic
)
foreach(target IN LISTS GRAPH_PROGRAMS)
  add_basic_common_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_graph/${target}.cpp"
    basic_common_graph_support
  )
endforeach()
