cmake_minimum_required(VERSION 3.30)

project(basic_common_cpp LANGUAGES CXX)

include(ExternalProject)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS serialization)

pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
find_library(GMP_LIBRARY NAMES gmp REQUIRED)

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
  Test_PrintStaticInfo
  Test_RealCubicField
  Test_QuadField
  Test_Rational
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
  SequenceApproximant
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
  SmithNormalFormIntegerMat
  Test_SubspaceCompletion
  RemoveFractionMatrix
  SolutionMatRepetitive
  Test_FindIsotropicMod
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
