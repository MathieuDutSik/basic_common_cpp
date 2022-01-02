cmake_minimum_required(VERSION 3.9)

project(bliss VERSION 0.73 DESCRIPTION "bliss library")


add_library(bliss SHARED
    bliss.cc
    heap.cc
    uintseqhash.cc
    bliss_C.cc
    orbit.cc
    utils.cc
    defs.cc
    partition.cc
    graph.cc
    timer.cc
)

set_target_properties(bliss PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(bliss PROPERTIES PUBLIC_HEADER defs.hh graph.hh partition.hh timer.hh utils.hh)

target_compile_features(bliss PRIVATE cxx_std_17)

# Install library
install(TARGETS bliss DESTINATION lib)

# Install library headers
file(GLOB HEADERS *.hh)
install(FILES ${HEADERS} DESTINATION include)


