cmake_minimum_required(VERSION 3.9)

project(nauty VERSION 1.0.0 DESCRIPTION "nauty library")

add_library(nauty SHARED
    gtools.c
    naurng.c
    nautinv.c
    nautil.c
    nauty.c
    naututil.c
    nausparse.c
    naugraph.c
    schreier.c
    traces.c
)

set_target_properties(nauty PROPERTIES VERSION ${PROJECT_VERSION})
#add_definitions(-g)

set(CMAKE_C_FLAGS " -O3")
#set(CMAKE_C_FLAGS " -O3 -g")
#set(CMAKE_C_FLAGS " -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -fsanitize=leak -Wstrict-aliasing")


# Install library
install(TARGETS nauty DESTINATION lib)

# Install library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include)


