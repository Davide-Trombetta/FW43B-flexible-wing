# stlInfo -- geometric summary of an STL surface (plain C++17, no OpenFOAM
# library: the same source builds with OpenFOAM.org and ESI). Built by
# Allrun.pre into tools/bin/stlInfo.
cmake_minimum_required(VERSION 3.12)
project(stlInfo CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(stlInfo stlInfo.cpp)
target_compile_options(stlInfo PRIVATE -Wall -Wextra)
set_target_properties(stlInfo PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin")
