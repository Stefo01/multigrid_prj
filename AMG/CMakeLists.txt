cmake_minimum_required(VERSION 3.12.0)
project(AlgebraicMultiGrid LANGUAGES CXX C)

include(../common/cmake-common.cmake)
include_directories(include)

file(GLOB source_files
  "src/*.cpp"
)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
add_executable(AMG ${source_files})

add_executable(test debugtest.cpp src/CSRMatrix.cpp src/Utilities.cpp)

