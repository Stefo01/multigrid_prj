cmake_minimum_required(VERSION 3.12.0)
project(GeometricMultiGrid LANGUAGES CXX C)

include(../common/cmake-common.cmake)
include_directories(include)

file(GLOB source_files
  "src/*.cpp"
)
add_executable(Multigrid ${source_files})



