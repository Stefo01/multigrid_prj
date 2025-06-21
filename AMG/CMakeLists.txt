cmake_minimum_required(VERSION 3.12.0)
project(AlgebraicMultiGrid LANGUAGES CXX C)

include(../common/cmake-common.cmake)
include_directories(include)

# Enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP. Adding flags.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

file(GLOB source_files
  "src/*.cpp"
)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# Main executable
add_executable(AMG ${source_files})
target_link_libraries(AMG PRIVATE OpenMP::OpenMP_CXX)

# Test executable
add_executable(test debugtest.cpp src/CSRMatrix.cpp src/Utilities.cpp src/FEM.cpp)
target_link_libraries(test PRIVATE OpenMP::OpenMP_CXX)
