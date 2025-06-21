# Minimum required CMake version
cmake_minimum_required(VERSION 3.10)

# Define the project name
project(Geometric Multigrid)

# Set the C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Ensure strict adherence to the standard

# Define directories
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/bin") # Output directory for executables/libraries

# Automatically find all .cpp source files in the src directory
file(GLOB_RECURSE SRC_FILES "${SOURCE_DIR}/*.cpp")

# Define common compilation flags
add_compile_options(-Wall -Wextra -O3)

# Option to enable OpenMP for parallel compilation
option(BUILD_PARALLEL "Build with OpenMP support for parallel execution" OFF)

if(BUILD_PARALLEL)
    # Find OpenMP and add its flags and libraries
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found. Adding OpenMP compilation flags and linking libraries.")
        add_compile_options(${OpenMP_CXX_FLAGS})
        add_link_options(${OpenMP_CXX_FLAGS}) # Add link options if needed for OpenMP
    else()
        message(FATAL_ERROR "OpenMP was requested but not found.")
    endif()
endif()

# Add the executable target
# The executable name will be 'Multigrid' as per your Makefile
add_executable(${PROJECT_NAME} ${SRC_FILES})

# Link include directories for the target
target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})

# Set the output directory for the executable
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BINARY_DIR})

# If OpenMP is enabled, link the OpenMP libraries to the executable
if(BUILD_PARALLEL AND OpenMP_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

# Define `make remake` behavior (optional, for convenience)
# This creates a custom target that cleans and then builds the default target.
# It's generally not needed as CMake's incremental builds are efficient.
add_custom_target(remake
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Cleans and then rebuilds the project."
)



