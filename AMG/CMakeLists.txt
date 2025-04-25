cmake_minimum_required(VERSION 3.10)
project(AMG_Project)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Explicitly set the main file
set(MAIN_SRC ${PROJECT_SOURCE_DIR}/src/main.cpp)

# Get all other .cpp files except main.cpp
file(GLOB_RECURSE OTHER_SOURCES
    ${PROJECT_SOURCE_DIR}/src/*.cpp
)

# Remove main.cpp from OTHER_SOURCES
list(REMOVE_ITEM OTHER_SOURCES ${MAIN_SRC})

# Combine sources
set(ALL_SOURCES ${MAIN_SRC} ${OTHER_SOURCES})

# Create executable
add_executable(${PROJECT_NAME} ${ALL_SOURCES})