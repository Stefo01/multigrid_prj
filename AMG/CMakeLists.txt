cmake_minimum_required(VERSION 3.10)
project(AMG_Project)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

# Source files from src/
file(GLOB_RECURSE SRC_FILES
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

# Add main file manually
set(MAIN_FILE ${PROJECT_SOURCE_DIR}/src/main.cpp)

# Add executable target
add_executable(AMGExecutable ${SRC_FILES} ${MAIN_FILE})