# Set minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(MyCUDAProject VERSION 1.0)

# Find CUDA package
find_package(CUDA REQUIRED)

set(CUDA_ARCHITECTURE "sm_86")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Specify the source directory
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

# Define source files
set(SRC_FILES ${SOURCE_DIR}/main.cu)

# Include CUDA directories
include_directories(${CUDA_INCLUDE_DIRS})

# Add executable for host and device code
cuda_add_executable(cuda_main ${SRC_FILES})

# Link CUDA libraries (if necessary)
target_link_libraries(cuda_main ${CUDA_LIBRARIES})
