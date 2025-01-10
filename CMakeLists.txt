# Set the minimum version of CMake
cmake_minimum_required(VERSION 3.18)
project(MyProject)

# Set the C++ standard and enable optional CUDA support
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

# Important Note: currently CUDA is unsupport
option(USE_CUDA "Enable CUDA support" OFF)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Include directories
include_directories(include)

# Find Boost (for unit testing)
# find_package(Boost REQUIRED COMPONENTS unit_test_framework)
# include_directories(${Boost_INCLUDE_DIRS})

# Add definitions based on options
add_definitions(-DCOMPUTING_IN_LOG_SPACE)

# Compiler flags for non-CUDA builds
if(NOT USE_CUDA)
    add_compile_options(-O3 -g -fopenmp -Werror -DDEBUG_INSIDE_ALGORITHM)
endif()

# Set CUDA flags if enabled
if(USE_CUDA)
    enable_language(CUDA)
    set(CUDA_SEPARABLE_COMPILATION ON)
    add_compile_options($<$<COMPILE_LANGUAGE:CUDA>: -arch=compute_86 -code=sm_86 -DUSE_CUDA>)
endif()

# Source Files - adjust filtering as necessary
file(GLOB_RECURSE SRC_FILES
    "src/*.cpp"
    "src/*.cu"
)
# Exclude specific files if CUDA is enabled
if(USE_CUDA)
    list(FILTER SRC_FILES EXCLUDE REGEX "src/kernels/.*\\.cpp")
    list(FILTER SRC_FILES EXCLUDE REGEX "src/algorithms/alg_inside_outside_main.cpp")
    list(FILTER SRC_FILES EXCLUDE REGEX "src/main.cpp")
    list(FILTER SRC_FILES EXCLUDE REGEX "src/statistics/.*\\.cpp")
endif()

# Further filtering, regardless of CUDA
list(FILTER SRC_FILES EXCLUDE REGEX "src/syntax_analysis.cpp")

# Object files library
add_library(shared_lib STATIC ${SRC_FILES})
# target_link_libraries(shared_lib ${Boost_LIBRARIES})

# Train PCFG executable targets
add_executable(train_pcfg src/main.cpp)
target_link_libraries(train_pcfg shared_lib yaml-cpp)

# Phase conversion executable
add_executable(phase_convert_executable src/phase_convert.cpp)
target_link_libraries(phase_convert_executable shared_lib yaml-cpp)

# Syntax analysis executable
add_executable(syntax_analysis src/syntax_analysis.cpp)
target_link_libraries(syntax_analysis shared_lib yaml-cpp)

add_executable(test_parsing src/test/test_parsing.cpp)
target_link_libraries(test_parsing shared_lib yaml-cpp)

# Conditional CUDA main executable
if(USE_CUDA)
    add_executable(main_cuda_executable src/main.cu)
    target_link_libraries(main_cuda_executable shared_lib yaml-cpp cuda cudart)
endif()

# Test executable target
file(GLOB TEST_SRCS "unit_test/*.cpp")
add_executable(unit_tests ${TEST_SRCS})
# target_link_libraries(unit_tests shared_lib ${Boost_LIBRARIES})

# Custom target for running tests
add_custom_target(run_tests
    COMMAND unit_tests
    DEPENDS unit_tests
    COMMENT "Running unit tests..."
)

# Output directories
set_target_properties(train_pcfg phase_convert_executable syntax_analysis test_parsing PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)


# Note: currently CUDA is unsupport
if(USE_CUDA)
    set_target_properties(main_cuda_executable PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
endif()

set_target_properties(shared_lib PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
