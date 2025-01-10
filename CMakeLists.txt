# Set the minimum version of CMake
cmake_minimum_required(VERSION 3.18)
project(MyProject)

# Set the C++ standard and enable optional CUDA support
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_C_COMPILER clang)
# set(CMAKE_CXX_COMPILER clang++)

# Important Note: currently CUDA is unsupport
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Include directories
include_directories(include)
# Add definitions based on options
add_definitions(-DCOMPUTING_IN_LOG_SPACE)

add_compile_options(-O3 -g -fopenmp -Werror -DDEBUG_INSIDE_ALGORITHM)

file(GLOB_RECURSE SRC_FILES
    "src/*.cpp"
)

list(FILTER SRC_FILES EXCLUDE REGEX "src/syntax_analysis.cpp")
list(FILTER SRC_FILES EXCLUDE REGEX "src/distribution/.*\\.cpp")

# Object files library
add_library(shared_lib STATIC ${SRC_FILES})

file(GLOB_RECURSE DISTRIBUTED_TRAINING_SRC_FILES
    "src/distribution/*.cpp"
)

add_library(distribution_lib STATIC ${DISTRIBUTED_TRAINING_SRC_FILES})
target_compile_options(distribution_lib PRIVATE -O0)  # Current share memory only run correctly in -O0
target_compile_options(shared_lib PRIVATE -O3)

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

add_executable(distributed_training_main src/distributed_training_main.cpp)
target_link_libraries(distributed_training_main shared_lib yaml-cpp distribution_lib)

add_executable(distributed_training_moderator src/distributed_training_moderator.cpp)
target_link_libraries(distributed_training_moderator shared_lib yaml-cpp distribution_lib)

target_compile_options(distributed_training_main PRIVATE -O0)
target_compile_options(distributed_training_moderator PRIVATE -O0)
target_compile_options(test_parsing PRIVATE -O3)
target_compile_options(syntax_analysis PRIVATE -O3)
target_compile_options(phase_convert_executable PRIVATE -O3)
target_compile_options(train_pcfg PRIVATE -O3)

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
set_target_properties(train_pcfg phase_convert_executable syntax_analysis test_parsing distributed_training_main distributed_training_moderator PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

set_target_properties(shared_lib PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)