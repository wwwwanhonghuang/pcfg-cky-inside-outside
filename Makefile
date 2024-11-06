# Include directories
INCLUDES = -Iinclude

# Compiler flags and libraries
TEST_LDFLAGS = -L/path/to/boost/lib -lboost_unit_test_framework

CXX = clang++ 

TARGET_COMPILER := $(CXX)
COMPILER_FLAGS := -std=c++17 -O3 -g -fopenmp -Werror -DDEBUG_INSIDE_ALGORITHM

# Directories
SRC_DIR := src
INCLUDE_DIR := include
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin
LIB_DIR := $(BUILD_DIR)/lib
UNIT_TEST_DIR := unit_test

# Source files and object files
SRC_FILES := $(wildcard $(SRC_DIR)/**/*.cpp $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

# Test files
TEST_SRCS := $(wildcard $(UNIT_TEST_DIR)/*.cpp)
TEST_OBJS := $(patsubst $(UNIT_TEST_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(TEST_SRCS))

# Executable names
TARGET_MAIN = $(BIN_DIR)/main_executable
TARGET_PHASE_CONVERT = $(BIN_DIR)/phase_convert_executable
TEST_EXE = $(BIN_DIR)/unit_tests
SHARED_LIB = $(LIB_DIR)/libshared.a

# Targets
.PHONY: all main phase_convert run_tests clean

# Default target to build main executable
all: main

# Build shared library
$(SHARED_LIB): $(OBJ_FILES) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(LIB_DIR):
	@mkdir -p $@

# Build main executable
main: $(TARGET_MAIN)

TARGET_MAIN_LD = -lyaml-cpp $(CUDA_LIBS) $(CUTENSOR_LIBS)
TARGET_PHASE_CONVERT_LD = -lyaml-cpp 
TARGET_MAIN_FLAGS = 

$(TARGET_MAIN): $(SRC_DIR)/main.cpp $(SHARED_LIB) | $(BIN_DIR)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(TARGET_MAIN_FLAGS) $(INCLUDES) -o $@ $(SRC_DIR)/main.cpp $(SHARED_LIB) $(TARGET_MAIN_LD) -latomic

# Build phase_convert executable
phase_convert: $(TARGET_PHASE_CONVERT)

$(TARGET_PHASE_CONVERT): $(SRC_DIR)/phase_convert.cpp $(SHARED_LIB) | $(BIN_DIR)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -o $@ $(SRC_DIR)/phase_convert.cpp $(SHARED_LIB) $(TARGET_PHASE_CONVERT_LD)

# Build unit tests executable
$(TEST_EXE): $(OBJ_FILES) $(TEST_OBJS) | $(BIN_DIR)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) $(TEST_LDFLAGS) -o $@ $(OBJ_FILES) $(TEST_OBJS)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(@D)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(UNIT_TEST_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(@D)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR):
	@mkdir -p $@

$(BIN_DIR):
	@mkdir -p $@

# Run tests
run_tests: $(TEST_EXE)
	./$(TEST_EXE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
