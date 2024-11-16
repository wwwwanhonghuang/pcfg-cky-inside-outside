
# Set USE_CUDA flag (0 or 1)
USE_CUDA := 0

# Include directories
INCLUDES = -Iinclude

# Compiler flags and libraries
TEST_LDFLAGS = -L/path/to/boost/lib -lboost_unit_test_framework
CXX = clang++ -DCOMPUTING_IN_LOG_SPACE
CUDA_COMPILER := nvcc

TARGET_COMPILER := $(CXX)

ifeq ($(USE_CUDA), 1)
COMPILER_FLAGS := -v -std=c++11 -arch compute_86 -code sm_86 -O3 -g -DUSE_CUDA -DCOMPUTING_IN_LOG_SPACE
else
COMPILER_FLAGS := -std=c++17 -O3 -g -fopenmp -Werror -DDEBUG_INSIDE_ALGORITHM -DCOMPUTING_IN_LOG_SPACE
endif



# Directories
SRC_DIR := src
INCLUDE_DIR := include
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin
LIB_DIR := $(BUILD_DIR)/lib
UNIT_TEST_DIR := unit_test

# Source files and object files

# If CUDA is enabled, filter out main.cpp and other unwanted files
ifneq ($(USE_CUDA), 0)
    SRC_FILES := $(wildcard $(SRC_DIR)/**/*.cpp $(SRC_DIR)/*.cpp)
    # Include both CUDA source files and regular source files
    CUDA_SRC_FILES := $(wildcard  $(SRC_DIR)/**/**/*.cu  $(SRC_DIR)/**/*.cu) $(wildcard $(SRC_DIR)/*.cu)
    SRC_FILES := $(SRC_FILES) $(CUDA_SRC_FILES)

    # Exclude unwanted .cpp files when CUDA is enabled
    SRC_FILES := $(filter-out $(SRC_DIR)/kernels/*.cpp, $(SRC_FILES))
    SRC_FILES := $(filter-out $(SRC_DIR)/algorithms/alg_inside_outside_main.cpp, $(SRC_FILES))
    SRC_FILES := $(filter-out $(SRC_DIR)/main.cpp, $(SRC_FILES))
    SRC_FILES := $(filter-out $(SRC_DIR)/statistics/*.cpp, $(SRC_FILES))
else
# Exclude unwanted .cpp files when CUDA is disabled
    SRC_FILES := $(wildcard $(SRC_DIR)/**/*.cpp $(SRC_DIR)/*.cpp)
    SRC_FILES := $(filter-out $(SRC_DIR)/kernels/cuda/*.cu, $(SRC_FILES))
endif

# Object files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(OBJ_FILES))

# Test files
TEST_SRCS := $(wildcard $(UNIT_TEST_DIR)/*.cpp)
TEST_OBJS := $(patsubst $(UNIT_TEST_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(TEST_SRCS))

# Executable names
TARGET_MAIN = $(BIN_DIR)/main_executable
TARGET_CUDA_MAIN = $(BIN_DIR)/main_cuda_executable
TARGET_PHASE_CONVERT = $(BIN_DIR)/phase_convert_executable
TEST_EXE = $(BIN_DIR)/unit_tests
SHARED_LIB = $(LIB_DIR)/libshared.a

# Targets
.PHONY: all main cuda_main phase_convert run_tests clean

# Default target to build main executable
all: main

# Build shared library
$(SHARED_LIB): $(OBJ_FILES) | $(LIB_DIR)
	$(AR) rcs $@ $^


# Build main executable
main: $(TARGET_MAIN)
cuda_main: $(TARGET_CUDA_MAIN)

TARGET_MAIN_LD = -lyaml-cpp
TARGET_CUDA_MAIN_LD =  -lyaml-cpp -lcuda -lcudart
TARGET_PHASE_CONVERT_LD = -lyaml-cpp 
TARGET_MAIN_FLAGS = 

$(TARGET_MAIN): $(SRC_DIR)/main.cpp $(SHARED_LIB) | $(BIN_DIR)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(TARGET_MAIN_FLAGS) $(INCLUDES) -o $@ $(SRC_DIR)/main.cpp $(SHARED_LIB) $(TARGET_MAIN_LD) -latomic

$(TARGET_CUDA_MAIN): $(SRC_DIR)/main.cu $(SHARED_LIB) | $(BIN_DIR)
	$(CUDA_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -o $@ $(SRC_DIR)/main.cu $(SHARED_LIB) $(TARGET_CUDA_MAIN_LD) 

# Build phase_convert executable
phase_convert: $(TARGET_PHASE_CONVERT)

$(TARGET_PHASE_CONVERT): $(SRC_DIR)/phase_convert.cpp $(SHARED_LIB) | $(BIN_DIR)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -o $@ $(SRC_DIR)/phase_convert.cpp $(SHARED_LIB) $(TARGET_PHASE_CONVERT_LD)

# Build unit tests executable
$(TEST_EXE): $(OBJ_FILES) $(TEST_OBJS) | $(BIN_DIR)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) $(TEST_LDFLAGS) -o $@ $(OBJ_FILES) $(TEST_OBJS)

# Compile source files to object files

ifeq ($(USE_CUDA), 1)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(@D)
	$(CUDA_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -c $< -o $@
else
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(@D)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -c $< -o $@
endif

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu | $(OBJ_DIR)
	@mkdir -p $(@D)
	$(CUDA_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(UNIT_TEST_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(@D)
	$(TARGET_COMPILER) $(COMPILER_FLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR):
	@mkdir -p $@

$(BIN_DIR):
	@mkdir -p $@

$(LIB_DIR):
	@mkdir -p $@

# Run tests
run_tests: $(TEST_EXE)
	./$(TEST_EXE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
