# Define paths and libraries

USE_CUDA=0
CUTENSOR_ROOT = ~/libcutensor
CUDA_INSTALL_PATH = /usr/local/cuda
INCLUDES = -I. -I$(CUDA_INSTALL_PATH)/include -I${CUTENSOR_ROOT}/include -I./utils -I./grammar

TEST_LDFLAGS = -L/path/to/boost/lib -lboost_unit_test_framework

FLAGS=-DDEBUG_INSIDE_ALGORITHM
# Define compilers
CXX = g++ -g -fopenmp $(FLAGS)
NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc -g -Xcompiler -fopenmp $(FLAGS)

ifeq ($(USE_CUDA), 1)
    TARGET_COMPILER := $(NVCC)
	CUDA_LIBS = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart
	CUTENSOR_LIBS = -L${CUTENSOR_ROOT}/lib/12/ -lcutensor
else
    TARGET_COMPILER := $(CXX)
    COMPILER_FLAGS := -std=c++11 $(FLAGS)
endif

# Source files
SRC_FILES := ./utils/tensor.cpp ./utils/data_structure.cpp ./utils/printer.cpp ./utils/data_encoding.cpp ./grammar/grammar.cpp ./grammar/grammar_parser.cpp utils/data_accessing.cpp utils/application_io.cpp dataset/dataset_helper.cpp utils/string_helper.cpp
HEADER_FILES := ./utils/tensor.hpp ./utils/data_structure.hpp ./utils/printer.hpp ./utils/data_encoding.h ./grammar/grammar.hpp ./grammar/grammar_parser.hpp utils/data_accessing.hpp utils/application_io.hpp dataset/dataset_helper.hpp utils/string_helper.hpp
MAIN_FILE = main.cpp
OBJ_FILES := $(SRC_FILES:.cpp=.o)
TEST_SRCS = ./unit_test/hashtable_test.cpp # $(wildcard ./unit_test/*.cpp)
# Target executable
TARGET = main
TEST_TARGET = test_executable


# Default target


# Build target
$(TARGET): $(MAIN_FILE) $(OBJ_FILES) $(HEADER_FILES)
	LD_LIBRARY_PATH=${CUTENSOR_ROOT}/lib/12/:${LD_LIBRARY_PATH} $(TARGET_COMPILER) $(COMPILER_FLAGS) $(MAIN_FILE) $(OBJ_FILES) $(INCLUDES) -o $(TARGET) $(CUDA_LIBS) $(CUTENSOR_LIBS)

all: $(TARGET)

# Clean up build files
clean:
	rm -f $(TARGET) $(OBJ_FILES) $(TEST_TARGET)

phase_convert: $(MAIN_FILE) $(OBJ_FILES) $(HEADER_FILES)
	LD_LIBRARY_PATH=${CUTENSOR_ROOT}/lib/12/:${LD_LIBRARY_PATH} $(TARGET_COMPILER) $(COMPILER_FLAGS) phase_convert.cpp $(OBJ_FILES) $(INCLUDES) -o phase_convert $(CUDA_LIBS) $(CUTENSOR_LIBS) 


%.o: %.cpp $(HEADER_FILES)
	$(CXX) $(INCLUDES) -c $< -o $@

%.d: %.cpp
	$(CXX) $(INCLUDES) -M $< > $@

-include $(SRC_FILES:.cpp=.d)