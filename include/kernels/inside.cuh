#ifndef CUH_INSIDE
#define CUH_INSIDE

#include <cstring>
#include <tuple>

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <unordered_set>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif


#include "utils/data_structure.hpp"
#include "macros.def"
#ifdef DEBUG_INSIDE_ALGORITHM
#include "grammar/grammar.hpp"
#include "utils/data_encoding.h"
#include "utils/data_accessing.hpp"
#endif

extern "C" {
#ifdef USE_CUDA
__global__ 
#endif
void kernel_inside_alpha_zerolization(double* alpha, int N, int MS);

#ifdef USE_CUDA
__global__ 
#endif
void kernel_inside_base_fill_alpha(  
        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        , uint32_t* symbol_A_vector
        #ifndef USE_CUDA
        , pcfg* grammar
        #endif
);


#ifdef USE_CUDA
__global__ 
#endif 
void kernel_inside_computeSpanKernel(const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
        #ifndef USE_CUDA
                std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
        #else
                uint32_t* inside_order_1_rule_iteration_path, uint32_t inside_order_1_rule_iteration_path_size
        #endif
        , uint32_t* symbol_A_vector
        #ifndef USE_CUDA
        , pcfg* grammar
        #endif
);
}
#endif