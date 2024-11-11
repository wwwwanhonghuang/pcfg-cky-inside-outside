#ifndef CUH_OUTSIDE
#define CUH_OUTSIDE
#include <omp.h>
#include <cstring>


#include <cstdint>
#include <vector>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif


#include "utils/data_accessing.hpp"
#include "grammar/grammar.hpp"
#include "macros.def"
#define DEBUG_OUTSIDE_CELL(x, y, X) if(i == x && j == y) { X }

extern "C" {
#ifdef USE_CUDA
__global__ 
#endif
void kernel_outside_main(double* mu, double* beta, 
                        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
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