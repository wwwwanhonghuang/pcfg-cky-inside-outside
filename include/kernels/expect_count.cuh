#ifndef CUH_EXPECT_COUNT
#define CUH_EXPECT_COUNT

#include <cstring>
#include <cstdint>
#include <omp.h>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#include "utils/data_accessing.hpp"
#include "grammar/grammar.hpp"
#include "macros.def"

#if DEBUG_MODE == 0
#define EXPERCT_COUNT_DEBUG_OUTPUT(CONDITION)
#else
#define EXPERCT_COUNT_DEBUG_OUTPUT(CONDITION) if(CONDITION > 0) \
                                            std::cout << "(" << i << ", " << j << ")" << \
                                            "expectation_count::" << \
                                            "gid = " << gid << " " << SYMBOL_STR(sym_A) << "->" << SYMBOL_STR(sym_B) \
                                            << " " <<  SYMBOL_STR(sym_C) << \
                                            " expectation count += " <<  mu_val << "(MU[" << gid << "," << i << ", " << j << "]" << std::endl;
#endif

extern "C" {
#ifdef USE_CUDA
__global__
#endif
void kernel_expect_count(double* count, double* mu, double* beta, 
    const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , pcfg* grammar
                        #endif
                        ,uint32_t* symbol_A_vector
);
}
#endif