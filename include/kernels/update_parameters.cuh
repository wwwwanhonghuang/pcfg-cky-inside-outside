#ifndef CUH_UPDATE_PARAMETERS
#define CUH_UPDATE_PARAMETERS
#include<omp.h>
#include "utils/data_structure.hpp"
#ifndef USE_CUDA
#include <cstring>
#endif
#ifdef USE_CUDA
__global__
#endif

#include <cstdint>
#include "grammar/grammar.hpp"

const float epsilon = 1e-12f;
inline float _calculate_new_possibility(float S, float f_gid);

void kernel_update_parameters(double* f, float* count, float* mu, float* beta, uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
    #ifdef USE_CUDA
        uint32_t* 
    #else
        common_32bit*
    #endif
        grammar_table, float* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            , pcfg* grammar
        #endif
);
#endif