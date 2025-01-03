#ifndef CUH_UPDATE_PARAMETERS
#define CUH_UPDATE_PARAMETERS
#include<omp.h>
#include "utils/data_structure.hpp"
#include <cstring>
#include <cstdint>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#include "grammar/grammar.hpp"

extern "C" {
#ifndef USE_CUDA
inline 
#endif
#ifdef USE_CUDA
__device__
#endif
double _calculate_new_possibility(double S, double f_gid);


#ifdef USE_CUDA
__global__ 
#endif
void kernel_update_parameters(double* f, double* count, double* mu, double* beta, const uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
        uint32_t*
        grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            , pcfg* grammar
        #endif
        , bool do_update
);
}
#endif