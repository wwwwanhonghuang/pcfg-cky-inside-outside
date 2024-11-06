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

inline double _calculate_new_possibility(double S, double f_gid);

void kernel_update_parameters(double* f, double* count, double* mu, double* beta, const uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
    
    #ifdef USE_CUDA
        uint32_t* 
    #else
        uint32_t*
    #endif
        grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            , pcfg* grammar
        #endif
        , bool do_update
);
#endif