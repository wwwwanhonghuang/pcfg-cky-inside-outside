#ifndef __CUH_VALIDATION
#define __CUH_VALIDATION
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
#include <cstdint>

#ifdef USE_CUDA
__global__
#endif
void validation_at_device(double *log_likelihoods, uint32_t valid_set_size, uint32_t* pretermination_lookuptable_device, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                    int sequence_length, int n_syms, int N, int T, int MS, int cnt_grammar,
                    uint32_t* inside_order_1_rule_iteration_path, int inside_order_1_rule_iteration_path_size);

#endif
