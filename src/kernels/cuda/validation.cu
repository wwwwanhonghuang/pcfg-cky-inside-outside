#ifdef USE_CUDA

#include "kernels/validation.cuh"
__global__
void validation_at_device(double *log_likelihoods, uint32_t valid_set_size, uint32_t* pretermination_lookuptable_device, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                    int sequence_length, int n_syms, int N, int T, int MS, int cnt_grammar,
                    uint32_t* inside_order_1_rule_iteration_path, int inside_order_1_rule_iteration_path_size){
    
}   
#endif