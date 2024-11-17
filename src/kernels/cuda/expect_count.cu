#ifdef USE_CUDA
#include "kernels/expect_count.cuh"
#include "utils/math.hpp"

#if DEBUG_MODE == 0
#define EXPERCT_COUNT_DEBUG_OUTPUT(CONDITION)
#else
#define EXPERCT_COUNT_DEBUG_OUTPUT(CONDITION) if(CONDITION) \
                                            std::cout << "(" << i << ", " << j << ")" << \
                                            "expectation_count::" << \
                                            "gid = " << gid << " " << SYMBOL_STR(sym_A) << "->" << SYMBOL_STR(sym_B) \
                                            << " " <<  SYMBOL_STR(sym_C) << \
                                            " expectation count += " <<  mu_val << "(MU[" << gid << "," << i << ", " << j << "]" << std::endl;
#endif


__device__ void lock(int* lockFlag) {
    while (atomicCAS(lockFlag, 0, 1) != 0) {
        // Spin-wait until the lock is acquired
    }
}

__device__ void unlock(int* lockFlag) {
    atomicExch(lockFlag, 0);  // Release the lock
}


__global__ void kernel_expect_count(double* count, double* mu, double* beta, const uint32_t* sequence, 
                        uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars, uint32_t* symbol_A_vector
    ){
    
    int thread_id_x = blockIdx.x * blockDim.x + threadIdx.x;
    int thread_id_y = blockIdx.y * blockDim.y + threadIdx.y;

    int total_threads_x = blockDim.x * gridDim.x;
    int total_threads_y = blockDim.y * gridDim.y;


    __shared__ int lockFlag;  // Shared lock flag
    double local_buffer_count[MAX_SEQUENCE_LENGTH]{0};
    __shared__ double Z;
    Z = ALPHA(0, 0, sequence_length - 1); 

    if(thread_id_x == 0 && thread_id_y == 0){
        for (int i = 0; i < n_grammars; i++) {
            count[i] = -INFINITY;
        }
        /* 0 is the id of S symbol. This expression assign alpha['S', 0, sequence_length - 1] to Z */
        lockFlag = 0;
    }
    
    __syncthreads();
    
    for(int span_length = 1; span_length <= sequence_length; span_length++){
        for(int i = thread_id_x; i < sequence_length - span_length + 1; i += total_threads_x){
            int j = i + span_length - 1;
            
            if(thread_id_x == 0 && thread_id_y == 0){
                for(int _i = 0; _i < n_grammars; _i++){
                    local_buffer_count[_i] = -INFINITY;
                }
            }

            __syncthreads();
            
            for(int gid = 0; gid < n_grammars; gid++){
                uint32_t sym_A = symbol_A_vector[gid];
                uint32_t sym_BC = *(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                uint32_t sym_B = (sym_BC >> 16) & 0xFFFF;
                uint32_t sym_C = (sym_BC) & 0xFFFF;
                double possibility = *(double*)(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1);

                double mu_val = MU(gid, i, j);
                LOG_SUM_EXP_SET(local_buffer_count[gid], mu_val);
            }
            
            __syncthreads();

            for(int gid = thread_id_y; gid < n_grammars; gid += total_threads_y){
                lock(&lockFlag);
                    LOG_SUM_EXP_SET(count[gid], local_buffer_count[gid]);
                unlock(&lockFlag);
            }
            __syncthreads();
        } 
    }
}
#endif