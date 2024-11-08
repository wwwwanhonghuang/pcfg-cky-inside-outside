

#include "kernels/expect_count.cuh"
#include "utils/math.hpp"

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


__device__ void lock(volatile int* lockFlag) {
    while (atomicCAS(lockFlag, 0, 1) != 0) {
        // Spin-wait until the lock is acquired
    }
}

__device__ void unlock(volatile int* lockFlag) {
    atomicExch(lockFlag, 0);  // Release the lock
}


__device__ void kernel_expect_count(double* count, double* mu, double* beta, const uint32_t* sequence, 
                        uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
    ){
    
    int thread_id_x = blockIdx.x * blockDim.x + threadIdx.x;
    int thread_id_y = blockIdx.y * blockDim.y + threadIdx.y;

    int total_threads_x = blockDim.x * gridDim.x;
    int total_threads_y = blockDim.y * gridDim.y;


    __shared__ int lockFlag;  // Shared lock flag
    __shared__ double local_buffer_count[n_grammars];
    
    __shared__ double Z = ALPHA(0, 0, sequence_length - 1); 

    if(thread_id_x == 0 && thread_id_y == 0){
        memset(count, 0, n_grammars * sizeof(double));

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
            
            for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                        PCFGItemIterator(N, grammar_index, grammar_table)){
                uint32_t sym_A = std::get<0>(item);
                uint32_t sym_B = std::get<1>(item);
                uint32_t sym_C = std::get<2>(item);
                double possibility = std::get<3>(item);
                uint32_t gid = std::get<4>(item);
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
