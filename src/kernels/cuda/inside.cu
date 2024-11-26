#ifdef USE_CUDA
#include <math.h>
#include <cstdio>
#include "kernels/inside.cuh"
#include "utils/math.hpp"
#include "constants.h"
#include "utils/data_encoding.h"

__global__  void kernel_inside_alpha_zerolization(double* alpha, int N, int MS){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_size = N * MS * MS;
    for (int i = idx; i < total_size; i += blockDim.x * gridDim.x) {
        alpha[i] = -INFINITY;
    }
}

__global__ void kernel_inside_base_fill_alpha(  
        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        , uint32_t* symbol_A_vector
        #ifndef USE_CUDA
        , pcfg* grammar
        #endif
) {
        int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
        int total_threads = blockDim.x * gridDim.x;
        bool changed = false;
        uint32_t nonterminate_remains[MAX_NONTERMINATES]{false};
        __shared__ int size;

        size = 0; // nonterminate_remains.size();
        if(thread_id == 0){
            /* thread i process sym_A = i, i + total_thread, i + 2 * total_thread...*/
            for(int sym_A = thread_id; sym_A < N; sym_A += total_threads){
                for(int i = 0; i < sequence_length; i++){
                    double p = -INFINITY;
                    uint64_t key = encode_key(sym_A, sequence[i]);
                    p = reverse_grammar_hashtable_get_value(pretermination_lookuptable, n_grammars * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS, key);
                    
                    if(p == -INFINITY) continue;
                            
                    ALPHA(sym_A, i, i) = p;
                }
            }
                        
            for(int sym_id = 0; sym_id < N; sym_id++){
                nonterminate_remains[sym_id] = true;
            }

            for(int _i = 0; _i < N; _i++) size += (nonterminate_remains[_i] ? 1 : 0);
        }
        

        __syncthreads();

        for(int _i = 0; _i < N; _i++){
            int _size = size;
            for(int i = thread_id; i < sequence_length; i += total_threads){
                    for(int gid = 0; gid < n_grammars; gid++){
                        uint32_t sym_A = symbol_A_vector[gid];
                        uint32_t sym_BC = *(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                        uint32_t sym_B = (sym_BC >> 16) & 0xFFFF;
                        uint32_t sym_C = (sym_BC) & 0xFFFF;
                        double possibility = *(double*)(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1);
                        
                        if(!IS_EPSILON(sym_C) || IS_TERMINATE(sym_B)) continue;
                        double alpha_B = ALPHA(sym_B, i, i);
                        
                        if(alpha_B + possibility > ALPHA(sym_A, i, i)){
                            ALPHA(sym_A, i, i) = alpha_B + possibility;
                            atomicExch(&nonterminate_remains[sym_A], false);
                            atomicSub(&size, 1);
                        }                    
                }
            }
            __syncthreads();
            if(_size == size) break;
        }
}

__global__ void kernel_inside_computeSpanKernel(
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
) {

        int thread_id_x = blockIdx.x * blockDim.x + threadIdx.x; // Thread ID for i
        int total_threads_x = blockDim.x * gridDim.x;
        double buffer[MAX_NONTERMINATES * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]; // record symbol. span length and i

        if (thread_id_x == 0) {
            for (int i = 0; i < N * MS * MS; ++i) {
                buffer[i] = -INFINITY;
            }
        }

        __syncthreads();

        for (int span_length = 2; span_length <= sequence_length; span_length++) {
                for (int i = thread_id_x; i <= sequence_length - span_length; i += total_threads_x) {
                    // the buffer need be access with lock, if parallel gid. we choose not parallel this axis.
                    for(int gid = 0; gid < n_grammars; gid ++){
                        int j = i + span_length - 1; // Ending index of the span
                        uint32_t sym_A = symbol_A_vector[gid];
                        uint32_t sym_BC = *(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                        uint32_t sym_B = (sym_BC >> 16) & 0xFFFF;
                        uint32_t sym_C = (sym_BC) & 0xFFFF;
                        double possibility = *(double*)(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1);

                        for (int k = i; k < j; k++) {
                            if(IS_EPSILON(sym_C)) continue;

                            // A->BC
                            double alpha_B = ALPHA_GET(sym_B, i, k);
                            double alpha_C = ALPHA_GET(sym_C, k + 1, j);
                            LOG_SUM_EXP_SET(buffer[sym_A * MS * MS + span_length * MS + i], alpha_B + alpha_C + possibility);
                            
                        }
                    }

                    __syncthreads();

                    // A->B
                    for(int idx_rule_id = thread_id_x; idx_rule_id < inside_order_1_rule_iteration_path_size; idx_rule_id += total_threads_x) {
                        int j = i + span_length - 1; // Ending index of the span
                        uint32_t gid = *(inside_order_1_rule_iteration_path + idx_rule_id * 2);
                        uint32_t sym_A = *(inside_order_1_rule_iteration_path + idx_rule_id * 2 + 1);
                        uint32_t* addr = (grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                        uint32_t sym_B = ((*addr) >> 16) & 0xFFFF;
                        double possibility = *(double*)(addr + 1);                        
                        double alpha_B = POSSIBILITY_ADD(buffer[sym_B * sequence_length * sequence_length + span_length * sequence_length + i], ALPHA_GET(sym_B, i, j)); // Cell [i, j] is newly calculated and has not been written back. We need access the buffer, rather than ALPHA(..)

                        LOG_SUM_EXP_SET(buffer[sym_A * MS * MS + span_length * MS + i], alpha_B + possibility);
                    }
                    
                    __syncthreads();

                    // write back.
                    for (int sym_A = thread_id_x; sym_A < N; sym_A += total_threads_x) {
                        int j = i + span_length - 1; // Ending index of the span
                        printf("set %d logsumexp ", ALPHA(sym_A, i, j), buffer[sym_A * MS * MS + span_length * MS + i]);
                         LOG_SUM_EXP_SET(ALPHA(sym_A, i, j), buffer[sym_A * MS * MS + span_length * MS + i]);  
                    }

                }
                __syncthreads();
        }
}
#endif