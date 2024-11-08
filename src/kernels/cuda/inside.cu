#include <math.h>
#include "kernels/inside.cuh"
#include "utils/math.hpp"
#include "constants.h"


__device__ void kernel_inside_alpha_zerolization(double* alpha, int N, int MS){
    memset(alpha, 0, N * MS * MS * sizeof(double));
    for (int i = 0; i < N * MS * MS; i++) {
        alpha[i] = -INFINITY;
    }
}

__device__ void kernel_inside_base_fill_alpha(  
        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
) {
        int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
        int total_thread = blockDim.x * gridDim.x;
        /* thread i process sym_A = i, i + total_thread, i + 2 * total_thread...*/
        for(int sym_A = thread_id; sym_A < N; sym_A += total_thread){
            for(int i = 0; i < sequence_length; i++){
                double p = 
                    #ifdef COMPUTING_IN_LOG_SPACE
                    -INFINITY
                    #else
                    0.0
                    #endif
                ;
                uint64_t key = encode_key(sym_A, sequence[i]);
                p = reverse_grammar_hashtable_get_value(pretermination_lookuptable, n_grammars * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS, key);
                
                #ifdef COMPUTING_IN_LOG_SPACE
                    if(p == -INFINITY) continue;
                #else
                    if(abs(p - 0) < grammar_minimal_possibility) continue;
                #endif           
                ALPHA(sym_A, i, i) = p;
            }
        }
        
        bool changed = false;
        bool nonterminate_remains[N] = {false};
        for(int sym_id = 0; sym_id < N; sym_id++){
            nonterminate_remains[sym_id] = true;
        }

        for(int i = 0; i < N; i++){
            int size = 0; // nonterminate_remains.size();
            for(int _i = 0; _i < N; _i++) size += (nonterminate_remains[_i] ? 1 : 0);
            for(int i = 0; i < sequence_length; i++){
                for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                        PCFGItemIterator(N, grammar_index, grammar_table)){
                    uint32_t sym_A = std::get<0>(item);
                    uint32_t sym_B = std::get<1>(item);
                    uint32_t sym_C = std::get<2>(item);
                    double possibility = std::get<3>(item);
                    
                    if(!IS_EPSILON(sym_C) || IS_TERMINATE(sym_B)) continue;
                    double alpha_B = ALPHA(sym_B, i, i);
                    
                    if(alpha_B + possibility > ALPHA(sym_A, i, i)){
                        ALPHA(sym_A, i, i) = alpha_B + possibility;
                        nonterminate_remains[sym_A] = false;
                    }                    
                }
                
            }
            int size_after = 0;
            for(int _i = 0; _i < N; _i++) size_after += (nonterminate_remains[_i] ? 1 : 0);
            if(size_after == size) break;
        }
}

__device__ void kernel_inside_computeSpanKernel(const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path,
        uint32_t* symbol_A_vector
        ) {

        int thread_id_x = blockIdx.x * blockDim.x + threadIdx.x; // Thread ID for i

        int total_threads_x = blockDim.x * gridDim.x;

        __shared__ double buffer[N * MS * MS]; // record symbol. span length and i
        if (thread_id_x == 0) {
            for (int i = 0; i < N * MS * MS; ++i) {
                buffer[i] = -INFINITY;
            }
        }
        __syncthreads();
        
        for (int span_length = 2; span_length <= sequence_length; span_length++) {
                for (int i = thread_id_x; i <= sequence_length - span_length; thread_id_x += total_threads_x) {
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
                    for(int idx_rule_id = thread_id_x; idx_rule_id < inside_order_1_rule_iteration_path.size(); idx_rule_id += thread_id_x) {
                        std::tuple<uint32_t, uint32_t>& rule_id = inside_order_1_rule_iteration_path[idx_rule_id];
                        int j = i + span_length - 1; // Ending index of the span
                        uint32_t gid = std::get<0>(rule_id);
                        uint32_t sym_A = std::get<1>(rule_id);
                        uint32_t* addr = (grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                        uint32_t sym_B = ((*addr) >> 16) & 0xFFFF;
                        double possibility = *(double*)(addr + 1);                        
                        double alpha_B = POSSIBILITY_ADD(buffer[sym_B * sequence_length * sequence_length + span_length * sequence_length + i], ALPHA_GET(sym_B, i, j)); // Cell [i, j] is newly calculated and has not been written back. We need access the buffer, rather than ALPHA(..)

                        LOG_SUM_EXP_SET(buffer[sym_A * MS * MS + span_length * MS + i], alpha_B + possibility);
                    }
                    
                    __syncthreads();

                    // write back.
                    for (int sym_A = thread_id_x; sym_A < N; sym_A += thread_id_x) {
                        int j = i + span_length - 1; // Ending index of the span
                         LOG_SUM_EXP_SET(ALPHA(sym_A, i, j), buffer[sym_A * MS * MS + span_length * MS + i]);  
                    }

                }
                __syncthreads();
        }
}
