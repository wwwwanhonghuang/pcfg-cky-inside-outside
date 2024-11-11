#ifdef USE_CUDA
#include "kernels/outside.cuh"
#include "utils/math.hpp"

__global__ void kernel_outside_main(double* mu, double* beta, const uint32_t* sequence, 
                        uint32_t* pretermination_lookuptable, 
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
){
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    int total_threads = blockDim.x * gridDim.x;
    double buffer_beta[256 * 512];
    
    if(thread_id == 0){
        memset(mu, 0, n_grammars * MS * MS * sizeof(double));
        memset(beta, 0, n_syms * MS * MS * sizeof(double)); // beta and mu should be initialized by the host.
        for (int i = 0; i < n_grammars * MS * MS; i++) {
            mu[i] = -INFINITY;
        }
        for (int i = 0; i < n_syms * MS * MS; i++) {
            beta[i] = -INFINITY;
        }
        for (int i = 0; i < n_syms * MS; i++){
            buffer_beta[i] = -INFINITY;
        }
        BETA(0, 0, sequence_length - 1) = 1.0;
    }

    __syncthreads();
    
    /* base case: S is the root of the whole sequence with possibility 1.0. */
    /* all cases:
                    1. B->AC      <1>
                    2. B->CA      <1>
                    3. B->A       <2>
                    4. B->C       (x)
                    5. B->w_A C   <3>
                    6. B->w_C A   <1>
                    7. B->A w_C   <1>
                    8. B->C w_A   <3>
                    9. B->w_A w_C <3>
                    10. B->w_Cw_A <3>
                    11. B->w_A    <4>
                    12. B->w_C    (x)
        , where w_A w_C are terminates.
    */


    // Case 1: A is non-terminate 
    // confirm these codes are correct.
    /* diagonal-order iteration */
    for(int span_length = sequence_length; span_length >= 1; span_length--){
        /* for one diagnal of the beta table, all cell can be parallelly computed. */
        if(thread_id == 0){
            for(int _i = 0; _i < n_syms * MS; _i ++){
                buffer_beta[_i] = -INFINITY;
            }
        }
        
        __syncthreads();
        
        for(int i = thread_id; i < sequence_length - span_length + 1; i += total_threads){
            int j = i + span_length - 1;
            
            // 1. 2. 6. 7.
            for(int gid = 0; gid < n_grammars; gid++){
                int j = i + span_length - 1; // Ending index of the span
                uint32_t _sym_A = symbol_A_vector[gid];
                uint32_t sym_BC = *(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                uint32_t _sym_B = (sym_BC >> 16) & 0xFFFF;
                uint32_t _sym_C = (sym_BC) & 0xFFFF;
                double possibility = *(double*)(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1);
        
                for(int k = 0; k < i; k++){       
                    uint32_t sym_B = _sym_A;
                    uint32_t sym_C = _sym_B;
                    uint32_t sym_A = _sym_C;
                    
                    // 1. 7. (A is nonterminate, and C is terminate or nonterminate. A at right side.)
                    if(IS_NONTERMINATE(sym_A) && !IS_EPSILON(sym_C)){
                        // C: [k, i - 1] part
                        double alpha_C = ALPHA_GET(sym_C, k, i - 1);
                        // B: [k, j] part
                        double beta_B = BETA(sym_B, k, j);

                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], possibility + alpha_C + beta_B);
                    }   
                }

                for(int k = j + 1; k < sequence_length; k++){       
                    uint32_t sym_B = _sym_A;
                    uint32_t sym_A = _sym_B;
                    uint32_t sym_C = _sym_C;

                    // 2. 6. (A is nonterminate, and C is terminate or nonterminate. A at left side.)
                    if(IS_NONTERMINATE(sym_A) && !IS_EPSILON(sym_C)){
                        // C: [k, i - 1] part
                        double alpha_C = ALPHA_GET(sym_C, j + 1, k);
                        // B: [k, j] part
                        double beta_B = BETA(sym_B, i, k);
                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i],  possibility + alpha_C + beta_B);
                        
                    }   
                }
            }
            
            // 3.
            for(int p = inside_order_1_rule_iteration_path_size - 1; p >= 0; p--){ 
                // B -> A, A is a nonterminate.t)
                uint32_t gid = *(inside_order_1_rule_iteration_path + p * 2);
                uint32_t sym_B = *(inside_order_1_rule_iteration_path + p * 2 + 1);
            
                uint32_t sym_A = grammar_table[(n_grammars + 1) * 1 + gid];
                double possibility = *(double*)(grammar_table + (n_grammars + 1) * 4 + gid * 2);
            

                double alpha_B = ALPHA_GET(sym_B, i, j);
                double beta_B = POSSIBILITY_ADD(BETA(sym_B, i, j), buffer_beta[sym_B * MS + i]);
                
                if(IS_TERMINATE(sym_A)) continue; 
                /* grammar become B -> A. In this condition, B -> A contributes possibility * beta_B
                    to A's outside possibility spanning i to j-th symbols in the sequence. 
                    We doesn't need to iterate split point k, as there only one symbol in the right side
                    of this rule. 'continue;' is uesed to skip k iterations. */
                
                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], possibility + beta_B);
                
            }


            // 5. 8. 9. 10.
            // Case 2: A is terminate 
            for(int gid = 0; gid < n_grammars; gid++){
                int j = i + span_length - 1; // Ending index of the span
                uint32_t _sym_A = symbol_A_vector[gid];
                uint32_t sym_BC = *(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                uint32_t _sym_B = (sym_BC >> 16) & 0xFFFF;
                uint32_t _sym_C = (sym_BC) & 0xFFFF;
                double possibility = *(double*)(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1);
            
                uint32_t sym_B = _sym_A;
                uint32_t sym_A = _sym_B;
                uint32_t sym_C = _sym_C;
                if(IS_NONTERMINATE(sym_A) || i != j) break;
                if(IS_EPSILON(sym_C)) continue;

                // 5. B -> w_A C
                if(IS_TERMINATE(sym_A) && IS_NONTERMINATE(sym_C) && (i == j)){
                    for(int k = j + 1; k < sequence_length; k++){
                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], BETA(sym_B, i, k) + ALPHA(sym_C, j + 1, k) + possibility);
                    }
                }

                // 9. B->w_C w_A 
                if(IS_TERMINATE(sym_A) && IS_TERMINATE(sym_C) && (i == j)){
                    if(j + 1 < sequence_length){
                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], ALPHA(sym_C, j + 1, j + 1) + BETA(sym_B, i, j + 1) + possibility);
                    }
                }
                
                sym_C = _sym_B;
                sym_A = _sym_C;

                // 8.
                if(IS_NONTERMINATE(sym_C) && IS_TERMINATE(sym_A) && (i == j)){
                    for(int k = 0; k <= i - 1; k++){
                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], ALPHA(sym_C, k, i - 1) + BETA(sym_B, k, j) + possibility);
                    }
                }

                // 10.
                if(IS_TERMINATE(sym_C) && IS_TERMINATE(sym_A) && (i == j)){
                    if(i - 1 >= 0){
                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], ALPHA(sym_C, i - 1, i - 1) + BETA(sym_B, i - 1, j) + possibility);
                    }                   
                }
            }


            // 11.
            for(int p = inside_order_1_rule_iteration_path_size - 1; p >= 0; p--){ 
                // B -> A
                uint32_t gid = *(inside_order_1_rule_iteration_path + p * 2);
                uint32_t sym_B = *(inside_order_1_rule_iteration_path + p * 2 + 1); 
                
                uint32_t sym_A = grammar_table[(n_grammars + 1) * 1 + gid];
                double possibility = *(double*)(grammar_table + (n_grammars + 1) * 4 + gid * 2);
                

                double alpha_B = ALPHA_GET(sym_B, i, j);
                double beta_B = POSSIBILITY_ADD(BETA(sym_B, i, j), buffer_beta[sym_B * MS + i]);
                if(i != j) break;
                if(IS_NONTERMINATE(sym_A)) continue; 
                // B->w_A
                
                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i],  possibility + beta_B);
            }


            // write back
            for(int sym_A = 0; sym_A < N; sym_A++){
                LOG_SUM_EXP_SET(BETA(sym_A, i, j), buffer_beta[sym_A * MS + i]);
            }
        } // parallel for end.
        __syncthreads();
    }

    __syncthreads();

    // fill mu[grammar_id, i, j]
    for (int span_length = 1; span_length < sequence_length + 1; span_length++) {
        double local_buffer_mu[512 * 256];

        if(thread_id == 0){
            for(int i = 0; i < MS * n_grammars; i++){
                local_buffer_mu[i] = -INFINITY;
            }
        }
        
        __syncthreads();

        for (int i = thread_id; i <= sequence_length - span_length; i += total_threads) {
            int j = i + span_length - 1; // Ending index of the spanx`
            for(int gid = 0; gid < n_grammars; gid++){
                uint32_t sym_A = symbol_A_vector[gid];
                uint32_t sym_BC = *(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                uint32_t sym_B = (sym_BC >> 16) & 0xFFFF;
                uint32_t sym_C = (sym_BC) & 0xFFFF;
                for (int k = i; k <= j; k++) { // TODO: k < j? k == j?
                
                    double possibility = *(double*)(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1);

                    if(k == j && !IS_EPSILON(sym_C)) continue;

                    double beta_A_i_j = BETA(sym_A, i, j);
                    double alpha_B_i_k = ALPHA_GET(sym_B, i, k);
                    double alpha_C_k_p1_j = ALPHA_GET(sym_C, k + 1, j);
                    
                    #ifdef COMPUTING_IN_LOG_SPACE
                        LOG_SUM_EXP_SET(local_buffer_mu[gid * MS + i], 
                                            possibility + beta_A_i_j + alpha_B_i_k + alpha_C_k_p1_j);                        
                    #else
                        local_buffer_mu[gid * MS + i] +=  possibility * beta_A_i_j * alpha_B_i_k * alpha_C_k_p1_j;
                    #endif
                }
            }
            
            // write back.
            for(int gid = 0; gid < n_grammars; gid++){
                LOG_SUM_EXP_SET(MU(gid, i, j), local_buffer_mu[gid * MS + i]);
            }
        } // end parallel for
        __syncthreads();
    }
}
#endif