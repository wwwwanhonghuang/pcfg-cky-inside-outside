
#include "kernels/outside.cuh"
#include "utils/math.hpp"

void kernel_outside_main(double* mu, double* beta, const uint32_t* sequence, 
                        uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , pcfg* grammar
                        #endif
){
    memset(mu, 0, n_grammars * MS * MS * sizeof(double));
    memset(beta, 0, n_syms * MS * MS * sizeof(double));

    #ifdef COMPUTING_IN_LOG_SPACE
    for (int i = 0; i < n_grammars * MS * MS; i++) {
        mu[i] = -INFINITY;
    }
    for (int i = 0; i < n_syms * MS * MS; i++) {
        beta[i] = -INFINITY;
    }
    #endif
    
    /* base case: S is the root of the whole sequence with possibility 1.0. */
    BETA(0, 0, sequence_length - 1) = 1.0;

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

    std::vector<double> buffer_beta(n_syms * MS,
            #ifdef COMPUTING_IN_LOG_SPACE
            -INFINITY
            #else
            0.0
            #endif
        );
    // Case 1: A is non-terminate 
    // confirm these codes are correct.
    /* diagonal-order iteration */
    for(int span_length = sequence_length; span_length >= 1; span_length--){
        /* for one diagnal of the beta table, all cell can be parallelly computed. */
        
        std::fill(buffer_beta.begin(), buffer_beta.end(),
            #ifdef COMPUTING_IN_LOG_SPACE
            -INFINITY
            #else
            0.0
            #endif
        );
        
        #pragma omp parallel
        {
            #pragma omp for
            for(int i = 0; i < sequence_length - span_length + 1; i++){
                int j = i + span_length - 1;
                
                // 1. 2. 6. 7.
                for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                    PCFGItemIterator(N, grammar_index, grammar_table)){                    
                    for(int k = 0; k < i; k++){       
                        uint32_t sym_B = std::get<0>(item);
                        uint32_t sym_C = std::get<1>(item);
                        uint32_t sym_A = std::get<2>(item);
                        double possibility = std::get<3>(item);
                        
                        // 1. 7. (A is nonterminate, and C is terminate or nonterminate. A at right side.)
                        if(IS_NONTERMINATE(sym_A) && !IS_EPSILON(sym_C)){
                            // C: [k, i - 1] part
                            double alpha_C = ALPHA_GET(sym_C, k, i - 1);
                            // B: [k, j] part
                            double beta_B = BETA(sym_B, k, j);

                            #ifdef COMPUTING_IN_LOG_SPACE
                                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], possibility + alpha_C + beta_B);                    
                            #else
                                // A: [i, j] part
                                buffer_beta[sym_A * MS + i] += possibility * alpha_C * beta_B;
                            #endif
                        }   
                    }

                    for(int k = j + 1; k < sequence_length; k++){       
                        uint32_t sym_B = std::get<0>(item);
                        uint32_t sym_A = std::get<1>(item);
                        uint32_t sym_C = std::get<2>(item);
                        double possibility = std::get<3>(item);

                        // 2. 6. (A is nonterminate, and C is terminate or nonterminate. A at left side.)
                        if(IS_NONTERMINATE(sym_A) && !IS_EPSILON(sym_C)){
                            // C: [k, i - 1] part
                            double alpha_C = ALPHA_GET(sym_C, j + 1, k);
                            // B: [k, j] part
                            double beta_B = BETA(sym_B, i, k);
                            #ifdef COMPUTING_IN_LOG_SPACE
                                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i],  possibility + alpha_C + beta_B);
                            #else
                                // A: [i, j] part
                                buffer_beta[sym_A * MS + i] += possibility * alpha_C * beta_B;
                            #endif
                        }   
                    }
                }

                // 3.
                for(std::vector<std::tuple<uint32_t, uint32_t>>::reverse_iterator it = inside_order_1_rule_iteration_path.rbegin(); 
                    // B -> A, A is a nonterminate.
                    it != inside_order_1_rule_iteration_path.rend(); ++it) {
                    std::tuple<uint32_t, uint32_t> rule_id = *it;
                    uint32_t gid = std::get<0>(rule_id);
                    uint32_t sym_B = std::get<1>(rule_id);

                    #ifndef ENABLE_GRAMMAR_VECTORIZATION_OPTIMIZATION
                        uint32_t* addr = (grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                        uint32_t sym_A = ((*addr) >> 16) & 0xFFFF;
                        double possibility = *(double*)(addr + 1);
                    #else
                        uint32_t sym_A = grammar_table[(n_grammars + 1) * 1 + gid];
                        double possibility = *(double*)(grammar_table + (n_grammars + 1) * 4 + gid * 2);
                    #endif

                    double alpha_B = ALPHA_GET(sym_B, i, j);
                    double beta_B = POSSIBILITY_ADD(BETA(sym_B, i, j), buffer_beta[sym_B * MS + i]);
                    
                    if(IS_TERMINATE(sym_A)) continue; 
                    /* grammar become B -> A. In this condition, B -> A contributes possibility * beta_B
                        to A's outside possibility spanning i to j-th symbols in the sequence. 
                        We doesn't need to iterate split point k, as there only one symbol in the right side
                        of this rule. 'continue;' is uesed to skip k iterations. */
                    
                    #ifdef COMPUTING_IN_LOG_SPACE
                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], possibility + beta_B);
                    #else
                        buffer_beta[sym_A * MS + i] += possibility * beta_B;
                    #endif
                }

                // 5. 8. 9. 10.
                // Case 2: A is terminate 
                for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                                                            PCFGItemIterator(N, grammar_index, grammar_table)){                    
                        
                    uint32_t sym_B = std::get<0>(item);
                    double possibility = std::get<3>(item);
                    
                    uint32_t sym_A = std::get<1>(item);
                    uint32_t sym_C = std::get<2>(item);
                    if(IS_NONTERMINATE(sym_A) || i != j) break;
                    if(IS_EPSILON(sym_C)) continue;

                    // 5. B -> w_A C
                    if(IS_TERMINATE(sym_A) && IS_NONTERMINATE(sym_C) && (i == j)){
                        for(int k = j + 1; k < sequence_length; k++){
                            #ifdef COMPUTING_IN_LOG_SPACE
                                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], BETA(sym_B, i, k) + ALPHA(sym_C, j + 1, k) + possibility);
                            #else
                                buffer_beta[sym_A * MS + i] += BETA(sym_B, i, k) * ALPHA(sym_C, j + 1, k) * possibility;
                            #endif
                        }
                    }

                    // 9. B->w_C w_A 
                    if(IS_TERMINATE(sym_A) && IS_TERMINATE(sym_C) && (i == j)){
                        if(j + 1 < sequence_length){
                            #ifdef COMPUTING_IN_LOG_SPACE
                                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], ALPHA(sym_C, j + 1, j + 1) + BETA(sym_B, i, j + 1) + possibility);
                            #else
                                buffer_beta[sym_A * MS + i] += ALPHA(sym_C, j + 1, j + 1) * BETA(sym_B, i, j + 1) * possibility;
                            #endif
                        }
                    }
                    
                    sym_C = std::get<1>(item);
                    sym_A = std::get<2>(item);

                    // 8.
                    if(IS_NONTERMINATE(sym_C) && IS_TERMINATE(sym_A) && (i == j)){
                        for(int k = 0; k <= i - 1; k++){
                            #ifdef COMPUTING_IN_LOG_SPACE
                                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], ALPHA(sym_C, k, i - 1) + BETA(sym_B, k, j) + possibility);
                            #else
                                buffer_beta[sym_A * MS + i] += ALPHA(sym_C, k, i - 1) * BETA(sym_B, k, j) * possibility;
                            #endif
                        }
                    }

                    // 10.
                    if(IS_TERMINATE(sym_C) && IS_TERMINATE(sym_A) && (i == j)){
                        if(i - 1 >= 0){
                            #ifdef COMPUTING_IN_LOG_SPACE
                                LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i], ALPHA(sym_C, i - 1, i - 1) + BETA(sym_B, i - 1, j) + possibility);
                            #else
                                buffer_beta[sym_A * MS + i] += ALPHA(sym_C, i - 1, i - 1) * BETA(sym_B, i - 1, j) * possibility;
                            #endif
                        }                   
                    }

                    
                
                }

                // 11.
                for(std::vector<std::tuple<uint32_t, uint32_t>>::reverse_iterator it = inside_order_1_rule_iteration_path.rbegin(); 
                    // B -> A
                    it != inside_order_1_rule_iteration_path.rend(); ++it) {
                    std::tuple<uint32_t, uint32_t> rule_id = *it;
                    uint32_t gid = std::get<0>(rule_id);
                    uint32_t sym_B = std::get<1>(rule_id);
                    #ifndef ENABLE_GRAMMAR_VECTORIZATION_OPTIMIZATION
                    uint32_t* addr = (grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                    uint32_t sym_A = ((*addr) >> 16) & 0xFFFF;
                    double possibility = *(double*)(addr + 1);
                    #else
                    uint32_t sym_A = grammar_table[(n_grammars + 1) * 1 + gid];
                    double possibility = *(double*)(grammar_table + (n_grammars + 1) * 4 + gid * 2);
                    #endif

                    double alpha_B = ALPHA_GET(sym_B, i, j);
                    double beta_B = POSSIBILITY_ADD(BETA(sym_B, i, j), buffer_beta[sym_B * MS + i]);
                    if(i != j) break;
                    if(IS_NONTERMINATE(sym_A)) continue; 
                    // B->w_A
                    
                    #ifdef COMPUTING_IN_LOG_SPACE
                        LOG_SUM_EXP_SET(buffer_beta[sym_A * MS + i],  possibility + beta_B);
                    #else
                        buffer_beta[sym_A * MS + i] += possibility * beta_B;
                    #endif
                }

                // write back
                for(int sym_A = 0; sym_A < N; sym_A++){
                    #ifdef COMPUTING_IN_LOG_SPACE
                        LOG_SUM_EXP_SET(BETA(sym_A, i, j), buffer_beta[sym_A * MS + i]);
                    #else
                        BETA_INCREASE(sym_A, i, j, buffer_beta[sym_A * MS + i]);
                    #endif
                }
            } // parallel for end.

        }
    }

    
    // fill mu[grammar_id, i, j]
    for (int span_length = 1; span_length < sequence_length + 1; span_length++) {
        std::vector<double> local_buffer_mu(MS * n_grammars,
            #ifdef COMPUTING_IN_LOG_SPACE
                -INFINITY
            #else
                0.0
            #endif
        );

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i <= sequence_length - span_length; i++) {
                int j = i + span_length - 1; // Ending index of the spanx`
                
                for (int k = i; k <= j; k++) { // TODO: k < j? k == j?
                    for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                        uint32_t sym_A = std::get<0>(item);
                        uint32_t sym_B = std::get<1>(item);
                        uint32_t sym_C = std::get<2>(item);
                        double possibility = std::get<3>(item);
                        uint32_t gid = std::get<4>(item);
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
                #ifdef COMPUTING_IN_LOG_SPACE
                    for(int gid = 0; gid < n_grammars; gid++){
                        LOG_SUM_EXP_SET(MU(gid, i, j), local_buffer_mu[gid * MS + i]);
                    }
                #else
                    for(int gid = 0; gid < n_grammars; gid++){
                        MU_INCREASE(gid, i, j, local_buffer_mu[gid * MS + i]);
                    }
                #endif
            
            } // end parallel for
        }

        
        
    }
}
