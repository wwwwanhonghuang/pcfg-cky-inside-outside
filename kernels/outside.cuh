#include <omp.h>
#ifndef USE_CUDA
#include <cstring>
#endif

#include "utils/data_accessing.hpp"
#include "macros.def"
#define DEBUG_OUTSIDE_CELL(x, y, X) if(i == x && j == y) { X }

#ifdef USE_CUDA
__global__
#endif
void kernel_outside_main(float* mu, float* beta, uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , pcfg* grammar
                        #endif
){
    #ifndef USE_CUDA
    memset(mu, 0, n_grammars * MS * MS * sizeof(float));
    memset(beta, 0, n_syms * MS * MS * sizeof(float));

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


    // Case 1: A is non-terminate 
    // confirm these codes are correct.
    /* diagonal-order iteration */
    for(int span_length = sequence_length; span_length >= 1; span_length--){
        /* for one diagnal of the beta table, all cell can be parallelly computed. */
        #pragma omp parallel for
        for(int i = 0; i < sequence_length - span_length + 1; i++){
            int j = i + span_length - 1;
            // 1. 2. 6. 7.
            for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                PCFGItemIterator(N, grammar_index, grammar_table)){                    
                for(int k = 0; k < i; k++){       
                    uint32_t sym_B = std::get<0>(item);
                    uint32_t sym_C = std::get<1>(item);
                    uint32_t sym_A = std::get<2>(item);
                    float possibility = std::get<3>(item);
                    
                    // 1. 7. (A is nonterminate, and C is terminate or nonterminate. A at right side.)
                    if(IS_NONTERMINATE(sym_A) && !IS_EPSILON(sym_C)){
                        // C: [k, i - 1] part
                        float alpha_C = ALPHA_GET(sym_C, k, i - 1);
                        // B: [k, j] part
                        float beta_B = BETA(sym_B, k, j);

                        #pragma omp atomic
                        // A: [i, j] part
                        BETA_INCREASE(sym_A, i, j, possibility * alpha_C * beta_B);   
                    }   
                }

                for(int k = j + 1; k < sequence_length; k++){       
                    uint32_t sym_B = std::get<0>(item);
                    uint32_t sym_A = std::get<1>(item);
                    uint32_t sym_C = std::get<2>(item);
                    float possibility = std::get<3>(item);

                    // 2. 6. (A is nonterminate, and C is terminate or nonterminate. A at left side.)
                    if(IS_NONTERMINATE(sym_A) && !IS_EPSILON(sym_C)){
                        // C: [k, i - 1] part
                        float alpha_C = ALPHA_GET(sym_C, j + 1, k);
                        // B: [k, j] part
                        float beta_B = BETA(sym_B, i, k);

                        #pragma omp atomic
                        // A: [i, j] part
                        BETA_INCREASE(sym_A, i, j, possibility * alpha_C * beta_B);   
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
                uint32_t* addr = (grammar_table + gid * 2);
                uint32_t sym_A = ((*addr) >> 16) & 0xFFFF;
                float alpha_B = ALPHA_GET(sym_B, i, j);
                float possibility = *(float*)(addr + 1);
                
                if(IS_TERMINATE(sym_A)) continue; 
                /* grammar become B -> A. In this condition, B -> A contributes possibility * beta_B
                    to A's outside possibility spanning i to j-th symbols in the sequence. 
                    We doesn't need to iterate split point k, as there only one symbol in the right side
                    of this rule. 'continue;' is uesed to skip k iterations. */
                
                #pragma omp atomic
                BETA_INCREASE(sym_A, i, j, possibility * BETA(sym_B, i, j));
            }

            // 5. 8. 9. 10.
            // Case 2: A is terminate 
            for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                                                        PCFGItemIterator(N, grammar_index, grammar_table)){                    
                    
                uint32_t sym_B = std::get<0>(item);
                float possibility = std::get<3>(item);
                
                uint32_t sym_A = std::get<1>(item);
                uint32_t sym_C = std::get<2>(item);
                if(IS_NONTERMINATE(sym_A) || i != j) break;
                if(IS_EPSILON(sym_C)) continue;

                // 5. B -> w_A C
                if(IS_TERMINATE(sym_A) && IS_NONTERMINATE(sym_C) && (i == j)){
                    for(int k = j + 1; k < sequence_length; k++){
                            BETA(sym_A, i, j) += BETA(sym_B, i, k) * ALPHA(sym_C, j + 1, k) * possibility;
                    }
                }

                // 9. B->w_C w_A 
                if(IS_TERMINATE(sym_A) && IS_TERMINATE(sym_C) && (i == j)){
                    if(j + 1 < sequence_length){
                        BETA(sym_A, i, j) += ALPHA(sym_C, j + 1, j + 1) * BETA(sym_B, i, j + 1) * possibility;
                    }
                }
                
                sym_C = std::get<1>(item);
                sym_A = std::get<2>(item);

                // 8.
                if(IS_NONTERMINATE(sym_C) && IS_TERMINATE(sym_A) && (i == j)){
                    for(int k = 0; k <= i - 1; k++){
                        BETA(sym_A, i, j) += ALPHA(sym_C, k, i - 1) * BETA(sym_B, k, j) * possibility;
                    }
                }

                // 10.
                if(IS_TERMINATE(sym_C) && IS_TERMINATE(sym_A) && (i == j)){
                    if(i - 1 >= 0){
                        BETA(sym_A, i, j) += ALPHA(sym_C, i - 1, i - 1) * BETA(sym_B, i - 1, j) * possibility;
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
                uint32_t* addr = (grammar_table + gid * 2);
                uint32_t sym_A = ((*addr) >> 16) & 0xFFFF;
                float alpha_B = ALPHA_GET(sym_B, i, j);
                float possibility = *(float*)(addr + 1);
                if(i != j) break;
                if(IS_NONTERMINATE(sym_A)) continue; 
                // B->w_A
                
                #pragma omp atomic
                BETA_INCREASE(sym_A, i, j, possibility * BETA(sym_B, i, j));
            }
        }
    }

    
    
    // fill mu[grammar_id, i, j]
    for (int span_length = 1; span_length < sequence_length + 1; span_length++) {
        #pragma omp parallel for
        for (int i = 0; i <= sequence_length - span_length; i++) {
            int j = i + span_length - 1; // Ending index of the spanx`
            for (int k = i; k <= j; k++) { // TODO: k < j? k == j?
                for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                    uint32_t sym_A = std::get<0>(item);
                    uint32_t sym_B = std::get<1>(item);
                    uint32_t sym_C = std::get<2>(item);
                    float possibility = std::get<3>(item);
                    uint32_t gid = std::get<4>(item);
                    if(k == j && !IS_EPSILON(sym_C)) continue;

                    float beta_A_i_j = BETA(sym_A, i, j);
                    float alpha_B_i_k = ALPHA_GET(sym_B, i, k);
                    float alpha_C_k_p1_j = ALPHA_GET(sym_C, k + 1, j);
                    #pragma omp atomic
                    MU_INCREASE(gid, i, j, possibility * beta_A_i_j * alpha_B_i_k * alpha_C_k_p1_j);
                }
            }
        }
    }

    
    #else
    std::err << "Error: CUDA Version Outside Algorithm is currently not be implemented." << std::endl;
    return nullptr; // NotImplemented
    #endif
}