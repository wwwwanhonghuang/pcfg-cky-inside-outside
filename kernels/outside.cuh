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
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                                #ifdef DEBUG_INSIDE_ALGORITHM
                                , pcfg* grammar
                                #endif
                        ){
    #ifndef USE_CUDA
    memset(mu, 0, n_grammars * MS * MS * sizeof(float));
    memset(beta, 0, N * MS * MS * sizeof(float));

    /* base case: S is the root of the whole sequence with possibility 1.0. */
    BETA(0, 0, sequence_length - 1) = 1.0;

    /* diagonal-order iteration */
    for(int span_length = 1; span_length < sequence_length + 1; span_length++){
        /* for one diagnal of the beta table, all cell can be parallelly computed. */
        #pragma omp parallel for
        for(int i = 0; i < sequence_length - span_length + 1; i++){
            int j = i + span_length - 1;
            for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                uint32_t sym_B = std::get<0>(item);
                uint32_t sym_C = std::get<1>(item);
                uint32_t sym_A = std::get<2>(item);
                float possibility = std::get<3>(item);

                /* the outside possibility of the symbol we are currently attempting to calculate (i.e., symbol A)
                 cannot be empty, and it must be a nonterminate. */
                if(IS_EPSILON(sym_A) || IS_TERMINATE(sym_A)){
                    continue; 
                }
                    
                for(int k = 0; k < i; k++){
                    /* Specifal condition 1: Sym_C is a terminate, we have form B -> 'w' A */
                    if(IS_TERMINATE(sym_C) && IS_NONTERMINATE(sym_A)){ 
                        /* In the situation that C is a terminate 'w', as B span k to j,
                        'w' mast span k to k, and A must span k + 1 (= i, i.e., i = k - 1) to j, so that we can have B -> CA. 
                        If the conditions above are not satisfied, this loop must contribute 0.0 to symbol A's outside
                        possibility.

                         Or it will consider contribute to the A's outside possibility p(B->CA) * outside['B', k, j].
                         This possibility can be comprehened as, when we know the possibility of 
                         parsing context (B -> * *)'s possibility, 
                         B -> 'w' A is one of the possibile condition, 
                         with possibity of p(B->'w' A) * outside['B', k, j].
                         A has possibility of p(B->'w' A) * outside['B', k, j] become a root for spanning i to j under 
                         grammar B -> 'w' A.
                         we need add this contributions into A's outside possibity spanning i, j. */
                        float condition = (sequence[k] == sym_C && (k == i - 1) ? 1.0f : 0.0f);
                        float delta = possibility * condition * BETA(sym_B, k, j);
                        BETA_INCREASE(sym_A, i, j, delta);
                        continue;
                    }
                    
                    // C: [k, i - 1] part
                    float alpha_C = ALPHA(sym_C, k, i - 1);
                    // B: [k, j] part
                    float beta_B = BETA(sym_B, k, j);

                    #pragma omp atomic
                    // A: [i, j] part
                    BETA_INCREASE(sym_A, i, j, possibility * alpha_C * beta_B);
                }
            }
            
            // B -> AC.
            for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                uint32_t sym_B = std::get<0>(item);
                uint32_t sym_A = std::get<1>(item);
                uint32_t sym_C = std::get<2>(item);
                float possibility = std::get<3>(item);
                if(IS_EPSILON(sym_C)){ // epsilon
                        /* grammar become B -> A. In this condition, B -> A contributes possibility * beta_B
                           to A's outside possibility spanning i to j-th symbols in the sequence. 
                           We doesn't need to iterate split point k, as there only one symbol in the right side
                           of this rule. 'continue;' is uesed to skip k iterations.
                        */
                        BETA_INCREASE(sym_A, i, j, possibility * BETA(sym_B, i, j));
                        continue;
                }

                /* it doesn't need to calculate a terminate's outside possibility. */
                if(IS_TERMINATE(sym_A)) continue; 

                for(int k = j + 1; k < sequence_length; k++){
                    /* in the condition of B -> A 'w'm and B span i-k.
                        'w' must span k-k, sequence[k] must equal to 'w', and A must span i-j, where j == k - 1. */
                    if(IS_NONTERMINATE(sym_A) && IS_TERMINATE(sym_C)){
                        float condition =  (sequence[k] == sym_C && k == j + 1 ? 1.0 : 0.0);
                        BETA_INCREASE(sym_A, i, j, possibility * condition * BETA(sym_B, i, k));
                        continue;
                    }
                        
                    /* C: [j + 1, k] part */
                    float alpha_C = ALPHA(sym_C, j + 1, k);
                    // B: [i, k] part
                    float beta_B = BETA(sym_B, i, k); 
                    
                    #pragma omp atomic
                    BETA_INCREASE(sym_A, i, j, possibility * alpha_C * beta_B);
                }                  
            }
        }
    }

    // fill mu[grammar_id, i, j]
    for (int span_length = 1; span_length < sequence_length + 1; span_length++) {
        #pragma omp parallel for
        for (int i = 0; i <= sequence_length - span_length; i++) {
            int j = i + span_length - 1; // Ending index of the spanx`
            for (int k = i; k < j; k++) { // TODO: k < j? k == j?
                for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                    uint32_t sym_A = std::get<0>(item);
                    uint32_t sym_B = std::get<1>(item);
                    uint32_t sym_C = std::get<2>(item);
                    float possibility = std::get<3>(item);
                    uint32_t gid = std::get<4>(item);

                    float beta_A_i_j = BETA(sym_A, i, j);

                    if(IS_TERMINATE(sym_B) && IS_TERMINATE(sym_C)){
                        if(IS_EPSILON(sym_C)){
                            // unreachable code.
                            #pragma omp atomic
                            MU_INCREASE(gid, i, j, possibility * beta_A_i_j * 
                                (sequence[i] == sym_B && i == j ? 1.0f : 0.0f));
                        }else{
                            float condition = (sequence[i] == sym_B && sequence[j] == sym_C &&  span_length == 2 ? 1.0f : 0.0f);
                            #pragma omp atomic
                            MU_INCREASE(gid, i, j, possibility * beta_A_i_j * condition);
                        }
                    }else if(IS_NONTERMINATE(sym_B) && IS_TERMINATE(sym_C)){
                        if(IS_EPSILON(sym_C)){
                            /* this condition may A -> B only be considered once in multiple k-axis loops. */
                            float limit_term = (i == k ? 1.0f : 0.0f); 
                            float alpha_B_i_j = ALPHA(sym_B, i, j);
                            #pragma omp atomic
                            MU_INCREASE(gid, i, j, possibility * beta_A_i_j * alpha_B_i_j * limit_term);
                        }else{
                            float condition = (sequence[j] == sym_C && k == j - 1 ? 1.0f : 0.0f);
                            float alpha_B_i_k = ALPHA(sym_B, i, k);
                            #pragma omp atomic
                            MU_INCREASE(gid, i, j, possibility * beta_A_i_j * alpha_B_i_k * condition);
                        }
                    }else if(sym_B >= N && sym_C < N){
                        float condition = (sequence[i] == sym_B && k == i ? 1.0f : 0.0f);
                        float alpha_C_k_p1_j = ALPHA(sym_C, k + 1, j);
                        #pragma omp atomic
                        MU_INCREASE(gid, i, j, possibility * beta_A_i_j * alpha_C_k_p1_j * condition);
                    }else{
                        float alpha_B_i_k = ALPHA(sym_B, i, k);
                        float alpha_C_k_p1_j = ALPHA(sym_C, k + 1, j);
                        #pragma omp atomic
                        MU_INCREASE(gid, i, j, possibility * beta_A_i_j * alpha_B_i_k * alpha_C_k_p1_j);
                    }
                }
            }
        }
    }

    
    #else
    std::err << "Error: CUDA Version Outside Algorithm is currently not be implemented." << std::endl;
    return nullptr; // NotImplemented
    #endif
}