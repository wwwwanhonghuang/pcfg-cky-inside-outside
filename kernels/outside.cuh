#include <omp.h>
#ifndef USE_CUDA
#include <cstring>
#endif

#include "utils/data_accessing.hpp"
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
    beta[0 * MS * MS + 0 * MS + sequence_length - 1] = 1.0; 

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
                if(sym_A == 0xFFFF || sym_A >= N){
                        continue; 
                }
                    
                for(int k = 0; k < i; k++){
                    /* Specifal condition 1: Sym_C is a terminate, we have form B -> 'w' A */
                    if(sym_C >= N && sym_A < N){ 
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
                        float beta_A_i_j = beta_get(beta, sym_A, i, j, MS);
                        float condition = (sequence[k] == sym_C && (k == i - 1) ? 1.0f : 0.0f);
                        float delta = possibility * condition * beta_get(beta, sym_B, k, j, MS);
                        beta_set(beta, beta_A_i_j + delta, sym_A, i, j, MS);                            
                        continue;
                    }
                    
                    // C: [k, i - 1] part
                    float alpha_C = alpha_get(alpha, sym_C, k, i - 1, MS);
                    // B: [k, j] part
                    float beta_B = beta_get(beta, sym_B, k, j, MS);

                    #pragma omp atomic
                    // A: [i, j] part
                    beta[sym_A * MS * MS + i * MS + j] += possibility * alpha_C * beta_B;
                }
            }
            
            // B -> AC.
            for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                uint32_t sym_B = std::get<0>(item);
                uint32_t sym_A = std::get<1>(item);
                uint32_t sym_C = std::get<2>(item);
                float possibility = std::get<3>(item);
                if(sym_C == 0xFFFF){ // epsilon
                        /* grammar become B -> A. In this condition, B -> A contributes possibility * beta_B
                           to A's outside possibility spanning i to j-th symbols in the sequence. 
                           We doesn't need to iterate split point k, as there only one symbol in the right side
                           of this rule. 'continue;' is uesed to skip k iterations.
                        */
                        beta[sym_A * MS * MS + i * MS + j] += possibility * beta[sym_B * MS * MS + i * MS + j];
                        continue;
                }

                /* it doesn't need to calculate a terminate's possibility. */
                if(sym_A >= N) continue; 

                for(int k = j + 1; k < sequence_length; k++){
                    /* in the condition of B -> A 'w'm and B span i-k.
                        'w' must span k-k, sequence[k] must equal to 'w', and A must span i-j, where j == k - 1. */
                    if(sym_A < N && sym_C >= N){
                        float condition =  (sequence[k] == sym_C && k == j + 1 ? 1.0 : 0.0);
                        beta[sym_A * MS * MS + i * MS + j] += 
                            possibility * condition * beta_get(beta, sym_B, i, k, MS);
                        continue;
                    }
                        
                    /* C: [j + 1, k] part */
                    float alpha_C = alpha_get(alpha, sym_C, j + 1, k, MS);
                    // B: [i, k] part
                    float beta_B = beta_get(beta, sym_B, i, k, MS); 
                    
                    #pragma omp atomic
                    beta[sym_A * MS * MS + i * MS + j] += possibility * alpha_C * beta_B;
                }                  
            }
        }
    }

    for (int span_length = 2; span_length <= sequence_length; span_length++) {
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

                    float beta_A_i_j = beta_get(beta, sym_A, i, j, MS);

                    if(sym_B >= N && sym_C >= N){
                        
                        if(sym_C == 0xFFFF){
                            // unreachable code.
                            #pragma omp atomic
                            mu[gid * MS * MS + i * MS + j] += possibility * beta_A_i_j * 
                                (sequence[i] == sym_B && i == j ? 1.0f : 0.0f);
                        }else{
                            float condition = (sequence[i] == sym_B && sequence[j] == sym_C &&  span_length == 2 ? 1.0f : 0.0f);
                            #pragma omp atomic
                            mu[gid * MS * MS + i * MS + j] += possibility * beta_A_i_j * condition;
                        }
                    }else if(sym_B < N && sym_C >= N){
                        if(sym_C == 0xFFFF && k == i){
                            float alpha_B_i_j = alpha_get(alpha, sym_B, i, j, MS);
                            #pragma omp atomic
                            mu[gid * MS * MS + i * MS + j] += possibility * beta_A_i_j * alpha_B_i_j;
                        }else{
                            float condition = (sequence[j] == sym_C && k == j - 1 ? 1.0f : 0.0f);
                            float alpha_B_i_k = alpha_get(alpha, sym_B, i, k, MS);
                            #pragma omp atomic
                            mu[gid * MS * MS + i * MS + j] += possibility * beta_A_i_j * alpha_B_i_k * condition;
                        }
                    }else if(sym_B >= N && sym_C < N){
                        float condition = (sequence[i] == sym_B && k == i ? 1.0f : 0.0f);
                        float alpha_C_k_p1_j = alpha_get(alpha, sym_B, k + 1, j, MS);
                        #pragma omp atomic
                        mu[gid * MS * MS + i * MS + j] += possibility * beta_A_i_j * alpha_C_k_p1_j * condition;
                    }else{
                        float alpha_B_i_k = alpha_get(alpha, sym_B, i, k, MS);
                        float alpha_C_k_p1_j = alpha_get(alpha, sym_B, k + 1, j, MS);
                        #pragma omp atomic
                        mu[gid * MS * MS + i * MS + j] += possibility * beta_A_i_j * alpha_B_i_k * alpha_C_k_p1_j;
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