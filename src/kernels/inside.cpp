#include <math.h>
#include "kernels/inside.cuh"
#include "utils/math.hpp"
#include "constants.h"

#ifdef USE_CUDA
__global__
#endif
void kernel_inside_alpha_zerolization(double* alpha, int N, int MS){
    #ifdef COMPUTING_IN_LOG_SPACE
        memset(alpha, 0, N * MS * MS * sizeof(double));
        for (int i = 0; i < N * MS * MS; i++) {
            alpha[i] = -INFINITY;
        }
    #else
        memset(alpha, 0, N * MS * MS * sizeof(double));
    #endif
}

#ifdef USE_CUDA
__global__ 
#endif
void kernel_inside_base_fill_alpha(  
        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
        , pcfg* grammar
        #endif
        ) {
        #pragma omp parallel for
        for(int sym_A = 0; sym_A < N; sym_A ++){
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
        std::unordered_set<uint32_t> nonterminate_remains; 
        
        for(int sym_id = 0; sym_id < N; sym_id++){
            nonterminate_remains.insert(sym_id);
        }

        for(int i = 0; i < N; i++){
            int size = nonterminate_remains.size();
            for(int i = 0; i < sequence_length; i++){
                for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                                                        PCFGItemIterator(N, grammar_index, grammar_table)){
                    uint32_t sym_A = std::get<0>(item);
                    uint32_t sym_B = std::get<1>(item);
                    uint32_t sym_C = std::get<2>(item);
                    double possibility = std::get<3>(item);
                    
                    if(!IS_EPSILON(sym_C) || IS_TERMINATE(sym_B)) continue;
                    double alpha_B = ALPHA(sym_B, i, i);
                    
                    #ifdef COMPUTING_IN_LOG_SPACE
                    if(alpha_B + possibility > alpha_get(alpha, sym_A, i, i, MS)){
                        ALPHA(sym_A, i, i) = alpha_B + possibility;
                        nonterminate_remains.erase(sym_A);
                    }
                    #else
                    if(alpha_B * possibility > alpha_get(alpha, sym_A, i, i, MS)){
                        ALPHA(sym_A, i, i) = alpha_B * possibility;
                        nonterminate_remains.erase(sym_A);
                    }
                    #endif
                }
                
            }
            if(nonterminate_remains.size() == size) break;
        }
}

void kernel_inside_computeSpanKernel(const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
        #ifdef DEBUG_INSIDE_ALGORITHM
        , pcfg* grammar
        #endif
        ) {
        std::vector<double> buffer(N * MS * MS, 
            #ifdef COMPUTING_IN_LOG_SPACE
                -INFINITY
            #else
                0.0
            #endif
        );
        // reduce the buffer relocations.
        for (int span_length = 2; span_length <= sequence_length; span_length++) {
            std::fill(buffer.begin(), buffer.end(),
                #ifdef COMPUTING_IN_LOG_SPACE
                    -INFINITY
                #else
                    0.0
                #endif
            );

            #pragma omp parallel for
            for (int i = 0; i <= sequence_length - span_length; i++) {
                int j = i + span_length - 1; // Ending index of the span
                // std::vector<double> local_buffer(N, 
                //     #ifdef COMPUTING_IN_LOG_SPACE
                //         -INFINITY
                //     #else
                //         0.0
                //     #endif
                // );
                for (int k = i; k < j; k++) {
                    // iterate all grammars
                    for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                            PCFGItemIterator(N, grammar_index, grammar_table)){
                        uint32_t sym_A = std::get<0>(item);
                        uint32_t sym_B = std::get<1>(item);
                        uint32_t sym_C = std::get<2>(item);
                        double possibility = std::get<3>(item);
                        uint32_t gid = std::get<4>(item);

                        if(IS_EPSILON(sym_C)) continue;

                        // A->BC
                        double alpha_B = ALPHA_GET(sym_B, i, k);
                        double alpha_C = ALPHA_GET(sym_C, k + 1, j);
                        #ifdef COMPUTING_IN_LOG_SPACE
                            LOG_SUM_EXP_SET(buffer[sym_A * MS + i], alpha_B + alpha_C + possibility);
                        #else
                            buffer[sym_A * MS + i] += alpha_B * alpha_C * possibility; 
                        #endif   
                    }

                    // A->B
                    for(std::tuple<uint32_t, uint32_t>& rule_id: inside_order_1_rule_iteration_path) {
                        uint32_t gid = std::get<0>(rule_id);
                        uint32_t sym_A = std::get<1>(rule_id);
                        uint32_t* addr = (grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                        uint32_t sym_B = ((*addr) >> 16) & 0xFFFF;
                        double alpha_B = POSSIBILITY_ADD(buffer[sym_B * MS + i], ALPHA_GET(sym_B, i, j)); // Cell [i, j] is newly calculated and has not been written back. We need access the buffer, rather than ALPHA(..)
                        double possibility = *(double*)(addr + 1);                        

                        #ifdef COMPUTING_IN_LOG_SPACE
                            LOG_SUM_EXP_SET(buffer[sym_A * MS + i], alpha_B + possibility);
                        #else
                            buffer[sym_A * MS + i] += alpha_B * possibility;
                        #endif
                    }
                }
                
                // write back.
                for (int sym_A = 0; sym_A < N; sym_A++) {
                    #ifdef COMPUTING_IN_LOG_SPACE
                        LOG_SUM_EXP_SET(ALPHA(sym_A, i, j), buffer[sym_A * MS + i]);
                    #else
                        ALPHA_INCREASE(sym_A, i, j, buffer[sym_A * MS + i]);
                    #endif
                }
            } // parallel for end.
        }
}
