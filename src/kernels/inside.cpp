#ifndef USE_CUDA
#include <math.h>
#include "kernels/inside.cuh"
#include "utils/math.hpp"
#include "constants.h"

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

void kernel_inside_base_fill_alpha(  
        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        , uint32_t* symbol_A_vector
        #ifdef DEBUG_INSIDE_ALGORITHM
        , pcfg* grammar
        #endif
        ) {
        #pragma omp parallel for
        for(int sym_A = 0; sym_A < N; sym_A ++){
            for(int i = 0; i < sequence_length; i++){
                double p = INIT_POSSIBILITY;
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

        for(int _i = 0; _i < N; _i++){
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
        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path, 
        uint32_t* symbol_A_vector,
        pcfg* grammar
        ) {
        std::vector<double> buffer(N * MS * MS, INIT_POSSIBILITY);
        // reduce the buffer relocations.
        for (int span_length = 2; span_length <= sequence_length; span_length++) {
            #pragma omp parallel
            {
                #pragma omp for
                for (int i = 0; i <= sequence_length - span_length; i++) {
                    for(int gid = 0; gid < n_grammars; gid++){
                        int j = i + span_length - 1; // Ending index of the span
                        uint32_t sym_A = grammar->symbol_A_vector[gid];
                        uint32_t sym_BC = *(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                        uint32_t sym_B = (sym_BC >> 16) & 0xFFFF;
                        uint32_t sym_C = (sym_BC) & 0xFFFF;
                        double possibility = *(double*)(grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1);

                        for (int k = i; k < j; k++) {
                            if(IS_EPSILON(sym_C)) continue;

                            // A->BC
                            double alpha_B = ALPHA_GET(sym_B, i, k);
                            double alpha_C = ALPHA_GET(sym_C, k + 1, j);
                            #ifdef COMPUTING_IN_LOG_SPACE
                                LOG_SUM_EXP_SET(buffer[sym_A * MS * MS + span_length * MS + i], alpha_B + alpha_C + possibility);
                            #else
                                buffer[sym_A * MS * MS + span_length * MS + i] += alpha_B * alpha_C * possibility; 
                            #endif   
                        }
                    }

                    //#pragma omp barrier

                    // A->B
                    for(std::tuple<uint32_t, uint32_t>& rule_id: inside_order_1_rule_iteration_path) {
                        int j = i + span_length - 1; // Ending index of the span
                        uint32_t gid = std::get<0>(rule_id);
                        uint32_t sym_A = std::get<1>(rule_id);

                        #ifndef ENABLE_GRAMMAR_VECTORIZATION_OPTIMIZATION
                            uint32_t* addr = (grammar_table + gid * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
                            uint32_t sym_B = ((*addr) >> 16) & 0xFFFF;
                            double possibility = *(double*)(addr + 1);
                        #else
                            uint32_t sym_B = grammar_table[(n_grammars + 1) * 1 + gid];
                            double possibility = *(double*)(grammar_table + (n_grammars + 1) * 4 + gid * 2);
                        #endif
                        
                        
                        double alpha_B = POSSIBILITY_ADD(buffer[sym_B * sequence_length * sequence_length + span_length * sequence_length + i], ALPHA_GET(sym_B, i, j)); // Cell [i, j] is newly calculated and has not been written back. We need access the buffer, rather than ALPHA(..)

                        #ifdef COMPUTING_IN_LOG_SPACE
                            LOG_SUM_EXP_SET(buffer[sym_A * MS * MS + span_length * MS + i], alpha_B + possibility);
                        #else
                            buffer[sym_A * MS * MS + span_length * MS + i] += alpha_B * possibility;
                        #endif
                    }
                    

                    //#pragma omp barrier

                    // write back.
                    for (int sym_A = 0; sym_A < N; sym_A++) {
                        int j = i + span_length - 1; // Ending index of the span
                        #ifdef COMPUTING_IN_LOG_SPACE
                            LOG_SUM_EXP_SET(ALPHA(sym_A, i, j), buffer[sym_A * MS * MS + span_length * MS + i]);
                        #else
                            ALPHA_INCREASE(sym_A, i, j, buffer[sym_A * MS * MS + span_length * MS + i]);
                        #endif
                    }
                }
            } // parallel for end.
        }
}
#endif
