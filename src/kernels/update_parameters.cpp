#ifndef USE_CUDA

#include "kernels/update_parameters.cuh"
#include "utils/data_encoding.h"
#include "constants.h"
#include "utils/math.hpp"


inline double _calculate_new_possibility(double S, double f_gid) {
    if(std::abs(f_gid) < std::log(grammar_minimal_possibility))
        f_gid = std::log(grammar_minimal_possibility);
    return f_gid - S;
}

#ifndef ENABLE_GRAMMAR_VECTORIZATION_OPTIMIZATION
    #define SYMBOL_AND_POSSIBILITY_EXTRACTION(B, C, P) \
                    uint32_t symbols = grammar_table[pt]; \
                    double P = *(double*)(grammar_table + pt + 1); \
                    uint32_t B = (symbols >> 16) & 0xFFFF; \
                    uint32_t C = symbols & 0xFFFF;
#else
    #define SYMBOL_AND_POSSIBILITY_EXTRACTION(B, C, P) \
                    uint32_t symbols = grammar_table[(n_grammars + 1) * 0 + pt]; \
                    double P = *(double*)(grammar_table + (n_grammars + 1) * 4 + pt * 2); \
                    uint32_t B = grammar_table[(n_grammars + 1) * 1 + pt]; \
                    uint32_t C = grammar_table[(n_grammars + 1) * 2 + pt];
#endif

#ifndef ENABLE_GRAMMAR_VECTORIZATION_OPTIMIZATION
    #define PT_INCREASE pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS
#else
    #define PT_INCREASE pt++
#endif

void kernel_update_parameters(double* f, double* count, double* mu, double* beta,
        const uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
        uint32_t*
        grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef USE_CUDA
            , pcfg* grammar
        #endif
        , bool do_update
){
        int gid = 0;

        // accumulate count.
        for(int sym_A = 0; sym_A < N; sym_A++){
            uint32_t grammar_pointer_current = *(grammar_index + sym_A);
            uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
            for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; PT_INCREASE){        
                LOG_SUM_EXP_SET(f[gid], count[gid]);
                gid++;
            }
        }

        // update the parameter if need.
        if(do_update){
            std::cout << "STATUS: parameter update." << std::endl;
            gid = 0;
            for(int sym_A = 0; sym_A < N; sym_A++){
                    double S = INIT_POSSIBILITY;

                    uint32_t grammar_pointer_current = *(grammar_index + sym_A);
                    uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
                    int gid_begin = gid;
                    for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; PT_INCREASE){
                        double f_gid = f[gid];
                        LOG_SUM_EXP_SET(S, 
                                (std::abs(f_gid - 0) < std::log(grammar_minimal_possibility) ? std::log(grammar_minimal_possibility) : f_gid));
                        gid ++;
                    }

                    grammar_pointer_current = *(grammar_index + sym_A);
                    grammar_pointer_next = *(grammar_index + sym_A + 1);
                    gid = gid_begin;

                    for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; PT_INCREASE){
                        SYMBOL_AND_POSSIBILITY_EXTRACTION(sym_B, sym_C, possibility);
                        double f_gid = f[gid];
                        double new_possibility = 
                            _calculate_new_possibility(S, 
                                (std::abs(f_gid - 0) < std::log(grammar_minimal_possibility) ?
                                 std::log(grammar_minimal_possibility) : f_gid));
                        
                        #ifndef ENABLE_GRAMMAR_VECTORIZATION_OPTIMIZATION
                            *(double*)(grammar_table + pt + 1) = new_possibility;
                        #else
                            *(double*)(grammar_table + (n_grammars + 1) * 4 + pt * 2) = new_possibility; 
                        #endif

                        if(IS_EPSILON(sym_C) && IS_TERMINATE(sym_B)){
                            uint64_t key = encode_key(sym_A, sym_B);
                            reverse_grammar_hashtable_set_value(
                                pretermination_lookuptable, n_grammars * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS, key, new_possibility);
                        }
                        gid++;
                    }
            }    
        }
        
    }
#endif