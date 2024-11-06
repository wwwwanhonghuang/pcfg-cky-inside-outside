
#include "kernels/update_parameters.cuh"
#include "utils/data_encoding.h"
#include "constants.h"
#ifdef USE_CUDA
__global__
#endif

inline double _calculate_new_possibility(double S, double f_gid) {
    #if COMPUTING_IN_LOG_SPACE
    if(std::abs(f_gid) < std::log(grammar_minimal_possibility))
        f_gid = std::log(grammar_minimal_possibility);
    return f_gid - S;
    #else
    if(std::abs(f_gid) < grammar_minimal_possibility)
        f_gid = grammar_minimal_possibility;
    return f_gid / S;
    #endif
}

void kernel_update_parameters(double* f, double* count, double* mu, double* beta,
        const uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
    #ifdef USE_CUDA
        uint32_t* 
    #else
        uint32_t*
    #endif
        grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            , pcfg* grammar
        #endif
        , bool do_update
){
        int gid = 0;

        for(int sym_A = 0; sym_A < N; sym_A++){
            uint32_t grammar_pointer_current = *(grammar_index + sym_A);
            uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
            
            for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS){
                uint32_t symbols = grammar_table[pt];
                double possibility = *(double*)(grammar_table + pt + 1);
                uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                uint32_t sym_C = symbols & 0xFFFF;
                
                #ifdef COMPUTING_IN_LOG_SPACE
                f[gid] = log_sum_exp(f[gid], count[gid]);
                #else
                f[gid] += count[gid];
                #endif

                gid++;
            }
        }

        if(do_update){
            std::cout << "STATUS: parameter update." << std::endl;
            gid = 0;
            for(int sym_A = 0; sym_A < N; sym_A++){
                    double S = 0.0;
                    uint32_t grammar_pointer_current = *(grammar_index + sym_A);
                    uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
                    int gid_begin = gid;
                    for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS){
                        uint32_t symbols = grammar_table[pt];
                        double possibility = *(double*)(grammar_table + pt + 1);
                        uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                        uint32_t sym_C = symbols & 0xFFFF;
                        double f_gid = f[gid];
                        #ifdef COMPUTING_IN_LOG_SPACE
                        S = log_sum_exp(S, 
                            (std::abs(f_gid - 0) < std::log(grammar_minimal_possibility) ? std::log(grammar_minimal_possibility) : f_gid));
                        #else
                        S = (std::abs(f_gid - 0) < grammar_minimal_possibility ? grammar_minimal_possibility : f_gid);
                        #endif
                        gid ++;
                    }

                    grammar_pointer_current = *(grammar_index + sym_A);
                    grammar_pointer_next = *(grammar_index + sym_A + 1);
                    gid = gid_begin;
                    for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS){
                        uint32_t symbols = grammar_table[pt];
                        double possibility = *(double*)(grammar_table + pt + 1);
                        uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                        uint32_t sym_C = symbols & 0xFFFF;
                        double f_gid = f[gid];
                        
                        double new_possibility = _calculate_new_possibility(S,  
                        #ifdef COMPUTING_IN_LOG_SPACE
                            (std::abs(f_gid - 0) < std::log(grammar_minimal_possibility) ? std::log(grammar_minimal_possibility) : f_gid)
                        #else
                            (std::abs(f_gid - 0) < grammar_minimal_possibility ? grammar_minimal_possibility : f_gid)
                        #endif
                        );
                        
                        // if (new_possibility < -grammar_minimal_possibility || new_possibility > 1.0L + grammar_minimal_possibility) {
                        //     std::cout << "Improper possibility updation, possibility = " << new_possibility 
                        //             << ", caused by " << f[gid] << "/" << S << std::endl;
                        //     assert(false);
                        // }

                        *(double*)(grammar_table + pt + 1) = new_possibility;
                        
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
