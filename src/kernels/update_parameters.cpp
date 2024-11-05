
#include "kernels/update_parameters.cuh"
#include "utils/data_encoding.h"
#include "constants.h"
#ifdef USE_CUDA
__global__
#endif

inline long double _calculate_new_possibility(long double S, long double f_gid) {
    if(std::abs(f_gid) < grammar_minimal_possibility)
        f_gid = grammar_minimal_possibility;
    return f_gid / S;
}

void kernel_update_parameters(long double* f, long double* count, long double* mu, long double* beta,
        const uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
    #ifdef USE_CUDA
        uint32_t* 
    #else
        uint32_t*
    #endif
        grammar_table, long double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            , pcfg* grammar
        #endif
){
        int gid = 0;

        for(int sym_A = 0; sym_A < N; sym_A++){
            uint32_t grammar_pointer_current = *(grammar_index + sym_A);
            uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
            
            for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS){
                uint32_t symbols = grammar_table[pt];
                long double possibility = *(long double*)(grammar_table + pt + 1);
                uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                uint32_t sym_C = symbols & 0xFFFF;
                
                f[gid] += count[gid];

                gid++;
            }
        }

        gid = 0;
        for(int sym_A = 0; sym_A < N; sym_A++){
                double S = 0.0;
                uint32_t grammar_pointer_current = *(grammar_index + sym_A);
                uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
                int gid_begin = gid;
                for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS){
                    uint32_t symbols = grammar_table[pt];
                    long double possibility = *(long double*)(grammar_table + pt + 1);
                    uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                    uint32_t sym_C = symbols & 0xFFFF;
                    long double f_gid = f[gid];
                    S += std::abs(f_gid - 0) < grammar_minimal_possibility ? grammar_minimal_possibility : f_gid;
                    gid ++;
                }

                grammar_pointer_current = *(grammar_index + sym_A);
                grammar_pointer_next = *(grammar_index + sym_A + 1);
                gid = gid_begin;
                for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS){
                    uint32_t symbols = grammar_table[pt];
                    long double possibility = *(long double*)(grammar_table + pt + 1);
                    uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                    uint32_t sym_C = symbols & 0xFFFF;
                    
                    long double new_possibility = _calculate_new_possibility(S, f[gid]);
                    
                    if (new_possibility < -grammar_minimal_possibility || new_possibility > 1.0L + grammar_minimal_possibility) {
                        std::cout << "Improper possibility updation, possibility = " << new_possibility 
                                << ", caused by " << f[gid] << "/" << S << std::endl;
                        assert(false);
                    }

                    *(long double*)(grammar_table + pt + 1) = new_possibility;
                    
                    if(IS_EPSILON(sym_C) && IS_TERMINATE(sym_B)){
                        uint64_t key = encode_key(sym_A, sym_B);
                        reverse_grammar_hashtable_set_value(
                            pretermination_lookuptable, n_grammars * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS, key, new_possibility);
                    }
                    
                    gid++;
                }
        }    
    }
