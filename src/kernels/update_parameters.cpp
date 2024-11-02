
#include "kernels/update_parameters.cuh"
#include "utils/data_encoding.h"
#ifdef USE_CUDA
__global__
#endif

inline float _calculate_new_possibility(float S, float f_gid) {
    if(abs(f_gid) < epsilon)
        f_gid = epsilon;
    return f_gid / S;
}

void kernel_update_parameters(double* f, float* count, float* mu, float* beta, uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
    #ifdef USE_CUDA
        uint32_t* 
    #else
        common_32bit*
    #endif
        grammar_table, float* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            , pcfg* grammar
        #endif
){
        int gid = 0;

        for(int sym_A = 0; sym_A < N; sym_A++){
            uint32_t grammar_pointer_current = *(grammar_index + sym_A);
            uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
            
            for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += 2){
                uint32_t symbols = grammar_table[pt].int32_value;
                float possibility = grammar_table[pt + 1].float32_value;
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
                for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += 2){
                    uint32_t symbols = grammar_table[pt].int32_value;
                    float possibility = grammar_table[pt + 1].float32_value;
                    uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                    uint32_t sym_C = symbols & 0xFFFF;
                    float f_gid = f[gid];
                    S += abs(f_gid - 0) < epsilon ? epsilon : f_gid;
                    gid ++;
                }

                grammar_pointer_current = *(grammar_index + sym_A);
                grammar_pointer_next = *(grammar_index + sym_A + 1);
                gid = gid_begin;
                for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += 2){
                    uint32_t symbols = grammar_table[pt].int32_value;
                    float possibility = grammar_table[pt + 1].float32_value;
                    uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                    uint32_t sym_C = symbols & 0xFFFF;
                    
                    float new_possibility = _calculate_new_possibility(S, f[gid]);
                    
                    if(!((new_possibility + epsilon) >= 0.0f && new_possibility - epsilon <= 1.0f)){
                        std::cout << "Improper possibility updation, possibility = " << new_possibility << ", caused by " <<
                        f[gid] << "/" << S 
                        << std::endl;
                        assert(false);
                    }
                    *(float*)(grammar_table + pt + 1) = new_possibility;
                    
                    if(IS_EPSILON(sym_C) && IS_TERMINATE(sym_B)){
                        uint64_t key = encode_key(sym_A, sym_B);
                        reverse_grammar_hashtable_set_value(
                            pretermination_lookuptable, n_grammars * 2, key, new_possibility);
                    }
                    
                    gid++;
                }
        }    
    }
