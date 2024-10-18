
#include<omp.h>
#include "data_structure.hpp"

#ifndef USE_CUDA
#include <cstring>
#endif
#ifdef USE_CUDA
__global__
#endif
void kernel_update_parameters(float* f, float* count, float* mu, float* beta, uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, 
                        #ifdef USE_CUDA
                        uint32_t* 
                        #else
                        common_32bit*
                        #endif
                        grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars){
    memset(f, 0, n_grammars * sizeof(float));
    
    int gid = 0;
    #pragma omp parallel for
    for(int sym_A = 0; sym_A < N; sym_A++){
        uint32_t grammar_pointer_current = *(grammar_index + sym_A);
        uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
        
        for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += 2){
            uint32_t symbols = grammar_table[pt].int32_value;
            float possibility = grammar_table[pt + 1].float32_value;
            uint32_t sym_B = (symbols >> 16) & 0xFFFF;
            uint32_t sym_C = symbols & 0xFFFF;
            #pragma omp atomic
            f[gid] += count[gid];
            #pragma omp atomic
            gid++;
        }
    }

    gid = 0;
    #pragma omp parallel for
    for(int sym_A = 0; sym_A < N; sym_A++){
            float S = 0.0;
            uint32_t grammar_pointer_current = *(grammar_index + sym_A);
            uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);

            for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += 2){
                uint32_t symbols = grammar_table[pt].int32_value;
                float possibility = grammar_table[pt + 1].float32_value;
                uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                uint32_t sym_C = symbols & 0xFFFF;
                #pragma omp atomic
                S += f[gid];
            }

            grammar_pointer_current = *(grammar_index + sym_A);
            grammar_pointer_next = *(grammar_index + sym_A + 1);
            for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; pt += 2){
                uint32_t symbols = grammar_table[pt].int32_value;
                float possibility = grammar_table[pt + 1].float32_value;
                uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                uint32_t sym_C = symbols & 0xFFFF;
                float new_possibility = (S != 0 ? f[gid] / S : 0);;
            
                #pragma omp atomic
                *(float*)(grammar_table + pt + 1) += new_possibility;
                gid++;
            }
    }    
}