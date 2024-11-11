
#include "kernels/update_parameters.cuh"
#include "utils/data_encoding.h"
#include "constants.h"
#include "utils/math.hpp"

__device__ double _calculate_new_possibility(double S, double f_gid) {
    if(std::abs(f_gid) < std::log(grammar_minimal_possibility))
        f_gid = std::log(grammar_minimal_possibility);
    return f_gid - S;
}

__global__ void kernel_update_parameters(double* f, double* count, double* mu, double* beta,
        const uint32_t* sequence, 
        uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, 
        uint32_t*
        grammar_table, double* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        , bool do_update
){
        int thread_id_x = blockIdx.x * blockDim.x + threadIdx.x;
        int total_threads_x = blockDim.x * gridDim.x;
        
        for(int sym_A = thread_id_x; sym_A < N; sym_A += total_threads_x){
            uint32_t grammar_pointer_current = *(grammar_index + sym_A);
            uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
            
            for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; 
                pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS
            ){
                uint32_t* addr = (grammar_table + pt);
                uint32_t symbols = *addr;
                double possibility = *(double*)(addr + 1);
                uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                uint32_t sym_C = (symbols) & 0xFFFF;
                int gid = pt / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
                LOG_SUM_EXP_SET(f[gid], count[gid]);
            }
        }

        __syncthreads();
        if(do_update){
            for(int sym_A = thread_id_x; sym_A < N; sym_A += total_threads_x){
                    double S = -INFINITY;
                    uint32_t grammar_pointer_current = *(grammar_index + sym_A);
                    uint32_t grammar_pointer_next = *(grammar_index + sym_A + 1);
                    for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; 
                        pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS
                    ){
                        uint32_t symbols = grammar_table[pt];
                        double possibility = *(double*)(grammar_table + pt + 1);
                        uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                        uint32_t sym_C = symbols & 0xFFFF;
                        int gid = pt / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
                        double f_gid = f[gid];
                        LOG_SUM_EXP_SET(S, 
                            (std::abs(f_gid - 0) < std::log(grammar_minimal_possibility) ? std::log(grammar_minimal_possibility) : f_gid));
                    
                    }

                    grammar_pointer_current = *(grammar_index + sym_A);
                    grammar_pointer_next = *(grammar_index + sym_A + 1);
                    for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; 
                        pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS
                    ){
                        uint32_t symbols = grammar_table[pt];
                        double possibility = *(double*)(grammar_table + pt + 1);
                        uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                        uint32_t sym_C = symbols & 0xFFFF;
                        int gid = pt / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
                        double f_gid = f[gid];

                        double new_possibility = _calculate_new_possibility(S,  
                        (std::abs(f_gid - 0) < std::log(grammar_minimal_possibility) ? std::log(grammar_minimal_possibility) : f_gid));
                        
                        *(double*)(grammar_table + pt + 1) = new_possibility;
                        
                        if(IS_EPSILON(sym_C) && IS_TERMINATE(sym_B)){
                            uint64_t key = encode_key(sym_A, sym_B);
                            reverse_grammar_hashtable_set_value(
                                pretermination_lookuptable, n_grammars * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS, key, new_possibility);
                        }
                    }
            }    
        }
        
    }
