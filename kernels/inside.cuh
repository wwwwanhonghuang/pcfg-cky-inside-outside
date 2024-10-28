#ifndef USE_CUDA
#include <cstring>
#include <tuple>
#endif
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#include "macros.def"
#ifdef DEBUG_INSIDE_ALGORITHM
#include "grammar/grammar.hpp"
#include "utils/data_encoding.h"
#include "utils/data_accessing.hpp"
#endif
/**
 * S: Max sequence length
 * length: Current input sequence length
 * alpha: inside possibility with shape (N, S, S)
**/
#ifdef USE_CUDA
__device__ 
#endif
__inline__ float reverse_grammar_hashtable_get_value(
    uint32_t* hashtable, int hashtable_size, uint64_t key){
    #ifndef USE_CUDA
        int left = (key >> 16) & 0xFFFF;
        int right1 = key & 0xFFFF;
        int cnt_hashtable_items = hashtable_size / 2;
        int pos = (key % cnt_hashtable_items) * 2;
        int cnt_trails = 0;
        for(int i = 0; i < cnt_hashtable_items; i++){
            if(hashtable[pos] == key){
                return ((float*) hashtable)[pos + 1];
            }else{
                pos = (pos + 2) % hashtable_size;
            }
        }
        return 0.0f;
    #else
    int pos = (key % cnt_hashtable_max_keys) * 2;
    int length = cnt_hashtable_max_keys * 2;
    for(int i = 0; i < cnt_hashtable_max_keys; i++){
        if(hashtable[pos] == key){
            return (float)hashtable[pos + 1];
        }else{
            pos = (pos + 2) % length;
        }
    }
    return 0.0f;
    #endif
};

#ifdef USE_CUDA
__device__ 
#endif
__inline__ void reverse_grammar_hashtable_set_value(
    uint32_t* hashtable, int hashtable_size, uint64_t key, float val){
    #ifndef USE_CUDA
        int left = (key >> 16) & 0xFFFF;
        int right1 = key & 0xFFFF;
        int cnt_hashtable_items = hashtable_size / 2;
        int pos = (key % cnt_hashtable_items) * 2;
        int cnt_trails = 0;
        for(int i = 0; i < cnt_hashtable_items; i++){
            if(hashtable[pos] == key){
                ((float*) hashtable)[pos + 1] = val;
            }else{
                pos = (pos + 2) % hashtable_size;
            }
        }
        std::cerr << "cannot find key " << key << " in reversed_grammar_hashtable" << std::endl;
        return;
    #else
    int pos = (key % cnt_hashtable_max_keys) * 2;
    int length = cnt_hashtable_max_keys * 2;
    for(int i = 0; i < cnt_hashtable_max_keys; i++){
        if(hashtable[pos] == key){
            return (float)hashtable[pos + 1];
        }else{
            pos = (pos + 2) % length;
        }
    }
    return 0.0f;
    #endif
};

#ifdef USE_CUDA
__global__
#endif
void kernel_inside_alpha_zerolization(float* alpha, int N, int MS){
    memset(alpha, 0, N * MS * MS * sizeof(float));
}

#ifdef USE_CUDA
__global__ 
#endif
void kernel_inside_base_fill_alpha(  
        uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
        , pcfg* grammar
        #endif
        ) {
    #ifndef USE_CUDA
        #pragma omp parallel for
        for(int sym_A = 0; sym_A < N; sym_A ++){
            for(int i = 0; i < sequence_length; i++){
                float p = 0.0;
                uint64_t key = encode_key(sym_A, sequence[i]);
                p = reverse_grammar_hashtable_get_value(pretermination_lookuptable, n_grammars * 2, key);
                
                if(abs(p - 0) < 1e-6) continue;               
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
                for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                                                        PCFGItemIterator(N, grammar_index, grammar_table)){
                    uint32_t sym_A = std::get<0>(item);
                    uint32_t sym_B = std::get<1>(item);
                    uint32_t sym_C = std::get<2>(item);
                    float possibility = std::get<3>(item);
                    
                    if(!IS_EPSILON(sym_C) || IS_TERMINATE(sym_B)) continue;
                    float alpha_B = ALPHA(sym_B, i, i);
                    
                    if(alpha_B * possibility > alpha_get(alpha, sym_A, i, i, MS)){
                        ALPHA(sym_A, i, i) = alpha_B * possibility;
                        nonterminate_remains.erase(sym_A);
                    }
                }
                
            }
            if(nonterminate_remains.size() == size) break;
        }
            
    #else
        int globalThreadIndex = blockIdx.x * blockDim.x + threadIdx.x;
        int threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
        int totalBlocks = gridDim.x * gridDim.y * gridDim.z;
        if (threadIdx.x == 0) {
            printf("Threads per Block: %d, Total Blocks: %d\n", threadsPerBlock, totalBlocks);
        }

        // computational compelxity in worst situation: O(S * N), compelxity in best situation O(S)
        // mapping strategy: one block one symbol, 
        //      i-th thread responsible for i, i + |Threads|, ..., i + 2|Threads|'s word in the sequence.
        int block_id = blockIdx.x;
        int thread_id = threadIdx.x;
        int B = ceilf(length / threadsPerBlock);
        for(int sym_id = block_id; sym_id < N; sym_id += totalBlocks){
            for(int i = thread_id; i < length && i < thread_id + B; i ++){
                float p = 0.0;
                if(preterminate_reversed_hashtable){
                    uint64_t key = ((uint64_t)sym_id << 32) | sequence[i];
                    p = reverse_grammar_hashtable_get_value(preterminate_reversed_hashtable, preterminate_reversed_hashtable_length, key);
                    alpha[sym_id * (S * S) + i * S + i] = p;
                }else{
                    _32bit_t* pt_begin = (_32bit_t*)&preterminate[sym_id];
                    _32bit_t* pt_end = (_32bit_t*)&preterminate[sym_id + 1];
                    
                    for(_32bit_t* pt = pt_begin; pt < pt_end; pt += 2){
                        if(*(pt) == sequence[i]){
                            p = *(pt + 1);
                            break;
                        }
                    }
                }
                
                // For computational effectiveness, we reused alpha space for each input. alpha is N x S x S tensor.
                alpha[sym_id * (S * S) + i * S + i] = p;
            }
        }
    #endif
}

#ifdef USE_CUDA
__global__ 
#endif
void kernel_inside_computeSpanKernel(uint32_t* sequence, uint32_t* pretermination_lookuptable, 
        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
        #ifdef DEBUG_INSIDE_ALGORITHM
        , pcfg* grammar
        #endif
        ) {
    #ifndef USE_CUDA                        
        for (int span_length = 2; span_length <= sequence_length; span_length++) {
            #pragma omp parallel for
            for (int i = 0; i <= sequence_length - span_length; i++) {
                int j = i + span_length - 1; // Ending index of the span
                for (int k = i; k < j; k++) {
                    // iterate all grammars
                    for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                        uint32_t sym_A = std::get<0>(item);
                        uint32_t sym_B = std::get<1>(item);
                        uint32_t sym_C = std::get<2>(item);
                        float possibility = std::get<3>(item);
                        uint32_t gid = std::get<4>(item);

                        if(IS_EPSILON(sym_C)) continue;

                        
                        // A->BC
                        float alpha_B = ALPHA_GET(sym_B, i, k);
                        float alpha_C = ALPHA_GET(sym_C, k + 1, j);
                        #pragma omp atomic
                        ALPHA_INCREASE(sym_A, i, j, alpha_B * alpha_C * possibility);   
                    }

                    // A->B
                    for(std::tuple<uint32_t, uint32_t>& rule_id: inside_order_1_rule_iteration_path) {
                        uint32_t gid = std::get<0>(rule_id);
                        uint32_t sym_A = std::get<1>(rule_id);
                        uint32_t* addr = (grammar_table + gid * 2);
                        uint32_t sym_B = ((*addr) >> 16) & 0xFFFF;
                        float alpha_B = ALPHA_GET(sym_B, i, j);
                        float possibility = *(float*)(addr + 1);
                        
                        #pragma omp atomic
                        ALPHA_INCREASE(sym_A, i, j, alpha_B * possibility);
                    }                        
                    
                }
            }
        }
    #else
        int globalThreadIndex = blockIdx.x * blockDim.x + threadIdx.x;
        
        // Loop over spans of length 2 to `length`
        // mapping strategy: for each span_length, conducting threads parallelly filling. 
        for (int span_length = 2; span_length <= length; span_length++) {
            for (int i = 0; i <= length - span_length; i++) {
                int j = i + span_length - 1; // Ending index of the span

                // Each thread handles multiple (i, j) pairs
                for (int k = i; k < j; k++) {
                    // iterate all grammars
                    float* grammar_pointer_current = *grammars;
                    float* grammar_pointer_next = *(grammars + 1);
                    int non_terminate_id = 0;
                    while(grammar_pointer_current != grammar_pointer_next){
                        for(float* pt = grammar_pointer_current; pt < grammar_pointer_next; pt += 2){
                            uint32_t symbols = *(uint32_t*)pt;
                            float possibility = *(float*)(pt + 1);
                            uint32_t sym_B = (symbols >> 16) & 0xFFFF;
                            uint32_t sym_C = symbols & 0xFFFF;
                            alpha[non_terminate_id * S * S + i * S + j] = 
                                alpha[sym_B * S * S + i * S + k] * alpha[sym_C * S * S + (k + 1) * S + j]
                                * possibility;
                        }
                        non_terminate_id += 1;
                    }
                }
            }
        }
    #endif
}
