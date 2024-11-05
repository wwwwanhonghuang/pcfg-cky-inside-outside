#ifndef H_DATA_STRUCTURE
#define H_DATA_STRUCTURE
#include <map>
#include <cstdint>
#include <iostream>
#include <bits/stdc++.h>
#include "macros.def"
template<typename TKey, typename TVal>
void map_insert(std::map<TKey, TVal>& map, TKey key, TVal val) {
    map.insert(std::make_pair(key, val));
}

template<typename TKey, typename TVal>
bool map_contains(const std::map<TKey, TVal>& map, TKey key) {
    return map.find(key) != map.end();
}



/**
 * S: Max sequence length
 * length: Current input sequence length
 * alpha: inside possibility with shape (N, S, S)
**/
#ifdef USE_CUDA
__device__ 
#endif
__inline__ long double reverse_grammar_hashtable_get_value(
    uint32_t* hashtable, int hashtable_size, uint64_t key){
    #ifndef USE_CUDA
        int left = (key >> 16) & 0xFFFF;
        int right1 = key & 0xFFFF;
        int cnt_hashtable_items = hashtable_size / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
        int pos = (key % cnt_hashtable_items) * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
        int cnt_trails = 0;
        for(int i = 0; i < cnt_hashtable_items; i++){
            if(hashtable[pos] == key){
                return *((long double*)(hashtable + pos + 1));
            }else{
                pos = (pos + BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS) % hashtable_size;
            }
        }
        return 0.0L;
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
    uint32_t* hashtable, int hashtable_size, uint64_t key, long double val){
    #ifndef USE_CUDA
        int left = (key >> 16) & 0xFFFF;
        int right1 = key & 0xFFFF;
        int cnt_hashtable_items = hashtable_size / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
        int pos = (key % cnt_hashtable_items) * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
        int cnt_trails = 0;
        for(int i = 0; i < cnt_hashtable_items; i++){
            if(hashtable[pos] == key){
                *((long double*)(hashtable + pos + 1)) = val;
                return;
            }else{
                pos = (pos + BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS) % hashtable_size;
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

#endif