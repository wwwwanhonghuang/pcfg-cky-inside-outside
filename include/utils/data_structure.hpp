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
inline double reverse_grammar_hashtable_get_value(uint32_t* hashtable, int hashtable_size, uint64_t key){
    
    int left = (key >> 16) & 0xFFFF;
    int right1 = key & 0xFFFF;
    int cnt_hashtable_items = hashtable_size / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
    int pos = (key % cnt_hashtable_items) * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
    int cnt_trails = 0;
    for(int i = 0; i < cnt_hashtable_items; i++){
        if(hashtable[pos] == key){
            return *((double*)(hashtable + pos + 1));
        }else{
            pos = (pos + BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS) % hashtable_size;
        }
    }
    return 
        #ifdef COMPUTING_IN_LOG_SPACE
        -INFINITY
        #else
        0.0
        #endif
    ;    
};

inline void reverse_grammar_hashtable_set_value(uint32_t* hashtable, int hashtable_size, uint64_t key, double val){
    int left = (key >> 16) & 0xFFFF;
    int right1 = key & 0xFFFF;
    int cnt_hashtable_items = hashtable_size / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
    int pos = (key % cnt_hashtable_items) * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
    int cnt_trails = 0;
    for(int i = 0; i < cnt_hashtable_items; i++){
        if(hashtable[pos] == key){
            *((double*)(hashtable + pos + 1)) = val;
            return;
        }else{
            pos = (pos + BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS) % hashtable_size;
        }
    }
    std::cerr << "cannot find key " << key << " in reversed_grammar_hashtable" << std::endl;
    return;
};
#endif