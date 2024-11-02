#ifndef H_DATA_STRUCTURE
#define H_DATA_STRUCTURE
#include <map>
#include <cstdint>
#include <iostream>
template<typename TKey, typename TVal>
void map_insert(std::map<TKey, TVal>& map, TKey key, TVal val) {
    map.insert(std::make_pair(key, val));
}

template<typename TKey, typename TVal>
bool map_contains(const std::map<TKey, TVal>& map, TKey key) {
    return map.find(key) != map.end();
}

union common_32bit {
    int32_t int32_value;
    float float32_value;

    common_32bit() : int32_value(0) {} // Default constructor
};


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
                return;
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

#endif