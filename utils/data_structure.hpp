#ifndef H_DATA_STRUCTURE
#define H_DATA_STRUCTURE
#include <map>
#include <cstdint>
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

#endif