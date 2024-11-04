#ifndef H_DEBUG
#define H_DEBUG
#include <string>
#include <iostream>
template<typename T, typename... Args>
void _debug_chunk(std::string description, Args... args) {
    std::cout << description << std::endl;
}
#endif
