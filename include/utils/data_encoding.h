#ifndef H_DATA_ENCODING
#define H_DATA_ENCODING
#include <stdint.h>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
#ifdef USE_CUDA
__device__
#endif
inline uint64_t encode_key(int symbol_1, int symbol_2){
    return (((uint64_t)symbol_1 << 16) & 0xFFFF0000) | symbol_2 & (0x0000FFFF);
}
#endif