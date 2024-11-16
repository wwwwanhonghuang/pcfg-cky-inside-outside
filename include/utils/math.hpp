#ifndef H_MATH
#define H_MATH
#include<algorithm>
#include<cmath>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
#ifdef USE_CUDA
template <typename T>
__device__ T device_max(T a, T b) {
    return (a > b) ? a : b;
}
#endif

#ifdef USE_CUDA
__device__
#endif
inline double log_sum_exp(double log_1, double log_2){
    if (log_1 == -INFINITY) return log_2;
    if (log_2 == -INFINITY) return log_1;
    #ifdef USE_CUDA
        return device_max(log_1, log_2) + std::log(1 + std::exp(-std::abs(log_1 - log_2)));
    #else
        return std::max(log_1, log_2) + std::log(1 + std::exp(-std::abs(log_1 - log_2)));
    #endif
}

inline double host_log_sum_exp(double log_1, double log_2){
    if (log_1 == -INFINITY) return log_2;
    if (log_2 == -INFINITY) return log_1;
    return std::max(log_1, log_2) + std::log(1 + std::exp(-std::abs(log_1 - log_2)));
}
#ifdef USE_CUDA
__device__
#endif
inline double add_logs(double log_1, double log_2){
    return log_1 + log_2;
}
#endif