#ifndef H_MATH
#define H_MATH
#include<algorithm>
#include<cmath>

inline double log_sum_exp(double log_1, double log_2){
    return std::max(log_1, log_2) + std::log(1 + std::exp(-std::abs(log_1 - log_2)));
}
inline double add_logs(double log_1, double log_2){
    return log_1 + log_2;
}
#endif