#ifndef H_ALG_INSIDE_OUTSIDE_MAIN
#define H_ALG_INSIDE_OUTSIDE_MAIN
#include <string>
#include <vector>
#include <tuple>
#include "kernels/inside.cuh"
#include "kernels/outside.cuh"
#include "kernels/update_parameters.cuh"
#include "kernels/expect_count.cuh"

long double* outside_algorithm(long double* mu, long double* beta, const uint32_t* sequence, 
    uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, long double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        ,pcfg* grammar
                        #endif
);


long double* em_algorithm_calculate_expection_count(long double* count, long double* mu, long double* beta, 
                        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, long double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        ,pcfg* grammar
                        #endif
);

long double* inside_algorithm(const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, long double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path, pcfg* grammar);

#endif