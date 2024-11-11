#ifndef H_ALG_INSIDE_OUTSIDE_MAIN
#define H_ALG_INSIDE_OUTSIDE_MAIN
#include <string>
#include <vector>
#include <tuple>
#include "kernels/inside.cuh"
#include "kernels/outside.cuh"
#include "kernels/update_parameters.cuh"
#include "kernels/expect_count.cuh"

double* outside_algorithm(double* mu, double* beta, const uint32_t* sequence, 
    uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        #ifndef USE_CUDA
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #else
                        uint32_t* inside_order_1_rule_iteration_path, uint32_t inside_order_1_rule_iteration_path_size
                        #endif
                        #ifndef USE_CUDA
                        ,pcfg* grammar
                        #endif
);


double* em_algorithm_calculate_expection_count(double* count, double* mu, double* beta, 
                        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        #ifndef USE_CUDA
                        ,pcfg* grammar
                        #endif
);

double* inside_algorithm(const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #ifndef USE_CUDA
                        , pcfg* grammar
                        #endif
    );

#endif