#include "algorithms/alg_inside_outside_main.h"

double* outside_algorithm(double* mu, double* beta, const uint32_t* sequence, 
                        uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, 
                        uint32_t* grammar_table, 
                        double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        ,pcfg* grammar
                        #endif
){
    kernel_outside_main(mu, beta, sequence, pretermination_lookuptable,
        grammar_index, grammar_table, alpha, sequence_length, n_syms, N, T, MS, n_grammars, inside_order_1_rule_iteration_path, grammar);
    return beta;
}

double* em_algorithm_calculate_expection_count(double* count, double* mu, double* beta, const uint32_t* sequence, 
                        uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, 
                        uint32_t* grammar_table, 
                        double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        ,pcfg* grammar
                        #endif
                        ){
    kernel_expect_count(count, mu, beta, sequence, pretermination_lookuptable,
        grammar_index, grammar_table, alpha, sequence_length, n_syms, N, T, MS, n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            ,  grammar
        #endif
    );
    return count;
}

double* inside_algorithm(const uint32_t* sequence, 
    uint32_t* pretermination_lookuptable, 
    uint32_t* grammar_index, 
    uint32_t* grammar_table, 
    double* alpha, 
    int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
    std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path, pcfg* grammar = nullptr){

    if(n_syms >= 65536) return nullptr;

    // 1. zerolization alpha.
    kernel_inside_alpha_zerolization(alpha, N, MS);

    // 2. fill alpha (base case).
    kernel_inside_base_fill_alpha(sequence, pretermination_lookuptable, grammar_index, grammar_table, alpha, 
                        sequence_length, n_syms, N, T, MS, n_grammars 
                        #ifdef DEBUG_INSIDE_ALGORITHM
                            ,grammar
                        #endif
    );
    

    // 3. fill alpha (recursive case).
    kernel_inside_computeSpanKernel(sequence, pretermination_lookuptable, grammar_index, grammar_table, alpha, 
        sequence_length, n_syms, N, T, MS, n_grammars,
        inside_order_1_rule_iteration_path
        #ifdef DEBUG_INSIDE_ALGORITHM
            , grammar
        #endif
        );

    return alpha;
};
