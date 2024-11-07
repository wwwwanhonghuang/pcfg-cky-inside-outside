#ifndef CUH_OUTSIDE
#define CUH_OUTSIDE
#include <omp.h>
#include <cstring>


#include <cstdint>
#include <vector>

#include "utils/data_accessing.hpp"
#include "grammar/grammar.hpp"
#include "macros.def"
#define DEBUG_OUTSIDE_CELL(x, y, X) if(i == x && j == y) { X }

void kernel_outside_main(double* mu, double* beta, 
                        const uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , pcfg* grammar
                        #endif
);
#endif