

#include "kernels/expect_count.cuh"


#if DEBUG_MODE == 0
#define EXPERCT_COUNT_DEBUG_OUTPUT(CONDITION)
#else
#define EXPERCT_COUNT_DEBUG_OUTPUT(CONDITION) if(CONDITION > 0) \
                                            std::cout << "(" << i << ", " << j << ")" << \
                                            "expectation_count::" << \
                                            "gid = " << gid << " " << SYMBOL_STR(sym_A) << "->" << SYMBOL_STR(sym_B) \
                                            << " " <<  SYMBOL_STR(sym_C) << \
                                            " expectation count += " <<  mu_val << "(MU[" << gid << "," << i << ", " << j << "]" << std::endl;
#endif


#ifdef USE_CUDA
__global__
#endif
void kernel_expect_count(long double* count, long double* mu, long double* beta, const uint32_t* sequence, 
                        uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, long double* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , pcfg* grammar
                        #endif

    ){
    memset(count, 0, n_grammars * sizeof(long double));

    /* 0 is the id of S symbol. This expression assign alpha['S', 0, sequence_length - 1] to Z */
    long double Z = ALPHA(0, 0, sequence_length - 1); 

    for(int span_length = 1; span_length <= sequence_length; span_length++){
        
        #pragma omp parallel for
        for(int i = 0; i < sequence_length - span_length + 1; i++){
            int j = i + span_length - 1;
            for(std::tuple<uint32_t, uint32_t, uint32_t, long double, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                uint32_t sym_A = std::get<0>(item);
                uint32_t sym_B = std::get<1>(item);
                uint32_t sym_C = std::get<2>(item);
                long double possibility = std::get<3>(item);
                uint32_t gid = std::get<4>(item);
                long double mu_val = MU(gid, i, j);
                
                #pragma omp atomic
                count[gid] += mu_val;
            }
        }
    }
    
    #if PRINT_EXPECTATION_COUNT == 1
        for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item :
                PCFGItemIterator(N, grammar_index, grammar_table)){
            uint32_t sym_A = std::get<0>(item);
            uint32_t sym_B = std::get<1>(item);
            uint32_t sym_C = std::get<2>(item);
            float possibility = std::get<3>(item);
            uint32_t gid = std::get<4>(item);
            std::cout << 
            "print::expectation_count:" <<
                "gid = " << gid << " " << SYMBOL_STR(sym_A) << "->" << SYMBOL_STR(sym_B)
                << " " <<  SYMBOL_STR(sym_C) <<
                " expectation count = " <<  count[gid] << std::endl;
        }
    #endif
}
