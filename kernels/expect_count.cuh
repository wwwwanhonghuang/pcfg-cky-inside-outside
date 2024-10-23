#ifndef USE_CUDA
#include <cstring>
#endif
#include "utils/data_accessing.hpp"
#include "grammar/grammar.hpp"
#include "macros.def"

#ifdef USE_CUDA
__global__
#endif

void kernel_expect_count(float* count, float* mu, float* beta, uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , pcfg* grammar
                        #endif

                        ){
    memset(count, 0, n_grammars * sizeof(float));
    float Z = ALPHA(0, 0, sequence_length - 1); // 0 is the id of S symbol. This expression assign alpha['S', 0, sequence_length - 1] to Z;

    for(int span_length = 1; span_length <= sequence_length; span_length++){
        //#pragma omp parallel for
        for(int i = 0; i < sequence_length - span_length + 1; i++){
            int j = i + span_length - 1;
            for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : PCFGItemIterator(N, grammar_index, grammar_table)){
                uint32_t sym_A = std::get<0>(item);
                uint32_t sym_B = std::get<1>(item);
                uint32_t sym_C = std::get<2>(item);
                float possibility = std::get<3>(item);
                uint32_t gid = std::get<4>(item);
                // #pragma omp atomic
                // if((gid >= 2 && gid <=6) || gid == 10 || gid == 11){
                //         std::cout << "count[" << gid << "] += MU(" << gid << "," << i << "," << j << " == " << MU(gid, i, j) << ")" << std::endl;
                // }
                count[gid] += MU(gid, i, j);
            }
        }
    }
}