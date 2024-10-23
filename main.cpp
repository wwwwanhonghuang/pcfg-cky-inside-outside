#ifdef USE_CUDA
#include <cuda_runtime.h>
// #include <cutensor.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <bits/stdc++.h>
#include <vector>
#include <iostream>

#include "macros.def"
#include "utils/tensor.hpp"
#include "kernels/inside.cuh"
#include "kernels/outside.cuh"
#include "kernels/expect_count.cuh"
#include "kernels/update_parameters.cuh"
#include "grammar/grammar.hpp"
#include "grammar/grammar_parser.hpp"
#include "utils/printer.hpp"
#include "utils/application_io.hpp"


#define PRINT_INSIDE 0
#define PRINT_OUTSIDE 0
#define PRINT_STEPS 0
#define PRINT_GRAMMAR_EACH_UPDATION_BEFORE 0
#define PRINT_GRAMMAR_EACH_UPDATION_AFTER 1

#define SANITARY_OUTPUT 1

#if SANITARY_OUTPUT == 1
#undef PRINT_INSIDE
#undef PRINT_OUTSIDE
#undef PRINT_STEPS
#undef PRINT_GRAMMAR_EACH_UPDATION_BEFORE
#undef PRINT_GRAMMAR_EACH_UPDATION_AFTER
#endif

void progress_bar(int progress, int total, int barWidth = 50) {
    float percentage = (float) progress / total;
    int pos = (int)(barWidth * percentage);
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) {
            std::cout << "=";
        } else if (i == pos) {
            std::cout << ">";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << int(percentage * 100.0) << " %  " << progress << "/" << total << "\r";
    std::cout.flush();  // Ensures the line is updated in place
}

float* outside_algorithm(float* mu, float* beta, uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        ,pcfg* grammar
                        #endif
                    ){
    #ifdef USE_CUDA
    <<<16, 16>>>
    #endif
    kernel_outside_main(mu, beta, sequence, pretermination_lookuptable,
        grammar_index, grammar_table, alpha, sequence_length, n_syms, N, T, MS, n_grammars, inside_order_1_rule_iteration_path, grammar);
    return beta;
}

float* em_algorithm_calculate_expection_count(float* count, float* mu, float* beta, uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        ,pcfg* grammar
                        #endif
                        ){
    #ifdef USE_CUDA
    <<<16, 16>>>
    #endif
    kernel_expect_count(count, mu, beta, sequence, pretermination_lookuptable,
        grammar_index, grammar_table, alpha, sequence_length, n_syms, N, T, MS, n_grammars
        #ifdef DEBUG_INSIDE_ALGORITHM
            ,  grammar
        #endif
    );
    return count;
}

float* inside_algorithm(uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars,
                        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path, pcfg* grammar = nullptr){
    
    if(n_syms >= 65536) return nullptr;

    // 1. zerolization alpha.
    kernel_inside_alpha_zerolization
        #ifdef USE_CUDA
            <<<16, 16>>>>
        #endif
    (alpha, N, sequence_length);

    // 2. fill alpha (base case).
    kernel_inside_base_fill_alpha
        #ifdef USE_CUDA
            <<<16, 16>>>>
        #endif
    (sequence, pretermination_lookuptable, grammar_index, grammar_table, alpha, 
                        sequence_length, n_syms, N, T, MS, n_grammars 
                        #ifdef DEBUG_INSIDE_ALGORITHM
                            ,grammar
                        #endif
    );
    // 3. fill alpha (recursive case).
    kernel_inside_computeSpanKernel
        #ifdef USE_CUDA
            <<<16, 16>>>>
        #endif
    (sequence, pretermination_lookuptable, grammar_index, grammar_table, alpha, 
                        sequence_length, n_syms, N, T, MS, n_grammars,
                        inside_order_1_rule_iteration_path
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , grammar
                        #endif
                        );
    
    return alpha;
};


int main(int argc, char* argv[])
{
    std::string grammar_filename = argc > 1 ? std::string(argv[1]) : "grammar_demo_2.pcfg";
    std::string input_filename = argc > 2 ? std::string(argv[2]) : "sequence.txt";
    
    pcfg* grammar = prepare_grammar(grammar_filename);

    float* alpha = new float[grammar->N() * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    float* beta = new float[(grammar->N() + grammar->T()) * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    float* mu = new float[grammar->cnt_grammar * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    float* count = new float[grammar->cnt_grammar]();
    float* f = new float[grammar->cnt_grammar]();

    std::cout << "Load sentences..." << std::endl;
    std::vector<std::vector<uint32_t>> sentences = parse_input_file(input_filename, grammar);
    std::cout << "Load sentences finished. Total instances:" << sentences.size() << std::endl;
    auto inside_order_1_rule_iteration_path = generate_inside_perterminate_iteration_paths(grammar);

    if(sentences.empty()) return 0;
    
    int sentence_id = 0;
    for(auto& sentence: sentences){
        progress_bar(sentence_id + 1, sentences.size());
        
        #if PRINT_GRAMMAR_EACH_UPDATION_BEFORE == 1
        std::cout << "grammar before iteration: " << sentence_id << " :" << std::endl;
        print_grammar(grammar);
        #endif
        int N = grammar->N();
        // std::cout << " -- proceed sentence: ";
        // for(auto&& word_id : sentence){
        //     std::cout << word_id << " ";
        // }
        // std::cout << std::endl;
        uint32_t* sequence = sentence.data();
        int sequence_length = sentence.size();

        #if PRINT_STEPS == 1
        std::cout << "1. Proceeding Inside Algorithm..." << std::endl;
        #endif
        
        inside_algorithm(sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (uint32_t*)(grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar,
            inside_order_1_rule_iteration_path
            #ifdef DEBUG_INSIDE_ALGORITHM
            , grammar
            #endif
        );

        #if PRINT_STEPS == 1
        std::cout << "Inside Algorithm Finished." << std::endl;
        #endif

        cky_printer printer;
        #if PRINT_INSIDE == 1
        printer.print_inside_outside_table(alpha,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);
        #endif

        #if PRINT_STEPS == 1
        std::cout << "2. Proceeding Outside Algorithm..." << std::endl;
        #endif

        outside_algorithm(mu, beta, sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (uint32_t*)(grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar,
            inside_order_1_rule_iteration_path
            #ifdef DEBUG_INSIDE_ALGORITHM
            ,grammar
            #endif
        );
        
        #if PRINT_STEPS == 1
        std::cout << "Outside Algorithm Finished." << std::endl;
        #endif

        #if PRINT_OUTSIDE == 1
        printer.print_inside_outside_table(beta,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);
        #endif

        #if PRINT_STEPS == 1
        std::cout << "3. Proceeding Calculate Expectation Count..." << std::endl;
        #endif

        kernel_expect_count(count, mu, beta, sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (uint32_t*)(grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, 
            grammar->cnt_grammar
            #ifdef DEBUG_INSIDE_ALGORITHM
            , grammar
            #endif
            );
        #if PRINT_STEPS == 1
        std::cout << "Calculate Expectation Count Finished." << std::endl;
        #endif
        #if PRINT_STEPS == 1
        std::cout << "4. Proceeding Update Parameters..." << std::endl;
        #endif

        kernel_update_parameters(f, count, mu, beta, sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(),
             MAX_SEQUENCE_LENGTH, grammar->cnt_grammar
             #ifdef DEBUG_INSIDE_ALGORITHM
            , grammar
            #endif
            );
        #if PRINT_STEPS == 1
        std::cout << "Update Parameter Finished." << std::endl;
        #endif

        #if PRINT_GRAMMAR_EACH_UPDATION_AFTER == 1
        std::cout << "grammar after iteration: " << sentence_id << " :" << std::endl;

        print_grammar(grammar);
        #endif
        sentence_id++;
    }

    std::cout << std::endl << "All finished" << std::endl;
    print_grammar(grammar);
    delete[] alpha;
    delete[] beta;
    delete[] mu;
    delete[] count;
    delete[] f;

    return 0;
}

