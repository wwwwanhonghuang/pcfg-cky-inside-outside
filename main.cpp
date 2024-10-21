#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#ifdef USE_CUDA
#include <cuda_runtime.h>
// #include <cutensor.h>
#endif
#include <map>

#include <unordered_map>
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

float* outside_algorithm(float* mu, float* beta, uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        ,pcfg* pcfg
                        #endif
                        ){
    
    #ifdef USE_CUDA
    <<<16, 16>>>
    #endif
    kernel_outside_main(mu, beta, sequence, pretermination_lookuptable,
        grammar_index, grammar_table, alpha, sequence_length, n_syms, N, T, MS, n_grammars, pcfg);
    return beta;
    
}

float* em_algorithm_calculate_expection_count(float* count, float* mu, float* beta, uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars){
    #ifdef USE_CUDA
    <<<16, 16>>>
    #endif
    kernel_expect_count(count, mu, beta, sequence, pretermination_lookuptable,
        grammar_index, grammar_table, alpha, sequence_length, n_syms, N, T, MS, n_grammars);
    return count;
}

float* inside_algorithm(uint32_t* sequence, uint32_t* pretermination_lookuptable, 
                        uint32_t* grammar_index, uint32_t* grammar_table, float* alpha, 
                        int sequence_length, int n_syms, int N, int T, int MS, int n_grammars, pcfg* pcfg = nullptr){
    std::cout << "n_syms = " << n_syms << " N = " << N << " T = " << T << std::endl;
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
                            ,pcfg
                        #endif
    );
    // 3. fill alpha (recursive case).
    kernel_inside_computeSpanKernel
        #ifdef USE_CUDA
            <<<16, 16>>>>
        #endif
    (sequence, pretermination_lookuptable, grammar_index, grammar_table, alpha, 
                        sequence_length, n_syms, N, T, MS, n_grammars
                        #ifdef DEBUG_INSIDE_ALGORITHM
                        , pcfg
                        #endif
                        );
    
    return alpha;
};

std::vector<std::vector<uint32_t>> parse_input_file(const std::string& file_path, pcfg* grammar){
    std::vector<std::vector<uint32_t>> sentences;
    std::string line;
	std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the input file at path: " << file_path << std::endl;
        return sentences;
    }
    int N = grammar->N();

    while (std::getline(file, line)) {
        if(line == "")
            continue;
        std::vector<uint32_t> input_words;

        std::string word;
        std::stringstream line_string_stream(line);
        while (getline(line_string_stream, word, ' ')) {
            input_words.push_back(grammar->terminate_map.find(std::string("\'") + word + std::string("\'"))->second + N);
        }
        sentences.push_back(input_words);
    }

    return sentences;
}
int main(int argc, char* argv[])
{
    std::string grammar_filename = argc > 1 ? std::string(argv[1]) : "grammar_demo_2.pcfg";
    std::string input_filename = argc > 2 ? std::string(argv[2]) : "sequence.txt";
    
    pcfg* grammar = prepare_grammar(grammar_filename);

    float* alpha = new float[grammar->N() * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH];
    float* beta = new float[grammar->N() * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH];
    float* mu = new float[grammar->cnt_grammar * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH];
    float* count = new float[grammar->cnt_grammar];
    float* f = new float[grammar->cnt_grammar];

    std::vector<std::vector<uint32_t>> sentences = parse_input_file(input_filename, grammar);
    if(sentences.empty()) return 0;
    
    for(auto& sentence: sentences){
        int N = grammar->N();
        uint32_t* sequence = sentence.data();
        int sequence_length = sentence.size();

        std::cout << "1. Proceeding Inside Algorithm..." << std::endl;
        inside_algorithm(sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (uint32_t*)(grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar
            #ifdef DEBUG_INSIDE_ALGORITHM
            , grammar
            #endif
        );
        
        std::cout << "Inside Algorithm Finished." << std::endl;
        cky_printer printer;
        printer.print_inside_outside_table(alpha,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);

        std::cout << "2. Proceeding Outside Algorithm..." << std::endl;
        outside_algorithm(mu, beta, sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (uint32_t*)(grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar
            #ifdef DEBUG_INSIDE_ALGORITHM
            ,grammar
            #endif
        );
        
        std::cout << "Outside Algorithm Finished." << std::endl;
        printer.print_inside_outside_table(beta,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);

        std::cout << "3. Proceeding Calculate Expectation Count..." << std::endl;
        kernel_expect_count(count, mu, beta, sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (uint32_t*)(grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar);
        std::cout << "Calculate Expectation Count Finished." << std::endl;

        std::cout << "4. Proceeding Update Parameters..." << std::endl;
        kernel_update_parameters(f, count, mu, beta, sequence, 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar);
        std::cout << "Update Parameter Finished." << std::endl;
    }

    return 0;
}