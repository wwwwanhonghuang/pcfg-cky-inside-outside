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
#include "dataset/dataset_helper.hpp"




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
    (alpha, N, MS);

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
    uint32_t log_itervals = argc > 3 ?  std::atoi(std::string(argv[3]).c_str()) : 10000; // 0xFFFFFFFF;
    
    pcfg* grammar = prepare_grammar(grammar_filename);
    auto inside_order_1_rule_iteration_path = generate_inside_perterminate_iteration_paths(grammar);

    float* alpha = new float[grammar->N() * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    float* beta = new float[(grammar->N() + grammar->T()) * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    float* mu = new float[grammar->cnt_grammar * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    float* count = new float[grammar->cnt_grammar]();
    double* f = new double[grammar->cnt_grammar]();

    std::cout << "Load sentences..." << std::endl;
    std::vector<std::vector<uint32_t>> sentences = parse_input_file(input_filename, grammar);
    std::cout << "Load sentences finished. Total instances:" << sentences.size() << std::endl;

    std::vector<std::vector<uint32_t>> train_set;
    std::vector<std::vector<uint32_t>> valid_set;
    double train_fraction = 0.8;
    split_dataset(sentences, train_set, valid_set, train_fraction);
    save_data_set_to_file("train_sentences.txt", train_set, grammar);
    save_data_set_to_file("validate_sentences.txt", valid_set, grammar);

    if(sentences.empty()) return 0;
    
    int sentence_id = 0;
    int n_sequences = sentences.size();
    int n_sequences_train = train_set.size();
    int n_sequences_val = valid_set.size();
    std::cout << "train set size = " << n_sequences_train << std::endl;
    std::cout << "val set size = " << n_sequences_val << std::endl;

    int T = 0;
    int n_epochs = 5;
    cky_printer printer;
    for(int epoch = 0; epoch < n_epochs; epoch++){
        // trainning
        for(int i = 0; i < n_sequences_train; i++){
            auto& sentence = train_set[i];
            progress_bar(i + 1, n_sequences_train);
            
            #if PRINT_GRAMMAR_EACH_UPDATION_BEFORE == 1
                std::cout << "grammar before iteration: " << sentence_id << " :" << std::endl;
                print_grammar(grammar);
            #endif
            int N = grammar->N();
            
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
            if((i + 1) == log_itervals){
                std::ofstream logfile_ostream = std::ofstream("./logs/log_" + std::to_string(i + 1) + "_epoch_id_" + std::to_string(epoch) + ".pcfg");
                if(!logfile_ostream){
                    std::cerr << "Error: Could not open log file for writing.\n";
                }
                print_grammar(grammar, logfile_ostream);
            }
            
        }

        // clear memory f at the end of an epoch.
        memset(f, 0, grammar->cnt_grammar * sizeof(double));
        std::ofstream logfile_ostream = std::ofstream("./logs/log_epoch_id_" + std::to_string(epoch) + ".pcfg");
        if(!logfile_ostream){
            std::cerr << "Error: Could not open log file for writing.\n";
        }
        print_grammar(grammar, logfile_ostream);
        
        
        // validation
        double log_likelihood = 0.0;
        for(int i = 0; i < n_sequences_val; i++){
            auto& sentence = valid_set[i];
            progress_bar(i + 1, n_sequences_val);
            
            int N = grammar->N();
            
            uint32_t* sequence = sentence.data();
            int sequence_length = sentence.size();
            
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
            if(alpha[0 + 0 + sequence_length - 1] >= 1){
                std::cout << alpha[0 + 0 + sequence_length - 1] << std::endl;
                print_grammar(grammar);
                for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                        PCFGItemIterator(grammar->N(), (uint32_t*)grammar->grammar_index, (uint32_t*)grammar->grammar_table)){
                    float p = std::get<3>(item);
                    assert(p <= 1);
                    std::cout << "gid: " << std::get<4>(item) << ":" << p << std::endl;
                }
                
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
                    
                printer.print_inside_outside_table(alpha,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);


                assert(alpha[0 + 0 + sequence_length - 1] <= 1);
            }
            
            
            log_likelihood += std::log(alpha[0 + 0 + sequence_length - 1]);            
        }
        double average_likelihood = (double)log_likelihood / (double)n_sequences_val;
        std::cout << "Average likelihood on validate set at epoch " << epoch << " = " ;
        std::cout << std::fixed << std::setprecision(9) << (double)(average_likelihood);
        std::cout << "  " << std::endl; 
        
    }
    

    std::cout << std::endl << "All finished" << std::endl;
    print_grammar(grammar);
    
    std::ofstream logfile_ostream = std::ofstream("./logs/log_final_" + std::to_string(sentences.size()) + ".pcfg");
    print_grammar(grammar, logfile_ostream);
    delete[] alpha;
    delete[] beta;
    delete[] mu;
    delete[] count;
    delete[] f;

    return 0;
}

