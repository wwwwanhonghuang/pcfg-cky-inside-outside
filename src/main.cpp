#ifndef USE_CUDA
#include "main.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <cmath>

#include "macros.def"
#include "utils/tensor.hpp"
#include "algorithms/alg_inside_outside_main.h"
#include "grammar/grammar.hpp"
#include "grammar/grammar_parser.hpp"
#include "utils/printer.hpp"
#include "utils/application_io.hpp"
#include "dataset/dataset_helper.hpp"
#include "constants.h"
#include "utils/math.hpp"

int main(int argc, char* argv[])
{
    std::string configuration_file_path = "config.yaml";
    if(argc > 1){
        configuration_file_path = std::string(argv[1]);
    }
    std::cout << "read configuration file " << configuration_file_path << std::endl;
    YAML::Node config = YAML::LoadFile(configuration_file_path);
    
    if (!config.IsDefined()) {
        std::cout << "Error: config.yaml could not be loaded!" << std::endl;
        return 1;
    }

    // 1. load configuration.
    std::string grammar_filename = config["main"]["grammar_file"].as<std::string>();
    std::string input_filename = config["main"]["input"].as<std::string>();
    uint32_t log_intervals = config["main"]["log_intervals"].as<int>();
    std::string log_path = config["main"]["log_path"].as<std::string>();
    int n_epochs = config["main"]["n_epochs"].as<int>();
    int batch_size_for_parameter_update = config["main"]["batch_size_for_parameter_update"].as<int>();
    bool save_f = config["main"]["log_f"]["enabled"].as<bool>();
    int log_f_intervals = -1;
    int limit_n_sentences = config["main"]["limit_n_sentences"].as<int>();
    if(save_f){
        log_f_intervals = config["main"]["log_f"]["intervals"].as<int>();
    }
    bool log_warning_in_training = config["main"]["log_warning_in_training"].as<bool>();

    // 2. parse grammar file.
    pcfg* grammar = nullptr;
    try {
        grammar = prepare_grammar(grammar_filename);
        if (grammar == nullptr) {
            throw std::runtime_error("Error: Failed to parse grammar file.");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    auto inside_order_1_rule_iteration_path = generate_inside_perterminate_iteration_paths(grammar);

    // 3. define matrices needed by the inside-outside algorithm.
    double* alpha = new double[grammar->N() * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    double* beta = new double[(grammar->N() + grammar->T()) * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    double* mu = new double[grammar->cnt_grammar * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();
    double* count = new double[grammar->cnt_grammar]();
    double* f = new double[grammar->cnt_grammar]();
    
    std::fill(f, f + grammar->cnt_grammar, -INFINITY);

    // 4. load corpus.
    std::cout << "Load sentences..." << std::endl;
    std::vector<std::vector<uint32_t>> sentences = parse_input_file(input_filename, grammar, limit_n_sentences);
    std::cout << "Load sentences finished. Total instances:" << sentences.size() << std::endl;
    if (sentences.empty()) {
        std::cerr << "Error: No sentences loaded." << std::endl;
        return 1;
    }

    // 5. split dataset if needed.
    bool is_split_dataset = config["main"]["split_data"]["enabled"].as<bool>();
    std::vector<std::vector<uint32_t>> train_set;
    std::vector<std::vector<uint32_t>> valid_set;
    if(is_split_dataset){
        double train_fraction = config["main"]["split_data"]["train_fraction"].as<double>();
        std::string train_set_file_save_path = config["main"]["split_data"]["train_dataset_path"].as<std::string>();
        std::string val_set_file_save_path = config["main"]["split_data"]["val_dataset_path"].as<std::string>();
        split_dataset(sentences, train_set, valid_set, train_fraction);
        save_data_set_to_file(train_set_file_save_path, train_set, grammar);
        save_data_set_to_file(val_set_file_save_path, valid_set, grammar);
    }else{
        train_set = std::move(sentences);
    }
    
    int n_sequences = sentences.size();
    int n_sequences_train = train_set.size();
    int n_sequences_val = valid_set.size();
    std::cout << "Train set size = " << n_sequences_train << std::endl;
    std::cout << "Validation set size = " << n_sequences_val << std::endl;

    // 6. training and validation.
    int T = 0;
    int MS = MAX_SEQUENCE_LENGTH;
    cky_printer printer;
    for(int epoch = 0; epoch < n_epochs; epoch++){
        for(int i = 0; i < n_sequences_train; i++){
            const std::vector<uint32_t>& sentence = train_set[i];
            
            progress_bar(i + 1, n_sequences_train);
            
            #if PRINT_GRAMMAR_EACH_UPDATION_BEFORE == 1
                std::cout << "grammar before iteration: " 
                        << i << " :" << std::endl;
                print_grammar(grammar);
            #endif
            int N = grammar->N();
            const uint32_t* sequence = sentence.data();
            int sequence_length = sentence.size();

            #if PRINT_STEPS == 1
            std::cout << "1. Proceeding Inside Algorithm..." << std::endl;
            #endif
            
            inside_algorithm(sequence, 
                (uint32_t*)(grammar->preterminate_rule_lookup_table),
                (uint32_t*)(grammar->grammar_index),
                (uint32_t*)(grammar->grammar_table),
                alpha,
                sequence_length, grammar->n_syms(), grammar->N(), 
                grammar->T(), MS, grammar->cnt_grammar,
                inside_order_1_rule_iteration_path
                    , grammar
            );
                
            if(ALPHA(0, 0, sequence_length - 1) > 1e-9){
                printer.print_inside_outside_table(alpha,  grammar->N(), grammar->T(), sequence_length, MS, grammar);
                std::cout << "assert failed: Log possibility should less than or equal to 0." << std::endl;
                assert(false);
            }            

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
                sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MS, grammar->cnt_grammar,
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
                sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MS, 
                grammar->cnt_grammar
                #ifdef DEBUG_INSIDE_ALGORITHM
                    , grammar
                #endif
                , grammar->symbol_A_vector
            );
            
            #if PRINT_STEPS == 1
                std::cout << "Calculate Expectation Count Finished." << std::endl;
            #endif

            #if PRINT_STEPS == 1
                std::cout << "4. Proceeding Update Parameters..." << std::endl;
            #endif
            bool update_parameter_immediately = 
            (i == n_sequences_train - 1) || 
            ((batch_size_for_parameter_update > 0) && (i % batch_size_for_parameter_update) == 0) ? true : false;

            kernel_update_parameters(f, count, mu, beta, sequence, 
                (uint32_t*)(grammar->preterminate_rule_lookup_table),
                (uint32_t*)(grammar->grammar_index),
                (uint32_t*)(grammar->grammar_table),
                alpha,
                sequence_length, grammar->n_syms(), grammar->N(), grammar->T(),
                MS, grammar->cnt_grammar
                    , grammar
                , update_parameter_immediately
            );

            #if PRINT_STEPS == 1
                std::cout << "Update Parameter Finished." << std::endl;
            #endif

            #if PRINT_GRAMMAR_EACH_UPDATION_AFTER == 1
                std::cout << "grammar after iteration: " << i << " :" << std::endl;
                print_grammar(grammar);
            #endif

            if(log_intervals > 0 && (i + 1) % log_intervals == 0){
                std::ofstream logfile_ostream = std::ofstream(log_path + "/log_" + 
                                        std::to_string(i + 1) + "_epoch_id_" + std::to_string(epoch) + ".pcfg");
                if(!logfile_ostream){
                    std::cerr << "Error: Could not open log file for writing.\n";
                }
                print_grammar(grammar, logfile_ostream);
            }
            if(save_f && log_f_intervals > 0 && (i + 1) % log_f_intervals == 0){
                log_f(log_path + "/log_" + std::to_string(i + 1) + "_epoch_id_" + std::to_string(epoch) + ".f", f, grammar);
            }
        }

        // clear memory f at the end of an epoch.
        std::fill(f, f + grammar->cnt_grammar, -INFINITY);
        std::ofstream logfile_ostream = std::ofstream(log_path + "/log_epoch_id_" + std::to_string(epoch) + ".pcfg");
        if(!logfile_ostream){
            std::cerr << "Error: Could not open log file for writing.\n";
        }
        print_grammar(grammar, logfile_ostream);

        // Validation
        double log_likelihood = -INFINITY;
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
                sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MS, grammar->cnt_grammar,
                inside_order_1_rule_iteration_path
                #ifndef USE_CUDA
                    , grammar
                #endif
            );
            
            double log_likelihood_this_sentence = alpha[0 + 0 + sequence_length - 1];
            if(log_likelihood_this_sentence >= 0 + 1e-9){
                printer.print_inside_outside_table(alpha,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);
                assert(false);
            }

            if (log_likelihood_this_sentence > -INFINITY) {
                LOG_SUM_EXP_SET(log_likelihood, log_likelihood_this_sentence);
            } else {
                #if SANITARY_OUTPUT == 0
                    std::cout << "Warning: ignore -inf log likelihood (alpha = -inf). Underflow may have happened.";
                #endif
            }
            
        }
        
        double average_likelihood = log_likelihood - std::log(n_sequences_val);

        std::cout << "Average log likelihood on validate set at epoch " << epoch << " = " ;
        std::cout << std::fixed << std::setprecision(56) << (double)(average_likelihood);
        std::cout << " = " << log_likelihood << "-" << std::log(n_sequences_val) << "  " << std::endl;
    }
    
    // 7. log results.
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


#endif