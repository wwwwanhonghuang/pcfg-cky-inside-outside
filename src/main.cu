#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
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
#include "kernels/expect_count.cuh"
#include "kernels/inside.cuh"
#include "kernels/outside.cuh"
#include "kernels/update_parameters.cuh"
#include "kernels/validation.cuh"

uint32_t* to_inside_order_1_rule_iteration_path_array(std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path){
    size_t size = inside_order_1_rule_iteration_path.size();
    uint32_t* inside_order_1_rule_iteration_path_array = new uint32_t[size * 2]();
    int i = 0;
    for(std::tuple<uint32_t, uint32_t>& rule_id: inside_order_1_rule_iteration_path) {
        uint32_t gid = std::get<0>(rule_id);
        uint32_t sym_A = std::get<1>(rule_id);
        inside_order_1_rule_iteration_path_array[i * 2] = gid;
        inside_order_1_rule_iteration_path_array[i * 2 + 1] = sym_A;
        i++;
    }
    return inside_order_1_rule_iteration_path_array;
}
int main(int argc, char* argv[])
{
    YAML::Node config = YAML::LoadFile("config.yaml");
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
    uint32_t* inside_order_1_rule_iteration_path_array = to_inside_order_1_rule_iteration_path_array(inside_order_1_rule_iteration_path);


    // 3. define matrices needed by the inside-outside algorithm.
    double* alpha, *beta, *mu, *count, *f;
    double init_val = -INFINITY;
    uint32_t* symbol_A_vector_device;

    cudaMalloc(&alpha, grammar->N() * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH * sizeof(double));
    cudaMalloc(&beta, (grammar->N() + grammar->T()) * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH * sizeof(double));
    cudaMalloc(&mu, grammar->cnt_grammar * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH * sizeof(double));
    cudaMalloc(&count, grammar->cnt_grammar * sizeof(double));
    cudaMalloc(&f, grammar->cnt_grammar * sizeof(double));
    cudaMalloc(&symbol_A_vector_device, grammar->cnt_grammar * sizeof(uint32_t));

    cudaMemcpy(f, &init_val, grammar->cnt_grammar * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(symbol_A_vector_device, grammar->symbol_A_vector, grammar->cnt_grammar * sizeof(uint32_t), cudaMemcpyHostToDevice);

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
    uint32_t *pretermination_lookuptable_device, *grammar_index_device, *grammar_table_device, *inside_order_1_rule_iteration_path_device;
    cudaMemcpy(&pretermination_lookuptable_device, grammar->preterminate_rule_lookup_table, 
        grammar->cnt_grammar * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(&grammar_index_device, grammar->grammar_index, (grammar->N() + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(&grammar_table_device, grammar->grammar_table, (grammar->cnt_grammar * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(&inside_order_1_rule_iteration_path_device, inside_order_1_rule_iteration_path_array, (inside_order_1_rule_iteration_path.size() * 2) * sizeof(uint32_t), cudaMemcpyHostToDevice);


    // 7. copy train set and validation set to the GPU device.
    std::vector<uint32_t*> train_set_device(train_set.size(), 0);
    std::vector<uint32_t*> valid_set_device(valid_set.size(), 0);
    for(int i = 0; i < train_set.size(); i++){
        cudaMalloc(&train_set_device[i], train_set[i].size() * sizeof(uint32_t));
        const uint32_t* sequence = train_set[i].data();
        cudaMemcpy(train_set_device[i], sequence, train_set[i].size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    }
    for(int i = 0; i < valid_set.size(); i++){
        cudaMalloc(&valid_set_device[i], valid_set[i].size() * sizeof(uint32_t));
        const uint32_t* sequence = train_set[i].data();
        cudaMemcpy(valid_set_device[i], sequence, valid_set[i].size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    }


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
            uint32_t* sequence_device = train_set_device[i];
            
            int sequence_length = sentence.size();

            #if PRINT_STEPS == 1
            std::cout << "1. Proceeding Inside Algorithm..." << std::endl;
            #endif

            // 1. zerolization alpha.
            kernel_inside_alpha_zerolization<<<16, 16>>>(alpha, N, MS);
            cudaDeviceSynchronize();

            // 2. fill alpha (base case).
            kernel_inside_base_fill_alpha<<<16, 16>>>
                (sequence_device, pretermination_lookuptable_device, grammar_index_device, grammar_table_device, alpha, 
                                sequence_length, grammar->n_syms(), N, T, MS, grammar->cnt_grammar, symbol_A_vector_device);
            
            cudaDeviceSynchronize();
            // 3. fill alpha (recursive case).
            kernel_inside_computeSpanKernel<<<16, 16>>>
                (sequence_device, pretermination_lookuptable_device, grammar_index_device, grammar_table_device, alpha, 
                    sequence_length, grammar->n_syms(), N, T, MS, grammar->cnt_grammar,
                    inside_order_1_rule_iteration_path_device, inside_order_1_rule_iteration_path.size(), 0
                #ifdef DEBUG_INSIDE_ALGORITHM
                    , grammar
                #endif
            );

            cudaDeviceSynchronize();
            #if PRINT_STEPS == 1
                std::cout << "Inside Algorithm Finished." << std::endl;
            #endif

            #if PRINT_INSIDE == 1
                printer.print_inside_outside_table(alpha,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);
            #endif


            #if PRINT_STEPS == 1
                std::cout << "2. Proceeding Outside Algorithm..." << std::endl;
            #endif


            kernel_outside_main<<<16, 16>>>(mu, beta, sequence_device, pretermination_lookuptable_device,
                grammar_index_device, grammar_table_device, alpha, sequence_length, 
                grammar->n_syms(), N, T, MS, grammar->cnt_grammar, inside_order_1_rule_iteration_path_device, inside_order_1_rule_iteration_path.size(),
                symbol_A_vector_device
                );

            cudaDeviceSynchronize();
            #if PRINT_STEPS == 1
                std::cout << "Outside Algorithm Finished." << std::endl;
            #endif

            #if PRINT_OUTSIDE == 1
                printer.print_inside_outside_table(beta,  grammar->N(), grammar->T(), sequence_length, MAX_SEQUENCE_LENGTH, grammar);
            #endif

            #if PRINT_STEPS == 1
                std::cout << "3. Proceeding Calculate Expectation Count..." << std::endl;
            #endif

            kernel_expect_count<<<16, 16>>>(count, mu, beta, sequence_device, 
                pretermination_lookuptable_device,
                grammar_index_device,
                grammar_table_device,
                alpha,
                sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MS, 
                grammar->cnt_grammar,
                symbol_A_vector_device
            );
            cudaDeviceSynchronize();

            
            #if PRINT_STEPS == 1
                std::cout << "Calculate Expectation Count Finished." << std::endl;
            #endif

            #if PRINT_STEPS == 1
                std::cout << "4. Proceeding Update Parameters..." << std::endl;
            #endif

            bool update_parameter_immediately = 
            (i == n_sequences_train - 1) || 
            ((batch_size_for_parameter_update > 0) && (i % batch_size_for_parameter_update) == 0) ? true : false;

            kernel_update_parameters<<<16, 16>>>(f, count, mu, beta, sequence_device, 
                pretermination_lookuptable_device,
                grammar_index_device,
                grammar_table_device,
                alpha,
                sequence_length, grammar->n_syms(), grammar->N(), 
                grammar->T(),
                MS, grammar->cnt_grammar, 
                update_parameter_immediately
            );
            cudaDeviceSynchronize();

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
                logfile_ostream.close();
            }
            if(save_f && log_f_intervals > 0 && (i + 1) % log_f_intervals == 0){
                double* host_f = new double[grammar->cnt_grammar]; // Allocate space on the host
                cudaMemcpy(host_f, f, sizeof(double) * grammar->cnt_grammar, cudaMemcpyDeviceToHost); // Copy from device to host
                std::ofstream logfile_ostream = std::ofstream(log_path + "/log_" + std::to_string(i + 1) + "_epoch_id_" + std::to_string(epoch) + ".f");
                log_f(host_f, grammar, logfile_ostream);
                logfile_ostream.close();
            }
            
        }

        // clear memory f at the end of an epoch.
        cudaMemcpy(f, &init_val, grammar->cnt_grammar * sizeof(double), cudaMemcpyHostToDevice);

        
        std::ofstream logfile_ostream = std::ofstream(log_path + "/log_epoch_id_" + std::to_string(epoch) + ".pcfg");
        if(!logfile_ostream){
            std::cerr << "Error: Could not open log file for writing.\n";
        }

        cudaMemcpy(grammar->grammar_table, f,  (grammar->cnt_grammar * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost); // Copy from device to host
        cudaMemcpy(grammar->preterminate_rule_lookup_table, f,  grammar->cnt_grammar * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS * sizeof(uint32_t), cudaMemcpyDeviceToHost); // Copy from device to host

        print_grammar(grammar, logfile_ostream);

        // // Validation
        // // we should do validation at device.
        // double log_likelihood = -INFINITY;
        // double* log_likelihoods_device;
        // cudaMalloc(&log_likelihoods_device, valid_set.size() * sizeof(double));
        
        // validation_at_device<<16, 16>>(
        //     log_likelihoods_device, valid_set.size(), 
        //     pretermination_lookuptable_device, grammar_index_device, grammar_table_device, alpha, 
        //             sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar,
        //             inside_order_1_rule_iteration_path_device, inside_order_1_rule_iteration_path.size());        
        // cudaDeviceSynchronize();

        // double log_likelihoods_host[valid_set.size()];
        // cudaMemcpy(log_likelihoods_device, log_likelihoods_host, valid_set.size() * sizeof(double), cudaMemcpyDeviceToHost); // Copy from device to host
        // for(int i = 0; i < valid_set.size(); i++){
        //     log_likelihood = log_sum_exp(log_likelihood, log_likelihoods_host[i]);
        // }
        // cudaFree(log_likelihoods_device);

        // double average_likelihood = log_likelihood - std::log(n_sequences_val);
        // std::cout << "Average log likelihood on validate set at epoch " << epoch << " = " ;
        // std::cout << std::fixed << std::setprecision(56) << (double)(average_likelihood);
        // std::cout << " = " << log_likelihood << 
        //     "-" << n_sequences_val << "  " << std::endl;
    }
    
    // 7. log results.
    std::cout << std::endl << "All finished" << std::endl;
    print_grammar(grammar);
    
    std::ofstream logfile_ostream = std::ofstream("./logs/log_final_" + std::to_string(sentences.size()) + ".pcfg");
    print_grammar(grammar, logfile_ostream);

    cudaFree(alpha);
    cudaFree(beta);
    cudaFree(mu);
    cudaFree(count);
    cudaFree(f);
    cudaFree(pretermination_lookuptable_device);
    cudaFree(grammar_index_device);
    cudaFree(grammar_table_device);

    for(int i = 0; i < train_set_device.size(); i++){
        cudaFree(train_set_device[i]);
    }
    for(int i = 0; i < valid_set_device.size(); i++){
        cudaFree(valid_set_device[i]);
    }

    free(inside_order_1_rule_iteration_path_array);

    return 0;
}

