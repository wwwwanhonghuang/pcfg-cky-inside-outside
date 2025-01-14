#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <cmath>
#include <filesystem>

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
#include <iostream>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstring>
#include <memory>
#include "distribution/shared_memory.hpp"
#include <string>
#include "kernels/update_parameters.cuh"

#define SHARED_MEMORY_NAME "/shared_mem"
int result = 0;

#define PT_INCREASE pt += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS

#define SYMBOL_AND_POSSIBILITY_EXTRACTION(B, C, P) \
                uint32_t symbols = grammar->grammar_table[pt]; \
                double P = *(double*)(grammar->grammar_table + pt + 1); \
                uint32_t B = (symbols >> 16) & 0xFFFF; \
                uint32_t C = symbols & 0xFFFF;

#define ACK(MSG_TYPE) (MSG_TYPE | (1 << 31))
#define CLEAR_MSG(MSG) MSG.status = EMPTY_SLOT;

void write_grammar_to_file(const std::string& grammar_file_path, pcfg* grammar) {
    std::filesystem::path file_path(grammar_file_path);
    std::filesystem::path directory_path = file_path.parent_path();
    
    if (!directory_path.empty() && !std::filesystem::exists(directory_path)) {
        try {
            std::filesystem::create_directories(directory_path);
            std::cout << "Directory created: " << directory_path << std::endl;
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error creating directory: " << e.what() << std::endl;
            return;
        }
    }

    try {
        std::ofstream logfile_ostream(grammar_file_path);
        if (!logfile_ostream.is_open()) {
            std::cerr << "Failed to open file: " << grammar_file_path << std::endl;
            return;
        }
        print_grammar(grammar, logfile_ostream);
        std::cout << "Grammar written to file: " << grammar_file_path << std::endl;
    } catch (const std::ios_base::failure& e) {
        std::cerr << "File operation error: " << e.what() << std::endl;
    }
}


void execution(int epoch, int partition_id){
    std::cout << "partition " << partition_id << 
            " begin execute epoch " << epoch << std::endl; 
    int begin_i = (partition_id - 1) * 50;
    int end_i = partition_id * 50;
    for(int i = begin_i; i < end_i; i++){
        result += i * partition_id;
    }
}

inline double _calculate_new_possibility(double S, double f_gid) {
    if (S == -INFINITY || f_gid == -INFINITY) 
        return -INFINITY;
    return f_gid - S;
}

int current_epoch = 0;

int main(int argc, char* argv[])
{
    std::string configuration_file_path = "config.yaml";

    if (argc < 2) {
        std::cerr << "Please provide the instance index (i).\n";
        return 1;
    }

    int partition_id = std::stoi(argv[1]);
    std::string program_name = std::string("pcfg-train-") +
            std::string(argv[1]);

    auto shared_memory = 
        std::make_shared<SharedMemory>(program_name.c_str(), CREATE_NEW, sizeof(MemoryStorage));
    
    std::cout << "shared memory created. Name = " << program_name << std::endl;
    auto storage = (MemoryStorage*)(shared_memory->get_data());
    
    int pos = 0;

    if(argc > 2){
        configuration_file_path = std::string(argv[2]);
    }
    /* END of Extend HEADER */

    std::cout << "configuration file path:: " << configuration_file_path << std::endl;
    std::cout << "partition id: " << partition_id << std::endl;

    std::cout << "read configuration file " << configuration_file_path << std::endl;
    YAML::Node config = YAML::LoadFile(configuration_file_path);
    YAML::Node cluster_config = YAML::LoadFile("cluster.yaml");

    if (!config.IsDefined() || !cluster_config.IsDefined()) {
        std::cout << "Error: config.yaml could not be loaded!" << std::endl;
        return 1;
    }

    // 1. load configurations.
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
    
    int sentence_from = cluster_config["pcfg-train"][std::string("pcfg-train-") + std::to_string(partition_id)].as<int>();
    int sentence_to = cluster_config["pcfg-train"][std::string("pcfg-train-") + std::to_string(partition_id)].as<int>();
    std::cout << "sentence_from = " << sentence_from << ", " << "sentence_to = " << sentence_to << std::endl;
    assert(false);

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
    std::vector<std::vector<uint32_t>> sentences = parse_input_file(input_filename, grammar, limit_n_sentences, MAX_SEQUENCE_LENGTH);
    std::cout << "Load sentences finished. Total instances:" << sentences.size() << std::endl;
    if (sentences.empty()) {
        std::cerr << "Error: No sentences loaded." << std::endl;
        return 1;
    }

    std::vector<std::vector<uint32_t>> train_set = std::move(sentences);
    
    int n_sequences = sentences.size();
    int n_sequences_train = train_set.size();
   
    std::cout << "Train set size of this partition " << n_sequences_train << std::endl;
    std::cout << "Total numbers of sentences = " << n_sequences << std::endl;
    std::cout << "train sentence from " << sentence_from << " to " << sentence_to << std::endl;
    
    int T = 0;
    int MS = MAX_SEQUENCE_LENGTH;
    cky_printer printer;

    while(true){
        auto& msg = 
            storage->application_messages[pos];
        pos = (pos + 1) % 1;
        if(msg.status == EMPTY_SLOT) continue;
        std::cout << "get msg, type = " << msg.msg_type << std::endl;
        switch (msg.msg_type)
        {
        case PARTITION_PREPARED:{
            std::cout << "response msg PARTITION_PREPARED" << std::endl;
            CLEAR_MSG(msg);
            const std::string& response = 
                std::string("Hi! This is a response for msg PARTITION_PREPARED. ")
            +   std::string("This partition responsible for sentence from ")
            +   std::to_string(sentence_from)
            +   std::string(" to ")
            +   std::to_string(sentence_to);
            
            Message response_msg;
            response_msg.status = MESSAGE_WAITING; 
            response_msg.msg_type = ACK(PARTITION_PREPARED);
            
            memcpy(response_msg.data, &grammar->cnt_grammar, sizeof(int));
            memcpy(response_msg.data + 4, response.c_str(), response.size() + 1);
            memcpy(&storage->network_communicator_messages[0], &response_msg, sizeof(response_msg));
            break;
        }
        case BEGIN_EPOCH:{
            CLEAR_MSG(msg);
            int epoch = -1;
            int partition_id = 0; 
            memcpy(&epoch, &msg.data[0], sizeof(int));
            
            /* Execution Block */
            for(int i = sentence_from; i < sentence_to; i++){
                const std::vector<uint32_t>& sentence = train_set[i];
                
                progress_bar(i + 1, sentence_to);
                
                int N = grammar->N();
                const uint32_t* sequence = sentence.data();
                int sequence_length = sentence.size();
                
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

                outside_algorithm(mu, beta, sequence, 
                    (uint32_t*)(grammar->preterminate_rule_lookup_table),
                    (uint32_t*)(grammar->grammar_index),
                    (uint32_t*)(grammar->grammar_table),
                    alpha,
                    sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MS, grammar->cnt_grammar,
                    inside_order_1_rule_iteration_path
                    ,grammar
                );

                kernel_expect_count(count, mu, beta, sequence, 
                    (uint32_t*)(grammar->preterminate_rule_lookup_table),
                    (uint32_t*)(grammar->grammar_index),
                    (uint32_t*)(grammar->grammar_table),
                    alpha,
                    sequence_length, grammar->n_syms(), grammar->N(), grammar->T(), MS, 
                    grammar->cnt_grammar
                    , grammar
                    , grammar->symbol_A_vector
                );
                
                kernel_update_parameters(f, count, mu, beta, sequence, 
                    (uint32_t*)(grammar->preterminate_rule_lookup_table),
                    (uint32_t*)(grammar->grammar_index),
                    (uint32_t*)(grammar->grammar_table),
                    alpha,
                    sequence_length, grammar->n_syms(), grammar->N(), grammar->T(),
                    MS, grammar->cnt_grammar
                        , grammar, false);
            }

            /* Execution Block */
            std::cout << "response msg BEGIN_EPOCH. result = " 
                << result << std::endl;

            Message response_msg;
            response_msg.status = MESSAGE_WAITING; 
            response_msg.msg_type = ACK(BEGIN_EPOCH);
            
            std::cout << "Client Results: " << std::endl;
            for(int gid = 0; gid < grammar->cnt_grammar; gid++){
                std::cout << "\tf[" << gid << "] = " << f[gid] << std::endl; 
            }

            memcpy(response_msg.data, 
                &grammar->cnt_grammar, sizeof(int));

            memcpy(response_msg.data + sizeof(int), 
                f, grammar->cnt_grammar * sizeof(double));

            memcpy(&storage->network_communicator_messages[0], &response_msg, 
                sizeof(response_msg));
            
            // clear memory f at the end of an epoch.
            std::fill(f, f + grammar->cnt_grammar, -INFINITY);
            break;
        }
        case NOTIFICATE_INTEGRATE_RESULT: {
            CLEAR_MSG(msg);
            double* integrated_f = new double[grammar->cnt_grammar]();
            int epoch = -1;
            memcpy(integrated_f, msg.data, sizeof(double) * grammar->cnt_grammar);
            memcpy(&epoch, msg.data + sizeof(double) * grammar->cnt_grammar, sizeof(int));

            std::cout << "Get integrated grammar results." << std::endl;
            for(int gid = 0; gid < grammar->cnt_grammar; gid++){
                std::cout << "\tintegrated_f[" << gid << "] = " << integrated_f[gid] << std::endl;
            }
            
            /* grammar parameter updation code here */
            std::cout << "STATUS: parameter update." << std::endl;
            int gid = 0;
            int N = grammar->N();
            for(int sym_A = 0; sym_A < N; sym_A++){
                double S = INIT_POSSIBILITY;

                uint32_t grammar_pointer_current = *(grammar->grammar_index + sym_A);
                uint32_t grammar_pointer_next = *(grammar->grammar_index + sym_A + 1);
                int gid_begin = gid;
                for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; PT_INCREASE){
                    double f_gid = integrated_f[gid];
                    
                    LOG_SUM_EXP_SET(S, f_gid);
                    gid ++;
                }

                grammar_pointer_current = *(grammar->grammar_index + sym_A);
                grammar_pointer_next = *(grammar->grammar_index + sym_A + 1);
                gid = gid_begin;

                for(uint32_t pt = grammar_pointer_current; pt < grammar_pointer_next; PT_INCREASE){
                    SYMBOL_AND_POSSIBILITY_EXTRACTION(sym_B, sym_C, possibility);
                    double f_gid = integrated_f[gid];
                    double new_possibility = _calculate_new_possibility(S, f_gid);
                    *(double*)(grammar->grammar_table + pt + 1) = new_possibility;

                    if(IS_EPSILON(sym_C) && IS_TERMINATE(sym_B)){
                        uint64_t key = encode_key(sym_A, sym_B);
                        reverse_grammar_hashtable_set_value(
                            grammar->preterminate_rule_lookup_table, 
                            grammar->cnt_grammar * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS, 
                            key, new_possibility);
                    }
                    gid++;
                }
            }
            /* grammar parameter updation code here */
            std::cout << "STATUS: parameter updated." << std::endl;
            std::string grammar_file_path = "./logs/grammar_log_partition_" +
                std::to_string(partition_id) + "_epoch_" + std::to_string(epoch) + ".pcfg";
            std::cout << "STATUS: save grammar file to " << grammar_file_path << std::endl;

            write_grammar_to_file(grammar_file_path, grammar);
            std::cout << "STATUS: response to communicator."  << std::endl;

            Message response_msg;
            response_msg.status = MESSAGE_WAITING; 
            response_msg.msg_type = ACK(NOTIFICATE_INTEGRATE_RESULT);
            memcpy(&storage->network_communicator_messages[0], &response_msg, sizeof(response_msg));
            break;
        }
        default:
            break;
        }
    }

    // 7. log results.
    std::cout << std::endl << "All finished" << std::endl;
    print_grammar(grammar);
    
    std::ofstream logfile_ostream = std::ofstream("./logs/log_final_" + 
                std::to_string(sentences.size()) + ".pcfg");
    print_grammar(grammar, logfile_ostream);

    delete[] alpha;
    delete[] beta;
    delete[] mu;
    delete[] count;
    delete[] f;

    return 0;
}




