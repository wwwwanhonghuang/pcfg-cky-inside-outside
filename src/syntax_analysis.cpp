#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <filesystem>

#include "macros.def"
#include "statistics/statistics.hpp"

#include "grammar/grammar.hpp"
#include "grammar/grammar_parser.hpp"
#include "utils/printer.hpp"
#include "utils/application_io.hpp"
#include "dataset/dataset_helper.hpp"
#include "algorithms/tree_parser.hpp"

namespace fs = std::filesystem;

void create_path_if_not_exists(const std::string& path){
    if (!fs::exists(path)) {
        try {
            // Create the directory if it doesn't exist
            fs::create_directories(path);  // create_directories will create intermediate directories if needed
            std::cout << "Directory created: " << path << std::endl;
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error creating directory: " << e.what() << std::endl;
        }
    }
}

int main(int argc, char* argv[])
{
    YAML::Node config = YAML::LoadFile("config.yaml");
    if (!config.IsDefined()) {
        std::cout << "Error: config.yaml could not be loaded!" << std::endl;
        return 1;
    }

    std::string grammar_filename = config["syntax_analysis"]["grammar_file"].as<std::string>();
    std::string input_filename = config["syntax_analysis"]["input"].as<std::string>();
    uint32_t log_itervals = config["syntax_analysis"]["log_intervals"].as<int>();
    std::string log_path = config["syntax_analysis"]["log_path"].as<std::string>();
    bool serialize_to_files = config["syntax_analysis"]["serialize_to_files"].as<bool>(); 
    std::string report_path = config["syntax_analysis"]["report_path"].as<std::string>();
    std::string tree_serialization_path = config["syntax_analysis"]["tree_serialization_path"].as<std::string>();

    create_path_if_not_exists(log_path);
    create_path_if_not_exists(tree_serialization_path);
    create_path_if_not_exists(report_path);

    pcfg* grammar = prepare_grammar(grammar_filename);
    auto inside_order_1_rule_iteration_path = generate_inside_perterminate_iteration_paths(grammar);

    double* alpha = new double[grammar->N() * MAX_SEQUENCE_LENGTH * MAX_SEQUENCE_LENGTH]();

    std::cout << "Load sentences..." << std::endl;
    std::vector<std::vector<uint32_t>> sentences = parse_input_file(input_filename, grammar);
    std::cout << "Load sentences finished. Total instances:" << sentences.size() << std::endl;

    if(sentences.empty()) return 0;

    cky_printer printer;
    size_t n_total_sentences = sentences.size();

    for(int i = 0; i < n_total_sentences; i++){
        auto& sentence = sentences[i];
        progress_bar(i + 1, n_total_sentences);
        if(sentence.size() > 500) {
            std::cout << "Warning: a sentence with length " << 
            sentence.size() << " is skipped." << std::endl;
            continue;
        }

        // std::cout << "sentences length = " << sentence.size() << std::endl;
        parsing::SyntaxTreeNode* root = parsing::SyntaxTreeNodeParser::parse(grammar, sentence, alpha, inside_order_1_rule_iteration_path);
        
        // std::cout << "parse finished" << std::endl;
        if(serialize_to_files){
            parsing::SyntaxTreeNodeSerializer::serialize_tree_to_file(tree_serialization_path + std::string("/sentence_") + 
                std::to_string(i + 1) + std::string(".txt"), root);
        }
        
        // std::cout << "serialize finished" << std::endl;

        std::string statistics_report = statistics::Statistician::report_all_statistics(root, alpha, sentence, grammar, 5);
        // std::cout << "report statistics finished" << std::endl;

        std::string report_filename = report_path + std::string("/sentence_") + std::to_string(i + 1) + std::string(".report");
        std::ofstream report_file_output_stream(report_filename);

        if(!report_file_output_stream){
            std::cerr << "Error: cannot open output file " << report_filename << std::endl;
            continue;
        }else{
            std::ostringstream sentence_serialization_stream;
            for(size_t word_id = 0; word_id < sentence.size(); word_id++){
                sentence_serialization_stream << sentence[word_id] << " ";
            }
            sentence_serialization_stream << std::endl;
            report_file_output_stream << sentence_serialization_stream.str() << 
            statistics_report << std::endl;
        }
    }
    std::cout << std::endl << "All finished" << std::endl;

    delete[] alpha;

    return 0;
}

