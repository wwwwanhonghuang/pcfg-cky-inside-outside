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
};
void extract_leaf_nodes(parsing::SyntaxTreeNode* node, std::vector<uint32_t>& leaf_nodes);
void test_leaf_order(parsing::SyntaxTreeNode* tree, const std::vector<uint32_t>& input_sequence);

void test_leaf_order(parsing::SyntaxTreeNode* tree, const std::vector<uint32_t>& input_sequence) {
    std::vector<uint32_t> leaf_nodes;
    extract_leaf_nodes(tree, leaf_nodes);

    assert(leaf_nodes == input_sequence && "Leaf nodes do not match input sequence!");
    std::cout << "Test passed: Leaf nodes match input sequence.\n";
};

void extract_leaf_nodes(parsing::SyntaxTreeNode* node, std::vector<uint32_t>& leaf_nodes) {
    if (!node) return;
    if (node->is_leaf()) {
        leaf_nodes.push_back(std::get<0>(node->value)); // Assuming `token_id` holds the terminal symbol
    } else {
        if(node->left) {
            extract_leaf_nodes(node->left, leaf_nodes);
        }

        if(node->right) {
            extract_leaf_nodes(node->right, leaf_nodes);
        }
    }
};

int main(int argc, char* argv[])
{
    YAML::Node config = YAML::LoadFile("config.yaml");
    if (!config.IsDefined()) {
        std::cout << "Error: config.yaml could not be loaded!" << std::endl;
        return 1;
    }

    std::string grammar_filename = config["test_parsing"]["grammar_file"].as<std::string>();
    std::string input_filename = config["test_parsing"]["input"].as<std::string>();

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

        parsing::SyntaxTreeNode* root = parsing::SyntaxTreeNodeParser::parse(grammar, sentence, alpha, inside_order_1_rule_iteration_path);
        if (!root) {
            std::cout << "Failed to parse sentence at index " << i << std::endl;
            continue;
        }

        try {
            test_leaf_order(root, sentence);
        } catch (const std::exception& e) {
            std::cerr << "Validation failed for sentence at index " << i << ": " << e.what() << std::endl;
            continue;
        }
    }

    std::cout << std::endl << "All finished" << std::endl;
    delete[] alpha;

    return 0;
}

