#ifndef USE_CUDA
#include <sstream>
#include "statistics/statistics.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>

namespace statistics{
    std::string Statistician::report_all_statistics(parsing::SyntaxTreeNode* node,
    double* alpha, std::vector<uint32_t> sentence, 
                                pcfg* grammar, int max_delays = 5, int max_span_length = 100){
        std::ostringstream oss;
        std::vector<uint32_t> derivations = to_derivations_by_preorder_iteration(node);
        uint32_t* sequence = sentence.data();
        double entropy_of_derivations = derivation_entropy(derivations);
        oss << "derivation_entropy: " << entropy_of_derivations << std::endl;
        double entropy_of_word = word_entropy(sentence);
        oss << "word_entropy: " << entropy_of_word << std::endl;

        for(int d = 1; d <= max_delays; d++){
            double word_delay_d_mutual_entropy = word_delay_L_mutual_entropy(sentence, d);
            oss << "word_delay_" << d << "_mutual_entropy: " << word_delay_d_mutual_entropy << std::endl;
        }

        for(int d = 1; d <= max_delays; d++){
            double derivation_delay_d_mutual_entropy = derivation_delay_L_mutual_entropy(derivations, d);
            oss << "derivation_delay_" << d << "_mutual_entropy: " << derivation_delay_d_mutual_entropy << std::endl;
        }
        
        for(int d = 1; d <= max_delays; d++){
            double word_delay_d_transitional_entropy = word_transitional_entropy_delay_L(sentence, d);
            oss << "word_delay_" << d << "_transitional_entropy: " << word_delay_d_transitional_entropy << std::endl;
        }

        for(int d = 1; d <= max_delays; d++){
            double derivation_delay_d_transitional_entropy = derivation_transitional_entropy_delay_L(derivations, d);
            oss << "derivation_delay_" << d << "_transitional_entropy: " << derivation_delay_d_transitional_entropy << std::endl;
        }

        for(int d = 1; d <= max_delays; d++){
            double word_delay_d_transfer_entropy = sequence_transfer_entropy_delay_L(sentence, d);
            oss << "word_delay_" << d << "_transfer_entropy: " << word_delay_d_transfer_entropy << std::endl;
        }

        for(int d = 1; d <= max_delays; d++){
            double derivation_delay_d_transfer_entropy = sequence_transfer_entropy_delay_L(derivations, d);
            oss << "derivation_delay_" << d << "_transfer_entropy: " << derivation_delay_d_transfer_entropy << std::endl;
        }

        int depth_of_tree = tree_depth(node);
        oss << "depth_of_tree: " << depth_of_tree << std::endl;

        for(int d = 1; d <= max_delays; d++){
            double d_layer_symbol_tree__entropy = L_layer_symbol_tree_mutual_entropy(grammar, node, d);
            oss << d << "_layer_symbol_tree__entropy: " << d_layer_symbol_tree__entropy << std::endl;
        }
        
        for(int d = 1; d <= max_delays; d++){
            double d_layer_derivation_tree__entropy = L_layer_derivation_tree_mutual_entropy(grammar, node, d);
            oss << d << "_layer_derivation_tree__entropy: " << d_layer_derivation_tree__entropy << std::endl;
        }

        for (int span_length = 2; span_length < std::min(max_span_length, (int) sentence.size()); span_length++){
            for (int k = sentence.size() - 1; k >= span_length; k--){
                double prefix_parse_entropy = prefix_L_parse_entropy(grammar, alpha, sentence.size(), k, span_length, sequence);
                oss << "pre_" << span_length << "_end_" << k << ": " << prefix_parse_entropy << std::endl;
            }
        }

        std::vector<std::vector<parsing::SyntaxTreeNode*>> paths = 
            get_paths<parsing::SyntaxTreeNode*>(node, [](parsing::SyntaxTreeNode* node)->parsing::SyntaxTreeNode*{return node;});

        double average_path_length = calculate_average_path_length(node, paths);
        int redundancy = calculate_redundancy(node);
        double symbol_entropy = calculate_tree_symbol_entropy(node);
        int traversal_steps = count_traversal_steps(node);
        double tree_skewness = calculate_skewness(node);
        double skewness = calculate_skewness(derivations);
        double density = calculate_density(node);
        oss << "average_path_length" << ": " << average_path_length << std::endl;
        oss << "redundancy" << ": " << redundancy << std::endl;
        oss << "traversal_steps" << ": " << traversal_steps << std::endl;
        oss << "tree_skewness" << ": " << tree_skewness << std::endl;
        oss << "skewness" << ": " << skewness << std::endl;
        oss << "density" << ": " << density << std::endl;
        oss << "symbol_entropy" << ": " << symbol_entropy << std::endl;

        std::vector<std::vector<parsing::SyntaxTreeNode*>> layers = 
            bfs_get_all_layers_value<parsing::SyntaxTreeNode*>(grammar, node, 
            [](parsing::SyntaxTreeNode* node)->parsing::SyntaxTreeNode*{return node;});
        
        double layer_average_derivation_entropy = calculate_layer_average_derivation_entropy(node, layers);
        double path_average_derivation_entropy = calculate_path_average_derivation_entropy(node);
        oss << "layer_average_derivation_entropy" << ": " << layer_average_derivation_entropy << std::endl;
        oss << "path_average_derivation_entropy" << ": " << path_average_derivation_entropy << std::endl;

        double layer_average_symbol_entropy = calculate_layer_average_symbol_entropy(node, layers);
        double path_average_symbol_entropy = calculate_path_average_symbol_entropy(node);
        oss << "layer_average_symbol_entropy" << ": " << layer_average_symbol_entropy << std::endl;
        oss << "path_average_symbol_entropy" << ": " << path_average_symbol_entropy << std::endl;

        std::vector<std::vector<int>> symbols_of_layers = 
            bfs_get_all_layers_value<int>(grammar, 
            node, [](parsing::SyntaxTreeNode* node)->int{return std::get<0>(node->value);});
        
        std::vector<std::vector<int>> derivations_of_layers = 
            bfs_get_all_layers_value<int>(grammar, 
            node, [](parsing::SyntaxTreeNode* node)->int{return std::get<5>(node->value);});

        double average_layer_symbol_skewness = calculate_average_skewness(symbols_of_layers);
        double average_layer_derivation_skewness = calculate_average_skewness(derivations_of_layers);
        oss << "average_layer_symbol_skewness" << ": " << average_layer_symbol_skewness << std::endl;
        oss << "average_layer_derivation_skewness" << ": " << average_layer_derivation_skewness << std::endl;

        double average_layer_symbol_KL_divergence = calculate_average_kl_divergence(symbols_of_layers);
        double average_layer_derivation_KL_divergence = calculate_average_kl_divergence(derivations_of_layers);
        oss << "average_layer_symbol_KL_divergence" << ": " << average_layer_symbol_KL_divergence << std::endl;
        oss << "average_layer_derivation_KL_divergence" << ": " << average_layer_derivation_KL_divergence << std::endl;
       
        return oss.str();        
    }
}

#endif