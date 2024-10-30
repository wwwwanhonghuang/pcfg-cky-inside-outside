#ifndef H_STATISTICS
#define H_STATISTICS

#include <string>
#include <vector>
#include <cmath>
#include <queue>
#include <functional>
#include <algorithm> 
#include "tree_parser.hpp"
#include "macros.def"
#include "grammar.hpp"
// entropy
void derivation_entropy(parse_tree* tree){

}

void word_entropy(parse_tree* tree){

}

void word_delay_L_mutual_entropy(parse_tree* tree, int L){

}
void derivation_delay_L_mutual_entropy(parse_tree* tree, int L){

}

void D_word_transitional_entropy(parse_tree* tree, int D){

}

void D_derivation_transitional_entropy(parse_tree* tree, int D){

}

// Function to calculate mutual information with a delay of L layers
template<typename T>
double calculate_delay_L_layer_mutual_information(std::vector<std::vector<T>> layers, int L){
    std::map<uint64_t, long> joint_counter;
    std::map<uint32_t, long> symbol_counter;

    // Counting occurances
    for(int pre_layer_id = 0; pre_layer_id < layers.size() - L; pre_layer_id ++){
        int delay_L_layer_id = pre_layer_id + L;
        
        for(auto&& prelayer_element: layers[pre_layer_id]){
            // Count occurrences in the joint distribution
            joint_counter[prelayer_element]++;

            for(auto&& postlayer_element: layers[delay_L_layer_id]){
                uint64_t key = ((prelayer_element << 16) & 0xFFFF0000) | (postlayer_element & 0xFFFF);
                joint_counter[key]++;
            }
        }
    }

    for(int layer_id = layers.size() - L; layer_id < layers.size(); layer_id++){
        for(auto&& layer_element: layers[layer_id]){
            symbol_counter[layer_element]++;
        }
    }

    // normalization to possibility
    std::map<uint64_t, double> joint_possibility;
    std::map<uint32_t, double> symbol_possibility;
    double Z = 0.0;
    for(auto& map_item : joint_counter){
        Z += map_item.second;
    }
    for(auto& map_item : joint_counter){
        joint_possibility[map_item.first] = map_item.second / Z;
    }
    Z = 0.0;
    for(auto& map_item : symbol_counter){
        Z += map_item.second;
    }
    for(auto& map_item : symbol_counter){
        symbol_possibility[map_item.first] = map_item.second / Z;
    }

    // mutual entropy
    double mutual_entropy = 0.0;
    for(uint32_t symbol1 = 0; symbol1 < grammar->n_syms(); symbol1++){
        for(uint32_t symbol2 = 0; symbol2 < grammar->n_syms(); symbol2++){
            double p_A = symbol_possibility.find(symbol1) == symbol_possibility.end() ? 0.0 : symbol_possibility.find(symbol1)->second;
            double p_B = symbol_possibility.find(symbol2) == symbol_possibility.end() ? 0.0 : symbol_possibility.find(symbol2)->second;
            double p_AB = 
                (joint_possibility.find(symbol1 << 16 | symbol2) == joint_possibility.end() ? 0.0 : joint_possibility.find(symbol1 << 16 | symbol2)->second) +
                (joint_possibility.find(symbol2 << 16 | symbol1) == joint_possibility.end() ? 0.0 : joint_possibility.find(symbol2 << 16 | symbol1)->second);
            bool is_zero = p_A == 0 || p_B == 0 || p_AB == 0;
            if(!is_zero)
                mutual_entropy +=  p_AB * log(p_AB / (p_A * p_B)); 
        }
    }
    return mutual_entropy;
}

template<typename T> std::vector<std::vector<T>> 
dfs_get_all_layers_value(pcfg* grammar, parse_tree* tree, std::function<T(const std::tuple<uint32_t, uint32_t, uint32_t, int, float, int>&)> value_selector){
    std::vector<std::vector<T>> layers;
    std::queue<parse_tree*> node_queue;
    
    // BFS get all layers' root symbols;
    node_queue.push(tree);
    while(!node_queue.empty()){
        std::vector<T> layer;

        int size = node_queue.size();
        for(int i = 0; i < size; i++){
            parse_tree* first_node = node_queue.front();
            node_queue.pop();
            layers.emplace_back(value_selector(tree->value)); // the A'id in A->BC/
            if(first_node->left != nullptr){
                node_queue.push(first_node->left);
            }else if(first_node->right != nullptr){
                node_queue.push(first_node->right);
            }
        }
        layers.emplace_back(layer);
    }
    return layers;
}

double L_layer_symbol_tree_transitional_entropy(pcfg* grammar, parse_tree* tree, int L){
    std::vector<std::vector<int>> layers = dfs_get_all_layers_value<int>(grammar, tree, [](const std::tuple<uint32_t, uint32_t, uint32_t, int, float, int>& value){
        return std::get<0>(value); 
    });
    return calculate_delay_L_layer_mutual_information(layers, L);
}

double L_layer_derivation_tree_transitional_entropy(pcfg* grammar, parse_tree* tree, int L){
    auto layers = dfs_get_all_layers_value<int>(grammar, tree, 
            [](const std::tuple<uint32_t, uint32_t, uint32_t, int, float, int>& value) {
        return std::get<5>(value);
    });
    return calculate_delay_L_layer_mutual_information(layers, L);
}

// !important
double prefix_L_parse_entropy(pcfg* grammar, float* alpha, int sequence_length, int end, int L){
    std::vector<float> p_s;
    int N = grammar->N();
    int MS = MAX_SEQUENCE_LENGTH;
    for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                PCFGItemIterator(N, grammar->grammar_index, grammar->grammar_table)){                    
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        float possibility = std::get<3>(item);
        if(IS_EPSILON(sym_C)){
            p_s.emplace_back(possibility * ALPHA(sym_B, std::max(end - L, 0), end));
        }else{
            for(int k = std::min(end - L, 0); k < end; k++){
                p_s.emplace_back(possibility * ALPHA(sym_B, k + 1, end));
            }
        }
    }
    
    // Normalize probabilities
    double total_probability = 0.0;
    for (auto&& p : p_s) {
        total_probability += p;
    }

    // Calculate entropy
    double entropy = 0.0;
    if (total_probability > 0) {
        for (auto&& p : p_s) {
            if (p > 0) {
                double normalized_p = p / total_probability; // Normalize p
                entropy += normalized_p * std::log(normalized_p); // log(p)
            }
        }
        entropy = -entropy; // Final entropy calculation
    }
    return entropy;
}

void statistics_main(parse_tree* tree, std::vector<std::string> statistics_items){
    
}
#endif