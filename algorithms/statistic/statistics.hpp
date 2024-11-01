#ifndef H_STATISTICS
#define H_STATISTICS

#include <string>
#include <vector>
#include <cmath>
#include <queue>
#include <functional>
#include <algorithm>
#include <deque>
#include <unordered_set>
#include "tree_parser.hpp"
#include "macros.def"
#include "grammar.hpp"



// TODO, mutual entropy in sliding window.
template<typename T>
double _sequence_entropy(std::vector<T> sequence){
    std::map<uint32_t, uint32_t> counter;
    for(auto&& element : sequence){
        counter[element]++;
    }
    uint32_t Z = 0;

    
    for(auto& map_item : counter){
        Z += map_item.second / Z;
    }

    double entropy = 0.0;
    for(auto& map_item : counter){
        double v = static_cast<double>(map_item.second) / Z;
        entropy += (v == 0 ? 0 : -v * std::log(v)); 
    }
    return entropy;
}


// entropy
double derivation_entropy(std::vector<uint32_t> derivations);

double word_entropy(std::vector<uint32_t> words);

double _sequence_delay_L_mutual_entropy(std::vector<uint32_t> words, int L);

double word_delay_L_mutual_entropy(std::vector<uint32_t> words, int L);


void derivation_delay_L_mutual_entropy(std::vector<uint32_t> derivations, int L);


std::vector<uint32_t> to_derivations_by_preorder_iteration(parse_tree* node);

// Function to calculate mutual information with a delay of L layers
template<typename T>
double calculate_delay_L_layer_mutual_information(pcfg* grammar, std::vector<std::vector<T>> layers, int L){
    std::map<uint64_t, long> joint_counter;
    std::map<uint32_t, long> symbol_counter;

    // Counting occurances
    for(int pre_layer_id = 0; pre_layer_id < layers.size() - L; pre_layer_id++){
        int delay_L_layer_id = pre_layer_id + L;
        
        for(auto&& prelayer_element: layers[pre_layer_id]){
            // Count occurrences in the joint distribution
            symbol_counter[prelayer_element]++;

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
        joint_possibility[map_item.first] = static_cast<double>(map_item.second) / Z;
    }
    Z = 0.0;
    for(auto& map_item : symbol_counter){
        Z += map_item.second;
    }
    for(auto& map_item : symbol_counter){
        symbol_possibility[map_item.first] = static_cast<double>(map_item.second) / Z;
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
    
    if (tree == nullptr) {
        return layers;
    }

    // BFS to select all layers' root's values;
    node_queue.push(tree);
    while(!node_queue.empty()){
        std::vector<T> layer;
        int size = node_queue.size();
        
        for(int i = 0; i < size; i++){
            parse_tree* first_node = node_queue.front();
            node_queue.pop();
            layers.emplace_back(value_selector(first_node->value)); // the A'id in A->BC/
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

int tree_depth(parse_tree* node);

double L_layer_symbol_tree_transitional_entropy(pcfg* grammar, parse_tree* tree, int L);
double L_layer_derivation_tree_transitional_entropy(pcfg* grammar, parse_tree* tree, int L);

// !important
double prefix_L_parse_entropy(pcfg* grammar, float* alpha, int sequence_length, int end, int L);
#endif