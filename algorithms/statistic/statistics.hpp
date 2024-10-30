#ifndef H_STATISTICS
#define H_STATISTICS

#include <string>
#include <vector>
#include <cmath>
#include <algorithm> 
#include "tree_parser.hpp"
#include "macros.def"
#include "grammar.hpp"
// entropy
void derivation_entropy(parse_tree* tree){
}

void word_entropy(parse_tree* tree){
}

void D_word_transitional_entropy(parse_tree* tree, int D){
}
void D_derivation_transitional_entropy(parse_tree* tree, int D){
}
void L_layer_tree_transitional_entropy(parse_tree* tree, int L){
}

double prefix_L_parse_entropy(pcfg* grammar, float* alpha, int sequence_length, int end, int L){
    std::vector<float> p_s;
    int N = grammar->N();
    for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                PCFGItemIterator(N, grammar->grammar_index, grammar->grammar_table)){                    
                    
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        float possibility = std::get<3>(item);
        if(IS_EPSILON(sym_C)){
            p_s.emplace_back(possibility * ALPHA(sym_B, max(end - L, 0), end));
        }else{
            for(int k = min(end - L, 0); k < end; k++){
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