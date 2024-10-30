#include "tree_parser.hpp"
#include "grammar.hpp"
#include "alg_inside_outside_main.h"


void serialize_helper(parse_tree* root, std::ostringstream& oss) {
    if (root == nullptr) {
        oss << "# ";
        return;
    }

    oss << std::get<0>(root->sym_id) << " "
        << std::get<1>(root->sym_id) << " "
        << std::get<2>(root->sym_id) << " "
        << std::get<3>(root->sym_id) << " "
        << std::get<4>(root->sym_id) << " ";

    serialize_helper(root->left, oss);
    serialize_helper(root->right, oss);
}


std::string serialize_tree(parse_tree* root){
    std::ostringstream oss;
    serialize_helper(root, oss);
    return oss.str();
}
parse_tree* deserialize_tree(std::string tree){
    std::string token;
    if (!(iss >> token) || token == "#") { // Check for null pointer
        return nullptr;
    }

    uint32_t id1 = std::stoul(token);
    uint32_t id2, id3;
    int id4;
    float id5;

    iss >> id2 >> id3 >> id4 >> id5;

    parse_tree* node = new parse_tree(std::make_tuple(id1, id2, id3, id4, id5));

    node->left = deserialize_helper(iss);
    node->right = deserialize_helper(iss);
    return node;
}
parse_tree* _parsing_helper(uint32_t symbol_id, int span_from, int span_to, pcfg* grammar){
    if(span_from > span_to){
        return nullptr;
    }
    
    float p = ALPHA(symbol_id, span_from, span_to);
    uint32_t best_symbol_B;
    uint32_t best_symbol_C;
    uint32_t best_k;
    float best_v = 0.0f;
    uint32_t sym_A = symbol_id;

    uint32_t current_offset = grammar->grammar_index[symbol_id].int32_value;
    uint32_t next_offset = grammar->grammar_index[symbol_id + 1].int32_value;
    while(current_offset < next_offset){
        uint32_t encode_symbol = (*(grammar->grammar_table + current_offset)).int32_value;
        uint32_t sym_B = (encode_symbol >> 16) & 0xFFFF;
        uint32_t sym_C = encode_symbol & 0xFFFF;
        float possibility = (*(grammar->grammar_table + current_offset + 1)).float32_value;

        if(!IS_EPSILON(sym_C)){
            for(int k = span_from; k < span_to; k++){
                float v = possibility * ALPHA(sym_B, span_from, k) * ALPHA(sym_C, k + 1, span_to);
                if(v > best_v){
                    best_v = v;
                    best_symbol_B = sym_B;
                    best_symbol_C = sym_C;
                    best_k = k;
                }
            }
            break;
        }else{
            float v = possibility * ALPHA(sym_B, span_from, span_to);
            if(v > best_v){
                best_v = v;
                best_symbol_B = sym_B;
                best_symbol_C = sym_C;
                best_k = -1;
            }
        }

        current_offset += 2;
    }
    
    parse_tree* tree_1 = _helper(best_symbol_B, span_from, best_k, grammar);
    parse_tree* tree_2 = _helper(best_symbol_C, best_k + 1, span_to, grammar);
    parse_tree* merged_parse_tree = merge_trees(0, best_symbol_B, best_symbol_C, best_k, best_v, tree_1, tree_2);

}
parse_tree merge_trees(uint32_t sym_A, uint32_t sym_B, uint32_t sym_C, int k, float p, parse_tree* left, parse_tree* right){
    parse_tree* result = new parse_tree();
    result->left = left;
    result->right = right;
    result->value = std::make_tuple(sym_A, sym_B, sym_C, k, p);
    return result;
}
parse_tree* parse(pcfg* grammar, std::vector<uint32_t> sequence, float* alpha, 
        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path){
    int sequence_length = sequence.size();
    inside_algorithm(sequence, 
                (uint32_t*)(grammar->preterminate_rule_lookup_table),
                (uint32_t*)(grammar->grammar_index),
                (uint32_t*)(grammar->grammar_table),
                alpha,
                sequence_length, grammar->n_syms(), grammar->N(), 
                grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar,
                inside_order_1_rule_iteration_path
                #ifdef DEBUG_INSIDE_ALGORITHM
                    , grammar
                #endif
    );
    int N = grammar->N();
    int argmax_nonterminate_id = 0;
    float max_inside_value = 0;
    
    return _parsing_helper(0, 0, sequence.size() - 1, grammar)



}