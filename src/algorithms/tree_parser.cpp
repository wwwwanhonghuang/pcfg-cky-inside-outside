#include "algorithms/tree_parser.hpp"
#include "grammar/grammar.hpp"
#include "algorithms/alg_inside_outside_main.h"
#include "macros.def"
namespace parsing
{
    void SyntaxTreeSerializer::_serialize_helper(SyntaxTreeNode* root, std::ostringstream& oss) {
        if (root == nullptr) {
            oss << "# ";
            return;
        }

        oss << std::get<0>(root->value) << " "
            << std::get<1>(root->value) << " "
            << std::get<2>(root->value) << " "
            << std::get<3>(root->value) << " "
            << std::get<4>(root->value) << " "
            << std::get<5>(root->value) << " ";

        _serialize_helper(root->left, oss);
        _serialize_helper(root->right, oss);
    };

    std::string SyntaxTreeSerializer::serialize_tree(SyntaxTreeNode* root){
        std::ostringstream oss;
        _serialize_helper(root, oss);
        return oss.str();
    };

    void SyntaxTreeSerializer::serialize_tree_to_file(const std::string& filepath, SyntaxTreeNode* root){
        std::string serialized_tree = serialize_tree(root);
        std::ofstream output_file_stream(filepath);
        if(!output_file_stream){
            std::cout << "Failed to open output file stream. Filepath = " << filepath << std::endl;
        }
        output_file_stream << serialized_tree;
    }

    SyntaxTreeNode* SyntaxTreeSerializer::deserialize_tree_from_file(const std::string& filepath){
        std::ifstream infile(filepath);
        if (!infile.is_open()) {
            throw std::runtime_error("Could not open file: " + filepath);
        }

        std::ostringstream oss;
        oss << infile.rdbuf();
        infile.close();

        return deserialize_tree(oss.str());
    }

    SyntaxTreeNode* SyntaxTreeSerializer::deserialize_tree(const std::string& tree){
        std::istringstream iss(tree);
        return deserialize_tree(iss);
    }

    SyntaxTreeNode* SyntaxTreeSerializer::deserialize_tree(std::istringstream& tree_iss){
        std::string token;
        if (!(tree_iss >> token) || token == "#") { // Check for null pointer
            return nullptr;
        }

        uint32_t id1 = std::stoul(token);
        uint32_t id2, id3;
        int id4;
        double id5;
        int id6;

        tree_iss >> id2 >> id3 >> id4 >> id5 >> id6;

        SyntaxTreeNode* node = new SyntaxTreeNode(std::make_tuple(id1, id2, id3, id4, id5, id6));

        node->left = deserialize_tree(tree_iss);
        node->right = deserialize_tree(tree_iss);
        return node;
    }


    std::shared_ptr<frfl::logger::Logger> SyntaxTreeParser::logger = 
            frfl::logger::Loggers::build_logger<frfl::logger::StdLogger>("SyntaxTreeParser");

    SyntaxTreeNode* SyntaxTreeParser::_parsing_helper(double* alpha, int MS, uint32_t symbol_id, int span_from, int span_to, pcfg* grammar, uint32_t* sequence){
        int N = grammar->N();

        if(span_from > span_to || IS_EPSILON(symbol_id)){
            return nullptr;
        }

        if (alpha == nullptr) {
            std::cerr << "Error: alpha is null!" << std::endl;
            return nullptr;
        }

        if (grammar == nullptr) {
            std::cerr << "Error: grammar is null!" << std::endl;
            return nullptr;
        }
        
        // terminate case
        if(IS_TERMINATE(symbol_id)){
            SyntaxTreeNode* node = new SyntaxTreeNode();
            node->value = std::make_tuple(symbol_id, 0xFFFF, 0xFFFF, span_from, 1.0f, 0xFFFF); 
            node->right = nullptr;
            node->left = nullptr;   
            return node;
        }
        double p = ALPHA_GET(symbol_id, span_from, span_to);
        uint32_t best_symbol_B = 0xFFFF;
        uint32_t best_symbol_C = 0xFFFF;
        uint32_t best_k = 0;
        uint32_t best_gid = 0;
        
        double best_v = -INFINITY;
        
        
        uint32_t sym_A = symbol_id;
        uint32_t current_offset = grammar->grammar_index[symbol_id];
        uint32_t next_offset = grammar->grammar_index[symbol_id + 1];

        // iterate all grammar rules that the left symbol is sym_A.
        while(current_offset < next_offset){
            uint32_t encode_symbol = (uint32_t)(*(grammar->grammar_table + current_offset));
            uint32_t sym_B = (encode_symbol >> 16) & 0xFFFF;
            uint32_t sym_C = encode_symbol & 0xFFFF;
            double possibility = *(double*)(grammar->grammar_table + current_offset + 1);
            uint32_t gid = (int)(current_offset / BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS);
            
            // An impossible case: span length = 1, however sym_C is not the empty (epsilon) in rule A -> BC.
            // Skip when encounter this condition.
            if(!IS_EPSILON(sym_C) && span_from == span_to) {
                current_offset += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
                continue;
            }

            // Find the best split point k and the value of possibility.
            if(!IS_EPSILON(sym_C)){
                for(int k = span_from; k < span_to; k++){
                    double v = possibility + ALPHA_GET(sym_B, span_from, k) + ALPHA_GET(sym_C, k + 1, span_to);

                    if(v > best_v){
                        best_v = v;
                        best_symbol_B = sym_B;
                        best_symbol_C = sym_C;
                        best_k = k;
                        best_gid = gid;
                    }
                }
            }else{
                double v = possibility + ALPHA_GET(sym_B, span_from, span_to);

                if(v > best_v){
                    best_v = v;
                    best_symbol_B = sym_B;
                    best_symbol_C = sym_C;
                    best_k = -1;
                    best_gid = gid;
                }
            }

            current_offset += BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
        }

        assert(best_symbol_B != 0xFFFF || best_symbol_C != 0xFFFF);

        if(span_from == span_to){
            SyntaxTreeNode* node = new SyntaxTreeNode();
            node->value = std::make_tuple(sym_A, best_symbol_B, best_symbol_C, span_from, best_v, best_gid); 
            node->right = nullptr;
            SyntaxTreeNode* tree_left = 
                _parsing_helper(alpha, MS, best_symbol_B, span_from, span_to, grammar, sequence);
            node->left = tree_left;        
            return node;
        }
        
        SyntaxTreeNode* tree_1 = _parsing_helper(alpha, MS, best_symbol_B, span_from, best_k, grammar, sequence);
        SyntaxTreeNode* tree_2 = _parsing_helper(alpha, MS, best_symbol_C, best_k + 1, span_to, grammar, sequence);
        SyntaxTreeNode* merged_SyntaxTreeNode = merge_trees(symbol_id, best_gid, best_symbol_B, best_symbol_C, best_k, best_v, tree_1, tree_2);
        return merged_SyntaxTreeNode;
    }

    SyntaxTreeNode* SyntaxTreeParser::merge_trees(uint32_t sym_A, int gid, uint32_t sym_B, uint32_t sym_C, int k, double p, SyntaxTreeNode* left, SyntaxTreeNode* right){
        assert((sym_B == 0xFFFF && !left) || std::get<0>(left->value) == sym_B);
        assert((sym_C == 0xFFFF && !right) || std::get<0>(right->value) == sym_C);
        SyntaxTreeNode* result = new SyntaxTreeNode();
        result->left = left;
        result->right = right;
        result->value = std::make_tuple(sym_A, sym_B, sym_C, k, p, gid);
        return result;
    }

    SyntaxTreeNode* SyntaxTreeParser::parse(pcfg* grammar, std::vector<uint32_t> sequence, double* alpha, 
            std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path){
        int sequence_length = sequence.size();
        inside_algorithm(sequence.data(), 
            (uint32_t*)(grammar->preterminate_rule_lookup_table),
            (uint32_t*)(grammar->grammar_index),
            (uint32_t*)(grammar->grammar_table),
            alpha,
            sequence_length, grammar->n_syms(), grammar->N(), 
            grammar->T(), MAX_SEQUENCE_LENGTH, grammar->cnt_grammar,
            inside_order_1_rule_iteration_path
            , grammar
        );
        int N = grammar->N();
        int argmax_nonterminate_id = 0;
        double max_inside_value = 0;
        assert(alpha[sequence.size() - 1] > -INFINITY);
        SyntaxTreeNode* node = _parsing_helper(alpha, MAX_SEQUENCE_LENGTH, 0, 0, sequence.size() - 1, grammar, sequence.data());
        return node;
    }
}

namespace frfl::logger {
    std::unordered_map<std::string, std::shared_ptr<Logger>> Loggers::loggers;

}