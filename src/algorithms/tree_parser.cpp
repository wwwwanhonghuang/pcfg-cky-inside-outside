#include "algorithms/tree_parser.hpp"
#include "grammar/grammar.hpp"
#include "algorithms/alg_inside_outside_main.h"
#include "macros.def"
namespace parsing
{
    void SyntaxTreeSerializer::_serialize_helper(SyntaxTree* root, std::ostringstream& oss) {
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

    std::string SyntaxTreeSerializer::serialize_tree(SyntaxTree* root){
        std::ostringstream oss;
        _serialize_helper(root, oss);
        return oss.str();
    };

    void SyntaxTreeSerializer::serialize_tree_to_file(const std::string& filepath, SyntaxTree* root){
        std::string serialized_tree = serialize_tree(root);
        std::ofstream output_file_stream(filepath);
        if(!output_file_stream){
            std::cout << "Failed to open output file stream. Filepath = " << filepath << std::endl;
        }
        output_file_stream << serialized_tree;
    }

    SyntaxTree* SyntaxTreeSerializer::deserialize_tree_from_file(const std::string& filepath){
        std::ifstream infile(filepath);
        if (!infile.is_open()) {
            throw std::runtime_error("Could not open file: " + filepath);
        }

        std::ostringstream oss;
        oss << infile.rdbuf();
        infile.close();

        return deserialize_tree(oss.str());
    }

    SyntaxTree* SyntaxTreeSerializer::deserialize_tree(const std::string& tree){
        std::istringstream iss(tree);
        return deserialize_tree(iss);
    }

    SyntaxTree* SyntaxTreeSerializer::deserialize_tree(std::istringstream& tree_iss){
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

        SyntaxTree* node = new SyntaxTree(std::make_tuple(id1, id2, id3, id4, id5, id6));

        node->left = deserialize_tree(tree_iss);
        node->right = deserialize_tree(tree_iss);
        return node;
    }


    std::shared_ptr<frfl::logger::Logger> SyntaxTreeParser::logger = 
            frfl::logger::Loggers::build_logger<frfl::logger::StdLogger>("SyntaxTreeParser");

    SyntaxTree* SyntaxTreeParser::_parsing_helper(double* alpha, int MS, uint32_t symbol_id, int span_from, int span_to, pcfg* grammar, uint32_t* sequence){
        int N = grammar->N();
        // std::cout << "in parser sym = " << symbol_id << " span_from = " << 
        //     span_from << " span_to = " << span_to << std::endl;

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
        
        // std::cout << "symbol_id = " << symbol_id << " span_from " << span_from << " span_to "
        //     << span_to << std::endl;
        
        // terminate case
        if(IS_TERMINATE(symbol_id)){
            // std::cout << " !!! - terminate: " << symbol_id << std::endl;
            SyntaxTree* node = new SyntaxTree();
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
            // std::cout << " - for " << sym_A << " -> " << sym_B << ", " << sym_C << std::endl;
            
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
                    
                    // std::cout << "     -- in k = " << k
                    //     << " alpha(" << sym_B << ", " << span_from << ", " << k << ")"
                    //     << " = " << ALPHA_GET(sym_B, span_from, k)
                    //     << " alpha(" << sym_C << ", " << k + 1 << ", " << span_to << ")"
                    //     << " = " <<  ALPHA_GET(sym_C, k + 1, span_to)
                    //     << " possibility = " << possibility << std::endl;

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

        // std::cout << "  best_symbol_B = " << best_symbol_B << std::endl;
        // std::cout << "  best_symbol_C = " << best_symbol_C << std::endl;
        assert(best_symbol_B != 0xFFFF || best_symbol_C != 0xFFFF);

        // no spliting.
        if(span_from == span_to){
            // std::cout << "  best_symbol_B = " << best_symbol_B << std::endl;
            SyntaxTree* node = new SyntaxTree();
            node->value = std::make_tuple(sym_A, best_symbol_B, best_symbol_C, span_from, best_v, best_gid); 
            node->right = nullptr;
            SyntaxTree* tree_left = _parsing_helper(alpha, MS, best_symbol_B, span_from, span_to, grammar, sequence);
            node->left = tree_left;        
            return node;
        }
        
        // split span at best k.
        SyntaxTree* tree_1 = _parsing_helper(alpha, MS, best_symbol_B, span_from, best_k, grammar, sequence);
        SyntaxTree* tree_2 = _parsing_helper(alpha, MS, best_symbol_C, best_k + 1, span_to, grammar, sequence);
        SyntaxTree* merged_SyntaxTree = merge_trees(symbol_id, best_gid, best_symbol_B, best_symbol_C, best_k, best_v, tree_1, tree_2);
        return merged_SyntaxTree;
    }

    SyntaxTree* SyntaxTreeParser::merge_trees(uint32_t sym_A, int gid, uint32_t sym_B, uint32_t sym_C, int k, double p, SyntaxTree* left, SyntaxTree* right){
        // std::cout << "merge " << sym_B << " " << sym_C << " -> " << sym_A << std::endl;
        assert((sym_B == 0xFFFF && !left) || std::get<0>(left->value) == sym_B);
        assert((sym_C == 0xFFFF && !right) || std::get<0>(right->value) == sym_C);
        SyntaxTree* result = new SyntaxTree();
        result->left = left;
        result->right = right;
        result->value = std::make_tuple(sym_A, sym_B, sym_C, k, p, gid);
        return result;
    }

    // parse a sequence into syntax tree
    SyntaxTree* SyntaxTreeParser::parse(pcfg* grammar, std::vector<uint32_t> sequence, double* alpha, 
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
        SyntaxTree* node = _parsing_helper(alpha, MAX_SEQUENCE_LENGTH, 0, 0, sequence.size() - 1, grammar, sequence.data());
        return node;
    }
} // namespace parsing

namespace frfl::logger {
    std::unordered_map<std::string, std::shared_ptr<Logger>> Loggers::loggers;

}