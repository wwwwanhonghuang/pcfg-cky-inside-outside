#ifndef H_TREE_PARSER
#define H_TREE_PARSER
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include <cstdint>
#include "grammar/grammar.hpp"
#include "utils/logger.hpp"

namespace parsing{    
    struct SyntaxTree{
    public:
        SyntaxTree* left;
        SyntaxTree* right;
        std::tuple<uint32_t, uint32_t, uint32_t, int, double, int> value;

        // sym_A sym_B sym_C k possibility gid
        SyntaxTree(std::tuple <uint32_t, uint32_t, uint32_t, int, double, int> value): 
            value(value), left(nullptr), right(nullptr) {}

        SyntaxTree():
            left(nullptr), right(nullptr) {}
    };

    class SyntaxTreeParser{
    public:
        static SyntaxTree* merge_trees(uint32_t sym_A, int gid, uint32_t sym_B, uint32_t sym_C, int k, 
            double p, SyntaxTree* left, SyntaxTree* right);
        static SyntaxTree* parse(pcfg* grammar, std::vector<uint32_t> sequence, double* alpha, 
            std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path);
        static void serialize_tree_to_file(std::string filepath, SyntaxTree* root);
        
        SyntaxTreeParser(){
        };
    private:
        static SyntaxTree* _parsing_helper(double* alpha, int MS, uint32_t symbol_id, int span_from, int span_to, pcfg* grammar, uint32_t* sequence);
        static std::shared_ptr<frfl::logger::Logger> logger;
    };
    
    class SyntaxTreeSerializer{
    public:
        static std::string serialize_tree(SyntaxTree* root);
        static void serialize_tree_to_file(const std::string& filepath, SyntaxTree* root);
        static SyntaxTree* deserialize_tree(const std::string& tree);
        static SyntaxTree* deserialize_tree(std::istringstream& tree);
        static SyntaxTree* deserialize_tree_from_file(const std::string& filepath);

    private:
        static void _serialize_helper(SyntaxTree* root, std::ostringstream& oss);

    };
}


#endif