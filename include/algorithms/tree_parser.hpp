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
struct parse_tree{
public:
    parse_tree* left;
    parse_tree* right;
    std::tuple<uint32_t, uint32_t, uint32_t, int, float, int> value;
    // sym_A sym_B sym_C k possibility gid
    parse_tree(std::tuple <uint32_t, uint32_t, uint32_t, int, float, int> value): 
        value(value), left(nullptr), right(nullptr) {}
    parse_tree(): 
        left(nullptr), right(nullptr) {}
};

std::string serialize_tree(parse_tree* root);
parse_tree* deserialize_tree(std::string tree);
parse_tree* deserialize_tree(std::istringstream& tree);

parse_tree* merge_trees(uint32_t sym_A, int gid, uint32_t sym_B, uint32_t sym_C, int k, float p, parse_tree* left, parse_tree* right);
parse_tree* parse(pcfg* grammar, std::vector<uint32_t> sequence, float* alpha, 
        std::vector<std::tuple<uint32_t, uint32_t>> inside_order_1_rule_iteration_path);
#endif