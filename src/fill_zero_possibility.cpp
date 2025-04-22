#include "main.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <cmath>

#include "macros.def"
#include "utils/tensor.hpp"
#include "algorithms/alg_inside_outside_main.h"
#include "grammar/grammar.hpp"
#include "grammar/grammar_parser.hpp"
#include "utils/printer.hpp"
#include "utils/application_io.hpp"
#include "dataset/dataset_helper.hpp"
#include "constants.h"
#include "utils/math.hpp"
#include "grammar/grammar_normalizer.hpp"

void print_grammar_tmp(pcfg* grammar, std::ostream& stream, double fill_zero = 0){
    int N = grammar->N();
    for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
        PCFGItemIterator(N, (uint32_t*) grammar->grammar_index, (uint32_t*) grammar->grammar_table)){
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        double possibility = std::get<3>(item);
        possibility = std::exp(possibility);
        if (possibility == 0){
            possibility = fill_zero;
        }
    
        uint32_t gid = std::get<4>(item);
        stream << "[" << gid << "] " << SYMBOL_STR(sym_A) << " -> " << SYMBOL_STR(sym_B) << " " <<
            SYMBOL_STR(sym_C)  << " [" << std::fixed << std::setprecision(56) <<
            possibility << "]" << std::endl;
    }
}
int main(int argc, char* argv[])
{
    std::string grammar_filename = "";
    if(argc > 1){
        grammar_filename = std::string(argv[1]);
    }else{
        return 1;
    }
    
    // 2. parse grammar file.
    pcfg* grammar = nullptr;

    try {
        grammar = prepare_grammar(grammar_filename);
        if (grammar == nullptr) {
            throw std::runtime_error("Error: Failed to parse grammar file.");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    
    print_grammar_tmp(grammar, std::cout, 1e-29);
}

