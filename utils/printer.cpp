#include "printer.hpp"


void print_grammar(pcfg* grammar, std::ostream& stream){
    int N = grammar->N();
    for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
        PCFGItemIterator(N, (uint32_t*) grammar->grammar_index, (uint32_t*) grammar->grammar_table)){
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        float possibility = std::get<3>(item);
        uint32_t gid = std::get<4>(item);
        stream << "[" << gid << "] " << SYMBOL_STR(sym_A) << " -> " << SYMBOL_STR(sym_B) << " " <<
            SYMBOL_STR(sym_C)  << " [" << possibility << "]" << std::endl;
    }
}