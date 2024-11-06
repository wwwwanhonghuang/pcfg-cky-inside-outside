#include <cmath>
#include "utils/printer.hpp"


void print_grammar(pcfg* grammar, std::ostream& stream){
    int N = grammar->N();
    for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
        PCFGItemIterator(N, (uint32_t*) grammar->grammar_index, (uint32_t*) grammar->grammar_table)){
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        double possibility = std::get<3>(item);
        #if COMPUTING_IN_LOG_SPACE
        possibility = std::exp(possibility);
        #endif
        uint32_t gid = std::get<4>(item);
        stream << "[" << gid << "] " << SYMBOL_STR(sym_A) << " -> " << SYMBOL_STR(sym_B) << " " <<
            SYMBOL_STR(sym_C)  << " [" << possibility << "]" << std::endl;
    }
}


void progress_bar(int progress, int total, int barWidth) {
    float percentage = (float) progress / total;
    int pos = (int)(barWidth * percentage);
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) {
            std::cout << "=";
        } else if (i == pos) {
            std::cout << ">";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << int(percentage * 100.0) << " %  " << progress << "/" << total << "\r";
    std::cout.flush();  // Ensures the line is updated in place
}
void log_f(std::string file_path, double* f, pcfg* grammar){
    auto stream = std::ofstream(file_path);
    int N = grammar->N();
    for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
        PCFGItemIterator(N, (uint32_t*) grammar->grammar_index, (uint32_t*) grammar->grammar_table)){
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        double possibility = std::get<3>(item);
        uint32_t gid = std::get<4>(item);
        stream << "[" << gid << "] " << SYMBOL_STR(sym_A) << " -> " << SYMBOL_STR(sym_B) << " " <<
            SYMBOL_STR(sym_C)  << " [" << f[gid] << "]" << std::endl;
    }
}