#ifndef H_PRINTER
#define H_PRINTER
#include <map>
#include <iostream>
#include <fstream>
#include "grammar/grammar.hpp"
#include "macros.def"

void log_f(std::string file_path, double* f, pcfg* grammar);
void log_f(double* f, pcfg* grammar, std::ostream& stream);
void progress_bar(int progress, int total, int barWidth = 50);

template<typename T1, typename T2>
void print_map(const std::map<T1, T2>& map){
    std::cout << "{" << std::endl;
    for(const auto& item : map){
        std::cout << item.first << ": " << item.second << ", " << std::endl;
    }
    std::cout << "}" << std::endl;
}

struct cky_printer{
private:
void _print_cell(double* alpha, int i, int j, int N, int T, int sequence_length, int MS, pcfg* grammar, const std::vector<int>& max_length_each_column){
    int base = i * MS + j;
    std::string cell = "";
    std::cout << "(" << i << "," << j << ")"<< "[";
    for(int sym_id = 0; sym_id < N; sym_id++){
        double p = alpha[sym_id * MS * MS + base];
        if(p == -INFINITY) continue;
        std::string item =  grammar->reversed_nonterminate_map[sym_id] + "[" + std::to_string(p) + "]";
        cell += item;
    }
    std::cout << cell.append(max_length_each_column[j - i] - cell.size() - 2, ' ');

    std::cout << "]";
};
void _pre_calculate_column_max_length(double* alpha, int i, int j, int N, int T, int sequence_length, int MS, pcfg* grammar, std::vector<int>& max_length_each_column){
    int base = i * MS + j;
    std::string cell = "";
    for(int sym_id = 0; sym_id < N; sym_id++){
        double p = alpha[sym_id * MS * MS + base];
        
        if(p == -INFINITY) continue;

        std::string item =  grammar->reversed_nonterminate_map[sym_id] + "[" + std::to_string(p) + "]";
        cell += item;
    }
    max_length_each_column[j - i] = cell.size() + 4 > max_length_each_column[j - i] ? cell.size() + 4 : max_length_each_column[j - i];
};
public:
    
    void print_inside_outside_table(double* alpha, int N, int T, int sequence_length, int MS, pcfg* grammar){
        std::vector<int> max_length_each_column(sequence_length); 
        for(int i = 0; i < sequence_length; i++){
            for(int j = i; j < sequence_length; j++){
                _pre_calculate_column_max_length(alpha, i, j, N, T, sequence_length, MS, grammar, max_length_each_column);
            }
        }
        for(int i = 0; i < sequence_length; i++){
            for(int j = i; j < sequence_length; j++){
                _print_cell(alpha, i, j, N, T, sequence_length, MS, grammar, max_length_each_column);
                std::cout << "\t";
            }
            std::cout << std::endl;
        }
    };
};

void print_grammar(pcfg* grammar, std::ostream& stream=std::cout);
#endif
