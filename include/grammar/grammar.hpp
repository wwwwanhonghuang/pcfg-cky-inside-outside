#ifndef H_GRAMMAR
#define H_GRAMMAR
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>
#include "utils/data_structure.hpp"
#include "macros.def"

struct pcfg_grammar_item{
public:
    std::string left;
    std::string right1;
    std::string right2;
    double possibility;
    
    pcfg_grammar_item(std::string left, std::string right1, std::string right2, double possibility): 
        left(left), right1(right1), right2(right2), possibility(possibility){};
};

struct pcfg{
public:
    // low-level grammar storage
    uint32_t* grammar_index;
    uint32_t* grammar_table;
    uint32_t* preterminate_rule_lookup_table;
    uint32_t* symbol_A_vector;

    // high-level grammar storage
    std::map<std::string, std::vector<pcfg_grammar_item>> grammar_items_map {};

    // symbol maps
    std::map<std::string, int> nonterminate_map{};
    std::map<std::string, int> terminate_map{};
    std::map<int, std::string> reversed_nonterminate_map{};
    std::map<int, std::string> reversed_terminate_map{};
    
    int cnt_grammar = 0;

    int N(){
        return this->nonterminate_map.size();
    }

    int T(){
        return this->terminate_map.size();
    }

    int n_syms(){
        return this->N() + this->T();
    }
    
    int get_sym_id(const std::string& symbol);
};


class PCFGItemIterator {
public:
    PCFGItemIterator(int N, uint32_t* grammar_index, uint32_t* grammar_table)
            : _N(N), _grammar_index(grammar_index), _grammar_table(grammar_table){
        this->_grammar_pointer_current = *grammar_index;
        this->_grammar_pointer_next = *(grammar_index + 1);
        this->pt = this->_grammar_pointer_current;
        this->_current_left_symbol_id = 0;
        this->_current_gid = 0;
    }

    std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> operator*() const;

    PCFGItemIterator& operator++();

    bool operator!=(const PCFGItemIterator& other) const;

    PCFGItemIterator begin() const;
    PCFGItemIterator end() const;

private:
    int _N;     
    uint32_t* _grammar_index;    
    uint32_t* _grammar_table; 
    uint32_t _grammar_pointer_current;
    uint32_t _grammar_pointer_next;
    uint32_t _current_left_symbol_id;
    uint32_t _current_gid;
    uint32_t pt;
};

pcfg_grammar_item parse_grammar_single_line(std::string line);
std::vector<std::tuple<uint32_t, uint32_t>> generate_inside_perterminate_iteration_paths(pcfg* grammar);
#endif

