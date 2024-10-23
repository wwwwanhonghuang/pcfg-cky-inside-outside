#include "grammar.hpp"
pcfg_grammar_item parse_grammar_single_line(std::string line){
    std::string left = "";
    std::string right1 = "";
    std::string right2 = "";
    std::string possibility_string = "";

    int state = 1;
    int length = line.size();
    int pos = 0;

    while(pos < length){
        char ch = line[pos];
        switch(state){
            case 1:
                if(ch == '-'){
                    state = 11;
                }else{
                    left += {ch};
                }
                pos += 1;
                break;
            case 11:
                if(ch == '>'){
                    state = 21;
                }else{
                    left += {ch};
                }
                pos += 1;
                break;
            case 21:
                if(ch == ' '){
                    state = 92;
                }else{
                    right1 += {ch};
                }
                pos += 1;
                break;
            case 92:
                if(ch == ' '){
                    pos++;
                }else if (ch == '[')
                {
                    state = 31;
                    pos++;
                }else{
                    state = 22;
                }
                break;
            case 22:
                if(ch == ' '){
                    state = 31;
                }else{
                    right2 += {ch};
                }
                pos += 1;
                break;
            case 31:
                if(ch != '[' && ch != ']' && ch != ' '){
                    possibility_string += {ch};
                }
                pos++;
                break;
        }
    }
    

    return pcfg_grammar_item(left, right1, right2, std::stof(possibility_string));
};

int pcfg::get_sym_id(const std::string& symbol){
    if(map_contains(this->nonterminate_map, symbol)){
        return this->nonterminate_map[symbol];
    }else if(map_contains(this->terminate_map, symbol)){
        return (this->terminate_map[symbol]) + this->N();
    }
    return -1;
}

std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> PCFGItemIterator::operator*() const {
        uint32_t sym_A = this->_current_left_symbol_id;
        uint32_t symbols = *(this->_grammar_table + this->pt);
        float possibility = *(float*)(this->_grammar_table + this->pt + 1);
        uint32_t sym_B = (symbols >> 16) & 0xFFFF;
        uint32_t sym_C = symbols & 0xFFFF;
        return std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t>(sym_A, sym_B, sym_C, possibility, this->_current_gid);
}

PCFGItemIterator& PCFGItemIterator::operator++() {
        if(this->pt + 2 < this->_grammar_pointer_next){
            this->pt += 2;
            this->_current_gid++;
        }else{
            this->_current_left_symbol_id++; 
            if(this->_current_left_symbol_id < this->_N){
                this->_grammar_pointer_current = *(this->_grammar_index + this->_current_left_symbol_id);
                this->_grammar_pointer_next = *(this->_grammar_index + this->_current_left_symbol_id + 1);
                this->pt = this->_grammar_pointer_current;
                this->_current_gid++;
            }
        }
        return *this;
}

bool PCFGItemIterator::operator!=(const PCFGItemIterator& other) const {
        return this->_current_left_symbol_id < other._current_left_symbol_id;
}

PCFGItemIterator PCFGItemIterator::begin() const{
        return PCFGItemIterator(this->_N, this->_grammar_index, this->_grammar_table);
}

PCFGItemIterator PCFGItemIterator::end() const{
        PCFGItemIterator iterator = PCFGItemIterator(this->_N, this->_grammar_index, this->_grammar_table);
        iterator._grammar_pointer_current = *(this->_grammar_index + this->_N);
        iterator._grammar_pointer_next = iterator._grammar_pointer_current;
        iterator.pt = this->_grammar_pointer_current;
        iterator._current_left_symbol_id = this->_N;
        iterator._current_gid =  *(this->_grammar_index + this->_N) / 2;
        return iterator;
}


std::vector<std::tuple<uint32_t, uint32_t>> generate_inside_perterminate_iteration_paths(pcfg* grammar){
    int n_syms = grammar->N() + grammar->T();
    int N = grammar->N();
    bool dependency_graph[n_syms * n_syms]{0};
    int edges = 0;
    std::vector<std::tuple<uint32_t, uint32_t>> rules = std::vector<std::tuple<uint32_t, uint32_t>>();
    std::vector<std::tuple<uint32_t, uint32_t>> results;
    
    for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                                PCFGItemIterator(N, (uint32_t*) grammar->grammar_index, (uint32_t*)grammar->grammar_table)){
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        float possibility = std::get<3>(item);
        uint32_t gid = std::get<4>(item);
        if(!IS_EPSILON(sym_C)) continue;
        dependency_graph[sym_A * n_syms + sym_B] = 1;
        edges++;
    }
    
    // topological sort
    for(int sym_B = 0; sym_B < n_syms; sym_B++){
        bool elimnatable = true;
        for(int sym_A = 0; sym_A < n_syms; sym_A++){
            assert(sym_A != sym_B || !dependency_graph[sym_B * n_syms + sym_A]);
            if(dependency_graph[sym_B * n_syms + sym_A]){
                elimnatable = false;
            }
        }
        if(!elimnatable) continue;
        for(int sym_A = 0; sym_A < n_syms; sym_A++){
            if(dependency_graph[sym_A * n_syms + sym_B]){
                dependency_graph[sym_A * n_syms + sym_B] = 0;
                rules.emplace_back(std::make_pair(sym_A, sym_B));
                edges--;
            }
        }
    }
    
    if(edges > 0){
        std::cout << "Cyclic preterminate production detected." << std::endl;
        assert(edges == 0);
    }
    
    for(auto&& syms : rules){
        uint32_t sym_A = std::get<0>(syms);
        uint32_t sym_B = std::get<1>(syms);

        for(std::tuple<uint32_t, uint32_t, uint32_t, float, uint32_t> item : 
                                PCFGItemIterator(N, (uint32_t*) grammar->grammar_index, (uint32_t*)grammar->grammar_table)){
        
            uint32_t _sym_A = std::get<0>(item);
            uint32_t _sym_B = std::get<1>(item);
            uint32_t sym_C = std::get<2>(item);
            float possibility = std::get<3>(item);
            uint32_t gid = std::get<4>(item);
            if(!IS_EPSILON(sym_C)) continue;
            if(sym_A == _sym_A && sym_B == _sym_B){
                results.emplace_back(std::make_pair(gid, sym_A));
            }
        }
    }

    // for(auto&& gid : results){
    //     uint32_t* addr = ((uint32_t*)grammar->grammar_table + gid * 2);
    //     uint32_t syms = *addr;
    //     uint32_t sym_A = (syms >> 16) & 0xFFFF;
    //     uint32_t sym_B = syms & 0xFFFF;

    //     std::cout << SYMBOL_STR(sym_A) << "->" << SYMBOL_STR(sym_B) << sym_B << std::endl;
    // }
    
    return results;
}
