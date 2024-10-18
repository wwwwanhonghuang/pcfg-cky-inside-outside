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