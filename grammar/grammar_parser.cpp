
#include <fstream>
#include <vector>
#include <iostream>

#include "grammar_parser.hpp"
#include "printer.hpp"
#include "data_structure.hpp"

#define _STRICK_CHECK
// #define DEBUG_PRINT_GRAMMAR_PARSING
// #define VERBOSE_GRAMMAR_STRUCTURE_BUILDING
#define DEBUG_PRINT_GRAMMAR_SYMBOL_MAP
std::ifstream open_grammar_file(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the grammar file at path: " << path << std::endl;
    }
    return file;
}

void _record_terminate_symbol(const std::string non_terminate_symbol, pcfg* grammar, 
        const pcfg_grammar_item& grammar_item, std::map<std::string, 
        std::vector<pcfg_grammar_item>>& grammar_items_map){
    if(grammar->terminate_map.find(non_terminate_symbol) == grammar->terminate_map.end()){
            int id = grammar->T();
            grammar->terminate_map.insert(std::make_pair(non_terminate_symbol, id));
            grammar->reversed_terminate_map.insert(std::make_pair(id, non_terminate_symbol));
    }
}

void _record_non_terminate_symbol(const std::string non_terminate_symbol, pcfg* grammar, 
        const pcfg_grammar_item& grammar_item, std::map<std::string, std::vector<pcfg_grammar_item>>& grammar_items_map){
    int id = grammar->N();
    if(grammar->nonterminate_map.find(non_terminate_symbol) == grammar->nonterminate_map.end()){
            grammar->nonterminate_map.insert(std::make_pair(non_terminate_symbol, id));
            grammar->reversed_nonterminate_map.insert(std::make_pair(id, non_terminate_symbol));
    }
}

pcfg* _parse_grammar_file(const std::string& path){
    pcfg* grammar = new pcfg();
    std::ifstream grammar_file = open_grammar_file(path);
    std::map<std::string, std::vector<pcfg_grammar_item>>& grammar_items_map = grammar->grammar_items_map;
    std::string line; 
    int cnt_recognized_grammars = 0;
    
    // Check if the file is successfully opened
    if (!grammar_file.is_open()) {
        std::cerr << "Error: Could not open the grammar file at path: " << path << std::endl;
        return nullptr;
    }

    while (std::getline(grammar_file, line)) {
        if(line == ""){
            continue;
        }
        pcfg_grammar_item rule = parse_grammar_single_line(line);
        std::cout << "[grammar: " << ++cnt_recognized_grammars << "] " << rule.left << "->" << rule.right1 << " " << rule.right2 << " " << rule.possibility << std::endl;
        
        _record_non_terminate_symbol(rule.left, grammar, rule, grammar_items_map);
        if(!map_contains(grammar_items_map, rule.left)){
            grammar_items_map.insert(std::make_pair(rule.left, std::vector<pcfg_grammar_item> {}));
        }

        grammar_items_map.find(rule.left)->second.push_back(rule);
        #ifdef DEBUG_PRINT_GRAMMAR_PARSING
            std::cout << "add rule: " << rule.left << " -> " << rule.right1 << " " << rule.right2 << std::endl;
        #endif
        if(rule.right1 != ""){
            if(rule.right1.size() > 2 && rule.right1[0] == '\'' && *rule.right1.rbegin() == '\''){
                std::string word = std::string(rule.right1.begin() + 1, rule.right1.begin() + rule.right1.size() - 1);
                if(grammar->terminate_map.find(word) == grammar->terminate_map.end()){
                    _record_terminate_symbol(rule.right1, grammar, rule, grammar_items_map);
                }
            }else{
                _record_non_terminate_symbol(rule.right1, grammar, rule, grammar_items_map);
            }
        }

        if(rule.right2 != ""){
            if(rule.right2.size() > 2 && rule.right2[0] == '\'' && *rule.right2.rbegin() == '\''){
                std::string word = std::string(rule.right2.begin() + 1, rule.right2.begin() + rule.right2.size() - 1);
                if(grammar->terminate_map.find(word) == grammar->terminate_map.end()){
                    _record_terminate_symbol(rule.right2, grammar, rule, grammar_items_map);
                }
            }else{
                _record_non_terminate_symbol(rule.right2, grammar, rule, grammar_items_map);
            }
        }

        
    }
    #ifdef DEBUG_PRINT_GRAMMAR_SYMBOL_MAP
        print_map(grammar->nonterminate_map);
        print_map(grammar->reversed_nonterminate_map);
        print_map(grammar->terminate_map);
        print_map(grammar->reversed_terminate_map);
    #endif

    #ifdef _STRICK_CHECK
        int s = 0;
        for(const auto& item : grammar_items_map){
            s += item.second.size();
        }
        if(s != cnt_recognized_grammars){
            std::cout << "grammar item count not match the numbers of valid lines in the grammar file." << std::endl;
        }
    #endif
    grammar->grammar_items_map = grammar_items_map;
    grammar->cnt_grammar = cnt_recognized_grammars;
    return grammar;
}

pcfg* _build_grammar_table(pcfg* grammar, common_32bit* non_terminate_grammars, common_32bit* grammar_vector){
    common_32bit offset;
    #ifdef VERBOSE_GRAMMAR_STRUCTURE_BUILDING
        std::cout << "non_terminate_grammars, addr: " << reinterpret_cast<uintptr_t>(non_terminate_grammars) << std::endl;
        std::cout << "grammar_vector, addr: " << reinterpret_cast<uintptr_t>(grammar_vector) << std::endl;
        std::cout << grammar->N() << std::endl;
    #endif
    for(int i = 0; i < grammar->N(); i++){
        #ifdef VERBOSE_GRAMMAR_STRUCTURE_BUILDING
            std::cout << "process non-terminate symbol id = " << i << "/" << grammar->N() << " symbol = " 
            << (grammar->reversed_nonterminate_map.at(i))  << std::endl;
        #endif

        if(!map_contains(grammar->reversed_nonterminate_map, i)){
            std::cerr << "Error: Non-terminal index " << i << " not found." << std::endl;
            delete[] non_terminate_grammars; // Clean up allocated memory
            delete[] grammar_vector;
            return nullptr;
        }

        if(!map_contains(grammar->grammar_items_map, grammar->reversed_nonterminate_map[i])){
            std::cout << "no rules for this non-terminate. non_terminate_grammars[" << i << "]=" << non_terminate_grammars[i].int32_value << std::endl;
            non_terminate_grammars[i] = offset;
        }else{
            auto& rules = grammar->grammar_items_map.at(grammar->reversed_nonterminate_map[i]);
            int item_count = rules.size();
            non_terminate_grammars[i] = offset;
            #ifdef VERBOSE_GRAMMAR_STRUCTURE_BUILDING
                std::cout << "non-terminate " << i << " item count = " << item_count << " non_terminate_grammars[" << i << "]=" << non_terminate_grammars[i].int32_value << std::endl;
            #endif
            for (int j = 0; j < item_count; j++) {
                pcfg_grammar_item item = rules[j];
                int right1_id = grammar->get_sym_id(item.right1);
                int right2_id = grammar->get_sym_id(item.right2);

                if(right1_id == -1){
                    std::cout << "Error: unrecognized symbol: " << item.right1 << std::endl;
                    return nullptr;
                }
                
                (*(grammar_vector + offset.int32_value)).int32_value =
                        ((right1_id << 16) & (0xFFFF0000) | right2_id & (0x0000FFFF));

                (*(grammar_vector + offset.int32_value + 1)).float32_value = 
                        item.possibility; 

                offset.int32_value += 2;
            }
        }

        if(i >= 1){
            long addr_diff = non_terminate_grammars[i].int32_value - non_terminate_grammars[i - 1].int32_value;
            // std::cout << "addr diff = " << addr_diff << std::endl;
            #ifdef _STRICK_CHECK
                if(addr_diff != 2 * grammar->grammar_items_map.at(grammar->reversed_nonterminate_map[i - 1]).size()){
                    std::cout << "check grammar item memory arrangement failed." << std::endl;
                    delete[] non_terminate_grammars;
                    delete[] grammar_vector;
                    return nullptr;
                }
            #endif
        }
    }
    non_terminate_grammars[grammar->N()].int32_value = offset.int32_value;
    grammar_vector[grammar->cnt_grammar * 2].int32_value = 0xFFFFFFFF; // End Marker

    
    #ifdef _STRICK_CHECK
        int processed_rule_count = 0;
        for(int nonterminate_id = 0; nonterminate_id < grammar->N(); nonterminate_id++) {
            int offset = non_terminate_grammars[nonterminate_id].int32_value;
            int next_offset = non_terminate_grammars[nonterminate_id + 1].int32_value;

            auto it = grammar->grammar_items_map.find(grammar->reversed_nonterminate_map[nonterminate_id]);
            if (it == grammar->grammar_items_map.end()) {
                std::cerr << "Error: No rules found for non-terminal ID " << nonterminate_id << std::endl;
                return nullptr;
            }
            auto& rules = it->second;

            int rule_index_inner = 0;
            while (offset < next_offset) {
                int encoded_symbols = (grammar_vector + offset)->int32_value;
                int symbol_id_right_1 = (encoded_symbols >> 16) & 0xFFFF;
                int symbol_id_right_2 = (encoded_symbols) & 0xFFFF;
                float possibility = (grammar_vector + offset + 1)->float32_value;
                offset += 2;

                // Decode the value
                pcfg_grammar_item rule = rules[rule_index_inner];
                int32_t decoded_value = (grammar->get_sym_id(rule.right1) << 16) & 0xFFFF0000 | 
                                        (grammar->get_sym_id(rule.right2)) & 0x0000FFFF;

                if (decoded_value != encoded_symbols) {
                    std::cout << "Error: check decoding grammars failed. Encoded value = " << encoded_symbols
                            << " decoded value = " << decoded_value << std::endl;
                    return nullptr;
                }
                if(abs(possibility - rule.possibility) > 1e-5){
                    std::cout << "Error: check possibility failed." << std::endl;
                    return nullptr;
                }

                rule_index_inner++;
                processed_rule_count++;
            }
        }

        if (processed_rule_count != grammar->cnt_grammar) {
            std::cout << "Error: check total grammar item failed. Expect " << grammar->cnt_grammar << ", got " << processed_rule_count << std::endl;
        }
    #endif
    grammar->grammar_index = non_terminate_grammars;
    grammar->grammar_table = grammar_vector;
    return grammar;
}


pcfg* _build_preterminate_grammar_lookup_table(pcfg* grammar, common_32bit* non_terminate_grammars, 
                common_32bit* grammar_vector){
    int cnt_grammar = grammar->cnt_grammar;
    int hashtable_length = cnt_grammar * 2;
    int cnt_hashtable_items = cnt_grammar;
    auto& grammar_items_map = grammar->grammar_items_map;
    common_32bit* hashtable = new common_32bit[hashtable_length]();
    for(int nonterminate_id = 0; nonterminate_id < grammar->N(); nonterminate_id++) {
        int offset = non_terminate_grammars[nonterminate_id].int32_value;
        int next_offset = non_terminate_grammars[nonterminate_id + 1].int32_value;

        auto it = grammar_items_map.find(grammar->reversed_nonterminate_map[nonterminate_id]);
        if (it == grammar_items_map.end()) {
                std::cerr << "Error: No rules found for non-terminal ID " << nonterminate_id << std::endl;
                return nullptr;
        }
        auto& rules = it->second;

        int rule_index_inner = 0;
        while (offset < next_offset) {
                int encoded_symbols = (grammar_vector + offset)->int32_value;
                int symbol_id_right_1 = (encoded_symbols >> 16) & 0xFFFF;
                int symbol_id_right_2 = (encoded_symbols) & 0xFFFF;
                pcfg_grammar_item rule = rules[rule_index_inner];

                if(symbol_id_right_1 >= grammar->N() && symbol_id_right_2 == 0xFFFF){
                    #ifdef DEBUG_PRINT_GRAMMAR_PARSING
                    std::cout << "Processing: " << nonterminate_id << " (" << rule.left << ") -> " 
                          << symbol_id_right_1 << " (" << rule.right1 << ") " 
                          << symbol_id_right_2 << " (" << rule.right2 << ")" << std::endl;
                    #endif
                    int key = ((nonterminate_id << 16) & (0xFFFF0000)) | ((symbol_id_right_1) & (0x0000FFFF));
                    int position = (key % cnt_hashtable_items) * 2;
                    int trail_count = 0;
                  
                    while (hashtable[position].int32_value != 0 && hashtable[position].int32_value != key) {
                        position = (position + 2) % hashtable_length;
                        if (++trail_count >= hashtable_length / 2) {
                            std::cout << "Error: Key cannot be found, and there is no available space for inserting a new element in the hash table." << std::endl; 
                            delete[] hashtable;
                            return nullptr;
                        }
                    }
                    hashtable[position].int32_value = key;
                    hashtable[position + 1].float32_value = rule.possibility;
                }
                offset += 2;
                rule_index_inner++;
        }
    }

    #ifdef _STRICK_CHECK

        // Iterate over the grammar rules
        for (int nonterminate_id = 0; nonterminate_id < grammar->N(); nonterminate_id++) {
            int offset = non_terminate_grammars[nonterminate_id].int32_value;
            int next_offset = non_terminate_grammars[nonterminate_id + 1].int32_value;

            auto it = grammar_items_map.find(grammar->reversed_nonterminate_map[nonterminate_id]);
            if (it == grammar_items_map.end()) {
                std::cerr << "Error: No rules found for non-terminal ID " << nonterminate_id << std::endl;
                return nullptr;
            }
            auto& rules = it->second;

            int rule_index_inner = 0;
            while (offset < next_offset) {
                int encoded_symbols = (grammar_vector + offset)->int32_value;
                int symbol_id_right_1 = (encoded_symbols >> 16) & 0xFFFF;
                int symbol_id_right_2 = (encoded_symbols) & 0xFFFF;
                pcfg_grammar_item rule = rules[rule_index_inner];
                // Check if the rule is in the form A -> w_i
                if(symbol_id_right_1 >= grammar->N() && symbol_id_right_2 == 0xFFFF){
                    int key = ((nonterminate_id << 16) & (0xFFFF0000)) | ((symbol_id_right_1) & (0x0000FFFF));
                    int position = (key % cnt_hashtable_items) * 2;
                    int trail_count = 0;

                    while (hashtable[position].int32_value != 0 && hashtable[position].int32_value != key) {
                        position = (position + 2) % hashtable_length;
                        if (++trail_count >= hashtable_length / 2) {
                            std::cerr << "Error: Key not found in hash table after full probing." << std::endl;
                            delete[] hashtable;
                            return nullptr;
                        }
                    }
                    // std::cout << "Found rule: " << rule.left << " -> " << rule.right1 << " in hash table." << "p=" << 
                    // hashtable[position + 1].float32_value << " in " << position << std::endl;
                }

                offset += 2; // Move to the next entry
                rule_index_inner++;
            }
        }

    #endif
    grammar->preterminate_rule_lookup_table = hashtable;
    return grammar;
}

pcfg* prepare_grammar(const std::string& path){
    pcfg* grammar = _parse_grammar_file(path);
    int cnt_grammar = grammar->cnt_grammar;
    auto& grammar_items_map = grammar->grammar_items_map;
    
    // build non_terminate_grammars
    common_32bit* non_terminate_grammars = new common_32bit[grammar->N() + 1]();
    common_32bit* grammar_vector = new common_32bit[cnt_grammar * 2 + 1]();
    
    if(_build_grammar_table(grammar, non_terminate_grammars, grammar_vector) == nullptr){
        delete[] non_terminate_grammars;
        delete[] grammar_vector;
        return nullptr;
    }

    if(_build_preterminate_grammar_lookup_table(grammar, non_terminate_grammars, grammar_vector) == nullptr){
        delete[] non_terminate_grammars;
        delete[] grammar_vector;
        return nullptr;
    }

    return grammar;
    
    
    // return nullptr;
};
