#include "grammar/grammar_normalizer.hpp"
#include "utils/data_structure.hpp"
#include "utils/data_encoding.h"
#include "utils/math.hpp"
#include "macros.def"
void normalize_grammar(pcfg* grammar){
    int N = grammar->N();
    
    // for each grammar set that each rule begin with a certain nonteriminate A, {A -> Rule_i}, normalize
    // each item
    for(auto& map_item : grammar->grammar_items_map){
        const std::string& nonterminate = map_item.first; // get the nonterminate symbol string
        uint32_t left_id = grammar->nonterminate_map[nonterminate]; // nonterminate symbol string to symbol id
        std::vector<pcfg_grammar_item>& grammar_items = map_item.second; // get rules that produce by the symbol (id = left_id).
        double logZ = 0; // normalization constant

        // calculate normalization constant
        for(auto& grammar_item : grammar_items){
            logZ = log_sum_exp(logZ, grammar_item.possibility);
        }
        
        for(auto& grammar_item : grammar_items){
            grammar_item.possibility -= logZ; // Grammar possibility in in log form. log (p/Z) = log p - (log Z).
        }
        


        // update grammar table, reverse_grammar_hashtable, and ternimate grammar lookup table.
        /* locate grammars that left is 'left_id' */
        uint32_t current_offset = (uint32_t)(grammar->grammar_index[left_id]);
        uint32_t end_offset = (uint32_t)(grammar->grammar_index[left_id + 1]);

        while(current_offset < end_offset){
            double possibility = *(double*)(grammar->grammar_table + current_offset + 1);
            (*(double*)(grammar->grammar_table + current_offset + 1)) = possibility / Z;
            uint32_t right_side_encode =  (uint32_t)(*(grammar->grammar_table + current_offset));
            uint32_t right_1 = (right_side_encode >> 16) & 0xFFFF;
            uint32_t right_2 = (right_side_encode) & 0xFFFF;
            uint64_t key = encode_key(left_id, right_1);
            
            /* if this is a pre-terminate grammar rule, we need update reverse_grammar_hashtable (key: encode(left_id, right_1))*/
            if(IS_EPSILON(right_2)){
                reverse_grammar_hashtable_set_value(
                    reinterpret_cast<uint32_t*>(grammar->preterminate_rule_lookup_table), 
                    grammar->cnt_grammar * BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS, key, possibility
                );
            }
            // move to next grammar rule item.
            current_offset +=  BYTE_4_CELL_PER_GRAMMAR_TABLE_ITEMS;
        }
    }
}