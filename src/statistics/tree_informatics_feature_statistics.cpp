#include "statistics/statistics.hpp"
namespace statistics{
    double Statistician::derivation_entropy(std::vector<uint32_t> derivations){
        return _sequence_entropy(derivations);
    }

    double Statistician::word_entropy(std::vector<uint32_t> words){
        return _sequence_entropy(words);
    }

    double Statistician::derivation_transitional_entropy_delay_L(std::vector<uint32_t> derivations, int L){
        return _sequence_delay_L_transitional_entropy(derivations, L);
    }

    double Statistician::word_transitional_entropy_delay_L(std::vector<uint32_t> words, int L){
        return _sequence_delay_L_transitional_entropy(words, L);
    }

    double Statistician::_sequence_delay_L_mutual_entropy(std::vector<uint32_t> words, int L){
        if (words.size() <= L) return 0;
        
        std::vector<uint32_t> pre_layer(words.begin(), words.begin() + words.size() - L);
        std::vector<uint32_t> post_layer(words.begin() + L, words.end());
        double entropy = mutual_information(pre_layer, post_layer);
        return entropy;
    }

    double Statistician::_sequence_delay_L_transitional_entropy(std::vector<uint32_t> sequence, int delay_L){
        std::map<std::pair<uint32_t, uint32_t>, uint32_t> transition_count;
        int length = sequence.size();
        for(int i = delay_L; i < length; i++){
            uint32_t l = sequence[i - delay_L];
            uint32_t r = sequence[i];
            transition_count[{l, r}]++;
        }

        double entropy = 0.0;
        double total_transitions = sequence.size() - delay_L;

        for(const auto& transition : transition_count){
            double p = transition.second / total_transitions;
            entropy += p * std::log(p);
        }
        return -entropy;
    }

    double Statistician::word_delay_L_mutual_entropy(std::vector<uint32_t> words, int L){
        return _sequence_delay_L_mutual_entropy(words, L);
    }

    double Statistician::derivation_delay_L_mutual_entropy(std::vector<uint32_t> derivations, int L){
        return _sequence_delay_L_mutual_entropy(derivations, L);
    }

    double Statistician::L_layer_symbol_tree_mutual_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L){
        if (tree == nullptr) {
            throw std::invalid_argument("Tree cannot be null");
        }
        std::vector<std::vector<int>> layers = bfs_get_all_layers_value<int>(
            grammar, 
            tree, 
            [](parsing::SyntaxTreeNode* node){
                return std::get<TREE_VALUE_INDEX_SYMBOL>(node->value); 
            }
        );
        
        return calculate_delay_L_layer_mutual_information(grammar, layers, L);
    }

    double Statistician::L_layer_derivation_tree_mutual_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L){
        if (tree == nullptr) {
            throw std::invalid_argument("Tree cannot be null");
        }
        
        auto layers = bfs_get_all_layers_value<int>(
            grammar, 
            tree, 
            [](parsing::SyntaxTreeNode* node) {
                return std::get<TREE_VALUE_INDEX_DERIVATION>(node->value);
            }
        );
        
        return calculate_delay_L_layer_mutual_information(grammar, layers, L);
    }

    double Statistician::prefix_L_parse_entropy(pcfg* grammar, double* alpha, int sequence_length, int end, int L, uint32_t* sequence){
        if (end < 0 || L < 0) {
            throw std::invalid_argument("end and L must be non-negative");
        }

        std::vector<double> p_s;
        int N = grammar->N();
        int MS = MAX_SEQUENCE_LENGTH;

        for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                    PCFGItemIterator(N, (uint32_t*) grammar->grammar_index, (uint32_t*) grammar->grammar_table)){                    
            uint32_t sym_A = std::get<0>(item);
            uint32_t sym_B = std::get<1>(item);
            uint32_t sym_C = std::get<2>(item);
            double possibility = std::get<3>(item);

            if(IS_EPSILON(sym_C)){
                p_s.emplace_back(possibility + ALPHA_GET(sym_B, std::max(end - L, 0), end));
            }else{
                for(int k = std::min(end - L, 0); k < end; k++){
                    p_s.emplace_back(possibility + ALPHA_GET(sym_B, k + 1, end));
                }
            }
        }
        
        double total_probability = 0.0;
        for (auto&& p : p_s) {
            total_probability += std::exp(p);
        }
        double entropy = 0.0;
        if (total_probability > 0) {
            for (double p : p_s) {
                
                p = std::exp(p);
                
                if (p > 0) {
                    double normalized_p = p / total_probability; // Normalize p
                    entropy += normalized_p * std::log(normalized_p); // log(p)
                }
            }
            entropy = -entropy;
        }
        return entropy;
    }
}