#ifndef USE_CUDA
#include "statistics/statistics.hpp"

double derivation_entropy(std::vector<uint32_t> derivations){  // vector of derivations (grammar IDs)
    return _sequence_entropy(derivations);
}

double word_entropy(std::vector<uint32_t> words){
    return _sequence_entropy(words);
}

double _sequence_delay_L_mutual_entropy(std::vector<uint32_t> words, int L){
    std::map<uint32_t, uint32_t> word_counter;
    std::map<uint64_t, uint32_t> word_joint_counter;
    std::map<uint32_t, double> word_possibility;
    std::map<uint64_t, double> word_joint_possibility;

    for(int i = 0; i < words.size() - L; i++){
        word_joint_counter[static_cast<uint64_t>(words[i]) << 32 | words[i + L]] ++;
    }
    for(int i = 0; i < words.size(); i++){
        word_counter[words[i]] ++;
    }
    long Z = 0;
    for(auto& map_item : word_counter){
        word_possibility[map_item.first] = static_cast<double>(map_item.second) / Z;
    }
    Z = 0;
    for(auto& map_item : word_joint_counter){
        word_joint_possibility[map_item.first] = static_cast<double>(map_item.second) / Z;
    }
    double entropy = 0.0;
    for(auto& map_item: word_joint_possibility){
        uint32_t symbol_1 = (map_item.first >> 32) & 0xFFFFFFFF;
        uint32_t symbol_2 = (map_item.first) & 0xFFFFFFFF;
        double p_ij = map_item.second;
        double p_i = word_possibility.count(symbol_1) > 0 ? word_possibility[symbol_1] : 0.0;
        double p_j = word_possibility.count(symbol_2) > 0 ? word_possibility[symbol_2] : 0.0;
        if (p_i > 0 && p_j > 0 && p_ij > 0) {
                entropy += p_ij * std::log(p_ij / (p_i * p_j));
        }
    }
    
    return entropy;
}

double word_delay_L_mutual_entropy(std::vector<uint32_t> words, int L){
    return _sequence_delay_L_mutual_entropy(words, L);
}


double derivation_delay_L_mutual_entropy(std::vector<uint32_t> derivations, int L){
    return _sequence_delay_L_mutual_entropy(derivations, L);
}

std::vector<uint32_t> to_derivations_by_preorder_iteration(parse_tree* node){
    std::deque<parse_tree*> stack;
    std::vector<uint32_t> result;
    if(node == nullptr)
        return result;
    stack.push_back(node);

    std::unordered_set<parse_tree*> visited;

    while(!stack.empty()){
        parse_tree* current = stack.back();
        
        if(visited.find(current) != visited.end()){
            result.push_back(std::get<5>(current->value));
            visited.erase(current);
            stack.pop_back();
        }else{
            visited.insert(current);
            if (current->right != nullptr) {
                stack.push_back(current->right);
            }
            if (current->left != nullptr) {
                stack.push_back(current->left);
            }
        }
        
    }
    return result;
}



int tree_depth(parse_tree* node){
    if(node == nullptr) 
        return 0;
    int left_right_max_depth = std::max(tree_depth(node->left), tree_depth(node->right));
    return left_right_max_depth + 1;
}

double L_layer_symbol_tree_transitional_entropy(pcfg* grammar, parse_tree* tree, int L){
    if (tree == nullptr) {
        throw std::invalid_argument("Tree cannot be null");
    }
    std::vector<std::vector<int>> layers = dfs_get_all_layers_value<int>(
        grammar, 
        tree, 
        [](const std::tuple<uint32_t, uint32_t, uint32_t, int, double, int>& value){
            return std::get<0>(value); 
        }
    );
    return calculate_delay_L_layer_mutual_information(grammar, layers, L);
}

double L_layer_derivation_tree_transitional_entropy(pcfg* grammar, parse_tree* tree, int L){
    if (tree == nullptr) {
        throw std::invalid_argument("Tree cannot be null");
    }
    auto layers = dfs_get_all_layers_value<int>(
        grammar, 
        tree, 
        [](const std::tuple<uint32_t, uint32_t, uint32_t, int, double, int>& value) {
            return std::get<5>(value);
        }
    );
    return calculate_delay_L_layer_mutual_information(grammar, layers, L);
}

// !important
double prefix_L_parse_entropy(pcfg* grammar, double* alpha, int sequence_length, int end, int L){
    if (end < 0 || L < 0) {
        throw std::invalid_argument("end and L must be non-negative");
    }

    std::vector<double> p_s;
    int N = grammar->N();
    int MS = MAX_SEQUENCE_LENGTH;

    for(std::tuple<uint32_t, uint32_t, uint32_t, double, uint32_t> item : 
                PCFGItemIterator(N, (uint32_t*)grammar->grammar_index, (uint32_t*)grammar->grammar_table)){                    
        uint32_t sym_A = std::get<0>(item);
        uint32_t sym_B = std::get<1>(item);
        uint32_t sym_C = std::get<2>(item);
        double possibility = std::get<3>(item);

        if(IS_EPSILON(sym_C)){
            p_s.emplace_back(possibility * ALPHA(sym_B, std::max(end - L, 0), end));
        }else{
            for(int k = std::min(end - L, 0); k < end; k++){
                p_s.emplace_back(possibility * ALPHA(sym_B, k + 1, end));
            }
        }
    }
    
    // Normalize probabilities
    double total_probability = 0.0;
    for (auto&& p : p_s) {
        total_probability += p;
    }

    // Calculate entropy
    double entropy = 0.0;
    if (total_probability > 0) {
        for (auto&& p : p_s) {
            if (p > 0) {
                double normalized_p = p / total_probability; // Normalize p
                entropy += normalized_p * std::log(normalized_p); // log(p)
            }
        }
        entropy = -entropy; // Final entropy calculation
    }
    return entropy;
}
#endif