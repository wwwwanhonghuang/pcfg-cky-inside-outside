#ifndef USE_CUDA
#include <sstream>
#include "statistics/statistics.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>

std::string report_all_statistics(parse_tree* node, double* alpha, std::vector<uint32_t> sentence, 
                            pcfg* grammar, int max_delays = 5){
    std::ostringstream oss;
    std::vector<uint32_t> derivations = to_derivations_by_preorder_iteration(node);
    uint32_t* sequence = sentence.data();
    double entropy_of_derivations = derivation_entropy(derivations);
    oss << "derivation_entropy: " << entropy_of_derivations << std::endl;
    // std::cout << "  - derivation_entropy outputted." << std::endl;
    double entropy_of_word = word_entropy(sentence);
    oss << "word_entropy: " << entropy_of_word << std::endl;
    // std::cout << "  - word_entropy outputted." << std::endl;

    for(int d = 1; d <= max_delays; d++){
        double word_delay_d_mutual_entropy = word_delay_L_mutual_entropy(sentence, d);
        oss << "word_delay_" << d << "_mutual_entropy: " << word_delay_d_mutual_entropy << std::endl;
    }
    // std::cout << "  - word_delay_L_mutual_entropy all outputted." << std::endl;


    for(int d = 1; d <= max_delays; d++){
        double derivation_delay_d_mutual_entropy = derivation_delay_L_mutual_entropy(sentence, d);
        oss << "derivation_delay_" << d << "_mutual_entropy: " << derivation_delay_d_mutual_entropy << std::endl;
    }
    // std::cout << "  - derivation_delay_L_mutual_entropy all outputted." << std::endl;


    int depth_of_tree = tree_depth(node);
    oss << "depth_of_tree: " << depth_of_tree << std::endl;
    // std::cout << "  - depth_of_tree outputted." << std::endl;

    for(int d = 1; d <= max_delays; d++){
        double d_layer_symbol_tree__entropy = L_layer_symbol_tree_mutual_entropy(grammar, node, d);
        oss << d << "_layer_symbol_tree__entropy: " << d_layer_symbol_tree__entropy << std::endl;
    }
    // std::cout << "  - L_layer_symbol_tree__entropy all outputted." << std::endl;

    for(int d = 1; d <= max_delays; d++){
        double d_layer_derivation_tree__entropy = L_layer_derivation_tree_mutual_entropy(grammar, node, d);
            oss << d << "_layer_derivation_tree__entropy: " << d_layer_derivation_tree__entropy << std::endl;
    }
    // std::cout << "  - L_layer_derivation_tree__entropy all outputted." << std::endl;

    for(int span_length = 2; span_length < sentence.size(); span_length++){
        for(int k = sentence.size() - 1; k >= span_length; k--){
            double prefix_parse_entropy = prefix_L_parse_entropy(grammar, alpha, sentence.size(), k, span_length, sequence);
            oss << "pre_" << span_length << "_end_" << k << ": " << prefix_parse_entropy << std::endl;
        }
    }
    // std::cout << "  - entropies of spans all outputted." << std::endl;


    return oss.str();
}

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
    for(int i = 0; i < (int)(words.size()) - L; i++){
        word_joint_counter[static_cast<uint64_t>(words[i]) << 32 | words[i + L]] ++;
    }
    for(int i = 0; i < (int)(words.size()); i++){
        word_counter[words[i]] ++;
    }
    long Z = 0;
    for(auto& map_item : word_counter){
        Z += static_cast<double>(map_item.second);
    }
    for(auto& map_item : word_counter){
        word_possibility[map_item.first] = static_cast<double>(map_item.second) / Z;
    }
    // std::cout << "  - in _sequence_delay_L_mutual_entropy:" << std::endl;
    // std::cout << "      - word_possibility updated:" << std::endl;

    Z = 0;
    for(auto& map_item : word_counter){
        Z += static_cast<double>(map_item.second);
    }
    for(auto& map_item : word_joint_counter){
        word_joint_possibility[map_item.first] = static_cast<double>(map_item.second) / Z;
    }
    // std::cout << "      - word_joint_possibility updated:" << std::endl;


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
    // std::cout << "      - entropy calculate finished:" << std::endl;
    // std::cout << "  - end _sequence_delay_L_mutual_entropy" << std::endl;

    
    return -entropy;
}

double _sequence_transitional_entropy(std::vector<uint32_t> sequence){
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> symbol_transition_counter;
    int length = sequence.size();
    for(int i = 1; i < length; i++){
        uint32_t l = sequence[i - 1];
        uint32_t r = sequence[i];
        symbol_transition_counter[l][r]++;
    }

    double entropy = 0.0;
    for(auto& counter_item: symbol_transition_counter){
        double Z = 0.0;
        for(auto& current_symbol_counter_item: counter_item.second){
            Z += current_symbol_counter_item.second;
        }
        
        if(Z > 0.0){
            for(auto& current_symbol_counter_item: counter_item.second){
                double p = current_symbol_counter_item.second / Z;
                entropy += p * std::log(p);
            }
        }
    }
    return -entropy;
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

double L_layer_symbol_tree_mutual_entropy(pcfg* grammar, parse_tree* tree, int L){
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
    //  for(auto& layer : layers){
    //     for(auto&& e : layer){
    //         std::cout << e << " " ;
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "layer size = " << layers.size() << std::endl;
    // std::ofstream ofs("1.log");
    // for(auto& layer: layers){
    //     for(int i = 0; i < layer.size(); i++){
    //         ofs << layer[i] << " ";
    //     }
    //     ofs << "\n";
    // }
    return calculate_delay_L_layer_mutual_information(grammar, layers, L);
}

double L_layer_derivation_tree_mutual_entropy(pcfg* grammar, parse_tree* tree, int L){
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
double prefix_L_parse_entropy(pcfg* grammar, double* alpha, int sequence_length, int end, int L, uint32_t* sequence){
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
    
    // Normalize probabilities
    double total_probability = 0.0;
    for (auto&& p : p_s) {
        total_probability += std::exp(p);
    }

    // Calculate entropy
    double entropy = 0.0;
    if (total_probability > 0) {
        for (double p : p_s) {
            
            p = std::exp(p);
            
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