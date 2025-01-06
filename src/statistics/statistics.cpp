#ifndef USE_CUDA
#include <sstream>
#include "statistics/statistics.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <deque>


namespace statistics{
    int Statistician::calculate_height(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + std::max(calculate_height(node->left), calculate_height(node->right));
    }

    double Statistician::calculate_layer_average_derivation_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers){
        return _calculate_layer_average_entropy<int>(node, layers, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<5>(_node->value);});
    }

    double Statistician::calculate_layer_average_symbol_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers) {
        return _calculate_layer_average_entropy<int>(node, layers, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<0>(_node->value);});
    }

    double Statistician::calculate_path_average_symbol_entropy(parsing::SyntaxTreeNode* node) {
        return _calculate_path_average_entropy<int>(node, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<0>(_node->value);});
    }

    double Statistician::calculate_path_average_derivation_entropy(parsing::SyntaxTreeNode* node) {
        return _calculate_path_average_entropy<int>(node, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<5>(_node -> value);});
    }


    double Statistician::calculate_tree_symbol_entropy(parsing::SyntaxTreeNode* node){
        std::vector<int> symbols;
        double entropy = 0.0;
        std::unordered_map<int, int> symbol_counter;
        std::queue<parsing::SyntaxTreeNode*> tree_node_queue;

        tree_node_queue.push(node);

        while(!tree_node_queue.empty()){
            parsing::SyntaxTreeNode* current_node = tree_node_queue.front();
            tree_node_queue.pop();
            int symbol_id = std::get<0>(current_node->value);

            if(symbol_id != 0xFFFF){
                symbols.push_back(symbol_id);
                symbol_counter[symbol_id]++;
            }

            if(current_node->left) tree_node_queue.push(current_node->left);
            if(current_node->right) tree_node_queue.push(current_node->right);
        }

        int Z = 0;
        for(const auto& symbol_record : symbol_counter){
            Z += symbol_record.second;
        }

        std::unordered_map<int, double> symbol_frequency;
        for(const auto& symbol_record : symbol_counter){
            symbol_frequency[symbol_record.first] = static_cast<double>(symbol_record.second) / Z;
        }

        for(const auto& frequency_record : symbol_frequency){
            double possibility = frequency_record.second;
            entropy += -std::log(possibility) * possibility;
        }

        return entropy;
    }

    int Statistician::calculate_total_depth(parsing::SyntaxTreeNode* node, std::vector<std::vector<parsing::SyntaxTreeNode*>>& paths) {
        if (paths.size() == 0) return 0;
        int total_length = 0;
        for(int i = 0; i < paths.size(); i++){
            total_length += paths.size();
        }
        return total_length;
    }

    int Statistician::count_nodes(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + count_nodes(node->left) + count_nodes(node->right);
    }

    double Statistician::calculate_average_path_length(parsing::SyntaxTreeNode* node, std::vector<std::vector<parsing::SyntaxTreeNode*>>& paths) {
        int total_depth = calculate_total_depth(node, paths);
        return static_cast<double>(total_depth) / paths.size();
    }

    // TODO: calculate the redunctancy
    bool Statistician::areSubtreesSimilar(parsing::SyntaxTreeNode* node1, parsing::SyntaxTreeNode* node2) {
        if (node1 != nullptr && node2 != nullptr) return true;
        if (node1 == nullptr || node2 == nullptr) return false;
        return (node1->value == node2->value) && 
            areSubtreesSimilar(node1->left, node2->left) &&
            areSubtreesSimilar(node1->right, node2->right);
    }

    int Statistician::calculateRedundancy(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        int redundancy = 0;
        if (areSubtreesSimilar(node->left, node->right)) redundancy = 1;
        return redundancy + calculateRedundancy(node->left) + calculateRedundancy(node->right);
    }

    // Metric : Traversal Time (Steps taken to traverse)
    int Statistician::countTraversalSteps(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + countTraversalSteps(node->left) + countTraversalSteps(node->right);
    }

    // Metric : Skewness (Imbalance of the tree)
    double Statistician::calculateSkewness(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0.0;
        int leftHeight = calculate_height(node->left);
        int rightHeight = calculate_height(node->right);
        return static_cast<double>(leftHeight - rightHeight);  // Positive skew means left-heavy
    }

    // Metric : Tree Density (How full the tree is)
    int Statistician::calculateMaxNodes(int height) {
        return std::pow(2, height) - 1;  // Full binary tree max nodes
    }

    double Statistician::calculateDensity(parsing::SyntaxTreeNode* node) {
        int node_count = count_nodes(node);
        int tree_height = calculate_height(node);
        int maxNodes = calculateMaxNodes(tree_height);
        return static_cast<double>(node_count) / maxNodes;
    }

    std::string Statistician::report_all_statistics(parsing::SyntaxTreeNode* node,
    double* alpha, std::vector<uint32_t> sentence, 
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
            double derivation_delay_d_mutual_entropy = derivation_delay_L_mutual_entropy(derivations, d);
            oss << "derivation_delay_" << d << "_mutual_entropy: " << derivation_delay_d_mutual_entropy << std::endl;
        }
        // std::cout << "  - derivation_delay_L_mutual_entropy all outputted." << std::endl;

        for(int d = 1; d <= max_delays; d++){
            double word_delay_d_transitional_entropy = word_transitional_entropy_delay_L(sentence, d);
            oss << "word_delay_" << d << "_transitional_entropy: " << word_delay_d_transitional_entropy << std::endl;
        }
        for(int d = 1; d <= max_delays; d++){
            double derivation_delay_d_transitional_entropy = derivation_transitional_entropy_delay_L(derivations, d);
            oss << "derivation_delay_" << d << "_transitional_entropy: " << derivation_delay_d_transitional_entropy << std::endl;
        }

        for(int d = 1; d <= max_delays; d++){
            double word_delay_d_transfer_entropy = sequence_transfer_entropy_delay_L(sentence, d);
            oss << "word_delay_" << d << "_transfer_entropy: " << word_delay_d_transfer_entropy << std::endl;
        }

        for(int d = 1; d <= max_delays; d++){
            double derivation_delay_d_transfer_entropy = sequence_transfer_entropy_delay_L(derivations, d);
            oss << "derivation_delay_" << d << "_transfer_entropy: " << derivation_delay_d_transfer_entropy << std::endl;
        }

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
        for (int span_length = 2; span_length < sentence.size(); span_length++){
            for (int k = sentence.size() - 1; k >= span_length; k--){
                double prefix_parse_entropy = prefix_L_parse_entropy(grammar, alpha, sentence.size(), k, span_length, sequence);
                oss << "pre_" << span_length << "_end_" << k << ": " << prefix_parse_entropy << std::endl;
            }
        }

        std::vector<std::vector<parsing::SyntaxTreeNode*>> paths = 
            get_paths<parsing::SyntaxTreeNode*>(node, [](parsing::SyntaxTreeNode* node)->parsing::SyntaxTreeNode*{return node;});

        double average_path_length = calculate_average_path_length(node, paths);
        int redundancy = calculateRedundancy(node);
        double symbol_entropy = calculate_tree_symbol_entropy(node);
        int traversal_steps = countTraversalSteps(node);
        double tree_skewness = calculateSkewness(node);
        double skewness = calculate_skewness(derivations);
        double density = calculateDensity(node);
        oss << "average_path_length" << ": " << average_path_length << std::endl;
        oss << "redundancy" << ": " << redundancy << std::endl;
        oss << "traversal_steps" << ": " << traversal_steps << std::endl;
        oss << "tree_skewness" << ": " << tree_skewness << std::endl;
        oss << "skewness" << ": " << skewness << std::endl;
        oss << "density" << ": " << density << std::endl;
        oss << "symbol_entropy" << ": " << symbol_entropy << std::endl;

        std::vector<std::vector<parsing::SyntaxTreeNode*>> layers = 
            bfs_get_all_layers_value<parsing::SyntaxTreeNode*>(grammar, node, 
            [](parsing::SyntaxTreeNode* node)->parsing::SyntaxTreeNode*{return node;});
        
        double layer_average_derivation_entropy = calculate_layer_average_derivation_entropy(node, layers);
        double path_average_derivation_entropy = calculate_path_average_derivation_entropy(node);
        oss << "layer_average_derivation_entropy" << ": " << layer_average_derivation_entropy << std::endl;
        oss << "path_average_derivation_entropy" << ": " << path_average_derivation_entropy << std::endl;
        std::cout << layers.size() << std::endl;
        std::cout << layer_average_derivation_entropy << std::endl;
        abort();
        double layer_average_symbol_entropy = calculate_layer_average_symbol_entropy(node, layers);
        double path_average_symbol_entropy = calculate_path_average_symbol_entropy(node);
        oss << "layer_average_symbol_entropy" << ": " << layer_average_symbol_entropy << std::endl;
        oss << "path_average_symbol_entropy" << ": " << path_average_symbol_entropy << std::endl;

        std::vector<std::vector<int>> symbols_of_layers = 
            bfs_get_all_layers_value<int>(grammar, 
            node, [](parsing::SyntaxTreeNode* node)->int{return std::get<0>(node->value);});
        
        std::vector<std::vector<int>> derivations_of_layers = 
            bfs_get_all_layers_value<int>(grammar, 
            node, [](parsing::SyntaxTreeNode* node)->int{return std::get<5>(node->value);});
            
       
        double average_layer_symbol_skewness = calculate_average_skewness(symbols_of_layers);
        double average_layer_derivation_skewness = calculate_average_skewness(derivations_of_layers);
        oss << "average_layer_symbol_skewness" << ": " << average_layer_symbol_skewness << std::endl;
        oss << "average_layer_derivation_skewness" << ": " << average_layer_derivation_skewness << std::endl;

        double average_layer_symbol_KL_divergence = calculate_average_kl_divergence(symbols_of_layers);
        double average_layer_derivation_KL_divergence = calculate_average_kl_divergence(derivations_of_layers);
        oss << "average_layer_symbol_KL_divergence" << ": " << average_layer_symbol_KL_divergence << std::endl;
        oss << "average_layer_derivation_KL_divergence" << ": " << average_layer_derivation_KL_divergence << std::endl;
       
        return oss.str();        
    }

    double Statistician::derivation_entropy(std::vector<uint32_t> derivations){  // vector of derivations (grammar IDs)
        return _sequence_entropy(derivations);
    }

    double Statistician::word_entropy(std::vector<uint32_t> words){
        return _sequence_entropy(words);
    }

    double Statistician::derivation_transitional_entropy_delay_L(std::vector<uint32_t> derivations, int L){  // vector of derivations (grammar IDs)
        return _sequence_delay_L_transitional_entropy(derivations, L);
    }

    double Statistician::word_transitional_entropy_delay_L(std::vector<uint32_t> words, int L){
        return _sequence_delay_L_transitional_entropy(words, L);
    }


    


    double Statistician::mutual_information(const std::vector<uint32_t>& seq_X, const std::vector<uint32_t>& seq_Y) {
        if(seq_X.size() != seq_Y.size()){
            std::cout << "Warning: X and Y has different length in the computation of mutual information. " << std::endl;
        }
        double H_X = calculate_entropy_from_element_sequence(seq_X);
        double H_Y = calculate_entropy_from_element_sequence(seq_Y);
        double H_XY = calculate_joint_entropy(seq_X, seq_Y);
        
        return H_X + H_Y - H_XY;
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

    std::vector<uint32_t> Statistician::to_derivations_by_preorder_iteration(parsing::SyntaxTreeNode* node){
        std::deque<parsing::SyntaxTreeNode*> stack;
        std::vector<uint32_t> result;
        if(node == nullptr)
            return result;
        stack.push_back(node);

        std::unordered_set<parsing::SyntaxTreeNode*> visited;

        while(!stack.empty()){
            parsing::SyntaxTreeNode* current = stack.back();
            
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

    int Statistician::tree_depth(parsing::SyntaxTreeNode* node){
        if(node == nullptr) 
            return 0;
        int left_right_max_depth = std::max(tree_depth(node->left), tree_depth(node->right));
        return left_right_max_depth + 1;
    }


    double Statistician::L_layer_symbol_tree_mutual_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L){
        if (tree == nullptr) {
            throw std::invalid_argument("Tree cannot be null");
        }
        std::vector<std::vector<int>> layers = bfs_get_all_layers_value<int>(
            grammar, 
            tree, 
            [](parsing::SyntaxTreeNode* node){
                return std::get<0>(node->value); 
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
                return std::get<5>(node->value);
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
        
        // normalize probabilities
        double total_probability = 0.0;
        for (auto&& p : p_s) {
            total_probability += std::exp(p);
        }

        // calculate entropy
        double entropy = 0.0;
        if (total_probability > 0) {
            for (double p : p_s) {
                
                p = std::exp(p);
                
                if (p > 0) {
                    double normalized_p = p / total_probability;
                    entropy += normalized_p * std::log(normalized_p);
                }
            }
            entropy = -entropy;
        }
        return entropy;
    }
}

#endif