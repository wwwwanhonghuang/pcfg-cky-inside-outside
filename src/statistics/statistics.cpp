#ifndef USE_CUDA
#include <sstream>
#include "statistics/statistics.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <deque>


namespace statistics{
    int Statistician::calculateHeight(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + std::max(calculateHeight(node->left), calculateHeight(node->right));
    }

    double Statistician::calculate_layer_average_derivation_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers){
        return _calculate_layer_average_entropy<int>(node, layers, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<0>(_node->value);});
    }

    double Statistician::calculate_layer_average_symbol_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers) {
        return _calculate_layer_average_entropy<int>(node, layers, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<5>(_node->value);});
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
            entropy += -std::log(frequency_record.second) * frequency_record.second;
        }

        return entropy;
    }

    int Statistician::calculateTotalDepth(parsing::SyntaxTreeNode* node, int depth = 0) {
        if (node == nullptr) return 0;
        return depth + calculateTotalDepth(node->left, depth + 1) + 
            calculateTotalDepth(node->right, depth + 1);
    }

    int Statistician::countNodes(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + countNodes(node->left) + countNodes(node->right);
    }

    double Statistician::averagePathLength(parsing::SyntaxTreeNode* node) {
        int totalDepth = calculateTotalDepth(node);
        int totalNodes = countNodes(node);
        return static_cast<double>(totalDepth) / totalNodes;
    }

    // Metric : Redundancy (Measuring Subtree Similarity)
    bool Statistician::areSubtreesSimilar(parsing::SyntaxTreeNode* node1, parsing::SyntaxTreeNode* node2) {
        if (node1 == nullptr && node2 == nullptr) return true;
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
        int leftHeight = calculateHeight(node->left);
        int rightHeight = calculateHeight(node->right);
        return static_cast<double>(leftHeight - rightHeight);  // Positive skew means left-heavy
    }

    // Metric : Tree Density (How full the tree is)
    int Statistician::calculateMaxNodes(int height) {
        return std::pow(2, height) - 1;  // Full binary tree max nodes
    }

    double Statistician::calculateDensity(parsing::SyntaxTreeNode* node) {
        int nodeCount = countNodes(node);
        int treeHeight = calculateHeight(node);
        int maxNodes = calculateMaxNodes(treeHeight);
        return static_cast<double>(nodeCount) / maxNodes;
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

        // std::cout << "  - entropies of spans all outputted." << std::endl;
        double average_path_length = averagePathLength(node);
        int redundancy = calculateRedundancy(node);
        double symbol_entropy = calculate_tree_symbol_entropy(node);
        int traversal_steps = countTraversalSteps(node);
        double skewness = calculateSkewness(node);
        double density = calculateDensity(node);
        oss << "average_path_length" << ": " << average_path_length << std::endl;
        oss << "redundancy" << ": " << redundancy << std::endl;
        oss << "traversal_steps" << ": " << traversal_steps << std::endl;
        oss << "skewness" << ": " << skewness << std::endl;
        oss << "density" << ": " << density << std::endl;
        oss << "symbol_entropy" << ": " << symbol_entropy << std::endl;
        std::vector<std::vector<parsing::SyntaxTreeNode*>> layers = 
            bfs_get_all_layers_value<parsing::SyntaxTreeNode*>(grammar, node, [](parsing::SyntaxTreeNode* node)->parsing::SyntaxTreeNode*{return node;});
        double layer_average_derivation_entropy = calculate_layer_average_derivation_entropy(node, layers);
        double path_average_derivation_entropy = calculate_path_average_derivation_entropy(node);
        oss << "layer_average_derivation_entropy" << ": " << layer_average_derivation_entropy << std::endl;
        oss << "path_average_derivation_entropy" << ": " << path_average_derivation_entropy << std::endl;

        double layer_average_symbol_entropy = calculate_layer_average_symbol_entropy(node, layers);
        double path_average_symbol_entropy = calculate_path_average_symbol_entropy(node);
        oss << "layer_average_symbol_entropy" << ": " << layer_average_symbol_entropy << std::endl;
        oss << "path_average_symbol_entropy" << ": " << path_average_symbol_entropy << std::endl;


        for(int d = 1; d <= max_delays; d++){
            double value_layer_symbol_transfer_entropy_delay_L = layer_symbol_transfer_entropy_delay_L(layers, d);
            oss << "layer_symbol_transfer_entropy_delay_" << d << ": " << value_layer_symbol_transfer_entropy_delay_L << std::endl;
        }

        for(int d = 1; d <= max_delays; d++){
            double value_layer_derivation_transfer_entropy_delay_L = layer_derivation_transfer_entropy_delay_L(layers, d);
            oss << "layer_derivation_transfer_entropy_delay_" << d << ": " << value_layer_derivation_transfer_entropy_delay_L << std::endl;
        }


        // metrics for layer symbol and derivation sequence distribution's skew..
            // Mean Divergence
            // Skewness
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
        oss << "average_layer_symbol_skewness" << ": " << average_layer_symbol_skewness << std::endl;
        oss << "average_layer_derivation_skewness" << ": " << average_layer_derivation_skewness << std::endl;


        

        return oss.str();

    }

    double Statistician::derivation_entropy(std::vector<uint32_t> derivations){  // vector of derivations (grammar IDs)
        return _sequence_entropy(derivations);
    }

    double Statistician::word_entropy(std::vector<uint32_t> words){
        return _sequence_entropy(words);
    }

    double Statistician::_sequence_delay_L_mutual_entropy(std::vector<uint32_t> words, int L){
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
        for(auto& map_item : word_joint_counter){
            Z += static_cast<double>(map_item.second);
        }

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

        return -entropy;
    }

    double Statistician::_sequence_transitional_entropy(std::vector<uint32_t> sequence){
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


    // static double L_layer_symbol_tree_transitional_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L) {

    // }

    // static double L_layer_derivation_tree_transitional_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L) {

    // }


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

    // ! important
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
        
        // Normalize probabilities
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
                    double normalized_p = p / total_probability; // Normalize p
                    entropy += normalized_p * std::log(normalized_p); // log(p)
                }
            }
            entropy = -entropy;
        }
        return entropy;
    }
}

#endif