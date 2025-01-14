#include "statistics/statistics.hpp"

namespace statistics{
    double Statistician::calculate_layer_average_derivation_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers){
        return _calculate_layer_average_entropy<int>(node, layers, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<TREE_VALUE_INDEX_DERIVATION>(_node->value);});
    }
    
    double Statistician::calculate_average_path_length(parsing::SyntaxTreeNode* node, std::vector<std::vector<parsing::SyntaxTreeNode*>>& paths) {
        int total_depth = calculate_total_depth(node, paths);
        return static_cast<double>(total_depth) / paths.size();
    }

    int Statistician::calculate_total_depth(parsing::SyntaxTreeNode* node, std::vector<std::vector<parsing::SyntaxTreeNode*>>& paths) {
        if (paths.size() == 0) return 0;
        int total_length = 0;
        for(int i = 0; i < paths.size(); i++){
            total_length += paths.size();
        }
        return total_length;
    }

    double Statistician::calculate_layer_average_symbol_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers) {
        return _calculate_layer_average_entropy<int>(node, layers, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<TREE_VALUE_INDEX_SYMBOL>(_node->value);});
    }

    double Statistician::calculate_path_average_symbol_entropy(parsing::SyntaxTreeNode* node) {
        return _calculate_path_average_entropy<int>(node, [&](parsing::SyntaxTreeNode* _node)->int{return std::get<TREE_VALUE_INDEX_SYMBOL>(_node->value);});
    }

    double Statistician::calculate_path_average_derivation_entropy(parsing::SyntaxTreeNode* node) {
        return _calculate_path_average_entropy<int>(node, [&](parsing::SyntaxTreeNode* _node)->
            int{return std::get<TREE_VALUE_INDEX_DERIVATION>(_node -> value);});
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
            int symbol_id = std::get<TREE_VALUE_INDEX_SYMBOL>(current_node->value);

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
}