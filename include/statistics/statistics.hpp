#ifndef USE_CUDA
#ifndef H_STATISTICS
#define H_STATISTICS

#include <string>
#include <vector>
#include <cmath>
#include <queue>
//#include <functional>
#include <algorithm>
#include <deque>
#include <unordered_set>
#include "algorithms/tree_parser.hpp"
#include "grammar/grammar.hpp"
#include "macros.def"

namespace statistics{
    class Statistician{
    public:
        template<typename T>
        static double _sequence_entropy(std::vector<T> sequence){
            std::map<uint32_t, uint32_t> counter;
            uint32_t Z = 0;
            
            for(auto&& element : sequence){
                counter[element]++;
            }

            for(auto& map_item : counter){
                Z += map_item.second;
            }

            double entropy = 0.0;
            for(auto& map_item : counter){
                double v = (Z == 0.0 ? 0.0 : static_cast<double>(map_item.second) / Z);
                entropy += (v == 0 ? 0 : -v * std::log(v));
            }
            return entropy;
        }

        // entropy
        static double derivation_entropy(std::vector<uint32_t> derivations);

        static double word_entropy(std::vector<uint32_t> words);

        static double _sequence_delay_L_mutual_entropy(std::vector<uint32_t> words, int L);

        static double derivation_delay_L_mutual_entropy(std::vector<uint32_t> derivations, int L);

        static std::vector<uint32_t> to_derivations_by_preorder_iteration(parsing::SyntaxTreeNode* node);

        // Function to calculate mutual information with a delay of L layers
        template<typename T>
        static double calculate_delay_L_layer_mutual_information(pcfg* grammar, std::vector<std::vector<T>> layers, int L){
            std::unordered_map<uint64_t, long> joint_counter;
            std::unordered_map<uint32_t, long> symbol_counter;

            int size_layers = layers.size();
            // counting occurances
            for(int pre_layer_id = 0; pre_layer_id < size_layers - L; pre_layer_id++){
                int delay_L_layer_id = pre_layer_id + L;
                // std::cout << " -- " << pre_layer_id << " " << std::endl;
                for(T prelayer_element: layers[pre_layer_id]){
                    // Count occurrences in the joint distribution
                    symbol_counter[prelayer_element]++;

                    for(auto&& postlayer_element: layers[delay_L_layer_id]){
                        uint64_t key = ((prelayer_element << 16) & 0xFFFF0000) | (postlayer_element & 0xFFFF);
                        // std::cout << " --- pre = " << prelayer_element << 
                        // "  --- post = " << postlayer_element << std::endl;
                        joint_counter[key]++;
                    }
                }
            }

            for(int layer_id = layers.size() - L; layer_id < layers.size(); layer_id++){
                for(auto&& layer_element: layers[layer_id]){
                    symbol_counter[layer_element]++;
                }
            }

            // Normalization to possibility
            std::map<uint64_t, double> joint_possibility;
            std::map<uint32_t, double> symbol_possibility;
            double Z = 0.0;
            for(auto& map_item : joint_counter){
                Z += map_item.second;
            }

            for(auto& map_item : joint_counter){
                joint_possibility[map_item.first] = (Z == 0.0 ? 0.0 : static_cast<double>(map_item.second) / Z);
            }

            Z = 0.0;
            
            for(auto& map_item : symbol_counter){
                Z += map_item.second;
            }

            for(auto& map_item : symbol_counter){
                symbol_possibility[map_item.first] = (Z == 0.0 ? 0.0 : static_cast<double>(map_item.second) / Z);
            }

            // lacks of epsilon
            // mutual entropy
            double mutual_entropy = 0.0;
            for(uint32_t symbol1 = 0; symbol1 <= grammar->n_syms(); symbol1++){
                for(uint32_t symbol2 = 0; symbol2 <= grammar->n_syms(); symbol2++){
                    if(symbol1 == grammar->n_syms()) symbol1 = 0xFFFF;
                    if(symbol2 == grammar->n_syms()) symbol2 = 0xFFFF;

                    double p_A = symbol_possibility.find(symbol1) == symbol_possibility.end() ? 0.0 : symbol_possibility.find(symbol1)->second;
                    double p_B = symbol_possibility.find(symbol2) == symbol_possibility.end() ? 0.0 : symbol_possibility.find(symbol2)->second;
                    double p_AB = 
                        (joint_possibility.find(symbol1 << 16 | symbol2) == joint_possibility.end() ? 0.0 : joint_possibility.find(symbol1 << 16 | symbol2)->second) +
                        (joint_possibility.find(symbol2 << 16 | symbol1) == joint_possibility.end() ? 0.0 : joint_possibility.find(symbol2 << 16 | symbol1)->second);
                    bool is_zero = p_A == 0 || p_B == 0 || p_AB == 0;
                    if(!is_zero)
                        mutual_entropy +=  p_AB * log(p_AB / (p_A * p_B)); 
                }
            }
            return mutual_entropy;
        }

        template<typename T> 
        static std::vector<std::vector<T>> 
        bfs_get_all_layers_value(pcfg* grammar, parsing::SyntaxTreeNode* root, std::function<T(const std::tuple<uint32_t, uint32_t, uint32_t, int, double, int>&)> value_selector){
            std::vector<std::vector<T>> layers;
            std::queue<parsing::SyntaxTreeNode*> node_queue;
            
            if (root == nullptr) {
                return layers;
            }

            // BFS to select all layers' root's values.
            node_queue.push(root);
            int layer_id = 0;
            while(!node_queue.empty()){
                std::vector<T> layer;
                int size = node_queue.size();
                // std::cout << "layer " << layer_id << ": " << std::endl;
                layer_id++;
                for(int i = 0; i < size; i++){
                    parsing::SyntaxTreeNode* first_node = node_queue.front();
                    node_queue.pop();

                    layer.emplace_back(value_selector(first_node->value)); // the A'id in A->BC/
                    // std::cout << "push - " << value_selector(first_node->value) << std::endl;
                    if(first_node->left != nullptr){
                        node_queue.push(first_node->left);
                    }
                    if(first_node->right != nullptr){
                        node_queue.push(first_node->right);
                    }
                }
                layers.emplace_back(layer);
            }
            return layers;
        }
        
        template<typename T>
        static double _calculate_layer_average_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers,
            std::function<T(parsing::SyntaxTreeNode*)> value_extractor) {
            double total_entropy = 0.0;
            int path_count = 0;
            for(auto& layer : layers){
                std::vector<int> values;
                values.resize(layer.size());
                std::transform(layer.begin(), layer.end(), values.begin(), value_extractor);

                total_entropy += _sequence_entropy(values);
            }
            return path_count > 0 ? total_entropy / path_count : 0.0;
        }

        
        template<typename T>
        double _calculate_path_average_entropy(parsing::SyntaxTreeNode* node, std::function<T(parsing::SyntaxTreeNode*)> value_extractor) {
            double total_entropy = 0.0;
            int path_count = 0;

            if (node == nullptr) return total_entropy;

            std::vector<parsing::SyntaxTreeNode*> stack;
            stack.push_back(node);
            std::unordered_set<parsing::SyntaxTreeNode*> black_set;

            auto is_black = [&black_set](parsing::SyntaxTreeNode* node) -> bool {
                return node == nullptr || black_set.find(node) != black_set.end();
            };

            while (!stack.empty()) {
                parsing::SyntaxTreeNode* peek_node = stack.back();

                // If it's a leaf node, calculate the entropy for the current path
                if (peek_node->is_leaf()) {
                    std::vector<int> values;

                    values.resize(stack.size());

                    std::transform(stack.begin(), stack.end(), values.begin(), value_extractor);

                    total_entropy += _sequence_entropy(values);
                    path_count++;
                    stack.pop_back();
                    black_set.emplace(peek_node);
                    continue;
                }

                // Traverse the left child if not visited
                if (!is_black(peek_node->left)) {
                    stack.push_back(peek_node->left);
                }
                // Otherwise, traverse the right child if not visited
                else if (!is_black(peek_node->right)) {
                    stack.push_back(peek_node->right);
                } 
                // If both children are visited, mark the current node as processed
                else {
                    black_set.emplace(peek_node);
                    stack.pop_back();
                }
            }

            // Prevent division by zero if no paths were processed
            return path_count > 0 ? total_entropy / path_count : 0.0;
        }
        

        static double calculate_layer_average_derivation_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers);
        static double calculate_path_average_derivation_entropy(parsing::SyntaxTreeNode* node);

        static int tree_depth(parsing::SyntaxTreeNode* node);
        static double calculate_path_average_symbol_entropy(parsing::SyntaxTreeNode* node);
        static double calculate_layer_average_symbol_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers);
        static double L_layer_symbol_tree_transitional_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L);
        static double L_layer_derivation_tree_transitional_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L);
        static double L_layer_derivation_tree_mutual_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L);
        static double L_layer_symbol_tree_mutual_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L);

        // !important
        static double prefix_L_parse_entropy(pcfg* grammar, double* alpha, int sequence_length, int end, int L, uint32_t* sequence);
        static double word_delay_L_mutual_entropy(std::vector<uint32_t> words, int L);

        static std::string report_all_statistics(parsing::SyntaxTreeNode* node, 
            double* alpha, std::vector<uint32_t> sentence, pcfg* grammar, int max_delays);

        // Metric : Average Path Length
        static int calculateHeight(parsing::SyntaxTreeNode* node);
        static int calculateTotalDepth(parsing::SyntaxTreeNode* node, int depth);

        static int countNodes(parsing::SyntaxTreeNode* node);

        static double averagePathLength(parsing::SyntaxTreeNode* node);

        // Metric : Redundancy (Measuring Subtree Similarity)
        static bool areSubtreesSimilar(parsing::SyntaxTreeNode* node1, parsing::SyntaxTreeNode* node2);

        static int calculateRedundancy(parsing::SyntaxTreeNode* node);

        // Metric : Entropy (Tree Disorder/Complexity)
        static double calculateEntropy(parsing::SyntaxTreeNode* node);

        // Metric : Traversal Time (Steps taken to traverse)
        static int countTraversalSteps(parsing::SyntaxTreeNode* node);

        // Metric : Skewness (Imbalance of the tree)
        static double calculateSkewness(parsing::SyntaxTreeNode* node);

        // Metric : Tree Density (How full the tree is)
        static int calculateMaxNodes(int height);

        static double calculateDensity(parsing::SyntaxTreeNode* node);

        static double calculate_tree_symbol_entropy(parsing::SyntaxTreeNode* node);

    private:
        static double _sequence_transitional_entropy(std::vector<uint32_t> sequence);
    };


}








#endif
#endif