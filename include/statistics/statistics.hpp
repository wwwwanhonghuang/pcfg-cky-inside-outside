#ifndef USE_CUDA
#ifndef H_STATISTICS
#define H_STATISTICS

#include <string>
#include <vector>
#include <cmath>
#include <queue>
#include <unordered_map>

//#include <functional>
#include <utility>

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
        static double sequence_transfer_entropy_delay_L(const std::vector<T>& sequence, int L){
            std::vector<double> Y_t;
            std::vector<double> X_t;
            std::vector<double> Z_t;
            double entropy = 0.0;
            if(sequence.size() <= L) return 0;
            for(int t = L; t < sequence.size() - 1; t++){
                Y_t.emplace_back(sequence[t]);
                X_t.emplace_back(sequence[t - L]);
                Z_t.emplace_back(sequence[t - L - 1]);
            }

            double H_Y_given_X = conditional_entropy(Y_t, X_t);
            double H_Y_given_XY = conditional_entropy_Z_given_XY(Y_t, X_t, Z_t, 1);
            entropy = H_Y_given_X - H_Y_given_XY;
            return entropy;

        }

        template<typename T>
        static double calculate_average_kl_divergence(const std::vector<std::vector<T>>& layers) {
            double total_kl = 0.0;
            int num_layers = layers.size();
            int n_samples = 0;
            // Calculate KL divergence for each layer and accumulate
            for (int i = 0; i < num_layers; ++i) {
                for(int j = i + 1; j < num_layers; ++j){
                    // Calculate distribution for the current layer
                    std::unordered_map<T, double> layer_distribution_i = calculate_distribution(layers[i]);
                    std::unordered_map<T, double> layer_distribution_j = calculate_distribution(layers[j]);

                    // Calculate KL divergence between the layer and the reference distribution
                    total_kl += calculate_kl_divergence(layer_distribution_i, layer_distribution_j);
                    n_samples++;
                }
            }

            // Return the average KL divergence across layers
            return total_kl / n_samples;
        }

        template<typename T>
        static double calculate_kl_divergence(const std::unordered_map<T, double>& p, const std::unordered_map<T, double>& q) {
            double kl_divergence = 0.0;

            for (const auto& pair : p) {
                T val = pair.first;
                double p_val = pair.second;
                double q_val = q.count(val) > 0 ? q.at(val) : 0.0;

                if (p_val > 0 && q_val > 0) {
                    kl_divergence += p_val * std::log(p_val / q_val);
                }
            }

            return kl_divergence;
        }

        template<typename T>
        static std::unordered_map<T, double> calculate_distribution(const std::vector<T>& values) {
            std::unordered_map<T, int> counts;
            for (T val : values) {
                counts[val]++;
            }

            std::unordered_map<T, double> distribution;
            int total = values.size();
            for (const auto& pair : counts) {
                distribution[pair.first] = static_cast<double>(pair.second) / total;
            }

            return distribution;
        }
        
        template<typename T>
        static double calculate_skewness(const std::vector<T>& values) {
            int n = values.size();
            if (n < 3) return 0.0;  // Skewness is not defined for n < 3

            // Calculate mean
            double mean = 0.0;
            for (T val : values) {
                mean += val;
            }
            mean /= n;

            // Calculate standard deviation
            double variance = 0.0;
            for (T val : values) {
                variance += (val - mean) * (val - mean);
            }
            variance /= n;
            double stddev = std::sqrt(variance);

            // Calculate skewness
            double skewness = 0.0;
            for (T val : values) {
                skewness += std::pow((val - mean) / stddev, 3);
            }
            skewness *= (n / ((n - 1) * (n - 2)));
            
            return skewness;
        }

        static double calculate_average_skewness(const std::vector<std::vector<int>>& layers) {
            double total_skewness = 0.0;
            int num_layers = layers.size();

            // Calculate skewness for each layer and accumulate
            for (int i = 0; i < num_layers; ++i) {
                total_skewness += calculate_skewness(layers[i]);
            }

            // Return the average skewness across layers
            return total_skewness / num_layers;
        }

        template <typename T>
        static double layer_symbol_transfer_entropy_delay_L(
            std::vector<std::vector<T>> layers, 
            int L
        ) {
            return _layer_average_transfer_entropy_delay_L(
                layers,
                L
            );
        }

        template <typename T>
        static double layer_derivation_transfer_entropy_delay_L(
            std::vector<std::vector<T>> layers, 
            int L
        ) {
            return _layer_average_transfer_entropy_delay_L(
                layers,
                L
            );
        }

        
        static double calculate_entropy_from_probabilities(const std::vector<double>& probabilities)
        {
            double entropy = 0.0;
            const double epsilon = 1e-10;  // Small value to avoid log(0)

            for (auto& probability : probabilities) {
                if (probability > epsilon) {
                    entropy -= probability * log(probability);
                }
            }

            return entropy;
        }
        template<typename T>
        static double calculate_entropy_from_element_sequence(const std::vector<T>& sequence)
        {
            std::unordered_map<T, int> count_map;

            for (int element : sequence) {
                count_map[element]++;
            }
            double entropy = 0.0;
            int total_count = sequence.size();
            
            for (auto& entry : count_map) {
                double probability = static_cast<double>(entry.second) / total_count;
                entropy -= probability * log(probability);
            }
            
            return entropy;
        }

        template<typename T>
        static double calculate_joint_entropy(const std::vector<T>& seq_X, const std::vector<T>& seq_Y)
        {
            if (seq_X.size() != seq_Y.size()) {
                throw std::invalid_argument("Sequences must have the same size.");
            }

            std::unordered_map<std::string, uint32_t> joint_count_map;
 
            // Count joint occurrences of (X_i, Y_i)
            for (size_t i = 0; i < seq_X.size(); ++i) {
                joint_count_map[std::to_string(seq_X[i]) + std::string(",") + std::to_string(seq_Y[i])]++;
            }

            double joint_entropy = 0.0;
            int total_count = seq_X.size();  // Total count is the number of pairs

            // Calculate joint entropy
            for (auto& entry : joint_count_map) {
                double probability = static_cast<double>(entry.second) / total_count;
                joint_entropy -= probability * log(probability);
            }

            return joint_entropy;
        }

        static double mutual_information(const std::vector<uint32_t>& seq_X, const std::vector<uint32_t>& seq_Y);
        static double word_transitional_entropy_delay_L(std::vector<uint32_t> words, int L);
        static double derivation_transitional_entropy_delay_L(std::vector<uint32_t> derivations, int L);
        




        template<typename T>
        static double conditional_entropy_Z_given_XY(const std::vector<T>& z, const std::vector<T>& x, const std::vector<T>& y, int num_bins) {
            std::unordered_map<std::string, int> joint_counts;
            if (!(z.size() == x.size() && z.size() == y.size() && y.size() == x.size())) {
                throw std::invalid_argument("Sequences must have the same size.");
            }

            for (size_t i = 0; i < z.size(); ++i) {
                T zz = z[i];
                T xx = x[i];
                T yy = y[i];

                std::string key = std::to_string(zz) + "," + std::to_string(xx) + "," + std::to_string(yy);
                joint_counts[key]++;
            }

            double total_count =  z.size();
            std::vector<double> joint_probs;
            for (const auto& pair : joint_counts) {
                joint_probs.push_back(pair.second / total_count);
            }

            // Calculate entropy H(Y_t+1 | Y_t-k, X_t-k)
            double joint_entropy_xyz  = calculate_entropy_from_probabilities(joint_probs);
            double joint_entropy_yz = calculate_joint_entropy(y, x);
            return joint_entropy_xyz - joint_entropy_yz;
        }

        template<typename T>
        static double _layer_average_transfer_entropy_delay_L(
            std::vector<std::vector<T>> layers, 
            int L
        ) {
            if (L <= 1) return 0.0;
            
            double entropy = 0.0;

            for (int layer_id = 0; layer_id < layers.size() - L - 1; layer_id++) {
                const auto& X = layers[layer_id];
                const auto& Y = layers[layer_id + L];
                const int lag = 1;
                std::vector<T> X_past(X.begin(), X.end() - lag);  // Past values of X (up to time t-k)
                std::vector<T> Y_past(Y.begin(), Y.end() - lag);  // Past values of Y (up to time t-k)
                std::vector<T> Y_future(Y.begin() + lag, Y.end()); // Future values of Y (at time t+1)


                // H(Y_{t + 1} | Y_t^{t-k})
                double H_Y_given_Y = conditional_entropy(Y_future, Y_past); 

                // H(Y_{t + 1} | X_t^{t-k}, Y_t^{t-k})
                double H_Y_given_XY = conditional_entropy_Z_given_XY(Y_future, Y_past, X_past, 1);
                entropy += H_Y_given_Y - H_Y_given_XY;
            }
            return entropy / (layers.size() - L);
        }


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

        static double derivation_entropy(std::vector<uint32_t> derivations);

        static double word_entropy(std::vector<uint32_t> words);

        static double _sequence_delay_L_mutual_entropy(std::vector<uint32_t> words, int L);

        static double derivation_delay_L_mutual_entropy(std::vector<uint32_t> derivations, int L);

        static std::vector<uint32_t> to_derivations_by_preorder_iteration(parsing::SyntaxTreeNode* node);

        template<typename T>
        static double conditional_entropy(const std::vector<T>& Y, const std::vector<T>& X){
            std::map<std::pair<int, int>, int> joint_freq;
            std::map<int, int> x_freq;
            std::map<int, int> y_freq;

            for (size_t i = 0; i < X.size(); ++i) {
                joint_freq[{X[i], Y[i]}]++;
                x_freq[X[i]]++;
                y_freq[Y[i]]++;
            }

            double joint_entropy = 0.0;
            double x_entropy = 0.0;
            double y_entropy = 0.0;
            double total = X.size();

            // Calculate joint entropy H(X, Y)
            for (const auto& entry : joint_freq) {
                double p = entry.second / total;
                joint_entropy -= p * std::log(p);
            }

            // Calculate entropy of X, H(X)
            for (const auto& entry : x_freq) {
                double p = entry.second / total;
                x_entropy -= p * std::log(p);
            }

            

            // Calculate entropy of Y, H(Y)
            for (const auto& entry : y_freq) {
                double p = entry.second / total;
                y_entropy -= p * std::log(p);
            }
            return joint_entropy - x_entropy;
        }


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
        bfs_get_all_layers_value(pcfg* grammar, parsing::SyntaxTreeNode* root, std::function<T(parsing::SyntaxTreeNode*)> value_selector){
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

                    layer.emplace_back(value_selector(first_node)); // the A'id in A->BC/
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
        static double _calculate_layer_average_entropy(parsing::SyntaxTreeNode* node, 
            const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers,
            std::function<T(parsing::SyntaxTreeNode*)> value_extractor) {
            double total_entropy = 0.0;
            int layer_id = 0;
            for(auto& layer : layers){
                std::vector<int> values;
                values.resize(layer.size());
                std::transform(layer.begin(), layer.end(), values.begin(), value_extractor);
                double layer_entropy = _sequence_entropy(values);
                total_entropy += layer_entropy;
                layer_id ++;
            }
            return layer_id > 0 ? total_entropy / layer_id : 0.0;
        }

        
        template<typename T>
        static double _calculate_path_average_entropy(parsing::SyntaxTreeNode* node, std::function<T(parsing::SyntaxTreeNode*)> value_extractor) {
            double total_entropy = 0.0;
            int path_count = 0;

            if (node == nullptr) return total_entropy;

            std::vector<parsing::SyntaxTreeNode*> stack;
            stack.push_back(node);
            std::unordered_set<parsing::SyntaxTreeNode*> black_set;

            auto is_black = [&black_set](parsing::SyntaxTreeNode* node) -> bool {
                return node == nullptr || black_set.find(node) != black_set.end();
            };
            std::vector<int> path_values;  


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

        
        template<typename T>
        static std::vector<std::vector<T>> get_paths(parsing::SyntaxTreeNode* node, std::function<T(parsing::SyntaxTreeNode*)> value_extractor) {
            std::vector<std::vector<T>> paths;

            if (node == nullptr) return paths;

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
                    std::vector<T> values;

                    values.resize(stack.size());
                    std::transform(stack.begin(), stack.end(), values.begin(), value_extractor);

                    paths.emplace_back(values);
                    stack.pop_back();
                    black_set.emplace(peek_node);
                    continue;
                }

                if (!is_black(peek_node->left)) {
                    stack.push_back(peek_node->left);
                }
                else if (!is_black(peek_node->right)) {
                    stack.push_back(peek_node->right);
                } 
                else {
                    black_set.emplace(peek_node);
                    stack.pop_back();
                }
            }
            return paths;
        }
        

        static double calculate_layer_average_derivation_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers);
        static double calculate_path_average_derivation_entropy(parsing::SyntaxTreeNode* node);

        static int tree_depth(parsing::SyntaxTreeNode* node);
        static double calculate_path_average_symbol_entropy(parsing::SyntaxTreeNode* node);
        static double calculate_layer_average_symbol_entropy(parsing::SyntaxTreeNode* node, const std::vector<std::vector<parsing::SyntaxTreeNode*>> layers);
        
        
        static double L_layer_derivation_tree_mutual_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L);
        static double L_layer_symbol_tree_mutual_entropy(pcfg* grammar, parsing::SyntaxTreeNode* tree, int L);

        // !important
        static double prefix_L_parse_entropy(pcfg* grammar, double* alpha, int sequence_length, int end, int L, uint32_t* sequence);
        static double word_delay_L_mutual_entropy(std::vector<uint32_t> words, int L);

        static std::string report_all_statistics(parsing::SyntaxTreeNode* node, 
            double* alpha, std::vector<uint32_t> sentence, pcfg* grammar, int max_delays);

        // Metric : Average Path Length
        static int calculate_height(parsing::SyntaxTreeNode* node);
        static int calculate_total_depth(parsing::SyntaxTreeNode* node, std::vector<std::vector<parsing::SyntaxTreeNode*>>& paths);

        static int count_nodes(parsing::SyntaxTreeNode* node);

        static double calculate_average_path_length(parsing::SyntaxTreeNode* node, std::vector<std::vector<parsing::SyntaxTreeNode*>>& paths);

        static bool areSubtreesSimilar(parsing::SyntaxTreeNode* node1, parsing::SyntaxTreeNode* node2);

        static int calculateRedundancy(parsing::SyntaxTreeNode* node);

        static double calculateEntropy(parsing::SyntaxTreeNode* node);

        static int countTraversalSteps(parsing::SyntaxTreeNode* node);

        static double calculateSkewness(parsing::SyntaxTreeNode* node);

        static int calculate_max_nodes(int height);

        static double calculateDensity(parsing::SyntaxTreeNode* node);

        static double calculate_tree_symbol_entropy(parsing::SyntaxTreeNode* node);

    private:
        static double _sequence_delay_L_transitional_entropy(std::vector<uint32_t> sequence, int delay_L);
    };


}








#endif
#endif