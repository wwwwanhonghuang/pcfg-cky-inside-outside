#include "statistics/statistics.hpp"
#include <vector>

namespace statistics{
    
    int Statistician::calculate_height(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + std::max(calculate_height(node->left), calculate_height(node->right));
    }
    int Statistician::count_nodes(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + count_nodes(node->left) + count_nodes(node->right);
    }

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

    int Statistician::countTraversalSteps(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0;
        return 1 + countTraversalSteps(node->left) + countTraversalSteps(node->right);
    }

    double Statistician::calculateSkewness(parsing::SyntaxTreeNode* node) {
        if (node == nullptr) return 0.0;
        int leftHeight = calculate_height(node->left);
        int rightHeight = calculate_height(node->right);
        return static_cast<double>(leftHeight - rightHeight);  // Positive skew means left-heavy
    }

    int Statistician::calculate_max_nodes(int height) {
        return std::pow(2, height) - 1;
    }

    double Statistician::calculate_density(parsing::SyntaxTreeNode* node) {
        int node_count = count_nodes(node);
        int tree_height = calculate_height(node);
        int max_nodes = calculate_max_nodes(tree_height);
        return static_cast<double>(node_count) / max_nodes;
    }
    int Statistician::tree_depth(parsing::SyntaxTreeNode* node){
        if(node == nullptr) 
            return 0;
        int left_right_max_depth = std::max(tree_depth(node->left), tree_depth(node->right));
        return left_right_max_depth + 1;
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
                result.push_back(std::get<TREE_VALUE_INDEX_DERIVATION>(current->value));
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


}