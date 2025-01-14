#include "statistics/statistics.hpp"

namespace statistics{
    double Statistician::mutual_information(const std::vector<uint32_t>& seq_X, const std::vector<uint32_t>& seq_Y) {
        if(seq_X.size() != seq_Y.size()){
            std::cout << "Warning: X and Y has different length in the computation of mutual information. " << std::endl;
        }
        double H_X = calculate_entropy_from_element_sequence(seq_X);
        double H_Y = calculate_entropy_from_element_sequence(seq_Y);
        double H_XY = calculate_joint_entropy(seq_X, seq_Y);
        
        return H_X + H_Y - H_XY;
    }

}