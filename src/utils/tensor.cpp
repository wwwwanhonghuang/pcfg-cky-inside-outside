#include "utils/tensor.hpp"

template<typename T>
T tensor<T>::operator[](int index){
    return this->tensor_pt[index];
};



template<typename T>
void tensor_printer::print(const tensor<T>& tensor){
    if (tensor.mode.empty()) {
            return;
    }
    std::cout << "[";
    T* pt = tensor.tensor_pt;
    
    if(tensor.mode.size() == 1){
        int length = tensor.mode[0];
        for(int i = 0; i < length; i++){
            std::cout << *(pt + i) << ",";
            if (i < length - 1) std::cout << ", ";
        }
    }else{
        int length_sum = 1;
        std::vector<int> submode_vector = std::vector<int>(tensor.mode + 1, tensor.mode.end());
        for(int i = 1; i < tensor.mode.size(); i++){
            length_sum *= tensor.mode[i];
        }
        T* pt_end = tensor.tensor_pt + length_sum * tensor.mode[0];

        for(T* current_pt = pt; pt < pt_end; pt += length_sum){
            this->print(tensor(current_pt, submode_vector));
            if (current_pt + length_sum < pt_end) std::cout << ", "; 
        }
    }
    std::cout << "]" << std::endl;
}


