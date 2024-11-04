#ifndef HPP_TENSOR
#define HPP_TENSOR

#include <iostream>
#include <vector>

template<typename T>
void print_vector(std::vector<T> vec){
    std::cout << "[";
    int size = vec.size();
    for (auto& e: vec) {
        std::cout << e << ",";
    }
    std::cout << "]" << std::endl;
};

template<class T>
struct tensor
{
    public:
        T* tensor_pt;
        std::vector<int> mode;
        T operator[](int offset);
        tensor(T* pt, const std::vector<int>& m) : tensor_pt(pt), mode(m) {}
};

struct tensor_printer{
    template<typename T>
    void print(const tensor<T>& tensor);
};
#endif