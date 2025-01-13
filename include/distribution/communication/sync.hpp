#ifndef SYNCHRONIZATION_HPP
#define SYNCHRONIZATION_HPP
#include <functional>
#include <mutex>

template<typename T>
struct MutexableVariable{
    explicit MutexableVariable(T initial_value = T{}) : value(initial_value) {}

    std::mutex mutex_variable;
    T value;
    // static MutexableVariable<T> pack(T value){
    //     MutexableVariable<T> instance;
    //     instance.value = value;
    //     return instance;
    // }
    std::unique_lock<std::mutex> lock(){
        return std::unique_lock<std::mutex>(mutex_variable);
    }
    void access_with_function(std::function<void(T&)> func){
        std::lock_guard<std::mutex> lock_guard(mutex_variable);
        func(value);
    }
    T get(){
        std::lock_guard<std::mutex> lock(mutex_variable);
        return value;
    }
    void set(const T& new_value) {
        std::lock_guard<std::mutex> lock(mutex_variable);
        value = new_value;
    }
};
#endif