#ifndef PACKAGEMANAGER_HPP
#define PACKAGEMANAGER_HPP
#include <unordered_map>
#include <mutex>
#include <iostream>
#include "distribution/communication/sync.hpp"
#include "distribution/communication/package.hpp"
#define RED     "\033[31m"      /* Red */
#define RESET   "\033[0m"
#define GREEN   "\033[32m"      /* Green */

class PackageManager {
    private:
        static PackageManager* manager;
        static std::mutex mtx;
        PackageManager(){}
        MutexableVariable<std::unordered_map<int, int>> next_seq_numbers;

    public:
        PackageManager(const PackageManager&) = delete; // Delete copy constructor
        PackageManager& operator=(const PackageManager&) = delete; // Delete copy assignment

        static PackageManager* get_instance() {
            static PackageManager instance;
            return &instance;
        }
    
        Package pack_msg_to_package(const int& sock, const Message& msg){
            Package package;
            next_seq_numbers.access_with_function([&package, &sock](auto& map)->void{
                package.sequence_number = map[sock];
                map[sock]++;
                std::cout << RED << "map[" << sock << "]" << map[sock] << RESET << std::endl;
            });
            package.msg = msg;
            return package;
        }
};
#endif
