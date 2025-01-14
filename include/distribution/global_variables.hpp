#ifndef GLOBAL_VARIABLES_HPP
#define GLOBAL_VARIABLES_HPP

#include <atomic>
#include <unordered_map>
#include <condition_variable>
#include <mutex>
#include "communication/sync.hpp"
#include "communication/client.hpp"
#include <cmath>
namespace GlobalState {
    inline std::atomic<bool> keep_running(true);

    inline MutexableVariable<std::unordered_map<int, Client>> client_map;  // Maps sock -> Client
    inline MutexableVariable<std::unordered_map<int, int>> partiton_id_to_sock;  // Maps sock -> Client
    inline MutexableVariable<std::unordered_map<int, int>> sock_to_partition_id;  // Maps sock -> Client

    inline MutexableVariable<std::vector<double>> global_result(std::vector<double>(256, -INFINITY));

    inline std::condition_variable partition_prepared_msg_cv;
    inline MutexableVariable<int> partition_prepared_msg_ack_count;

    inline std::condition_variable begin_epoch_msg_cv;
    inline MutexableVariable<std::unordered_map<int, int>> begin_epoch_msg_ack_count;

    inline std::condition_variable epoch_completed_msg_cv;
    inline MutexableVariable<std::unordered_map<int, int>> epoch_completed_ack_count;

    inline std::condition_variable application_cv;
    inline std::mutex application_mutex;

    inline MutexableVariable<std::vector<double>> integrated_result(std::vector<double>(256, -INFINITY));

    inline std::condition_variable integrated_result_confirmation_cv;
    inline MutexableVariable<std::unordered_map<int, int>> integrated_result_confirmation_ack_count;

    
    inline std::condition_variable integrated_result_prepared_cv;
    inline MutexableVariable<std::unordered_map<int, int>> integrated_result_prepared_ack_count;
    
    inline int cnt_grammar = -1;
}

#endif
