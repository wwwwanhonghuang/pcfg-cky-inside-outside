#include "distribution/services.hpp"
#include "distribution/common/message.h"
#include "distribution/global_variables.hpp"
#include "distribution/communication/package.hpp"
#include <iostream>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/mman.h>
#include <memory>
#include <thread>
#include <cstring>
#include <fcntl.h>
#include <queue>
#include <unordered_map>
#include <cassert>
using namespace GlobalState;

#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define RESET   "\033[0m"

struct ComparePackage {
    bool operator()(const Package& p1, const Package& p2) {
        return p1.sequence_number > p2.sequence_number;
    }
};

MutexableVariable<std::priority_queue<Package, std::vector<Package>, ComparePackage>> saved_packages;

MutexableVariable<int> expect_seq_number;

inline int get_expect_seq_number(){
    return expect_seq_number.get();
}

inline void increase_package_seq(){
    expect_seq_number.access_with_function([](auto& v)->void{
        v++;
        // std::cout << GREEN << "next expect seq = " << v << std::endl;
    });
}



void save_package(const Package& package) {
    saved_packages.value.push(package);
}

const Package& most_early_one() {
    return saved_packages.value.top();
}

void remove_most_early_one() {
    saved_packages.value.pop();
}
bool has_saved_packages() {
    return !saved_packages.value.empty();
}

void process(const Package& package){
    Message msg_receive = package.msg;
    std::cout << BLUE << "[Process] seq = " << package.sequence_number 
              << " msg_type = " 
              << msg_receive.msg_type << "\n" << RESET;

    switch (msg_receive.msg_type)
    {
    case PARTITION_PREPARED: {
        {
            partition_prepared_msg_ack_count.access_with_function([](auto& v)->void{ v++;});
            std::cout << "ACK PARTITION_PREPARED" << " partition_prepared_msg_ack_count = " << 
                        partition_prepared_msg_ack_count.value << std::endl;
        }
        partition_prepared_msg_cv.notify_one(); 
        break;
    }
    case EPOCH_COMPLETE: {
        int epoch = -1;
        int client_partition_id = -1;
        int client_result = -1;
        memcpy(&client_partition_id, msg_receive.data, sizeof(int));
        memcpy(&epoch, &msg_receive.data[4], sizeof(int));
        memcpy(&client_result, &msg_receive.data[8], sizeof(int));

        std::cout << "[Client Service] client " <<  client_partition_id << " has completed epoch " << epoch <<
                    " reported result = " << client_result << std::endl;
        {
            global_result.access_with_function([&client_result](auto& v)->void{v += client_result;});
        }

        {
            epoch_completed_ack_count.access_with_function([&epoch](auto& map)->void{
                map[epoch]++;
            });
            std::cout << "ACK EPOCH_COMPLETE" << " epoch_completed_ack_count[" << 
                epoch << "] = " << epoch_completed_ack_count.get()[epoch] << std::endl;
        }
        
        std::cout << "[Client Service] Notify one thread that waiting for epoch_completed_msg_cv" << std::endl;
        epoch_completed_msg_cv.notify_one();
        break;
    }
    case BEGIN_EPOCH:{
        int epoch = -1;
        int client_partition_id = -1;
        memcpy(&client_partition_id, &msg_receive.data[4], sizeof(int));
        memcpy(&epoch, &msg_receive.data[0], sizeof(int));
        std::cout << "[Client Service] client " <<  client_partition_id << " prepare proceed epoch " << epoch << std::endl;
        {
            begin_epoch_msg_ack_count.access_with_function([&epoch](auto& map)->void{map[epoch] ++;} );

            std::cout << "ACK BEGIN_EPOCH" << " begin_epoch_msg_ack_count[" << epoch <<
                "] = " << begin_epoch_msg_ack_count.get()[epoch] << std::endl;
        }
        
        std::cout << "[Client Service] Notify one thread that waiting for begin_epoch_msg_cv" << std::endl;
        begin_epoch_msg_cv.notify_one();
        break;
    }
    case NOTIFICATE_INTEGRATE_RESULT: {
        int client_result = 0;
        int result_this_partition = integrated_result.get();
        memcpy(&client_result, &msg_receive.data[0], sizeof(int));
        if(client_result != result_this_partition){
            std::cout << "Error: client reported result " << client_result << " not equal to "
                    << " this partition's result " << result_this_partition << std::endl;
            assert(false);
        }
        std::cout << "[Client Service] Notify one thread that waiting for integrated_result_msg_cv" << std::endl;
        begin_epoch_msg_cv.notify_one();
        break;
    }
    case INTEGRATED_RESULT_PREPARED: {
        {
            int epoch = -1;
            int client_partition_id = -1;
            memcpy(&epoch, &msg_receive.data[0], sizeof(int));
            memcpy(&client_partition_id, &msg_receive.data[4], sizeof(int));
            integrated_result_prepared_ack_count.access_with_function([&epoch](auto& map)->void{
                map[epoch]++;
                std::cout << "ACK INTEGRATED_RESULT_PREPARED" << " integrated_result_prepared_cv = " << 
                        map[epoch] << std::endl;
            });
        
        }
        integrated_result_prepared_cv.notify_one(); 
        break;
    }
    default:
        break;
    }
}
#define TOTAL_CLIENTS 3
bool should_ignore(int partition_id, int sock, int seq_num){
    return (partition_id  - 1)== (seq_num % TOTAL_CLIENTS);
}
void handle_client(int client_sock, int partition_id) {
    std::cout << "Handling client in thread: " << std::this_thread::get_id() << "\n";
    

    fcntl(client_sock, F_SETFL, O_NONBLOCK);

    while (true) {
        // Fetch the last processed package sequence number for this client
        int expect_seq_number = get_expect_seq_number();
        if(should_ignore(partition_id, client_sock, expect_seq_number)){
            increase_package_seq();
            std::cout << YELLOW << "partition " << partition_id 
            << " ignore package with seq = "
            << expect_seq_number << RESET << std::endl;
            continue;
        }
        // Process saved packages in order if they match the expected sequence
        if (has_saved_packages()){
            const Package& earliest_package = most_early_one();
            if(earliest_package.sequence_number ==  expect_seq_number) {
                process(most_early_one());
                remove_most_early_one(); 
                increase_package_seq();
            }
        }

        Package package_receive;

        // Check if a new package has arrived
        ssize_t bytes_read = read(client_sock, &package_receive, sizeof(Package));
        if (bytes_read == 0) {
            std::cout << "Client disconnected: sock " << client_sock << "\n";
            close(client_sock);
            break; 
        } 

        if (bytes_read > 0) {
            std::cout << "Received: package seq = " << package_receive.sequence_number
                << " from sock " << client_sock  
                << " msg_type = " << package_receive.msg.msg_type << "\n";

            int seq_num = package_receive.sequence_number;
            expect_seq_number = get_expect_seq_number();

            if (seq_num == expect_seq_number) {
                // Process the package if it's the next in sequence
                process(package_receive);
                increase_package_seq();
            } else {
                // Save the out-of-order package for later processing
                save_package(package_receive);
            }
        }else if(bytes_read < 0) continue;
    }
    



}
