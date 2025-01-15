#include <iostream>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstring>
#include <memory>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <atomic>
#include <mutex>
#include <thread>
#include <yaml-cpp/yaml.h>
#include <condition_variable>
#include <functional>
#include "distribution/common/message.h"
#include "distribution/shared_memory.hpp"
#include "distribution/communication/client.hpp"
#include "distribution/communication/sync.hpp"
#include "distribution/global_variables.hpp"
#include "distribution/services.hpp"
#include "distribution/message_impl.h"
#include "distribution/communication/package.hpp"
#include "distribution/package_manager.hpp"
#include <cassert>
#include "utils/math.hpp"

#define RED     "\033[31m"      /* Red */
#define RESET   "\033[0m"
#define GREEN   "\033[32m"      /* Green */
#define BOLD_BLUE   "\033[1;34m"      /* Bold Blue */


using namespace GlobalState;


void push_msg_to_shared_memory_rr(Message& msg, std::shared_ptr<SharedMemory> shared_memory){
    static int pos = 0;
    auto storage = (MemoryStorage*) shared_memory->get_data();
   
    while(storage->application_messages[pos].status != EMPTY_SLOT){
        pos = (pos + 1) % MSG_COUNT;
    }

    std::cout << "find msg position " << pos << std::endl;
    memcpy(&storage->application_messages[pos], &msg, sizeof(msg));
}

void serverLoop(uint32_t port, int partition_id) {
    // Create a server socket
    int server_sock = socket(AF_INET, SOCK_STREAM, 0);
    if (server_sock < 0) {
        perror("Server socket creation failed");
        return;
    }

    // Configure the server address
    struct sockaddr_in server_addr = {};
    server_addr.sin_family = AF_INET;
    server_addr.sin_addr.s_addr = INADDR_ANY;
    server_addr.sin_port = htons(port);

    // Bind the socket to the specified port
    if (bind(server_sock, (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0) {
        perror("Bind failed");
        close(server_sock);
        return;
    }

    // Listen for incoming connections
    if (listen(server_sock, 5) < 0) {
        perror("Listen failed");
        close(server_sock);
        return;
    }

    std::cout << "Server is listening on port " << port << "\n";

    while (keep_running) {
        // Set the socket to non-blocking
        fd_set read_fds;
        struct timeval timeout = {1, 0}; // 1-second timeout
        FD_ZERO(&read_fds);
        FD_SET(server_sock, &read_fds);

        int activity = select(server_sock + 1, &read_fds, NULL, NULL, &timeout);
        if (activity < 0 && errno != EINTR) {
            perror("Select error");
            break;
        }

        if (activity > 0 && FD_ISSET(server_sock, &read_fds)) {
            // Accept a connection
            struct sockaddr_in client_addr = {};
            socklen_t client_len = sizeof(client_addr);
            int client_sock = accept(server_sock, (struct sockaddr*)&client_addr, &client_len);
            if (client_sock < 0) {
                perror("Accept failed");
                continue;
            }

            std::cout << "Connection accepted" << "client socket = " << client_sock << ", creating thread to handle client.\n";
            std::thread client_thread(handle_client, client_sock, partition_id);
            client_thread.detach();
        }
    }
}

void broadcast_message(const Message& message) {
    client_map.access_with_function([&message](auto& map)->void{
        for (const auto& [sock, client] : map) {
            Package package = PackageManager::get_instance()->pack_msg_to_package(client.sock, message);
            ssize_t bytes_sent = send(client.sock, &package, sizeof(package), 0);
            if (bytes_sent < 0) {
                perror(("Failed to send message to " + client.name).c_str());
            } else {
                std::cout << "\t[Broadcast package] seq = " 
                    << package.sequence_number
                    << ", message type = " << message.msg_type 
                    << ", broad cast to " << sock << std::endl;
            }
        }
    });
}

void broadcast_message(int base_seq_number, int partition_id, const Message& message) {
    client_map.access_with_function([&message, &base_seq_number, &partition_id](auto& map)->void{
        for (const auto& [sock, client] : map) {
            Package package;
            package.sequence_number = base_seq_number + partition_id - 1;
            package.msg = message;
            ssize_t bytes_sent = send(client.sock, &package, sizeof(package), 0);
            if (bytes_sent < 0) {
                perror(("Failed to send message to " + client.name).c_str());
            } else{
                std::cout << "\t[Broadcast package] seq = " 
                << "[b"
                << base_seq_number << " + " << partition_id
                << " - 1"
                << "] "
                << package.sequence_number
                << ", message type = " << message.msg_type 
                << ", broad cast to " << sock << std::endl;
            }
        }
    });
}

void connect_to_other_partitions(int& total_clients, int& connected_client, 
        int& client_index, const YAML::Node& clients, 
        int partition_id, const std::string& program_name){
    while(connected_client < total_clients) {
        const YAML::Node& client = clients[client_index];
        std::string name = client["name"].as<std::string>();
        int current_client_index = client_index;
        client_index = (client_index + 1) % total_clients;

        if (name == program_name) {
            continue;
        }
        
        client_map.lock();
        partiton_id_to_sock.lock();
        // Key is not ID but sock. Follow code not work. 
        if(partiton_id_to_sock.value.find(current_client_index) != partiton_id_to_sock.value.end()) continue;

        std::string ip = client["ip"].as<std::string>();
        uint32_t port = client["port"].as<uint32_t>();
        // std::cout << ip << " " << port << " " << name << std::endl;
        int sock = socket(AF_INET, SOCK_STREAM, 0);
        if (sock < 0) {
            perror("Socket creation failed");
            continue;
        }
        struct sockaddr_in client_addr = {};
        client_addr.sin_family = AF_INET;
        client_addr.sin_port = htons(port);
        inet_pton(AF_INET, ip.c_str(), &client_addr.sin_addr);

        if (connect(sock, (struct sockaddr*)&client_addr, sizeof(client_addr)) == 0) {
            std::cout << "\t- connect " << ip << ":" << port << " success." << " sock ="
                << sock << " \n";
            connected_client ++;
            Client client;
            client.state = CONNETED;
            client.sock = sock;
            {
                // Lock the mutex to safely modify the shared client_map
                client_map.value[sock] = client;
                client_map.value[sock].partition_id = current_client_index;
                partiton_id_to_sock.access_with_function([&current_client_index, &sock](auto& map){
                    map[current_client_index] = sock;
                });
                sock_to_partition_id.access_with_function([&partition_id, &sock, &current_client_index](auto& map){
                    std::cout << "set sock " << sock << " partition id = " << current_client_index + 1 << std::endl;
                    map[sock] = current_client_index + 1;
                });
            }
        } else {
            // perror("\t- connect failed");
            close(sock);
        }
    }
    
    while(true){
        client_map.lock();
        if(client_map.value.size() < total_clients - 1) continue;
        break;
    }

    std::cout << "All clients connected. " << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Please provide the instance index (i).\n";
        return 1;
    }
    

    int partition_id = std::stoi(argv[1]);
    YAML::Node config = YAML::LoadFile("cluster.yaml");
    std::string program_name = std::string("pcfg-train-") + std::to_string(partition_id);
    uint32_t server_port = 9239 + partition_id;
    std::cout << "Creating server for self at port " << server_port << "\n";
    std::thread server_thread(serverLoop, server_port, partition_id);
    const YAML::Node& clients = config["cluster"]["clients"];
    int total_clients = clients.size();
    int connected_client = 1;
    int client_index = 0;
    
    

    std::cout << "Open share memory. " << std::endl;

    int size = sizeof(MemoryStorage);
    auto shared_memory = std::make_shared<SharedMemory>(program_name.c_str(), CREATE_NEW, size);
    auto storage = (MemoryStorage*)shared_memory->get_data();
   
    /* 1. Notify the application network topology are prepared and connected. */
    Message network_component_prepared_msg = gen_network_component_prepared_msg(partition_id);
    push_msg_to_shared_memory_rr(network_component_prepared_msg, shared_memory);

 
    /* 2. Wait ACK from the application */
    while(storage->network_communicator_messages[0].status == EMPTY_SLOT){}
    std::cout << "application repied: " <<  
        storage->network_communicator_messages[0].data + sizeof(int) << std::endl;
    storage->network_communicator_messages[0].status = EMPTY_SLOT;
    
    memcpy(&cnt_grammar, storage->network_communicator_messages[0].data, sizeof(int));
    std::cout << "n_grammar = " << cnt_grammar << std::endl;
    std::cout << "partition_id = " << partition_id << std::endl; 
    connect_to_other_partitions(total_clients, connected_client, client_index, clients, 
        partition_id, program_name);

    /* 3. Broadcast partition prepared message to all other partitions. */
    std::cout << "Broadcast prepared message." << std::endl;
    Message partition_prepared_msg = gen_partition_prepared_msg(partition_id);
    broadcast_message(0, partition_id, partition_prepared_msg);

    /* 4. Wait ACKs of PARTITION_PREPARED . */
    std::cout << total_clients << std::endl;
    {
        std::unique_lock<std::mutex> lock = partition_prepared_msg_ack_count.lock();
        partition_prepared_msg_cv.wait(lock, [&total_clients] { return partition_prepared_msg_ack_count.value == total_clients - 1; });
    }
    
    std::cout << "[barrier passed] All partition prepared!" << std::endl;
    std::cout << RED << "[!Important] Barrier 1: All partition arrive front pre-epoch-0." << RESET << std::endl;    
 

    int epoch = 0;
    const int MAX_EPOCHS = 10;
    const int package_per_epoch = total_clients * 4;
    
    while(epoch < MAX_EPOCHS){
        std::cout << std::endl;
        std::cout << "[Main Loop] partition " << program_name 
                  << " at the beginning of epoch " << epoch << std::endl;
        
        // 5.1 Notify Application begin a new epoch.
        Message epoch_begin_msg = gen_epoch_begin_message(epoch, partition_id, cnt_grammar);
        push_msg_to_shared_memory_rr(epoch_begin_msg, shared_memory);
        broadcast_message(package_per_epoch * epoch + total_clients * 1, partition_id, epoch_begin_msg);
        {
            std::unique_lock<std::mutex> lock = begin_epoch_msg_ack_count.lock();
            begin_epoch_msg_cv.wait(lock, [&total_clients, &epoch] { return begin_epoch_msg_ack_count.value[epoch] == total_clients - 1; });
        }
        std::cout << RED << "[Main Loop] [barrier passed] All partition prepare to proceed epoch "
                  << RESET << epoch << "!" << std::endl;
        
            

        /* Application Execution */
        // 5.2 Wait application finished.
        sleep_time = 0; //30 * 60;
        std::cout << "[Main Loop] wait application execution. " << std::endl;
        {
            std::unique_lock<std::mutex> lock(application_mutex);
            while(storage->network_communicator_messages[0].status == EMPTY_SLOT){
                sleep(sleep_time);
            }
            storage->network_communicator_messages[0].status = EMPTY_SLOT;
            int client_cnt_grammars = -1;
            memcpy(&client_cnt_grammars, storage->network_communicator_messages[0].data, sizeof(int));
            assert(client_cnt_grammars == cnt_grammar);
        }
        sleep_time = 0;
        double* this_partition_f = new double[cnt_grammar]();
        memcpy(this_partition_f, storage->network_communicator_messages[0].data + sizeof(int), sizeof(double) * cnt_grammar);

        /* Application Execution */
        std::cout << "[Main Loop] application finish and reply result " << std::endl;
        
        std::cout << GREEN << "Finish Epoch " << epoch << std::endl;
        for(int i = 0; i < cnt_grammar; i++){
            std::cout << "f[" << i << "] = " << this_partition_f[i] << std::endl;
        }
        std::cout << RESET << std::endl;

        storage->network_communicator_messages[0].status = EMPTY_SLOT;
       

        // 5.3 Broadcast epoch finished message.
        Message epoch_finished_msg = gen_epoch_finished_msg(partition_id, epoch, this_partition_f, cnt_grammar);
        std::cout << "[Main Loop] Prepare and broadcast epoch " << epoch << "finished message." << std::endl;
        broadcast_message(package_per_epoch * epoch + total_clients * 2, partition_id,
            epoch_finished_msg);

        // 5.4 Wait clients until all clients finish current epoch.
        {
            auto lock = epoch_completed_ack_count.lock();
            epoch_completed_msg_cv.wait(lock, [&total_clients, &epoch] { return epoch_completed_ack_count.value[epoch] == total_clients - 1; });
        }


        double* local_integrated_result = new double[cnt_grammar];
        std::cout << "[Main Loop] [barrier passed] All partition Completed Epoch " << epoch << "!" << std::endl;
        std::cout <<  BOLD_BLUE << "[Main Loop] Integrated Result = ";
        for(int gid = 0; gid < cnt_grammar; gid++){
            local_integrated_result[gid] =  log_sum_exp(global_result.get()[gid], this_partition_f[gid]);
            std::cout  << "\t - " << global_result.get()[gid] << " + " << this_partition_f[gid] 
                    << " = " << local_integrated_result[gid] << std::endl;
        }
        std::cout << RESET << std::endl;

        // 5.5 Notify the application the integrated results.
        integrated_result.access_with_function([local_integrated_result](auto& v)->void{
            memcpy(v.data(), local_integrated_result, sizeof(double) * cnt_grammar);
        });

        std::cout << RED << "[!Important] Inner Epoch" << epoch << 
            " Partition calculate integration result finished." << RESET << std::endl;


        Message msg_integrated_result_notification = gen_notificate_integrate_result_msg(local_integrated_result, epoch, cnt_grammar);
        push_msg_to_shared_memory_rr(msg_integrated_result_notification, shared_memory);

        // 5.6 Wait for ack from application
        {
            std::unique_lock<std::mutex> lock(application_mutex);
            while(storage->network_communicator_messages[0].status == EMPTY_SLOT){
            }
            storage->network_communicator_messages[0].status = EMPTY_SLOT;
        }
        std::cout << "[Main Loop] Integrated result processed by application." << std::endl;
        std::cout << RED << "[!Important] Inner Epoch " << epoch << 
            " Partition application processed integration results." << RESET << std::endl;


        Message msg_integrated_result_prepared = gen_integrated_result_prepared_msg(partition_id, epoch);
        broadcast_message(package_per_epoch * epoch + total_clients * 3, partition_id , msg_integrated_result_prepared);
        {
            auto lock = integrated_result_prepared_ack_count.lock();
            integrated_result_prepared_cv.wait(lock, [&total_clients, &epoch] { return integrated_result_prepared_ack_count.value[epoch] == total_clients - 1; });
        }
        std::cout << RED << "[!Important] Inner Epoch " << epoch << 
            " Barrier 2: All partition prepared integrated results in epoch " << epoch << "." 
            << RESET << std::endl;    

        broadcast_message(package_per_epoch * epoch + total_clients * 4, 
                        partition_id, msg_integrated_result_notification);
        {
            auto lock = integrated_result_confirmation_ack_count.lock();
            integrated_result_confirmation_cv.wait(lock, [&total_clients, &epoch] { 
                return integrated_result_confirmation_ack_count.value[epoch] == total_clients - 1; 
            });
        }

        std::cout << RED << "[!Important] Inner Epoch" << epoch << 
            " Barrier 3: All partition confirmed integrated results of epoch " << epoch << "." 
            << RESET << std::endl;

        epoch ++;

        {
            global_result.access_with_function([](auto& v)->void{
                for(int grammar_id = 0; grammar_id < cnt_grammar; grammar_id++){
                    GlobalState::global_result.value[grammar_id] = -INFINITY;
                    GlobalState::integrated_result.value[grammar_id] = -INFINITY;
                }
            });
        }

        std::cout << std::endl;
        std::cout <<  GREEN << "========================= END OF EPOCH " 
                << epoch 
                << "=========================" << RESET << std::endl;
    }

    std::cin.get();
    abort();

    return 0;
}
