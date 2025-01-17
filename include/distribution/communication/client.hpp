#ifndef CLIENT_HPP
#define CLIENT_HPP
#include <string>
typedef enum {
    UNCONNECTED,
    CONNETED
} ClientStates;

struct Client
{
    ClientStates state;
    int sock;
    std::string name;
    int partition_id;
};
#endif