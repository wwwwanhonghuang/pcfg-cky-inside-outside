#ifndef PACKAGE_HPP
#define PACKAGE_HPP
#include "distribution/common/message.h"

typedef struct {
    int sequence_number;
    int acknowlege_number;
    Message msg;
} Package;
#endif