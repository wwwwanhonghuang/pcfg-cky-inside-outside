#ifndef MESSAGE_IMPL_H
#define MESSAGE_IMPL_H
#include <cstdint>
#include <cstring>
#include "common/message.h"

Message gen_integrated_result_prepared_msg(int partition_id, int epoch){
    Message msg;
    msg.status = MESSAGE_WAITING; 
    msg.msg_type = INTEGRATED_RESULT_PREPARED;
    memcpy(&msg.data[0], &epoch, sizeof(int));
    memcpy(&msg.data[4], &partition_id, sizeof(int));
    return msg;
};

Message gen_network_component_prepared_msg(int partition_id){
    Message msg;
    msg.status = MESSAGE_WAITING; 
    msg.msg_type = PARTITION_PREPARED;
    memcpy(msg.data, &partition_id, sizeof(int));
    return msg;
};

Message gen_partition_prepared_msg(int partition_id){
    Message prepared_msg;
    prepared_msg.status = MESSAGE_WAITING;
    prepared_msg.msg_type = PARTITION_PREPARED;
    memcpy(prepared_msg.data, &partition_id, sizeof(int));
    return prepared_msg;
};

Message gen_epoch_finished_msg(int partition_id, int epoch, int result){
    Message epoch_finished_msg;
    epoch_finished_msg.status = MESSAGE_WAITING;
    epoch_finished_msg.msg_type = EPOCH_COMPLETE;
    memcpy(epoch_finished_msg.data, &partition_id, sizeof(int));
    memcpy(&epoch_finished_msg.data[4], &epoch, sizeof(int));
    memcpy(&epoch_finished_msg.data[8], &result, sizeof(int));
    return epoch_finished_msg;
}

Message gen_notificate_integrate_result_msg(int integrated_result, int epoch){
    Message msg_integrated_result_notification;
    msg_integrated_result_notification.status = MESSAGE_WAITING;
    msg_integrated_result_notification.msg_type = NOTIFICATE_INTEGRATE_RESULT;
    memcpy(msg_integrated_result_notification.data, &integrated_result, sizeof(int));
    memcpy(&msg_integrated_result_notification.data[4], &epoch, sizeof(int));

    return msg_integrated_result_notification;
}

Message gen_epoch_begin_message(int epoch, int partition_id){
    Message epoch_begin_msg;
    epoch_begin_msg.msg_type = BEGIN_EPOCH;
    epoch_begin_msg.status = MESSAGE_WAITING;
    memcpy(&epoch_begin_msg.data[0], &epoch, sizeof(int));
    memcpy(&epoch_begin_msg.data[4], &partition_id, sizeof(int));
    return epoch_begin_msg;
}


#endif
