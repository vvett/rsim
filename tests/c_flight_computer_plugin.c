#include "FlightComputerAbi.h"

#include <stdlib.h>

static const RsimFlightComputerDescriptor descriptor = {
    sizeof(RsimFlightComputerDescriptor),
    RSIM_FLIGHT_COMPUTER_ABI_VERSION,
    "c_smoke",
    0, NULL,
    0, NULL,
    0, NULL,
    0, NULL
};

const RsimFlightComputerDescriptor* rsim_flight_computer_get_descriptor(void) {
    return &descriptor;
}

void* rsim_flight_computer_create(
    const RsimFlightComputerConfiguration* configuration,
    size_t configuration_count,
    char* error_message,
    size_t error_message_size) {
    (void) configuration;
    (void) error_message;
    (void) error_message_size;
    return configuration_count == 0 ? malloc(1) : NULL;
}

int rsim_flight_computer_step(
    void* state, double time, double dt,
    const double* sensors, size_t sensor_count,
    double* actuators, size_t actuator_count,
    double* telemetry, size_t telemetry_count,
    char* error_message, size_t error_message_size) {
    (void) time;
    (void) dt;
    (void) sensors;
    (void) actuators;
    (void) telemetry;
    (void) error_message;
    (void) error_message_size;
    return state && sensor_count == 0 && actuator_count == 0 &&
           telemetry_count == 0 ? 0 : 1;
}

void rsim_flight_computer_destroy(void* state) {
    free(state);
}
