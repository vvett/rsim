#ifndef RSIM_FLIGHT_COMPUTER_ABI_H
#define RSIM_FLIGHT_COMPUTER_ABI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RSIM_FLIGHT_COMPUTER_ABI_VERSION 2u

typedef struct RsimFlightComputerDescriptor {
    uint32_t struct_size;
    uint32_t abi_version;
    const char* plugin_name;
    size_t sensor_count;
    const char* const* sensor_names;
    size_t actuator_count;
    const char* const* actuator_names;
    size_t configuration_count;
    const char* const* configuration_names;
    size_t telemetry_count;
    const char* const* telemetry_names;
} RsimFlightComputerDescriptor;

typedef struct RsimFlightComputerConfiguration {
    uint32_t struct_size;
    const char* name;
    double value;
} RsimFlightComputerConfiguration;

typedef const RsimFlightComputerDescriptor*
    (*RsimFlightComputerGetDescriptor)(void);
typedef void* (*RsimFlightComputerCreate)(
    const RsimFlightComputerConfiguration* configuration,
    size_t configuration_count,
    char* error_message,
    size_t error_message_size);
typedef int (*RsimFlightComputerStep)(
    void* state,
    double simulation_time,
    double time_step,
    const double* sensor_values,
    size_t sensor_count,
    double* actuator_commands,
    size_t actuator_count,
    double* telemetry_values,
    size_t telemetry_count,
    char* error_message,
    size_t error_message_size);
typedef void (*RsimFlightComputerDestroy)(void* state);

#ifdef __cplusplus
}
#endif

#endif
