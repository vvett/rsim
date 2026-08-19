#include "FlightComputerAbi.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <new>

namespace {

const char* Sensors[] = {"altitude", "verticalVelocity"};
const char* Actuators[] = {"drogueDeploy", "mainDeploy"};
const char* Configuration[] = {
    "mainDeployAltitudeM", "apogeeVelocityDeadbandMps"};
const char* Telemetry[] = {
    "ascending", "drogueCommanded", "mainCommanded"};

const RsimFlightComputerDescriptor Descriptor{
    sizeof(RsimFlightComputerDescriptor),
    RSIM_FLIGHT_COMPUTER_ABI_VERSION,
    "vespula",
    2U, Sensors,
    2U, Actuators,
    2U, Configuration,
    3U, Telemetry};

struct VespulaState {
    double main_deploy_altitude = 609.6;
    double apogee_velocity_deadband = 0.1;
    bool ascending = false;
    bool drogue = false;
    bool main = false;
};

void setError(char* destination, std::size_t size, const char* message) {
    if (destination && size > 0U) {
        std::snprintf(destination, size, "%s", message);
    }
}

}  // namespace

extern "C" const RsimFlightComputerDescriptor*
rsim_flight_computer_get_descriptor() {
    return &Descriptor;
}

extern "C" void* rsim_flight_computer_create(
    const RsimFlightComputerConfiguration* configuration,
    std::size_t configuration_count,
    char* error_message,
    std::size_t error_message_size) {
    if (!configuration || configuration_count != 2U) {
        setError(error_message, error_message_size,
                 "Vespula requires exactly two configuration values");
        return nullptr;
    }
    for (std::size_t index = 0; index < configuration_count; ++index) {
        if (configuration[index].struct_size !=
            sizeof(RsimFlightComputerConfiguration)) {
            setError(error_message, error_message_size,
                     "configuration structure size mismatch");
            return nullptr;
        }
    }
    auto* state = new (std::nothrow) VespulaState{};
    if (!state) {
        setError(error_message, error_message_size, "allocation failed");
        return nullptr;
    }
    state->main_deploy_altitude = configuration[0].value;
    state->apogee_velocity_deadband = configuration[1].value;
    return state;
}

extern "C" int rsim_flight_computer_step(
    void* opaque,
    double,
    double,
    const double* inputs,
    std::size_t input_count,
    double* outputs,
    std::size_t output_count,
    double* telemetry,
    std::size_t telemetry_count,
    char* error_message,
    std::size_t error_message_size) {
    if (!opaque || !inputs || !outputs || !telemetry ||
        input_count != 2U || output_count != 2U || telemetry_count != 3U) {
        setError(error_message, error_message_size,
                 "Vespula runtime buffer size mismatch");
        return 1;
    }

    auto& state = *static_cast<VespulaState*>(opaque);
    const double altitude = inputs[0];
    const double vertical_velocity = inputs[1];
    if (vertical_velocity > state.apogee_velocity_deadband) {
        state.ascending = true;
    }
    if (state.ascending &&
        vertical_velocity < -state.apogee_velocity_deadband) {
        state.drogue = true;
    }
    if (state.drogue && altitude <= state.main_deploy_altitude &&
        vertical_velocity < 0.0) {
        state.main = true;
    }

    outputs[0] = state.drogue ? 1.0 : 0.0;
    outputs[1] = state.main ? 1.0 : 0.0;
    telemetry[0] = state.ascending ? 1.0 : 0.0;
    telemetry[1] = state.drogue ? 1.0 : 0.0;
    telemetry[2] = state.main ? 1.0 : 0.0;
    return 0;
}

extern "C" void rsim_flight_computer_destroy(void* opaque) {
    delete static_cast<VespulaState*>(opaque);
}
