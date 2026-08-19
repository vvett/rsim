#include "FlightComputer.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <array>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace rsim {
namespace {

void* openLibrary(const std::string& path) {
#ifdef _WIN32
    return reinterpret_cast<void*>(LoadLibraryA(path.c_str()));
#else
    return dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
}

void closeLibrary(void* library) noexcept {
#ifdef _WIN32
    if (library) {
        FreeLibrary(reinterpret_cast<HMODULE>(library));
    }
#else
    if (library) {
        dlclose(library);
    }
#endif
}

std::string libraryError() {
#ifdef _WIN32
    return "Windows loader error " + std::to_string(GetLastError());
#else
    const char* error = dlerror();
    return error ? error : "unknown dynamic-loader error";
#endif
}

template<class Function>
Function loadSymbol(void* library, const char* name) {
#ifdef _WIN32
    void* address = reinterpret_cast<void*>(GetProcAddress(
        reinterpret_cast<HMODULE>(library), name));
#else
    void* address = dlsym(library, name);
#endif
    if (!address) {
        throw std::runtime_error(
            std::string("flight-computer plugin is missing symbol: ") + name);
    }
    return reinterpret_cast<Function>(address);
}

template<class Value>
std::vector<std::shared_ptr<Value>> resolveNamedObjects(
    const char* const* required_names,
    std::size_t required_count,
    const std::unordered_map<std::string, std::shared_ptr<Value>>& supplied,
    const char* category) {
    if (supplied.size() != required_count) {
        throw std::invalid_argument(
            std::string("flight-computer ") + category +
            " count does not match plugin schema");
    }
    std::vector<std::shared_ptr<Value>> resolved;
    resolved.reserve(required_count);
    for (std::size_t index = 0; index < required_count; ++index) {
        if (!required_names[index]) {
            throw std::runtime_error(
                std::string("plugin has a null ") + category + " name");
        }
        const auto found = supplied.find(required_names[index]);
        if (found == supplied.end() || !found->second) {
            throw std::invalid_argument(
                std::string("missing flight-computer ") + category +
                " binding: " + required_names[index]);
        }
        resolved.push_back(found->second);
    }
    return resolved;
}

}  // namespace

FlightComputer::FlightComputer(
    std::string name,
    std::string plugin_path,
    std::unordered_map<std::string, std::shared_ptr<ScalarSensor>> sensors,
    std::unordered_map<std::string, std::shared_ptr<IdealActuator>> actuators,
    std::unordered_map<std::string, double> configuration)
    : Model(std::move(name)), plugin_path_(std::move(plugin_path)) {
    setPhase(ModelPhase::FlightControl);
    library_ = openLibrary(plugin_path_);
    if (!library_) {
        throw std::runtime_error(
            std::string("could not load flight-computer plugin: ") +
            libraryError());
    }

    try {
        const auto get_descriptor = loadSymbol<RsimFlightComputerGetDescriptor>(
            library_, "rsim_flight_computer_get_descriptor");
        const RsimFlightComputerDescriptor* descriptor = get_descriptor();
        if (!descriptor ||
            descriptor->struct_size != sizeof(RsimFlightComputerDescriptor) ||
            descriptor->abi_version != RSIM_FLIGHT_COMPUTER_ABI_VERSION) {
            throw std::runtime_error(
                "flight-computer descriptor size or ABI version mismatch");
        }

        sensors_ = resolveNamedObjects(
            descriptor->sensor_names, descriptor->sensor_count,
            sensors, "sensor");
        actuators_ = resolveNamedObjects(
            descriptor->actuator_names, descriptor->actuator_count,
            actuators, "actuator");

        std::vector<RsimFlightComputerConfiguration> native_configuration;
        native_configuration.reserve(descriptor->configuration_count);
        if (configuration.size() != descriptor->configuration_count) {
            throw std::invalid_argument(
                "flight-computer configuration count does not match plugin schema");
        }
        for (std::size_t index = 0;
             index < descriptor->configuration_count; ++index) {
            const char* required_name = descriptor->configuration_names[index];
            const auto found = configuration.find(required_name);
            if (found == configuration.end()) {
                throw std::invalid_argument(
                    std::string("missing flight-computer configuration: ") +
                    required_name);
            }
            native_configuration.push_back({
                sizeof(RsimFlightComputerConfiguration),
                required_name,
                found->second});
        }

        inputs_.resize(descriptor->sensor_count);
        outputs_.resize(descriptor->actuator_count);
        telemetry_.resize(descriptor->telemetry_count);
        for (std::size_t index = 0;
             index < descriptor->telemetry_count; ++index) {
            const char* telemetry_name = descriptor->telemetry_names[index];
            if (!telemetry_name) {
                throw std::runtime_error("plugin has a null telemetry name");
            }
            addComputedStateVariable(
                [this, index]() -> StateValue { return telemetry_[index]; },
                telemetry_name);
        }

        const auto create = loadSymbol<RsimFlightComputerCreate>(
            library_, "rsim_flight_computer_create");
        step_ = loadSymbol<RsimFlightComputerStep>(
            library_, "rsim_flight_computer_step");
        destroy_ = loadSymbol<RsimFlightComputerDestroy>(
            library_, "rsim_flight_computer_destroy");
        std::array<char, 256> error{};
        plugin_state_ = create(
            native_configuration.data(), native_configuration.size(),
            error.data(), error.size());
        if (!plugin_state_) {
            throw std::runtime_error(
                std::string("flight-computer plugin creation failed: ") +
                error.data());
        }
    } catch (...) {
        closeLibrary(library_);
        library_ = nullptr;
        throw;
    }
}

FlightComputer::~FlightComputer() {
    if (destroy_ && plugin_state_) {
        destroy_(plugin_state_);
    }
    if (library_) {
        closeLibrary(library_);
    }
}

void FlightComputer::update() {
    for (std::size_t index = 0; index < sensors_.size(); ++index) {
        inputs_[index] = sensors_[index]->valid()
            ? sensors_[index]->value() : 0.0;
    }

    std::array<char, 256> error{};
    int status = 0;
    try {
        status = step_(
            plugin_state_, simulationTime(), timeStep(),
            inputs_.data(), inputs_.size(),
            outputs_.data(), outputs_.size(),
            telemetry_.data(), telemetry_.size(),
            error.data(), error.size());
    } catch (...) {
        throw std::runtime_error(
            "flight-computer plugin threw across its C ABI boundary");
    }
    if (status != 0) {
        throw std::runtime_error(
            "flight-computer plugin step failed with status " +
            std::to_string(status) + ": " + error.data());
    }
    for (std::size_t index = 0; index < actuators_.size(); ++index) {
        actuators_[index]->command(outputs_[index]);
    }
}

const std::string& FlightComputer::pluginPath() const noexcept {
    return plugin_path_;
}

std::vector<const void*> FlightComputer::actuatorIdentities() const {
    std::vector<const void*> identities;
    identities.reserve(actuators_.size());
    for (const auto& actuator : actuators_) {
        identities.push_back(actuator.get());
    }
    return identities;
}

}  // namespace rsim
