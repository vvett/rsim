#pragma once

#include "Actuators.h"
#include "FlightComputerAbi.h"
#include "Sensors.h"

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

namespace rsim {

// Runs separately compiled flight software using only fixed native buffers.
class FlightComputer final : public Model {
public:
    FlightComputer(std::string name, std::string plugin_path,
                   std::unordered_map<std::string,
                                      std::shared_ptr<ScalarSensor>> sensors,
                   std::unordered_map<std::string,
                                      std::shared_ptr<IdealActuator>> actuators,
                   std::unordered_map<std::string, double> configuration = {});
    ~FlightComputer() override;
    void update() override;
    [[nodiscard]] const std::string& pluginPath() const noexcept;
    [[nodiscard]] std::vector<const void*> actuatorIdentities() const;
private:
    std::string plugin_path_; void* library_ = nullptr; void* plugin_state_ = nullptr;
    RsimFlightComputerStep step_ = nullptr; RsimFlightComputerDestroy destroy_ = nullptr;
    std::vector<std::shared_ptr<ScalarSensor>> sensors_;
    std::vector<std::shared_ptr<IdealActuator>> actuators_;
    std::vector<double> inputs_, outputs_, telemetry_;
};
}  // namespace rsim
