#pragma once

#include "Model.h"
#include "Propulsion.h"

#include <memory>
#include <string>

namespace rsim {

// Commands one pressurant-branch valve from downstream tank pressure using
// symmetric hysteresis. Schedule this model before its Propulsion network.
class BangBangPressureController final : public Model {
public:
    BangBangPressureController(
        std::string name,
        std::shared_ptr<Tank> controlled_tank,
        std::shared_ptr<Valve> branch_valve,
        double pressure_setpoint,
        double deadband_half_width,
        bool initially_open = false,
        UpdateRate update_rate = UpdateRate::everyStep());

    void update() override;

    [[nodiscard]] const std::shared_ptr<Tank>& controlledTank() const noexcept;
    [[nodiscard]] const std::shared_ptr<Valve>& branchValve() const noexcept;
    [[nodiscard]] double pressureSetpoint() const noexcept;
    [[nodiscard]] double deadbandHalfWidth() const noexcept;
    [[nodiscard]] double lowerThreshold() const noexcept;
    [[nodiscard]] double upperThreshold() const noexcept;
    [[nodiscard]] bool valveOpen() const noexcept;

private:
    std::shared_ptr<Tank> controlled_tank_;
    std::shared_ptr<Valve> branch_valve_;
    double pressure_setpoint_;
    double deadband_half_width_;
    bool valve_open_;
};

}  // namespace rsim
