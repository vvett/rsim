#include "PressureController.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {

BangBangPressureController::BangBangPressureController(
    std::string name,
    std::shared_ptr<Tank> controlled_tank,
    std::shared_ptr<Valve> branch_valve,
    double pressure_setpoint,
    double deadband_half_width,
    bool initially_open,
    UpdateRate update_rate)
    : Model(std::move(name), update_rate),
      controlled_tank_(std::move(controlled_tank)),
      branch_valve_(std::move(branch_valve)),
      pressure_setpoint_(pressure_setpoint),
      deadband_half_width_(deadband_half_width),
      valve_open_(initially_open) {
    if (!controlled_tank_ || !branch_valve_) {
        throw std::invalid_argument(
            "bang-bang controller tank and valve must not be null");
    }
    if (!std::isfinite(pressure_setpoint_) ||
        !std::isfinite(deadband_half_width_) ||
        !(pressure_setpoint_ > deadband_half_width_) ||
        deadband_half_width_ < 0.0) {
        throw std::invalid_argument(
            "bang-bang pressure setpoint must exceed a finite nonnegative "
            "deadband half-width");
    }
    branch_valve_->setOpening(valve_open_ ? 1.0 : 0.0);
    addComputedStateVariable([this]() -> StateValue {
        return controlled_tank_->pressure(); }, "tankPressure", "Pa");
    addStateVariable(valve_open_, "valveOpen");
}

void BangBangPressureController::update() {
    const double pressure = controlled_tank_->pressure();
    if (pressure < lowerThreshold()) {
        valve_open_ = true;
    } else if (pressure > upperThreshold()) {
        valve_open_ = false;
    }
    branch_valve_->setOpening(valve_open_ ? 1.0 : 0.0);
}

const std::shared_ptr<Tank>&
BangBangPressureController::controlledTank() const noexcept {
    return controlled_tank_;
}

const std::shared_ptr<Valve>&
BangBangPressureController::branchValve() const noexcept {
    return branch_valve_;
}

double BangBangPressureController::pressureSetpoint() const noexcept {
    return pressure_setpoint_;
}

double BangBangPressureController::deadbandHalfWidth() const noexcept {
    return deadband_half_width_;
}

double BangBangPressureController::lowerThreshold() const noexcept {
    return pressure_setpoint_ - deadband_half_width_;
}

double BangBangPressureController::upperThreshold() const noexcept {
    return pressure_setpoint_ + deadband_half_width_;
}

bool BangBangPressureController::valveOpen() const noexcept {
    return valve_open_;
}

}  // namespace rsim
