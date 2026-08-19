#include "Sensors.h"
#include "ForceTorque.h"
#include "Vehicle.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {

ScalarSensor::ScalarSensor(std::string name, const Model& source,
                           std::string alias, UpdateRate update_rate)
    : Model(std::move(name), update_rate) {
    for (const StateVariable& variable : source.stateVariables()) {
        if (variable.alias == alias) { read_ = variable.read; break; }
    }
    if (!read_) throw std::invalid_argument("source does not publish state variable: " + alias);
    setPhase(ModelPhase::Sensors);
    addStateVariable(value_, "value");
    addStateVariable(valid_, "valid");
    addStateVariable(sample_time_, "sampleTime", "s");
    addStateVariable(sequence_, "sequence");
}

void ScalarSensor::update() {
    value_ = std::visit([](auto value) { return static_cast<double>(value); }, read_());
    valid_ = std::isfinite(value_);
    sample_time_ = simulationTime() + timeStep();
    ++sequence_;
}
void ScalarSensor::initialize() { update(); }
double ScalarSensor::value() const noexcept { return value_; }
bool ScalarSensor::valid() const noexcept { return valid_; }
double ScalarSensor::sampleTime() const noexcept { return sample_time_; }
std::uint64_t ScalarSensor::sequence() const noexcept { return static_cast<std::uint64_t>(sequence_); }

SpecificForceSensor::SpecificForceSensor(
    std::string name, Vehicle& vehicle, UpdateRate update_rate)
    : Model(std::move(name), update_rate), vehicle_(vehicle) {
    setPhase(ModelPhase::Sensors);
    addComputedStateVariable(
        [this]() -> StateValue { return specific_force_in_body_.x(); },
        "specificForceX", "m/s^2");
    addComputedStateVariable(
        [this]() -> StateValue { return specific_force_in_body_.y(); },
        "specificForceY", "m/s^2");
    addComputedStateVariable(
        [this]() -> StateValue { return specific_force_in_body_.z(); },
        "specificForceZ", "m/s^2");
    addStateVariable(valid_, "valid");
    addStateVariable(sample_time_, "sampleTime", "s");
    addStateVariable(sequence_, "sequence");
}

void SpecificForceSensor::update() {
    const double mass = vehicle_.rigidBody().mass();
    Vector3d applied_force = Vector3d::Zero();
    const ForceTorqueRegistry& registry = vehicle_.forces();
    for (std::size_t index = 0U; index < registry.size(); ++index) {
        const ForceTorque& force = registry.at(index);
        if (force.enabled() &&
            force.classification() == ForceClassification::Applied) {
            applied_force += force.force();
        }
    }
    valid_ = std::isfinite(mass) && mass > 0.0 && applied_force.allFinite();
    if (valid_) {
        specific_force_in_body_ = applied_force / mass;
    } else {
        specific_force_in_body_.setZero();
    }
    sample_time_ = simulationTime() + timeStep();
    ++sequence_;
}

void SpecificForceSensor::initialize() {
    update();
}

const Vector3d& SpecificForceSensor::specificForceInBody() const noexcept {
    return specific_force_in_body_;
}

bool SpecificForceSensor::valid() const noexcept {
    return valid_;
}

double SpecificForceSensor::sampleTime() const noexcept {
    return sample_time_;
}

std::uint64_t SpecificForceSensor::sequence() const noexcept {
    return static_cast<std::uint64_t>(sequence_);
}

}  // namespace rsim
