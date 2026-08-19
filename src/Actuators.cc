#include "Actuators.h"
#include "Parachute.h"
#include "Propulsion.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {
IdealActuator::IdealActuator(std::string name) : Model(std::move(name)) {
    setPhase(ModelPhase::Actuators);
}

ValveActuator::ValveActuator(std::string name, Valve& valve)
    : IdealActuator(std::move(name)), valve_(valve),
      commanded_(valve.opening()), actual_(valve.opening()) {
    addStateVariable(commanded_, "commanded"); addStateVariable(actual_, "actual");
}
void ValveActuator::command(double opening) {
    if (!std::isfinite(opening) || opening < 0.0 || opening > 1.0) {
        throw std::invalid_argument(
            "valve command must be between zero and one");
    }
    commanded_ = opening;
}
void ValveActuator::update() { valve_.setOpening(commanded_); actual_ = valve_.opening(); }
double ValveActuator::commanded() const noexcept { return commanded_; }
double ValveActuator::actual() const noexcept { return actual_; }
double ValveActuator::actualValue() const noexcept { return actual_; }

ParachuteActuator::ParachuteActuator(std::string name, Parachute& parachute)
    : IdealActuator(std::move(name)), parachute_(parachute) {
    addStateVariable(commanded_, "commanded"); addStateVariable(actual_, "actual");
}
void ParachuteActuator::commandDeploy() noexcept { commanded_ = true; }
void ParachuteActuator::command(double value) {
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw std::invalid_argument(
            "parachute command must be between zero and one");
    }
    if (value > 0.5) {
        commandDeploy();
    }
}
void ParachuteActuator::update() { if (commanded_) parachute_.deploy(); actual_ = parachute_.deployed(); }
bool ParachuteActuator::commanded() const noexcept { return commanded_; }
bool ParachuteActuator::actual() const noexcept { return actual_; }
double ParachuteActuator::actualValue() const noexcept { return actual_ ? 1.0 : 0.0; }

ModelEnableActuator::ModelEnableActuator(std::string name, Model& target)
    : IdealActuator(std::move(name)), target_(target),
      commanded_(target.enabled()), actual_(target.enabled()) {
    addStateVariable(commanded_, "commanded");
    addStateVariable(actual_, "actual");
}
void ModelEnableActuator::command(double value) {
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw std::invalid_argument(
            "model-enable command must be between zero and one");
    }
    commanded_ = value > 0.5;
}
void ModelEnableActuator::update() {
    target_.setEnabled(commanded_);
    actual_ = target_.enabled();
}
bool ModelEnableActuator::commanded() const noexcept { return commanded_; }
bool ModelEnableActuator::actual() const noexcept { return actual_; }
double ModelEnableActuator::actualValue() const noexcept {
    return actual_ ? 1.0 : 0.0;
}
}  // namespace rsim
