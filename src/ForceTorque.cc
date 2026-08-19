#include "ForceTorque.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace rsim {

ForceTorque::ForceTorque(std::string name)
    : Model(std::move(name), UpdateRate::everyStep()) {
    addComputedStateVariable(
        [this]() -> StateValue { return force_.x(); }, "forceX", "N");
    addComputedStateVariable(
        [this]() -> StateValue { return force_.y(); }, "forceY", "N");
    addComputedStateVariable(
        [this]() -> StateValue { return force_.z(); }, "forceZ", "N");
    addComputedStateVariable(
        [this]() -> StateValue { return torque_.x(); }, "torqueX", "N*m");
    addComputedStateVariable(
        [this]() -> StateValue { return torque_.y(); }, "torqueY", "N*m");
    addComputedStateVariable(
        [this]() -> StateValue { return torque_.z(); }, "torqueZ", "N*m");
}

const Vector3d& ForceTorque::force() const noexcept { return force_; }

Vector3d ForceTorque::torqueAbout(
    const Vector3d& reference_point) const noexcept {
    return torque_ + (application_point_ - reference_point).cross(force_);
}

Vector6d ForceTorque::forceTorqueAbout(
    const Vector3d& reference_point) const noexcept {
    Vector6d result;
    result << force_, torqueAbout(reference_point);
    return result;
}

const Vector3d& ForceTorque::applicationPoint() const noexcept {
    return application_point_;
}
ForceClassification ForceTorque::classification() const noexcept {
    return ForceClassification::Applied;
}

void ForceTorque::setLoad(Vector3d force,
                          Vector3d torque,
                          Vector3d application_point) noexcept {
    force_ = std::move(force);
    torque_ = std::move(torque);
    application_point_ = std::move(application_point);
}

ConstantForceTorque::ConstantForceTorque(std::string name)
    : ForceTorque(std::move(name)) {}

void ConstantForceTorque::setLoad(Vector3d force,
                                  Vector3d torque,
                                  Vector3d application_point) noexcept {
    ForceTorque::setLoad(
        std::move(force), std::move(torque), std::move(application_point));
}

void ConstantForceTorque::update() {}

ForceTorqueRegistry::ForceTorqueRegistry(
    std::string name, UpdateRate update_rate)
    : Model(std::move(name), update_rate) { setPhase(ModelPhase::Forces); }

void ForceTorqueRegistry::add(ForceTorque& force_model) {
    if (contains(force_model)) {
        throw std::invalid_argument(
            "force model is already registered: " + force_model.name());
    }
    force_models_.push_back(&force_model);
}

void ForceTorqueRegistry::remove(ForceTorque& force_model) {
    const auto found = std::find(
        force_models_.begin(), force_models_.end(), &force_model);
    if (found == force_models_.end()) {
        throw std::invalid_argument(
            "force model is not registered: " + force_model.name());
    }
    force_models_.erase(found);
}

bool ForceTorqueRegistry::contains(
    const ForceTorque& force_model) const noexcept {
    return std::find(force_models_.begin(), force_models_.end(),
                     &force_model) != force_models_.end();
}

std::size_t ForceTorqueRegistry::size() const noexcept {
    return force_models_.size();
}

ForceTorque& ForceTorqueRegistry::at(std::size_t index) {
    return *force_models_.at(index);
}

const ForceTorque& ForceTorqueRegistry::at(std::size_t index) const {
    return *force_models_.at(index);
}

void ForceTorqueRegistry::update() {
    for (ForceTorque* force_model : force_models_) {
        if (!force_model->enabled()) {
            continue;
        }
        // Children receive the registry's native C++ timing context.
        force_model->setUpdateContext({
            simulationTime(), timeStep(), globalStep(), substep()});
        force_model->update();
    }
}

Vector3d ForceTorqueRegistry::getForce() const noexcept {
    Vector3d total = Vector3d::Zero();
    for (const ForceTorque* force_model : force_models_) {
        if (force_model->enabled()) {
            total += force_model->force();
        }
    }
    return total;
}

Vector3d ForceTorqueRegistry::getTorque(
    const Vector3d& reference_point) const noexcept {
    Vector3d total = Vector3d::Zero();
    for (const ForceTorque* force_model : force_models_) {
        if (force_model->enabled()) {
            total += force_model->torqueAbout(reference_point);
        }
    }
    return total;
}

Vector6d ForceTorqueRegistry::getForceTorque(
    const Vector3d& reference_point) const noexcept {
    Vector6d result;
    result << getForce(), getTorque(reference_point);
    return result;
}

}  // namespace rsim
