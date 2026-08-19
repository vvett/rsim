#include "Vehicle.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace rsim {

Vehicle::Vehicle(std::string name,
                 std::shared_ptr<Frame> parent_frame,
                 double mass,
                 Matrix3d moment_of_inertia,
                 Vector3d center_of_gravity,
                 UpdateRate vehicle_rate,
                 UpdateRate rigid_body_rate,
                 UpdateRate force_rate,
                 UpdateRate kinematics_rate,
                 std::vector<std::shared_ptr<Model>> models,
                 std::vector<std::shared_ptr<ForceTorque>> force_torques)
    : Model(name, vehicle_rate),
      // make_unique constructs each owned model and returns its unique_ptr.
      rigid_body_(std::make_unique<RigidBody>(
          mass,
          std::move(moment_of_inertia),
          std::move(center_of_gravity),
          name + "-rigid-body",
          rigid_body_rate)),
      forces_(std::make_unique<ForceTorqueRegistry>(
          name + "-forces", force_rate)),
      // Dereferencing unique_ptr passes the owned objects to reference parameters.
      kinematics_(std::make_unique<Kinematics>(
          name + "-body",
          std::move(parent_frame),
          *forces_,
          *rigid_body_,
          kinematics_rate)) {
    try {
        for (std::shared_ptr<Model>& model : models) {
            addModel(std::move(model));
        }
        for (std::shared_ptr<ForceTorque>& force_torque : force_torques) {
            addSharedForceTorque(std::move(force_torque));
        }
    } catch (...) {
        // A failed constructor must not leave caller-owned models attached.
        for (const std::shared_ptr<Model>& model : models_) {
            model->detachFromVehicle();
        }
        for (std::size_t index = 0; index < forces_->size(); ++index) {
            forces_->at(index).detachFromVehicle();
        }
        throw;
    }
}

Vehicle::Vehicle(std::string name,
                 FrameGraph& frame_graph,
                 const std::string& integration_frame_name,
                 double mass,
                 Matrix3d moment_of_inertia,
                 Vector3d center_of_gravity,
                 UpdateRate vehicle_rate,
                 UpdateRate rigid_body_rate,
                 UpdateRate force_rate,
                 UpdateRate kinematics_rate,
                 std::vector<std::shared_ptr<Model>> models,
                 std::vector<std::shared_ptr<ForceTorque>> force_torques)
    : Vehicle(std::move(name),
              frame_graph.frame(integration_frame_name),
              mass,
              std::move(moment_of_inertia),
              std::move(center_of_gravity),
              vehicle_rate,
              rigid_body_rate,
              force_rate,
              kinematics_rate,
              std::move(models),
              std::move(force_torques)) {}

Vehicle::~Vehicle() {
    for (const std::shared_ptr<Model>& model : models_) {
        model->detachFromVehicle();
    }
    for (std::size_t index = 0; index < forces_->size(); ++index) {
        forces_->at(index).detachFromVehicle();
    }
}

void Vehicle::update() {
    // FlightStateModel is optional and runs through the generic model list.
}

VehicleState Vehicle::state() const noexcept { return state_; }

void Vehicle::setState(VehicleState state) noexcept { state_ = state; }

RigidBody& Vehicle::rigidBody() noexcept { return *rigid_body_; }

const RigidBody& Vehicle::rigidBody() const noexcept { return *rigid_body_; }

ForceTorqueRegistry& Vehicle::forces() noexcept { return *forces_; }

const ForceTorqueRegistry& Vehicle::forces() const noexcept { return *forces_; }

void Vehicle::addForceTorque(ForceTorque& force_torque) {
    if (forces_->contains(force_torque)) {
        throw std::invalid_argument(
            "force model is already registered: " + force_torque.name());
    }
    force_torque.attachToVehicle(*this);
    try {
        forces_->add(force_torque);
    } catch (...) {
        force_torque.detachFromVehicle();
        throw;
    }
}

void Vehicle::addSharedForceTorque(
    std::shared_ptr<ForceTorque> force_torque) {
    if (!force_torque) {
        throw std::invalid_argument("force model must not be null");
    }
    if (forces_->contains(*force_torque)) {
        throw std::invalid_argument(
            "force model is already registered: " + force_torque->name());
    }
    force_torque->attachToVehicle(*this);
    try {
        forces_->add(*force_torque);
    } catch (...) {
        force_torque->detachFromVehicle();
        throw;
    }
    try {
        shared_force_torques_.push_back(std::move(force_torque));
    } catch (...) {
        ForceTorque& added = forces_->at(forces_->size() - 1U);
        forces_->remove(added);
        added.detachFromVehicle();
        throw;
    }
}

void Vehicle::removeForceTorque(ForceTorque& force_torque) {
    forces_->remove(force_torque);
    force_torque.detachFromVehicle();
    shared_force_torques_.erase(
        std::remove_if(
            shared_force_torques_.begin(), shared_force_torques_.end(),
            [&force_torque](const std::shared_ptr<ForceTorque>& stored) {
                return stored.get() == &force_torque;
            }),
        shared_force_torques_.end());
}

void Vehicle::addModel(std::shared_ptr<Model> model) {
    if (!model) {
        throw std::invalid_argument("additional model must not be null");
    }
    const auto found = std::find_if(
        models_.begin(), models_.end(),
        [&model](const std::shared_ptr<Model>& stored) {
            return stored.get() == model.get();
        });
    if (found != models_.end()) {
        throw std::invalid_argument(
            "additional model is already registered: " + model->name());
    }
    model->attachToVehicle(*this);
    try {
        models_.push_back(std::move(model));
    } catch (...) {
        model->detachFromVehicle();
        throw;
    }
}

void Vehicle::removeModel(Model& model) {
    const auto found = std::find_if(
        models_.begin(), models_.end(),
        [&model](const std::shared_ptr<Model>& stored) {
            return stored.get() == &model;
        });
    if (found == models_.end()) {
        throw std::invalid_argument(
            "additional model is not registered: " + model.name());
    }
    (*found)->detachFromVehicle();
    models_.erase(found);
}

std::size_t Vehicle::modelCount() const noexcept { return models_.size(); }

Model& Vehicle::modelAt(std::size_t index) { return *models_.at(index); }

const Model& Vehicle::modelAt(std::size_t index) const {
    return *models_.at(index);
}

Kinematics& Vehicle::kinematics() noexcept { return *kinematics_; }

const Kinematics& Vehicle::kinematics() const noexcept { return *kinematics_; }

const std::shared_ptr<Frame>& Vehicle::bodyFrame() const noexcept {
    return kinematics_->bodyFrame;
}

}  // namespace rsim
