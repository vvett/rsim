#pragma once

#include "ForceTorque.h"
#include "Frames.h"
#include "Kinematics.h"
#include "MassProperties.h"
#include "Model.h"

#include <memory>
#include <string>
#include <vector>

namespace rsim {

// Identifies the vehicle's current high-level flight phase.
enum class VehicleState {
    OnRail,
    Thrusting,
    Coasting,
    Descending,
    DrogueDeployed,
    MainDeployed,
    Landed
};

// Owns the closely related models that together describe one vehicle.
class Vehicle : public Model {
public:
    // Construct the vehicle and all of its required child models together.
    Vehicle(std::string name,
            std::shared_ptr<Frame> parent_frame,
            double mass,
            Matrix3d moment_of_inertia,
            Vector3d center_of_gravity = Vector3d::Zero(),
            UpdateRate vehicle_rate = UpdateRate::everyStep(),
            UpdateRate rigid_body_rate = UpdateRate::everyStep(),
            UpdateRate force_rate = UpdateRate::everyStep(),
            UpdateRate kinematics_rate = UpdateRate::everyStep(),
            std::vector<std::shared_ptr<Model>> models = {},
            std::vector<std::shared_ptr<ForceTorque>> force_torques = {});

    // Disconnect borrowed force-model contexts before owned state is destroyed.
    ~Vehicle() override;

    // Construct the vehicle using an integration frame selected by graph name.
    Vehicle(std::string name,
            FrameGraph& frame_graph,
            const std::string& integration_frame_name,
            double mass,
            Matrix3d moment_of_inertia,
            Vector3d center_of_gravity = Vector3d::Zero(),
            UpdateRate vehicle_rate = UpdateRate::everyStep(),
            UpdateRate rigid_body_rate = UpdateRate::everyStep(),
            UpdateRate force_rate = UpdateRate::everyStep(),
            UpdateRate kinematics_rate = UpdateRate::everyStep(),
            std::vector<std::shared_ptr<Model>> models = {},
            std::vector<std::shared_ptr<ForceTorque>> force_torques = {});

    // Leave orchestration to separately scheduled child models.
    void update() override;

    // Return or replace the current flight phase.
    [[nodiscard]] VehicleState state() const noexcept;
    void setState(VehicleState state) noexcept;

    // Return the vehicle-owned mass-properties model.
    [[nodiscard]] RigidBody& rigidBody() noexcept;
    [[nodiscard]] const RigidBody& rigidBody() const noexcept;

    // Return the vehicle-owned ordered force-model registry.
    [[nodiscard]] ForceTorqueRegistry& forces() noexcept;
    [[nodiscard]] const ForceTorqueRegistry& forces() const noexcept;

    // Register a caller-owned C++ force model without taking ownership.
    void addForceTorque(ForceTorque& force_torque);

    // Retain and register a shared force model, as used by the Python frontend.
    void addSharedForceTorque(std::shared_ptr<ForceTorque> force_torque);

    // Unregister a force model and release any shared ownership held by Vehicle.
    void removeForceTorque(ForceTorque& force_torque);

    // Retain an additional model in deterministic insertion order.
    void addModel(std::shared_ptr<Model> model);
    void removeModel(Model& model);
    [[nodiscard]] std::size_t modelCount() const noexcept;
    [[nodiscard]] Model& modelAt(std::size_t index);
    [[nodiscard]] const Model& modelAt(std::size_t index) const;

    // Return the vehicle-owned motion model.
    [[nodiscard]] Kinematics& kinematics() noexcept;
    [[nodiscard]] const Kinematics& kinematics() const noexcept;

    // Return the body frame owned indirectly by Kinematics.
    [[nodiscard]] const std::shared_ptr<Frame>& bodyFrame() const noexcept;

private:
    VehicleState state_ = VehicleState::OnRail;

    // unique_ptr makes Vehicle the sole owner while keeping stable object addresses.
    std::unique_ptr<RigidBody> rigid_body_;
    std::unique_ptr<ForceTorqueRegistry> forces_;
    std::vector<std::shared_ptr<Model>> models_;
    std::vector<std::shared_ptr<ForceTorque>> shared_force_torques_;
    std::unique_ptr<Kinematics> kinematics_;
};

}  // namespace rsim
