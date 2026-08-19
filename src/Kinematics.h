#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <memory>
#include <string>
#include <vector>

#include "ForceTorque.h"
#include "Frames.h"
#include "MassProperties.h"

namespace rsim {

class BodyFrameDefinition;

// Stores and updates the translational and rotational state of a rigid body.
class Kinematics : public Model {
public:
    // Construct the body frame and its definition beneath parent_frame.
    Kinematics(std::string body_frame_name,
               std::shared_ptr<Frame> parent_frame,
               ForceTorqueRegistry& force_torque_registry,
               RigidBody& rigid_body,
               UpdateRate update_rate = UpdateRate::everyStep());

    // Resolve the integration frame by name from a registered frame graph.
    Kinematics(std::string body_frame_name,
               FrameGraph& frame_graph,
               const std::string& integration_frame_name,
               ForceTorqueRegistry& force_torque_registry,
               RigidBody& rigid_body,
               UpdateRate update_rate = UpdateRate::everyStep());

    Kinematics(const Kinematics&) = delete;
    Kinematics& operator=(const Kinematics&) = delete;
    Kinematics(Kinematics&&) = delete;
    Kinematics& operator=(Kinematics&&) = delete;

    // Solve the current accelerations and publish the updated body-frame motion.
    void update() override;

    // Replace the complete state that RK4 advances during each update.
    void setState(Vector3d position_in_parent,
                  Vector3d velocity_in_parent,
                  Eigen::Quaterniond parent_from_body_rotation,
                  Vector3d angular_velocity_in_body);

    // Return the body origin's position expressed in the parent frame.
    [[nodiscard]] const Vector3d& positionInParent() const noexcept;

    // Return the body origin's velocity expressed in the parent frame.
    [[nodiscard]] const Vector3d& velocityInParent() const noexcept;

    // Return the rotation that converts body components into parent components.
    [[nodiscard]] const Eigen::Quaterniond& rotationIntoParent() const noexcept;

    // Return the body's angular velocity expressed in body coordinates.
    [[nodiscard]] const Vector3d& angularVelocityInBody() const noexcept;

    // Return the frame in which position and translational velocity are integrated.
    [[nodiscard]] const std::shared_ptr<Frame>& integrationFrame() const noexcept;

    // Owns the body frame and, through it, the body frame definition.
    std::shared_ptr<Frame> bodyFrame;

    // Non-owning alias to the definition owned by bodyFrame.
    BodyFrameDefinition* bodyFrameDefinition = nullptr;

    // Non-owning read-only alias to the motion owned by bodyFrameDefinition.
    const FrameMotion* bodyFrameMotion = nullptr;

private:
    friend class BodyFrameDefinition;

    // Return body-frame translational acceleration followed by angular acceleration.
    [[nodiscard]] Vector6d solveAccels(
        const Vector3d& angular_velocity_in_body) const;

    ForceTorqueRegistry& bodyForceTorques;
    RigidBody& bodyMassProperties;
    std::shared_ptr<Frame> integration_frame_;

    std::vector<double> time;

    // State used by the body-frame equations of motion.
    Vector3d bodyAngularVelocity = Vector3d::Zero();
    Vector3d bodyAngularAcceleration = Vector3d::Zero();
    Vector3d bodyVelocity = Vector3d::Zero();
    Vector3d bodyAcceleration = Vector3d::Zero();

    // Body pose and rates expressed in the parent frame for FrameMotion.
    Vector3d bodyPositionInParent = Vector3d::Zero();
    Eigen::Quaterniond parentFromBodyRotation = Eigen::Quaterniond::Identity();
    Vector3d bodyOriginVelocityInParent = Vector3d::Zero();
    Vector3d bodyOriginAccelerationInParent = Vector3d::Zero();
    Vector3d bodyAngularVelocityInParent = Vector3d::Zero();
    Vector3d bodyAngularAccelerationInParent = Vector3d::Zero();
};

// Calculates the body frame's current motion relative to its parent.
class BodyFrameDefinition final : public FrameDefinition {
public:
    // Retain a live reference to the Kinematics object that supplies the state.
    explicit BodyFrameDefinition(const Kinematics& body_kinematics);

    // Copy the latest parent-expressed kinematic state into motion_.
    void update() override;

private:
    const Kinematics& kinematics_;
};

}  // namespace rsim
