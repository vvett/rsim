#include "Kinematics.h"

#include "Integration.h"

#include <stdexcept>
#include <iostream>
#include <utility>

namespace rsim {

namespace {

// Keeps all values advanced by RK4 in one coupled state object.
struct KinematicState {
    Vector3d position_in_parent = Vector3d::Zero();
    Vector3d velocity_in_parent = Vector3d::Zero();
    Eigen::Vector4d quaternion_wxyz{1.0, 0.0, 0.0, 0.0};
    Vector3d angular_velocity_in_body = Vector3d::Zero();
};

// Add two states component by component for the RK4 weighted sums.
KinematicState operator+(const KinematicState& left,
                         const KinematicState& right) {
    return {
        left.position_in_parent + right.position_in_parent,
        left.velocity_in_parent + right.velocity_in_parent,
        left.quaternion_wxyz + right.quaternion_wxyz,
        left.angular_velocity_in_body + right.angular_velocity_in_body
    };
}

// Scale every state component by the same RK4 coefficient.
KinematicState operator*(double scale, const KinematicState& state) {
    return {
        scale * state.position_in_parent,
        scale * state.velocity_in_parent,
        scale * state.quaternion_wxyz,
        scale * state.angular_velocity_in_body
    };
}

// Convert stored w-x-y-z coefficients into Eigen's quaternion type.
Eigen::Quaterniond quaternionFromWxyz(const Eigen::Vector4d& coefficients) {
    Eigen::Quaterniond rotation{
        coefficients(0), coefficients(1), coefficients(2), coefficients(3)};
    if (rotation.norm() == 0.0) {
        throw std::runtime_error("kinematic orientation cannot be zero");
    }
    return rotation.normalized();
}

// Return q-dot for a body-to-parent quaternion and body-expressed angular rate.
Eigen::Vector4d quaternionDerivative(
    const Eigen::Quaterniond& parent_from_body,
    const Vector3d& angular_velocity_in_body) {
    // q_dot = 1/2 q times omega_quaternion for body-expressed angular velocity.
    const Eigen::Quaterniond angular_rate{
        0.0,
        angular_velocity_in_body.x(),
        angular_velocity_in_body.y(),
        angular_velocity_in_body.z()};
    const Eigen::Quaterniond derivative = parent_from_body * angular_rate;
    return 0.5 * Eigen::Vector4d{
        derivative.w(), derivative.x(), derivative.y(), derivative.z()};
}

}  // namespace

Kinematics::Kinematics(std::string body_frame_name,
                       std::shared_ptr<Frame> parent_frame,
                       ForceTorqueRegistry& force_torque_registry,
                       RigidBody& rigid_body,
                       UpdateRate update_rate)
    : Model(body_frame_name + "-kinematics", update_rate),
      bodyForceTorques(force_torque_registry),
      bodyMassProperties(rigid_body),
      integration_frame_(parent_frame) {
    setPhase(ModelPhase::Kinematics);
    if (!parent_frame) {
        throw std::invalid_argument("body frame must have a parent frame");
    }
    if (!parent_frame->isInertial()) {
        std::cerr << "rsim warning: integration frame '"
                  << parent_frame->name()
                  << "' is non-inertial; transport accelerations are not yet "
                     "included\n";
    }
    addComputedStateVariable([this]() -> StateValue {
        return bodyPositionInParent.x(); }, "positionX", "m");
    addComputedStateVariable([this]() -> StateValue {
        return bodyPositionInParent.y(); }, "positionY", "m");
    addComputedStateVariable([this]() -> StateValue {
        return bodyPositionInParent.z(); }, "positionZ", "m");
    addComputedStateVariable([this]() -> StateValue {
        return bodyOriginVelocityInParent.x(); }, "velocityX", "m/s");
    addComputedStateVariable([this]() -> StateValue {
        return bodyOriginVelocityInParent.y(); }, "velocityY", "m/s");
    addComputedStateVariable([this]() -> StateValue {
        return bodyOriginVelocityInParent.z(); }, "velocityZ", "m/s");
    addComputedStateVariable([this]() -> StateValue {
        return bodyVelocity.x(); }, "bodyVelocityX", "m/s");
    addComputedStateVariable([this]() -> StateValue {
        return bodyVelocity.y(); }, "bodyVelocityY", "m/s");
    addComputedStateVariable([this]() -> StateValue {
        return bodyVelocity.z(); }, "bodyVelocityZ", "m/s");
    addComputedStateVariable([this]() -> StateValue {
        return bodyOriginAccelerationInParent.x(); }, "accelerationX", "m/s^2");
    addComputedStateVariable([this]() -> StateValue {
        return bodyOriginAccelerationInParent.y(); }, "accelerationY", "m/s^2");
    addComputedStateVariable([this]() -> StateValue {
        return bodyOriginAccelerationInParent.z(); }, "accelerationZ", "m/s^2");

    auto definition = std::make_unique<BodyFrameDefinition>(*this);
    bodyFrameDefinition = definition.get();
    bodyFrameMotion = &definition->motionIntoParent();

    bodyFrame = std::make_shared<Frame>(
        std::move(body_frame_name),
        InertialStatus::NonInertial,
        std::move(parent_frame),
        std::move(definition));

    bodyFrame->update();
}

Kinematics::Kinematics(std::string body_frame_name,
                       FrameGraph& frame_graph,
                       const std::string& integration_frame_name,
                       ForceTorqueRegistry& force_torque_registry,
                       RigidBody& rigid_body,
                       UpdateRate update_rate)
    : Kinematics(std::move(body_frame_name),
                 frame_graph.frame(integration_frame_name),
                 force_torque_registry,
                 rigid_body,
                 update_rate) {}

Vector6d Kinematics::solveAccels(
    const Vector3d& angular_velocity_in_body) const {
    // The result contains body-expressed linear then angular acceleration.
    double mass = bodyMassProperties.mass();
    Vector3d cg_loc = bodyMassProperties.cg();
    Matrix3d moi_cg = bodyMassProperties.momentOfInertia();

    // The spatial inertia below is assembled about the body-frame origin, so
    // its input wrench must be expressed about that same point. Passing torque
    // about the CG here would account for the CG lever arm inconsistently.
    Vector6d ft_body = bodyForceTorques.getForceTorque(Vector3d::Zero());

    Matrix3d I3 = Matrix3d::Identity();

    Matrix3d cskew;
    cskew << 0, -cg_loc(2, 0), cg_loc(1, 0),
                cg_loc(2, 0), 0, -cg_loc(0, 0),
                -cg_loc(1, 0), cg_loc(0, 0), 0;

    Matrix3d wskew;
    wskew << 0, -angular_velocity_in_body(2), angular_velocity_in_body(1),
                angular_velocity_in_body(2), 0, -angular_velocity_in_body(0),
                -angular_velocity_in_body(1), angular_velocity_in_body(0), 0;

    Matrix6d inertia_tensor;
    Matrix3d IQ1 = mass * I3;
    Matrix3d IQ2 = -mass * cskew;
    Matrix3d IQ3 = mass * cskew;
    Matrix3d IQ4 = moi_cg - mass * cskew * cskew;
    inertia_tensor << IQ1, IQ2,
                      IQ3, IQ4;

    Vector6d coriolis_vector;
    Vector3d CQ1 = mass * wskew * wskew * cg_loc;
    Vector3d CQ2 = wskew * (moi_cg - mass * cskew * cskew) *
                   angular_velocity_in_body;
    coriolis_vector << CQ1,
                       CQ2;

    const auto inertia_decomposition = inertia_tensor.ldlt();
    if (inertia_decomposition.info() != Eigen::Success ||
        !inertia_decomposition.isPositive()) {
        throw std::runtime_error(
            "generalized inertia tensor must be positive definite");
    }

    const Vector6d accel_vector =
        inertia_decomposition.solve(ft_body - coriolis_vector);
    if (inertia_decomposition.info() != Eigen::Success ||
        !accel_vector.allFinite()) {
        throw std::runtime_error(
            "failed to solve generalized rigid-body acceleration");
    }

    return accel_vector;
}

void Kinematics::update() {
    KinematicState initial_state{
        bodyPositionInParent,
        bodyOriginVelocityInParent,
        Eigen::Vector4d{parentFromBodyRotation.w(),
                        parentFromBodyRotation.x(),
                        parentFromBodyRotation.y(),
                        parentFromBodyRotation.z()},
        bodyAngularVelocity
    };

    // Each stage solves forces and moments using that stage's coupled state.
    const auto derivative = [this](double, const KinematicState& state) {
        const Eigen::Quaterniond rotation =
            quaternionFromWxyz(state.quaternion_wxyz);
        const Vector6d accelerations =
            solveAccels(state.angular_velocity_in_body);

        KinematicState rate;
        rate.position_in_parent = state.velocity_in_parent;
        rate.velocity_in_parent =
            rotation * accelerations.head<3>();
        rate.quaternion_wxyz = quaternionDerivative(
            rotation, state.angular_velocity_in_body);
        rate.angular_velocity_in_body = accelerations.tail<3>();
        return rate;
    };

    const KinematicState integrated = rungeKutta4(
        initial_state, simulationTime(), timeStep(), derivative);

    bodyPositionInParent = integrated.position_in_parent;
    bodyOriginVelocityInParent = integrated.velocity_in_parent;
    parentFromBodyRotation =
        quaternionFromWxyz(integrated.quaternion_wxyz);
    bodyAngularVelocity = integrated.angular_velocity_in_body;

    // Re-evaluate acceleration at the final state for the published snapshot.
    const Vector6d final_accelerations = solveAccels(bodyAngularVelocity);
    bodyAcceleration = final_accelerations.head<3>();
    bodyAngularAcceleration = final_accelerations.tail<3>();
    bodyVelocity = parentFromBodyRotation.conjugate() *
                   bodyOriginVelocityInParent;
    bodyOriginAccelerationInParent =
        parentFromBodyRotation * bodyAcceleration;
    bodyAngularVelocityInParent =
        parentFromBodyRotation * bodyAngularVelocity;
    bodyAngularAccelerationInParent =
        parentFromBodyRotation * bodyAngularAcceleration;

    bodyFrame->update();
}

void Kinematics::setState(
    Vector3d position_in_parent,
    Vector3d velocity_in_parent,
    Eigen::Quaterniond parent_from_body_rotation,
    Vector3d angular_velocity_in_body) {
    if (parent_from_body_rotation.norm() == 0.0) {
        throw std::invalid_argument("kinematic orientation cannot be zero");
    }
    bodyPositionInParent = std::move(position_in_parent);
    bodyOriginVelocityInParent = std::move(velocity_in_parent);
    parentFromBodyRotation = parent_from_body_rotation.normalized();
    bodyAngularVelocity = std::move(angular_velocity_in_body);
    bodyVelocity = parentFromBodyRotation.conjugate() *
                   bodyOriginVelocityInParent;
    bodyFrame->update();
}

const Vector3d& Kinematics::positionInParent() const noexcept {
    return bodyPositionInParent;
}

const Vector3d& Kinematics::velocityInParent() const noexcept {
    return bodyOriginVelocityInParent;
}

const Eigen::Quaterniond& Kinematics::rotationIntoParent() const noexcept {
    return parentFromBodyRotation;
}

const Vector3d& Kinematics::angularVelocityInBody() const noexcept {
    return bodyAngularVelocity;
}

const std::shared_ptr<Frame>& Kinematics::integrationFrame() const noexcept {
    return integration_frame_;
}

BodyFrameDefinition::BodyFrameDefinition(const Kinematics& body_kinematics)
    : kinematics_(body_kinematics) {}

void BodyFrameDefinition::update() {
    motion_ = FrameMotion{
        Transform{kinematics_.parentFromBodyRotation,
                  kinematics_.bodyPositionInParent},
        kinematics_.bodyOriginVelocityInParent,
        kinematics_.bodyOriginAccelerationInParent,
        kinematics_.bodyAngularVelocityInParent,
        kinematics_.bodyAngularAccelerationInParent
    };
}

}  // namespace rsim
