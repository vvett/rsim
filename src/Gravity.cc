#include "Gravity.h"

#include "Frames.h"
#include "Kinematics.h"
#include "MassProperties.h"
#include "Vehicle.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {

Gravity::Gravity(std::string name,
                 double gravitational_parameter,
                 Vector3d center_in_integration_frame)
    : ForceTorque(std::move(name)),
      gravitational_parameter_(gravitational_parameter),
      center_in_integration_frame_(
          std::move(center_in_integration_frame)) {
    if (!std::isfinite(gravitational_parameter_) ||
        gravitational_parameter_ <= 0.0) {
        throw std::invalid_argument(
            "gravity gravitational_parameter must be finite and positive");
    }
    if (!center_in_integration_frame_.allFinite()) {
        throw std::invalid_argument(
            "gravity center_in_integration_frame must be finite");
    }
}

void Gravity::update() {
    const Vehicle& vehicle = attachedVehicle();
    const RigidBody& rigid_body = vehicle.rigidBody();
    const Kinematics& kinematics = vehicle.kinematics();
    const std::shared_ptr<Frame>& integration_frame =
        kinematics.integrationFrame();

    if (!integration_frame) {
        throw std::logic_error("gravity requires an integration frame");
    }
    if (!integration_frame->isInertial()) {
        throw std::logic_error(
            "gravity requires an inertial integration frame; '" +
            integration_frame->name() + "' is non-inertial");
    }

    const Vector3d cg_in_body = rigid_body.cg();
    const Eigen::Quaterniond& integration_from_body =
        kinematics.rotationIntoParent();

    // Position includes the body-origin offset plus the rotated CG offset.
    const Vector3d cg_in_integration =
        kinematics.positionInParent() + integration_from_body * cg_in_body;
    const Vector3d radial_position =
        cg_in_integration - center_in_integration_frame_;
    const double radius_squared = radial_position.squaredNorm();

    if (!std::isfinite(radius_squared) || radius_squared <= 0.0) {
        throw std::runtime_error(
            "gravity is undefined at the central-body center");
    }

    // Newtonian central gravity is a = -mu r / |r|^3.
    const double inverse_radius_cubed =
        1.0 / (radius_squared * std::sqrt(radius_squared));
    const Vector3d force_in_integration =
        -rigid_body.mass() * gravitational_parameter_ *
        inverse_radius_cubed * radial_position;

    if (!force_in_integration.allFinite()) {
        throw std::runtime_error(
            "gravity produced a non-finite force; check radius and mass");
    }

    // ForceTorque loads use body axes, while the orbital position uses integration axes.
    const Vector3d force_in_body =
        integration_from_body.conjugate() * force_in_integration;
    setLoad(force_in_body, Vector3d::Zero(), cg_in_body);
}

ForceClassification Gravity::classification() const noexcept {
    return ForceClassification::Gravity;
}

double Gravity::gravitationalParameter() const noexcept {
    return gravitational_parameter_;
}

const Vector3d& Gravity::centerInIntegrationFrame() const noexcept {
    return center_in_integration_frame_;
}

}  // namespace rsim
