#include "Parachute.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {

namespace {

bool stateReached(VehicleState current, VehicleState target) noexcept {
    return static_cast<int>(current) >= static_cast<int>(target);
}

void validateDeploymentState(VehicleState state) {
    if (state != VehicleState::DrogueDeployed &&
        state != VehicleState::MainDeployed) {
        throw std::invalid_argument(
            "parachute deployment state must be DrogueDeployed or "
            "MainDeployed");
    }
}

}  // namespace

Parachute::Parachute(
    std::string name,
    double drag_coefficient,
    double reference_area,
    VehicleState deployment_state,
    double air_density,
    Vector3d wind_velocity_in_integration_frame,
    Vector3d application_point_in_body)
    : ForceTorque(std::move(name)),
      drag_coefficient_(drag_coefficient),
      reference_area_(reference_area),
      deployment_state_(deployment_state),
      air_density_(air_density),
      wind_velocity_in_integration_frame_(
          std::move(wind_velocity_in_integration_frame)),
      application_point_in_body_(std::move(application_point_in_body)) {
    if (!std::isfinite(drag_coefficient_) || drag_coefficient_ < 0.0 ||
        !std::isfinite(reference_area_) || reference_area_ <= 0.0 ||
        !std::isfinite(air_density_) || air_density_ < 0.0 ||
        !wind_velocity_in_integration_frame_.allFinite() ||
        !application_point_in_body_.allFinite()) {
        throw std::invalid_argument("parachute configuration is invalid");
    }
    validateDeploymentState(deployment_state_);
}

std::shared_ptr<Parachute> Parachute::fromDiameter(
    std::string name,
    double drag_coefficient,
    double diameter,
    VehicleState deployment_state,
    double air_density,
    Vector3d wind_velocity_in_integration_frame,
    Vector3d application_point_in_body) {
    if (!std::isfinite(diameter) || diameter <= 0.0) {
        throw std::invalid_argument("parachute diameter must be positive");
    }
    constexpr double pi = 3.14159265358979323846;
    return std::make_shared<Parachute>(
        std::move(name), drag_coefficient,
        0.25 * pi * diameter * diameter, deployment_state, air_density,
        std::move(wind_velocity_in_integration_frame),
        std::move(application_point_in_body));
}

std::shared_ptr<Parachute> Parachute::fromTerminalDescentSpeed(
    std::string name,
    double drag_coefficient,
    double terminal_descent_speed,
    double design_mass,
    double design_air_density,
    VehicleState deployment_state,
    Vector3d wind_velocity_in_integration_frame,
    Vector3d application_point_in_body) {
    if (!std::isfinite(drag_coefficient) || drag_coefficient <= 0.0 ||
        !std::isfinite(terminal_descent_speed) ||
        terminal_descent_speed <= 0.0 ||
        !std::isfinite(design_mass) || design_mass <= 0.0 ||
        !std::isfinite(design_air_density) || design_air_density <= 0.0) {
        throw std::invalid_argument(
            "terminal-descent sizing inputs must be finite and positive");
    }
    // Balance weight and quadratic drag at the setup-time design condition.
    const double area =
        2.0 * design_mass * StandardGravity /
        (design_air_density * drag_coefficient *
         terminal_descent_speed * terminal_descent_speed);
    return std::make_shared<Parachute>(
        std::move(name), drag_coefficient, area, deployment_state,
        design_air_density, std::move(wind_velocity_in_integration_frame),
        std::move(application_point_in_body));
}

void Parachute::setAtmosphere(
    double air_density,
    Vector3d wind_velocity_in_integration_frame) {
    if (!std::isfinite(air_density) || air_density < 0.0 ||
        !wind_velocity_in_integration_frame.allFinite()) {
        throw std::invalid_argument(
            "parachute atmosphere must contain finite nonnegative density");
    }
    air_density_ = air_density;
    wind_velocity_in_integration_frame_ =
        std::move(wind_velocity_in_integration_frame);
}

void Parachute::deploy() noexcept { deployed_ = true; }

double Parachute::dragCoefficient() const noexcept {
    return drag_coefficient_;
}

double Parachute::referenceArea() const noexcept { return reference_area_; }

double Parachute::equivalentDiameter() const noexcept {
    constexpr double pi = 3.14159265358979323846;
    return std::sqrt(4.0 * reference_area_ / pi);
}

VehicleState Parachute::deploymentState() const noexcept {
    return deployment_state_;
}

double Parachute::airDensity() const noexcept { return air_density_; }

const Vector3d& Parachute::windVelocityInIntegrationFrame() const noexcept {
    return wind_velocity_in_integration_frame_;
}

bool Parachute::deployed() const noexcept { return deployed_; }

void Parachute::onVehicleAttached(Vehicle& vehicle) {
    vehicle_ = &vehicle;
    kinematics_ = &vehicle.kinematics();
}

void Parachute::onVehicleDetached() noexcept {
    vehicle_ = nullptr;
    kinematics_ = nullptr;
}

void Parachute::update() {
    if (!vehicle_ || !kinematics_) {
        throw std::logic_error(
            "parachute must be added to a vehicle before update");
    }
    if (!deployed_ &&
        stateReached(vehicle_->state(), deployment_state_)) {
        deployed_ = true;
    }
    if (!deployed_) {
        setLoad(Vector3d::Zero(), Vector3d::Zero(),
                application_point_in_body_);
        return;
    }

    // Relative air velocity is vehicle velocity minus wind, in the integration frame.
    const Vector3d relative_velocity =
        kinematics_->velocityInParent() -
        wind_velocity_in_integration_frame_;
    const double speed = relative_velocity.norm();
    Vector3d drag_in_integration_frame = Vector3d::Zero();
    if (speed > 0.0) {
        drag_in_integration_frame =
            -0.5 * air_density_ * drag_coefficient_ * reference_area_ *
            speed * relative_velocity;
    }

    // ForceTorque loads are body-expressed, so rotate drag back from the parent.
    const Vector3d drag_in_body =
        kinematics_->rotationIntoParent().conjugate() *
        drag_in_integration_frame;
    setLoad(drag_in_body, Vector3d::Zero(), application_point_in_body_);
}

}  // namespace rsim
