#include "Frames.h"
#include "Gravity.h"
#include "Simulator.h"
#include "Vehicle.h"

#include <Eigen/Geometry>

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace {

constexpr double earth_radius = 6378137.0;

bool near(double left, double right, double relative_tolerance = 1.0e-12) {
    return std::abs(left - right) <=
        relative_tolerance * std::max(std::abs(left), std::abs(right));
}

std::shared_ptr<rsim::Frame> inertialFrame() {
    return std::make_shared<rsim::Frame>(
        "J2000", rsim::InertialStatus::Inertial);
}

}  // namespace

int main() {
    const auto frame = inertialFrame();
    auto gravity = std::make_shared<rsim::Gravity>();
    rsim::Vehicle vehicle(
        "rocket", frame, 12.0, rsim::Matrix3d::Identity(),
        rsim::Vector3d::Zero(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(), {}, {gravity});

    vehicle.kinematics().setState(
        {earth_radius, 0.0, 0.0}, rsim::Vector3d::Zero(),
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    vehicle.forces().update();

    const double expected_surface_force =
        12.0 * rsim::Gravity::EarthGravitationalParameter /
        (earth_radius * earth_radius);
    assert(near(gravity->force().x(), -expected_surface_force));
    assert(near(gravity->force().y(), 0.0));
    assert(gravity->torqueAbout(vehicle.rigidBody().cg()).isZero(1.0e-12));

    // A 90-degree body rotation maps integration -X into body +Y.
    const Eigen::Quaterniond integration_from_body(
        Eigen::AngleAxisd(0.5 * std::acos(-1.0), rsim::Vector3d::UnitZ()));
    vehicle.kinematics().setState(
        {earth_radius, 0.0, 0.0}, rsim::Vector3d::Zero(),
        integration_from_body, rsim::Vector3d::Zero());
    vehicle.forces().update();
    assert(std::abs(gravity->force().x()) < 1.0e-10);
    assert(near(gravity->force().y(), expected_surface_force));

    // Doubling radius reduces inverse-square gravity to one quarter.
    vehicle.kinematics().setState(
        {2.0 * earth_radius, 0.0, 0.0}, rsim::Vector3d::Zero(),
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    vehicle.forces().update();
    assert(near(gravity->force().norm(), 0.25 * expected_surface_force));

    bool invalid_mu_rejected = false;
    try {
        rsim::Gravity invalid("invalid", 0.0);
    } catch (const std::invalid_argument&) {
        invalid_mu_rejected = true;
    }
    assert(invalid_mu_rejected);

    bool missing_vehicle_rejected = false;
    try {
        rsim::Gravity unattached;
        unattached.update();
    } catch (const std::logic_error&) {
        missing_vehicle_rejected = true;
    }
    assert(missing_vehicle_rejected);

    bool zero_radius_rejected = false;
    vehicle.kinematics().setState(
        rsim::Vector3d::Zero(), rsim::Vector3d::Zero(),
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    try {
        vehicle.forces().update();
    } catch (const std::runtime_error&) {
        zero_radius_rejected = true;
    }
    assert(zero_radius_rejected);

    const auto rotating_frame = std::make_shared<rsim::Frame>(
        "ECEF", rsim::InertialStatus::NonInertial);
    auto ambiguous_gravity = std::make_shared<rsim::Gravity>();
    rsim::Vehicle ambiguous_vehicle(
        "ambiguous", rotating_frame, 1.0, rsim::Matrix3d::Identity(),
        rsim::Vector3d::Zero(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(), {}, {ambiguous_gravity});
    ambiguous_vehicle.kinematics().setState(
        {earth_radius, 0.0, 0.0}, rsim::Vector3d::Zero(),
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    bool noninertial_rejected = false;
    try {
        ambiguous_vehicle.forces().update();
    } catch (const std::logic_error&) {
        noninertial_rejected = true;
    }
    assert(noninertial_rejected);

    // Simulator schedules the registry, which forwards native timing to gravity.
    vehicle.kinematics().setState(
        {earth_radius, 0.0, 0.0}, rsim::Vector3d::Zero(),
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    rsim::Simulator simulator(0.01);
    simulator.addVehicle(vehicle);
    simulator.step();
    assert(gravity->globalStep() == 0U);
    assert(near(gravity->timeStep(), 0.01));
    assert(vehicle.kinematics().positionInParent().x() < earth_radius);

    return 0;
}
