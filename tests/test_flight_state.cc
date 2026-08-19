#include "FlightState.h"
#include "Parachute.h"
#include "Simulator.h"

#include <Eigen/Dense>

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

bool near(double left, double right, double tolerance = 1e-9) {
    return std::abs(left - right) <= tolerance;
}

}  // namespace

int main() {
    auto world = std::make_shared<rsim::Frame>(
        "world", rsim::InertialStatus::Inertial);
    auto engine = std::make_shared<rsim::ConstantForceTorque>("engine");
    engine->setLoad(rsim::Vector3d{100.0, 0.0, 0.0},
                    rsim::Vector3d::Zero());

    rsim::FlightStateConfiguration configuration;
    configuration.legacy_automatic_recovery = true;
    configuration.rail_length = 10.0;
    configuration.thrust_on_threshold = 20.0;
    configuration.thrust_off_threshold = 10.0;
    configuration.drogue_trigger.maximum_altitude = 2000.0;
    configuration.drogue_trigger.maximum_vertical_velocity = 0.0;
    configuration.main_trigger.maximum_altitude = 300.0;

    auto flight_state = std::make_shared<rsim::FlightStateModel>(
        "flight-state", configuration,
        std::vector<std::shared_ptr<rsim::ForceTorque>>{engine});
    auto drogue = std::make_shared<rsim::Parachute>(
        "drogue", 1.0, 1.0, rsim::VehicleState::DrogueDeployed, 1.0);
    auto main = std::make_shared<rsim::Parachute>(
        "main", 2.0, 2.0, rsim::VehicleState::MainDeployed, 1.0);

    rsim::Vehicle vehicle(
        "rocket", world, 10.0, rsim::Matrix3d::Identity(),
        rsim::Vector3d::Zero(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(),
        std::vector<std::shared_ptr<rsim::Model>>{flight_state},
        std::vector<std::shared_ptr<rsim::ForceTorque>>{
            engine, drogue, main});
    rsim::Simulator simulator(0.01);
    simulator.addVehicle(vehicle);

    // Duplicate registration is rejected without detaching the original model.
    bool duplicate_rejected = false;
    try {
        vehicle.addForceTorque(*drogue);
    } catch (const std::invalid_argument&) {
        duplicate_rejected = true;
    }
    assert(duplicate_rejected);

    vehicle.kinematics().setState(
        rsim::Vector3d{5.0, 0.0, 0.0}, rsim::Vector3d::Zero(),
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::OnRail);
    assert(!drogue->deployed());
    assert(drogue->force().isZero());

    vehicle.kinematics().setState(
        rsim::Vector3d{11.0, 0.0, 0.0}, rsim::Vector3d{5.0, 0.0, 0.0},
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::Thrusting);
    assert(near(flight_state->railDistance(), 11.0));
    assert(near(flight_state->thrustMagnitude(), 100.0));

    engine->setLoad(rsim::Vector3d::Zero(), rsim::Vector3d::Zero());
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::Coasting);

    // A later thrust increase cannot move the monotonic state machine backward.
    engine->setLoad(rsim::Vector3d{100.0, 0.0, 0.0},
                    rsim::Vector3d::Zero());
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::Coasting);
    engine->setLoad(rsim::Vector3d::Zero(), rsim::Vector3d::Zero());

    vehicle.kinematics().setState(
        rsim::Vector3d{1000.0, 0.0, 0.0},
        rsim::Vector3d{-10.0, 0.0, 0.0},
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::Descending);
    assert(!drogue->deployed());

    // Flight-state evaluation precedes forces, so the drogue loads this step.
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::DrogueDeployed);
    assert(drogue->deployed());
    assert(!main->deployed());

    vehicle.kinematics().setState(
        rsim::Vector3d{250.0, 0.0, 0.0},
        rsim::Vector3d{-10.0, 0.0, 0.0},
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::MainDeployed);
    assert(main->deployed());
    // Drogue contributes 50 N and main contributes 200 N opposite descent.
    assert(vehicle.forces().getForce().isApprox(
        rsim::Vector3d{250.0, 0.0, 0.0}, 1e-9));

    vehicle.kinematics().setState(
        rsim::Vector3d::Zero(), rsim::Vector3d{-1.0, 0.0, 0.0},
        Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    simulator.step();
    assert(vehicle.state() == rsim::VehicleState::Landed);

    const auto circular = rsim::Parachute::fromDiameter(
        "round", 1.5, 2.0, rsim::VehicleState::MainDeployed);
    constexpr double pi = 3.14159265358979323846;
    assert(near(circular->referenceArea(), pi));
    assert(near(circular->equivalentDiameter(), 2.0));

    // Terminal-speed sizing uses a fixed setup mass and density.
    const double design_mass = 100.0;
    const double design_density = 0.9;
    const double design_speed = 10.0;
    const double design_cd = 1.5;
    const auto sized = rsim::Parachute::fromTerminalDescentSpeed(
        "sized", design_cd, design_speed, design_mass, design_density,
        rsim::VehicleState::MainDeployed);
    const double expected_area =
        2.0 * design_mass * rsim::Parachute::StandardGravity /
        (design_density * design_cd * design_speed * design_speed);
    assert(near(sized->referenceArea(), expected_area));
    assert(near(sized->equivalentDiameter(),
                std::sqrt(4.0 * expected_area / pi)));

    // A second vehicle checks hysteretic apogee and 2000-ft AGL crossings.
    rsim::FlightStateConfiguration recovery_configuration;
    recovery_configuration.legacy_automatic_recovery = true;
    recovery_configuration.rail_length = 0.0;
    recovery_configuration.apogee_velocity_deadband = 0.1;
    recovery_configuration.drogue_trigger.deploy_at_apogee = true;
    recovery_configuration.main_trigger.maximum_altitude = 609.6;
    recovery_configuration.main_trigger
        .require_descending_altitude_crossing = true;
    auto recovery_state = std::make_shared<rsim::FlightStateModel>(
        "recovery-state", recovery_configuration);
    auto recovery_drogue = std::make_shared<rsim::Parachute>(
        "recovery-drogue", 0.8, 0.5,
        rsim::VehicleState::DrogueDeployed, 1.0);
    auto recovery_main = std::make_shared<rsim::Parachute>(
        "recovery-main", 1.7, 4.0,
        rsim::VehicleState::MainDeployed, 1.0);
    rsim::Vehicle recovery_vehicle(
        "recovery-rocket", world, 100.0, rsim::Matrix3d::Identity(),
        rsim::Vector3d::Zero(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(), rsim::UpdateRate::everyStep(),
        rsim::UpdateRate::everyStep(),
        std::vector<std::shared_ptr<rsim::Model>>{recovery_state},
        std::vector<std::shared_ptr<rsim::ForceTorque>>{
            recovery_drogue, recovery_main});
    rsim::Simulator recovery_simulator(0.01);
    recovery_simulator.addVehicle(recovery_vehicle);

    const auto setRecoveryState = [&recovery_vehicle](double altitude,
                                                       double velocity) {
        recovery_vehicle.kinematics().setState(
            rsim::Vector3d{altitude, 0.0, 0.0},
            rsim::Vector3d{velocity, 0.0, 0.0},
            Eigen::Quaterniond::Identity(), rsim::Vector3d::Zero());
    };
    setRecoveryState(1000.0, 20.0);
    recovery_simulator.step();
    assert(recovery_vehicle.state() == rsim::VehicleState::Coasting);

    // Zero and small noisy samples do not count as a downward crossing.
    setRecoveryState(1010.0, 0.05);
    recovery_simulator.step();
    setRecoveryState(1010.0, -0.05);
    recovery_simulator.step();
    assert(!recovery_state->apogeeCrossed());
    assert(!recovery_drogue->deployed());

    // The first sample below the negative deadband latches apogee deployment.
    setRecoveryState(1009.0, -0.2);
    recovery_simulator.step();
    assert(recovery_state->apogeeCrossed());
    assert(recovery_vehicle.state() == rsim::VehicleState::DrogueDeployed);
    assert(recovery_drogue->deployed());
    assert(!recovery_main->deployed());

    setRecoveryState(700.0, -20.0);
    recovery_simulator.step();
    assert(recovery_vehicle.state() == rsim::VehicleState::DrogueDeployed);
    setRecoveryState(609.6, -20.0);
    recovery_simulator.step();
    assert(recovery_vehicle.state() == rsim::VehicleState::MainDeployed);
    assert(recovery_main->deployed());

    // Deployment and phase remain latched even if later samples rise again.
    setRecoveryState(800.0, 10.0);
    recovery_simulator.step();
    assert(recovery_vehicle.state() == rsim::VehicleState::MainDeployed);
    assert(recovery_drogue->deployed());
    assert(recovery_main->deployed());

    bool bad_thresholds_rejected = false;
    try {
        rsim::FlightStateConfiguration invalid;
        invalid.thrust_on_threshold = 1.0;
        invalid.thrust_off_threshold = 2.0;
        rsim::FlightStateModel rejected("invalid", invalid);
    } catch (const std::invalid_argument&) {
        bad_thresholds_rejected = true;
    }
    assert(bad_thresholds_rejected);

    bool bad_deployment_state_rejected = false;
    try {
        rsim::Parachute rejected(
            "invalid", 1.0, 1.0, rsim::VehicleState::Coasting);
    } catch (const std::invalid_argument&) {
        bad_deployment_state_rejected = true;
    }
    assert(bad_deployment_state_rejected);

    bool bad_terminal_sizing_rejected = false;
    try {
        static_cast<void>(rsim::Parachute::fromTerminalDescentSpeed(
            "invalid-sized", 1.0, 0.0, 10.0, 1.0,
            rsim::VehicleState::MainDeployed));
    } catch (const std::invalid_argument&) {
        bad_terminal_sizing_rejected = true;
    }
    assert(bad_terminal_sizing_rejected);

    return 0;
}
