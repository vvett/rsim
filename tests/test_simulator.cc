#include "Integration.h"
#include "Simulator.h"
#include "Vehicle.h"

#include <Eigen/Dense>

#include <cassert>
#include <cmath>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Records each scheduler context so timing behavior can be checked directly.
class RecordingModel final : public rsim::Model {
public:
    RecordingModel(std::string name, rsim::UpdateRate rate)
        : Model(std::move(name), rate) {}

    void update() override {
        times.push_back(simulationTime());
        time_steps.push_back(timeStep());
        substeps.push_back(substep());
    }

    std::vector<double> times;
    std::vector<double> time_steps;
    std::vector<std::uint32_t> substeps;
};

// Records force-model update order while publishing a fixed test load.
class RecordingForce final : public rsim::ForceTorque {
public:
    RecordingForce(std::string name, int marker, std::vector<int>& order,
                   rsim::Vector3d force)
        : ForceTorque(std::move(name)), marker_(marker), order_(order),
          force_(std::move(force)) {}

    void update() override {
        order_.push_back(marker_);
        setLoad(force_, rsim::Vector3d::Zero());
    }

private:
    int marker_;
    std::vector<int>& order_;
    rsim::Vector3d force_;
};

// Changes one vehicle state during a step to verify next-step timestep selection.
class StateChangingModel final : public rsim::Model {
public:
    StateChangingModel(rsim::Vehicle& vehicle, rsim::VehicleState next_state)
        : Model("state-change"), vehicle_(vehicle), next_state_(next_state) {}

    void update() override { vehicle_.setState(next_state_); }

private:
    rsim::Vehicle& vehicle_;
    rsim::VehicleState next_state_;
};

bool near(double left, double right, double tolerance = 1e-12) {
    return std::abs(left - right) <= tolerance;
}

}  // namespace

int main() {
    RecordingModel fast("fast", rsim::UpdateRate::timesPerStep(2));
    RecordingModel slow("slow", rsim::UpdateRate::everyNSteps(3));
    rsim::Simulator simulator(0.1);
    simulator.addModel(fast);
    simulator.addModel(slow);
    simulator.runSteps(3);

    assert(fast.times.size() == 6U);
    assert(near(fast.times[0], 0.0));
    assert(near(fast.times[1], 0.05));
    assert(near(fast.times[5], 0.25));
    assert(near(fast.time_steps[0], 0.05));
    assert(fast.substeps[1] == 1U);

    assert(slow.times.size() == 1U);
    assert(near(slow.times[0], 0.0));
    assert(near(slow.time_steps[0], 0.3));
    assert(near(simulator.time(), 0.3));

    bool duplicate_rejected = false;
    try {
        simulator.addModel(fast);
    } catch (const std::invalid_argument&) {
        duplicate_rejected = true;
    }
    assert(duplicate_rejected);

    // RK4 should closely reproduce one short step of dy/dt = y.
    const double integrated = rsim::rungeKutta4(
        1.0, 0.0, 0.1,
        [](double, double state) { return state; });
    assert(near(integrated, std::exp(0.1), 1e-6));

    // addVehicle registers the owner and its three children exactly once.
    auto world = std::make_shared<rsim::Frame>(
        "world", rsim::InertialStatus::Inertial);
    rsim::Vehicle vehicle("rocket", world, 10.0,
                          Eigen::Matrix3d::Identity());
    rsim::Simulator vehicle_simulator(0.01);
    vehicle_simulator.addVehicle(vehicle);
    vehicle_simulator.step();
    assert(vehicle.bodyFrame()->parent() == world);
    assert(vehicle.kinematics().globalStep() == 0U);
    assert(near(vehicle_simulator.time(), 0.01));

    // The registry updates and accumulates models in insertion order.
    std::vector<int> force_update_order;
    RecordingForce first_force("first-force", 1, force_update_order,
                               Eigen::Vector3d{1.0, 0.0, 0.0});
    RecordingForce second_force("second-force", 2, force_update_order,
                                Eigen::Vector3d{0.0, 2.0, 0.0});
    rsim::ForceTorqueRegistry combined;
    combined.add(first_force);
    combined.add(second_force);
    rsim::Simulator force_simulator(0.1);
    force_simulator.addModel(combined);
    force_simulator.step();
    assert((force_update_order == std::vector<int>{1, 2}));
    assert(combined.getForce().isApprox(Eigen::Vector3d{1.0, 2.0, 0.0}));
    assert(first_force.simulationTime() == combined.simulationTime());
    second_force.setEnabled(false);
    assert(combined.getForce().isApprox(Eigen::Vector3d{1.0, 0.0, 0.0}));

    // Each model's application point is retained when torque is accumulated.
    rsim::ConstantForceTorque offset_force("offset-force");
    offset_force.setLoad(Eigen::Vector3d{1.0, 0.0, 0.0},
                         Eigen::Vector3d::Zero(),
                         Eigen::Vector3d{0.0, 1.0, 0.0});
    combined.add(offset_force);
    assert(combined.getTorque(Eigen::Vector3d::Zero()).isApprox(
        Eigen::Vector3d{0.0, 0.0, -1.0}));
    combined.remove(offset_force);
    assert(!combined.contains(offset_force));

    // A state transition during one step changes dt only on the following step.
    rsim::Vehicle state_vehicle(
        "state-timestep", world, 1.0, Eigen::Matrix3d::Identity());
    auto state_changer = std::make_shared<StateChangingModel>(
        state_vehicle, rsim::VehicleState::Thrusting);
    state_vehicle.addModel(state_changer);
    rsim::Simulator state_simulator(0.1);
    state_simulator.setVehicleStateTimeSteps(
        state_vehicle,
        {{rsim::VehicleState::OnRail, 0.04},
         {rsim::VehicleState::Thrusting, 0.01},
         {rsim::VehicleState::MainDeployed, 0.2}});
    state_simulator.addVehicle(state_vehicle);
    state_simulator.step();
    assert(near(state_simulator.time(), 0.04));
    assert(near(state_simulator.lastTimeStep(), 0.04));
    assert(state_vehicle.state() == rsim::VehicleState::Thrusting);
    state_simulator.step();
    assert(near(state_simulator.time(), 0.05));
    assert(near(state_simulator.lastTimeStep(), 0.01));

    // Unmapped states and cleared policies use the configured base timestep.
    state_vehicle.setState(rsim::VehicleState::Coasting);
    state_simulator.step();
    assert(near(state_simulator.time(), 0.15));
    state_vehicle.setState(rsim::VehicleState::MainDeployed);
    state_simulator.clearVehicleStateTimeSteps(state_vehicle);
    state_simulator.step();
    assert(near(state_simulator.time(), 0.25));
    state_simulator.setGlobalTimeStep(0.12);
    assert(near(state_simulator.globalTimeStep(), 0.12));
    state_simulator.step();
    assert(near(state_simulator.time(), 0.37));

    // Slow models receive the true sum of differently sized global steps.
    rsim::Vehicle variable_vehicle(
        "variable-slow", world, 1.0, Eigen::Matrix3d::Identity());
    RecordingModel variable_slow(
        "variable-slow-recorder", rsim::UpdateRate::everyNSteps(3));
    rsim::Simulator variable_simulator(0.2);
    variable_simulator.addVehicle(variable_vehicle);
    variable_simulator.addModel(variable_slow);
    variable_simulator.setVehicleStateTimeSteps(
        variable_vehicle,
        {{rsim::VehicleState::OnRail, 0.1},
         {rsim::VehicleState::Thrusting, 0.04},
         {rsim::VehicleState::Coasting, 0.02}});
    variable_simulator.step();
    variable_vehicle.setState(rsim::VehicleState::Thrusting);
    variable_simulator.step();
    variable_vehicle.setState(rsim::VehicleState::Coasting);
    variable_simulator.step();
    assert(near(variable_simulator.time(), 0.16));
    assert(variable_slow.times.size() == 1U);
    assert(near(variable_slow.times.front(), 0.0));
    assert(near(variable_slow.time_steps.front(), 0.16));

    // Fast models subdivide only the selected timestep for the current step.
    rsim::Vehicle fast_vehicle(
        "variable-fast", world, 1.0, Eigen::Matrix3d::Identity());
    RecordingModel variable_fast(
        "variable-fast-recorder", rsim::UpdateRate::timesPerStep(4));
    rsim::Simulator fast_simulator(0.2);
    fast_simulator.addVehicle(fast_vehicle);
    fast_simulator.addModel(variable_fast);
    fast_simulator.setVehicleStateTimeSteps(
        fast_vehicle, {{rsim::VehicleState::OnRail, 0.08}});
    fast_simulator.step();
    assert(variable_fast.times.size() == 4U);
    assert(near(variable_fast.time_steps.front(), 0.02));
    assert(near(variable_fast.times.back(), 0.06));

    // The smallest request wins when multiple vehicle policies are active.
    rsim::Vehicle first_vehicle(
        "first-policy", world, 1.0, Eigen::Matrix3d::Identity());
    rsim::Vehicle second_vehicle(
        "second-policy", world, 1.0, Eigen::Matrix3d::Identity());
    rsim::Simulator multi_vehicle_simulator(0.1);
    multi_vehicle_simulator.addVehicle(first_vehicle);
    multi_vehicle_simulator.addVehicle(second_vehicle);
    multi_vehicle_simulator.setVehicleStateTimeSteps(
        first_vehicle, {{rsim::VehicleState::OnRail, 0.06}});
    multi_vehicle_simulator.setVehicleStateTimeSteps(
        second_vehicle, {{rsim::VehicleState::OnRail, 0.025}});
    multi_vehicle_simulator.step();
    assert(near(multi_vehicle_simulator.time(), 0.025));
    assert(near(multi_vehicle_simulator.lastTimeStep(), 0.025));

    // Enabled and rate changes discard partial slow-model intervals.
    RecordingModel reset_model(
        "reset-model", rsim::UpdateRate::everyNSteps(3));
    rsim::Simulator reset_simulator(0.1);
    reset_simulator.addModel(reset_model);
    reset_simulator.step();
    reset_model.setEnabled(false);
    reset_model.setEnabled(true);
    reset_simulator.runSteps(2);
    assert(reset_model.times.empty());
    reset_simulator.step();
    assert(reset_model.times.size() == 1U);
    assert(near(reset_model.times.front(), 0.1));
    assert(near(reset_model.time_steps.front(), 0.3));
    reset_simulator.step();
    reset_model.setUpdateRate(rsim::UpdateRate::everyNSteps(2));
    reset_simulator.step();
    assert(reset_model.times.size() == 1U);
    reset_simulator.step();
    assert(reset_model.times.size() == 2U);
    assert(near(reset_model.time_steps.back(), 0.2));

    bool invalid_base_rejected = false;
    try {
        rsim::Simulator invalid(std::numeric_limits<double>::infinity());
    } catch (const std::invalid_argument&) {
        invalid_base_rejected = true;
    }
    assert(invalid_base_rejected);

    bool invalid_policy_rejected = false;
    try {
        state_simulator.setVehicleStateTimeSteps(
            state_vehicle, {{rsim::VehicleState::OnRail, 0.0}});
    } catch (const std::invalid_argument&) {
        invalid_policy_rejected = true;
    }
    assert(invalid_policy_rejected);

    return 0;
}
