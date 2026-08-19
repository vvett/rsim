#pragma once

#include "ForceTorque.h"
#include "Vehicle.h"

#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace rsim {
class Parachute;
class FlightConditions;
class Aerodynamics;

// Defines when a descending vehicle may enter one parachute flight phase.
struct DeploymentTrigger {
    double maximum_altitude = std::numeric_limits<double>::infinity();
    double maximum_vertical_velocity = std::numeric_limits<double>::infinity();
    double earliest_simulation_time = 0.0;

    // Require the first downward velocity crossing after confirmed ascent.
    bool deploy_at_apogee = false;

    // Require and latch a descending crossing of maximum_altitude from above.
    bool require_descending_altitude_crossing = false;
};

// Contains the frame assumptions and thresholds for monotonic flight phases.
struct FlightStateConfiguration {
    // Position used as zero altitude and as the start of rail displacement.
    Vector3d origin_in_integration_frame = Vector3d::Zero();

    // Direction used to measure altitude and vertical velocity.
    Vector3d up_direction_in_integration_frame = Vector3d::UnitX();

    // Direction used to measure progress along the launch rail.
    Vector3d rail_direction_in_integration_frame = Vector3d::UnitX();
    double rail_length = 0.0;

    // Separate thresholds prevent thrusting/coasting chatter near cutoff.
    double thrust_on_threshold = 1.0;
    double thrust_off_threshold = 0.5;
    double descending_velocity_threshold = -0.1;

    // Velocities inside this symmetric band do not declare an apogee crossing.
    double apogee_velocity_deadband = 0.1;

    bool drogue_stage_enabled = true;
    DeploymentTrigger drogue_trigger{
        std::numeric_limits<double>::infinity(), 0.0, 0.0};
    bool main_stage_enabled = true;
    DeploymentTrigger main_trigger{
        300.0, std::numeric_limits<double>::infinity(), 0.0};

    double ground_altitude = 0.0;
    double maximum_landing_vertical_velocity = 0.5;

    // Preserve threshold-driven recovery only for explicitly legacy setups.
    bool legacy_automatic_recovery = false;
};

// Advances a vehicle through deterministic native C++ flight phases.
class TruthState final : public Model {
public:
    explicit TruthState(
        std::string name = "flight-state",
        FlightStateConfiguration configuration = {},
        std::vector<std::shared_ptr<ForceTorque>> thrust_sources = {},
        UpdateRate update_rate = UpdateRate::everyStep());

    // Evaluate at most one forward-only phase transition per update.
    void update() override;
    void initialize() override;

    // Add a native force model whose resultant is treated as engine thrust.
    void addThrustSource(std::shared_ptr<ForceTorque> thrust_source);

    // Derive recovery phases from actual parachute hardware instead of thresholds.
    void setRecoveryParachutes(const Parachute* drogue, const Parachute* main) noexcept;

    // Return the configuration fixed during setup.
    [[nodiscard]] const FlightStateConfiguration& configuration() const noexcept;

    // Return the most recently evaluated scalar flight conditions.
    [[nodiscard]] double altitude() const noexcept;
    [[nodiscard]] double verticalVelocity() const noexcept;
    [[nodiscard]] double railDistance() const noexcept;
    [[nodiscard]] double thrustMagnitude() const noexcept;

    // Return whether a noise-resistant downward apogee crossing has occurred.
    [[nodiscard]] bool apogeeCrossed() const noexcept;

    // Print native truth telemetry at the requested simulation-time interval.
    void enableConsoleOutput(const FlightConditions& flight_conditions,
                             const Aerodynamics& aerodynamics,
                             double period = 1.0);

protected:
    // Cache native vehicle state references when Vehicle owns this model.
    void onVehicleAttached(Vehicle& vehicle) override;
    void onVehicleDetached() noexcept override;

private:
    [[nodiscard]] bool triggerReached(
        const DeploymentTrigger& trigger,
        bool altitude_crossing_latched) const noexcept;

    FlightStateConfiguration configuration_;
    std::vector<std::shared_ptr<ForceTorque>> thrust_sources_;
    Vehicle* vehicle_ = nullptr;
    Kinematics* kinematics_ = nullptr;
    double altitude_ = 0.0;
    double vertical_velocity_ = 0.0;
    double rail_distance_ = 0.0;
    double thrust_magnitude_ = 0.0;
    double previous_altitude_ = 0.0;
    double previous_vertical_velocity_ = 0.0;
    bool has_previous_sample_ = false;
    bool ascent_observed_ = false;
    bool apogee_crossed_ = false;
    bool drogue_altitude_crossed_ = false;
    bool main_altitude_crossed_ = false;
    const Parachute* drogue_parachute_ = nullptr;
    const Parachute* main_parachute_ = nullptr;
    const FlightConditions* console_flight_conditions_ = nullptr;
    const Aerodynamics* console_aerodynamics_ = nullptr;
    double console_period_ = 0.0;
    double next_console_time_ = 0.0;
};

// Preserve source compatibility while directing new code to TruthState.
using FlightStateModel [[deprecated("use TruthState")]] = TruthState;

}  // namespace rsim
