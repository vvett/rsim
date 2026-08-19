#include "FlightState.h"
#include "Environment.h"
#include "Parachute.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace rsim {

namespace {

void validateTrigger(const DeploymentTrigger& trigger,
                     const char* trigger_name) {
    if (std::isnan(trigger.maximum_altitude) ||
        std::isnan(trigger.maximum_vertical_velocity) ||
        !std::isfinite(trigger.earliest_simulation_time) ||
        trigger.earliest_simulation_time < 0.0) {
        throw std::invalid_argument(
            std::string(trigger_name) + " deployment trigger is invalid");
    }
    if (trigger.require_descending_altitude_crossing &&
        !std::isfinite(trigger.maximum_altitude)) {
        throw std::invalid_argument(
            std::string(trigger_name) +
            " altitude crossing requires a finite maximum altitude");
    }
}

}  // namespace

TruthState::TruthState(
    std::string name,
    FlightStateConfiguration configuration,
    std::vector<std::shared_ptr<ForceTorque>> thrust_sources,
    UpdateRate update_rate)
    : Model(std::move(name), update_rate),
      configuration_(std::move(configuration)) {
    if (!configuration_.origin_in_integration_frame.allFinite() ||
        !configuration_.up_direction_in_integration_frame.allFinite() ||
        !configuration_.rail_direction_in_integration_frame.allFinite()) {
        throw std::invalid_argument("flight-state vectors must be finite");
    }
    if (configuration_.up_direction_in_integration_frame.norm() == 0.0 ||
        configuration_.rail_direction_in_integration_frame.norm() == 0.0) {
        throw std::invalid_argument(
            "up and rail directions must be nonzero");
    }
    if (!std::isfinite(configuration_.rail_length) ||
        configuration_.rail_length < 0.0 ||
        !std::isfinite(configuration_.thrust_on_threshold) ||
        !std::isfinite(configuration_.thrust_off_threshold) ||
        configuration_.thrust_off_threshold < 0.0 ||
        configuration_.thrust_on_threshold <
            configuration_.thrust_off_threshold ||
        !std::isfinite(configuration_.descending_velocity_threshold) ||
        configuration_.descending_velocity_threshold > 0.0 ||
        !std::isfinite(configuration_.apogee_velocity_deadband) ||
        configuration_.apogee_velocity_deadband <= 0.0 ||
        !std::isfinite(configuration_.ground_altitude) ||
        !std::isfinite(
            configuration_.maximum_landing_vertical_velocity)) {
        throw std::invalid_argument("flight-state thresholds are invalid");
    }
    validateTrigger(configuration_.drogue_trigger, "drogue");
    validateTrigger(configuration_.main_trigger, "main");

    // Normalization makes all projections actual distances or velocities.
    configuration_.up_direction_in_integration_frame.normalize();
    configuration_.rail_direction_in_integration_frame.normalize();
    for (std::shared_ptr<ForceTorque>& source : thrust_sources) {
        addThrustSource(std::move(source));
    }
    addStateVariable(altitude_, "truthAltitude", "m");
    addStateVariable(vertical_velocity_, "verticalVelocity", "m/s");
    addStateVariable(rail_distance_, "railDistance", "m");
    addStateVariable(thrust_magnitude_, "thrustMagnitude", "N");
    addStateVariable(apogee_crossed_, "apogeeCrossed");
    addComputedStateVariable([this]() -> StateValue {
        return static_cast<std::int64_t>(
            vehicle_ ? vehicle_->state() : VehicleState::OnRail);
    }, "vehicleState");
}

void TruthState::addThrustSource(
    std::shared_ptr<ForceTorque> thrust_source) {
    if (!thrust_source) {
        throw std::invalid_argument("thrust source must not be null");
    }
    for (const std::shared_ptr<ForceTorque>& stored : thrust_sources_) {
        if (stored.get() == thrust_source.get()) {
            throw std::invalid_argument(
                "thrust source is already registered: " +
                thrust_source->name());
        }
    }
    thrust_sources_.push_back(std::move(thrust_source));
}

void TruthState::setRecoveryParachutes(
    const Parachute* drogue, const Parachute* main) noexcept {
    drogue_parachute_ = drogue;
    main_parachute_ = main;
}

const FlightStateConfiguration& TruthState::configuration() const noexcept {
    return configuration_;
}

double TruthState::altitude() const noexcept { return altitude_; }

double TruthState::verticalVelocity() const noexcept {
    return vertical_velocity_;
}

double TruthState::railDistance() const noexcept {
    return rail_distance_;
}

double TruthState::thrustMagnitude() const noexcept {
    return thrust_magnitude_;
}

bool TruthState::apogeeCrossed() const noexcept {
    return apogee_crossed_;
}

void TruthState::enableConsoleOutput(
    const FlightConditions& flight_conditions,
    const Aerodynamics& aerodynamics, double period) {
    if (!std::isfinite(period) || period <= 0.0) {
        throw std::invalid_argument("console output period must be positive and finite");
    }
    console_flight_conditions_ = &flight_conditions;
    console_aerodynamics_ = &aerodynamics;
    console_period_ = period;
    next_console_time_ = 0.0;
}

void TruthState::onVehicleAttached(Vehicle& vehicle) {
    vehicle_ = &vehicle;
    kinematics_ = &vehicle.kinematics();
    has_previous_sample_ = false;
    ascent_observed_ = false;
    apogee_crossed_ = false;
    drogue_altitude_crossed_ = false;
    main_altitude_crossed_ = false;
}

void TruthState::onVehicleDetached() noexcept {
    vehicle_ = nullptr;
    kinematics_ = nullptr;
}

bool TruthState::triggerReached(
    const DeploymentTrigger& trigger,
    bool altitude_crossing_latched) const noexcept {
    const bool altitude_reached =
        trigger.require_descending_altitude_crossing
            ? altitude_crossing_latched
            : altitude_ <= trigger.maximum_altitude;
    return altitude_reached &&
           vertical_velocity_ <= trigger.maximum_vertical_velocity &&
           simulationTime() >= trigger.earliest_simulation_time &&
           (!trigger.deploy_at_apogee || apogee_crossed_);
}

void TruthState::update() {
    if (!vehicle_ || !kinematics_) {
        throw std::logic_error(
            "flight-state model must be added to a vehicle before update");
    }

    const Vector3d displacement =
        kinematics_->positionInParent() -
        configuration_.origin_in_integration_frame;
    altitude_ = displacement.dot(
        configuration_.up_direction_in_integration_frame);
    vertical_velocity_ = kinematics_->velocityInParent().dot(
        configuration_.up_direction_in_integration_frame);
    rail_distance_ = displacement.dot(
        configuration_.rail_direction_in_integration_frame);

    Vector3d total_thrust = Vector3d::Zero();
    for (const std::shared_ptr<ForceTorque>& source : thrust_sources_) {
        if (source->enabled()) {
            total_thrust += source->force();
        }
    }
    thrust_magnitude_ = total_thrust.norm();

    if (console_flight_conditions_ &&
        simulationTime() + 1.0e-9 >= next_console_time_) {
        const Vector3d& velocity = kinematics_->velocityInParent();
        const double speed = velocity.norm();
        double angle_from_vertical = 0.0;
        if (speed > 1.0e-12) {
            const double cosine = std::clamp(
                velocity.dot(configuration_.up_direction_in_integration_frame) /
                    speed,
                -1.0, 1.0);
            angle_from_vertical = std::acos(cosine) * 180.0 / std::acos(-1.0);
        }
        std::cout << std::fixed
                  << "t=" << std::setw(7) << std::setprecision(2)
                  << simulationTime() << " s | altitude=" << std::setw(9)
                  << altitude_ << " m | mass=" << std::setw(8)
                  << vehicle_->rigidBody().mass() << " kg | thrust="
                  << std::setw(9) << thrust_magnitude_
                  << " N | velocity_inertial=[" << std::setw(8)
                  << velocity.x() << ", " << std::setw(8) << velocity.y()
                  << ", " << std::setw(8) << velocity.z()
                  << "] m/s | Mach=" << std::setw(6) << std::setprecision(3)
                  << console_flight_conditions_->mach()
                  << " | angle_from_vertical=" << std::setw(6)
                  << std::setprecision(2) << angle_from_vertical
                  << " deg | static_margin=" << std::setw(7)
                  << std::setprecision(3)
                  << console_aerodynamics_->staticStabilityMargin()
                  << " cal\n"
                  << std::flush;
        do {
            next_console_time_ += console_period_;
        } while (next_console_time_ <= simulationTime() + 1.0e-9);
    }

    const double deadband = configuration_.apogee_velocity_deadband;
    if (vertical_velocity_ > deadband) {
        ascent_observed_ = true;
    }
    if (!apogee_crossed_ && ascent_observed_ && has_previous_sample_ &&
        previous_vertical_velocity_ >= -deadband &&
        vertical_velocity_ < -deadband) {
        // The deadband ignores zero and small noisy reversals around apogee.
        apogee_crossed_ = true;
    }

    const auto updateAltitudeCrossing = [this](
        const DeploymentTrigger& trigger, bool& crossing_latched) {
        if (!crossing_latched &&
            trigger.require_descending_altitude_crossing &&
            has_previous_sample_ &&
            previous_altitude_ > trigger.maximum_altitude &&
            altitude_ <= trigger.maximum_altitude &&
            vertical_velocity_ < 0.0) {
            crossing_latched = true;
        }
    };
    updateAltitudeCrossing(
        configuration_.drogue_trigger, drogue_altitude_crossed_);
    updateAltitudeCrossing(
        configuration_.main_trigger, main_altitude_crossed_);

    const VehicleState state = vehicle_->state();
    VehicleState next_state = state;
    const bool hardware_driven_recovery = drogue_parachute_ || main_parachute_;
    const bool automatic_recovery =
        configuration_.legacy_automatic_recovery &&
        !hardware_driven_recovery;
    if ((state == VehicleState::Descending ||
         state == VehicleState::DrogueDeployed ||
         state == VehicleState::MainDeployed) &&
        altitude_ <= configuration_.ground_altitude &&
        vertical_velocity_ <=
            configuration_.maximum_landing_vertical_velocity) {
        next_state = VehicleState::Landed;
    } else {
        // Each case advances only forward, so noisy measurements cannot oscillate.
        switch (state) {
            case VehicleState::OnRail:
                if (rail_distance_ >= configuration_.rail_length) {
                    next_state =
                        thrust_magnitude_ >=
                                configuration_.thrust_on_threshold
                            ? VehicleState::Thrusting
                            : VehicleState::Coasting;
                }
                break;
            case VehicleState::Thrusting:
                if (thrust_magnitude_ <=
                    configuration_.thrust_off_threshold) {
                    next_state = VehicleState::Coasting;
                }
                break;
            case VehicleState::Coasting:
                if (vertical_velocity_ <=
                    configuration_.descending_velocity_threshold) {
                    const bool deploy_drogue_now =
                        automatic_recovery &&
                        configuration_.drogue_stage_enabled &&
                        configuration_.drogue_trigger.deploy_at_apogee &&
                        triggerReached(configuration_.drogue_trigger,
                                       drogue_altitude_crossed_);
                    // An apogee event can skip the transient Descending label so
                    // the force registry sees deployment in this same timestep.
                    next_state = deploy_drogue_now
                        ? VehicleState::DrogueDeployed
                        : VehicleState::Descending;
                }
                break;
            case VehicleState::Descending:
                if (main_parachute_ && main_parachute_->deployed()) {
                    next_state = VehicleState::MainDeployed;
                } else if (drogue_parachute_ && drogue_parachute_->deployed()) {
                    next_state = VehicleState::DrogueDeployed;
                } else if (automatic_recovery &&
                           configuration_.drogue_stage_enabled &&
                    triggerReached(configuration_.drogue_trigger,
                                   drogue_altitude_crossed_)) {
                    next_state = VehicleState::DrogueDeployed;
                } else if (!configuration_.drogue_stage_enabled &&
                           configuration_.main_stage_enabled &&
                           triggerReached(configuration_.main_trigger,
                                          main_altitude_crossed_)) {
                    // Disabling the drogue stage explicitly permits direct main deployment.
                    next_state = VehicleState::MainDeployed;
                }
                break;
            case VehicleState::DrogueDeployed:
                if (main_parachute_ && main_parachute_->deployed()) {
                    next_state = VehicleState::MainDeployed;
                } else if (automatic_recovery &&
                           configuration_.main_stage_enabled &&
                    triggerReached(configuration_.main_trigger,
                                   main_altitude_crossed_)) {
                    next_state = VehicleState::MainDeployed;
                }
                break;
            case VehicleState::MainDeployed:
            case VehicleState::Landed:
                break;
        }
    }

    if (next_state != state) {
        vehicle_->setState(next_state);
    }
    previous_altitude_ = altitude_;
    previous_vertical_velocity_ = vertical_velocity_;
    has_previous_sample_ = true;
}
void TruthState::initialize() { update(); }

}  // namespace rsim
