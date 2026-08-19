#pragma once

#include "ForceTorque.h"
#include "Vehicle.h"

#include <string>

namespace rsim {

// Applies native aerodynamic drag after its configured flight phase is reached.
class Parachute final : public ForceTorque {
public:
    static constexpr double StandardGravity = 9.80665;

    Parachute(std::string name,
              double drag_coefficient,
              double reference_area,
              VehicleState deployment_state,
              double air_density = 1.225,
              Vector3d wind_velocity_in_integration_frame = Vector3d::Zero(),
              Vector3d application_point_in_body = Vector3d::Zero());

    // Construct a circular parachute using area = pi*diameter^2/4.
    [[nodiscard]] static std::shared_ptr<Parachute> fromDiameter(
        std::string name,
        double drag_coefficient,
        double diameter,
        VehicleState deployment_state,
        double air_density = 1.225,
        Vector3d wind_velocity_in_integration_frame = Vector3d::Zero(),
        Vector3d application_point_in_body = Vector3d::Zero());

    // Size a fixed area for a target steady descent speed at design conditions.
    [[nodiscard]] static std::shared_ptr<Parachute> fromTerminalDescentSpeed(
        std::string name,
        double drag_coefficient,
        double terminal_descent_speed,
        double design_mass,
        double design_air_density,
        VehicleState deployment_state,
        Vector3d wind_velocity_in_integration_frame = Vector3d::Zero(),
        Vector3d application_point_in_body = Vector3d::Zero());

    // Publish zero load before deployment and quadratic relative-wind drag after.
    void update() override;

    // Replace the fixed native atmosphere used by this first implementation.
    void setAtmosphere(
        double air_density,
        Vector3d wind_velocity_in_integration_frame = Vector3d::Zero());

    // Permanently deploy this parachute regardless of the current vehicle state.
    void deploy() noexcept;

    [[nodiscard]] double dragCoefficient() const noexcept;
    [[nodiscard]] double referenceArea() const noexcept;
    [[nodiscard]] double equivalentDiameter() const noexcept;
    [[nodiscard]] VehicleState deploymentState() const noexcept;
    [[nodiscard]] double airDensity() const noexcept;
    [[nodiscard]] const Vector3d& windVelocityInIntegrationFrame() const noexcept;
    [[nodiscard]] bool deployed() const noexcept;

protected:
    void onVehicleAttached(Vehicle& vehicle) override;
    void onVehicleDetached() noexcept override;

private:
    double drag_coefficient_;
    double reference_area_;
    VehicleState deployment_state_;
    double air_density_;
    Vector3d wind_velocity_in_integration_frame_;
    Vector3d application_point_in_body_;
    Vehicle* vehicle_ = nullptr;
    Kinematics* kinematics_ = nullptr;
    bool deployed_ = false;
};

}  // namespace rsim
