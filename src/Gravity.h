#pragma once

#include "ForceTorque.h"

#include <string>

namespace rsim {

// Calculates central-body gravity from the attached vehicle's current state.
class Gravity final : public ForceTorque {
public:
    // IERS conventional Earth gravitational parameter in cubic metres per second squared.
    static constexpr double EarthGravitationalParameter = 3.986004418e14;

    // Construct inverse-square gravity about a fixed point in the integration frame.
    explicit Gravity(
        std::string name = "gravity",
        double gravitational_parameter = EarthGravitationalParameter,
        Vector3d center_in_integration_frame = Vector3d::Zero());

    // Recalculate body-coordinate gravity at the current vehicle center of gravity.
    void update() override;

    [[nodiscard]] ForceClassification classification() const noexcept override;

    // Return the central body's gravitational parameter in m^3/s^2.
    [[nodiscard]] double gravitationalParameter() const noexcept;

    // Return the central-body center expressed in the integration frame.
    [[nodiscard]] const Vector3d& centerInIntegrationFrame() const noexcept;

private:
    double gravitational_parameter_;
    Vector3d center_in_integration_frame_;
};

}  // namespace rsim
