#pragma once

#include "ForceTorque.h"
#include "Model.h"
#include "defs.h"
#include "utility/RegularGridInterpolator.h"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace rsim {

// Holds the six interpolated body-axis aerodynamic coefficients.
struct AerodynamicCoefficients {
    double axial = 0.0;
    double side_y = 0.0;
    double side_z = 0.0;
    double roll_moment = 0.0;
    double pitch_moment = 0.0;
    double yaw_moment = 0.0;
};

// Supplies the independent variables and dimensional dynamic pressure.
struct AerodynamicConditions {
    double mach = 0.0;
    double alpha_deg = 0.0;
    double phi_aero_deg = 0.0;
    double reynolds = 0.0;
    double dynamic_pressure = 0.0;
    // Air-relative vehicle velocity expressed in body coordinates.
    Vector3d relative_air_velocity_in_body = Vector3d::Zero();
};

// Stores one force/torque result expressed in body coordinates.
struct AerodynamicLoad {
    Vector3d force = Vector3d::Zero();
    Vector3d torque = Vector3d::Zero();
};

// Reads and interpolates the restricted MAT format written by rsim.aero.
class AerodynamicDatabase {
public:
    explicit AerodynamicDatabase(const std::string& path);

    // Interpolate all six coefficients at one aerodynamic operating point.
    [[nodiscard]] AerodynamicCoefficients evaluate(
        double mach,
        double alpha_deg,
        double phi_aero_deg,
        double reynolds) const;

    // Interpolate CP in the database's nose-to-tail axial coordinate.
    [[nodiscard]] double evaluateCenterPressureX(
        double mach, double alpha_deg, double phi_aero_deg,
        double reynolds) const;

    [[nodiscard]] double referenceArea() const noexcept;
    [[nodiscard]] double referenceLength() const noexcept;
    [[nodiscard]] double momentReferenceX() const noexcept;
    [[nodiscard]] double minimumAlphaDeg() const noexcept;
    [[nodiscard]] double maximumAlphaDeg() const noexcept;
    [[nodiscard]] double stabilityReferenceLength() const noexcept;

private:
    // Interpolate one stored coefficient with shared four-axis conventions.
    [[nodiscard]] double interpolate(
        const std::vector<double>& values,
        double mach,
        double alpha_deg,
        double phi_aero_deg,
        double reynolds) const;

    using Interpolator = utility::RegularGridInterpolator<4>;
    std::optional<Interpolator> interpolator_;
    std::array<std::vector<double>, 6> coefficients_;
    std::vector<double> center_pressure_x_;
    double reference_area_ = 0.0;
    double reference_length_ = 0.0;
    double moment_reference_x_ = 0.0;
    double stability_reference_length_ = 0.0;
    double phi_period_deg_ = 360.0;
};

// Calculates one aerodynamic load as a native C++ ForceTorque model.
class Aerodynamics : public ForceTorque {
public:
    Aerodynamics(std::string name,
                 AerodynamicDatabase database);

    // Replace the flight conditions used by the next model update.
    void setConditions(AerodynamicConditions conditions);

    // Return the current conditions, coefficients, and dimensional load.
    [[nodiscard]] const AerodynamicConditions& conditions() const noexcept;
    [[nodiscard]] const AerodynamicCoefficients& coefficients() const noexcept;
    [[nodiscard]] const AerodynamicLoad& load() const noexcept;
    [[nodiscard]] double centerPressureX() const noexcept;
    [[nodiscard]] double staticStabilityMargin() const noexcept;

    // Evaluate a load without changing or publishing model state.
    [[nodiscard]] AerodynamicLoad calculateLoad(
        const AerodynamicConditions& conditions) const;

    // Interpolate, dimensionalize, and publish the aerodynamic load.
    void update() override;

private:
    AerodynamicDatabase database_;
    AerodynamicConditions conditions_;
    AerodynamicCoefficients coefficients_;
    AerodynamicLoad load_;
    double center_pressure_x_ = 0.0;
    double static_stability_margin_ = 0.0;
};

}  // namespace rsim
