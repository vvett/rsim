#pragma once

#include "Aerodynamics.h"
#include "Model.h"
#include "defs.h"

#include <string>
#include <vector>

namespace rsim {
class Vehicle;
class Parachute;

// Owns an in-memory altitude table copied from Python during setup.
class AtmosphereTable {
public:
    AtmosphereTable(std::vector<double> altitude, std::vector<double> density,
                    std::vector<double> speed_of_sound,
                    std::vector<double> dynamic_viscosity,
                    std::vector<Vector3d> wind = {});
    [[nodiscard]] double density(double altitude) const;
    [[nodiscard]] double speedOfSound(double altitude) const;
    [[nodiscard]] double dynamicViscosity(double altitude) const;
    [[nodiscard]] Vector3d wind(double altitude) const;
    [[nodiscard]] bool outsideAltitudeRange(double altitude) const noexcept;
private:
    double interpolate(const std::vector<double>& values, double altitude) const;
    std::vector<double> altitude_, density_, speed_of_sound_, dynamic_viscosity_;
    std::vector<Vector3d> wind_;
};

// Computes atmosphere and aerodynamic conditions entirely inside native updates.
class FlightConditions final : public Model {
public:
    FlightConditions(std::string name, Vehicle& vehicle, AtmosphereTable atmosphere,
                     Vector3d up_direction = Vector3d::UnitX(),
                     double reference_length = 1.0,
                     Vector3d altitude_origin = Vector3d::Zero());
    void addAerodynamics(Aerodynamics& aerodynamics);
    void addParachute(Parachute& parachute);
    void update() override;
    void initialize() override;
    [[nodiscard]] double altitude() const noexcept;
    [[nodiscard]] double mach() const noexcept;
private:
    Vehicle& vehicle_; AtmosphereTable atmosphere_; Vector3d up_; double reference_length_; Vector3d altitude_origin_;
    std::vector<Aerodynamics*> aerodynamics_; std::vector<Parachute*> parachutes_;
    double altitude_ = 0.0;
    double density_ = 1.225;
    double mach_ = 0.0;
    double reynolds_ = 0.0;
    double dynamic_pressure_ = 0.0;
    double alpha_deg_ = 0.0;
    double phi_aero_deg_ = 0.0;
    double velocity_angle_from_vertical_deg_ = 0.0;
    Vector3d relative_air_velocity_body_ = Vector3d::Zero();
    bool atmosphere_clamped_ = false;
};
}  // namespace rsim
