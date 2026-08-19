#include "Environment.h"

#include "Parachute.h"
#include "Vehicle.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {

AtmosphereTable::AtmosphereTable(
    std::vector<double> altitude,
    std::vector<double> density,
    std::vector<double> speed_of_sound,
    std::vector<double> dynamic_viscosity,
    std::vector<Vector3d> wind)
    : altitude_(std::move(altitude)),
      density_(std::move(density)),
      speed_of_sound_(std::move(speed_of_sound)),
      dynamic_viscosity_(std::move(dynamic_viscosity)),
      wind_(std::move(wind)) {
    const std::size_t size = altitude_.size();
    if (size < 2U || density_.size() != size ||
        speed_of_sound_.size() != size ||
        dynamic_viscosity_.size() != size) {
        throw std::invalid_argument(
            "atmosphere arrays must have equal length of at least two");
    }
    for (std::size_t index = 0U; index < size; ++index) {
        const bool altitude_increases =
            index == 0U || altitude_[index] > altitude_[index - 1U];
        if (!std::isfinite(altitude_[index]) ||
            !std::isfinite(density_[index]) || density_[index] <= 0.0 ||
            !std::isfinite(speed_of_sound_[index]) ||
            speed_of_sound_[index] <= 0.0 ||
            !std::isfinite(dynamic_viscosity_[index]) ||
            dynamic_viscosity_[index] <= 0.0 || !altitude_increases) {
            throw std::invalid_argument("atmosphere table values are invalid");
        }
    }

    if (wind_.empty()) {
        wind_.assign(size, Vector3d::Zero());
    }
    if (wind_.size() != size) {
        throw std::invalid_argument(
            "atmosphere wind must have one vector per altitude");
    }
    for (const Vector3d& wind_value : wind_) {
        if (!wind_value.allFinite()) {
            throw std::invalid_argument("atmosphere wind values must be finite");
        }
    }
}

double AtmosphereTable::interpolate(
    const std::vector<double>& values, double altitude) const {
    if (altitude <= altitude_.front()) {
        return values.front();
    }
    if (altitude >= altitude_.back()) {
        return values.back();
    }
    const auto upper =
        std::upper_bound(altitude_.begin(), altitude_.end(), altitude);
    const std::size_t lower_index = static_cast<std::size_t>(
        upper - altitude_.begin() - 1);
    const double fraction =
        (altitude - altitude_[lower_index]) /
        (altitude_[lower_index + 1U] - altitude_[lower_index]);
    return values[lower_index] +
           fraction * (values[lower_index + 1U] - values[lower_index]);
}

double AtmosphereTable::density(double altitude) const {
    return interpolate(density_, altitude);
}

double AtmosphereTable::speedOfSound(double altitude) const {
    return interpolate(speed_of_sound_, altitude);
}

double AtmosphereTable::dynamicViscosity(double altitude) const {
    return interpolate(dynamic_viscosity_, altitude);
}

Vector3d AtmosphereTable::wind(double altitude) const {
    Vector3d result;
    for (Eigen::Index axis = 0; axis < 3; ++axis) {
        std::vector<double> values;
        values.reserve(wind_.size());
        for (const Vector3d& item : wind_) {
            values.push_back(item(axis));
        }
        result(axis) = interpolate(values, altitude);
    }
    return result;
}

bool AtmosphereTable::outsideAltitudeRange(double altitude) const noexcept {
    return altitude < altitude_.front() || altitude > altitude_.back();
}

FlightConditions::FlightConditions(
    std::string name,
    Vehicle& vehicle,
    AtmosphereTable atmosphere,
    Vector3d up_direction,
    double reference_length,
    Vector3d altitude_origin)
    : Model(std::move(name)),
      vehicle_(vehicle),
      atmosphere_(std::move(atmosphere)),
      up_(std::move(up_direction)),
      reference_length_(reference_length),
      altitude_origin_(std::move(altitude_origin)) {
    setPhase(ModelPhase::Environment);
    if (!up_.allFinite() || up_.norm() == 0.0 ||
        !std::isfinite(reference_length_) || reference_length_ <= 0.0) {
        throw std::invalid_argument("flight-condition geometry is invalid");
    }
    up_.normalize();

    addStateVariable(altitude_, "altitude", "m");
    addStateVariable(density_, "density", "kg/m^3");
    addStateVariable(mach_, "mach");
    addStateVariable(reynolds_, "reynolds");
    addStateVariable(dynamic_pressure_, "dynamicPressure", "Pa");
    addStateVariable(alpha_deg_, "alphaDeg", "deg");
    addStateVariable(phi_aero_deg_, "phiAeroDeg", "deg");
    addStateVariable(velocity_angle_from_vertical_deg_,
                     "velocityAngleFromVerticalDeg", "deg");
    addComputedStateVariable(
        [this]() -> StateValue { return relative_air_velocity_body_.x(); },
        "relativeAirVelocityBodyX", "m/s");
    addComputedStateVariable(
        [this]() -> StateValue { return relative_air_velocity_body_.y(); },
        "relativeAirVelocityBodyY", "m/s");
    addComputedStateVariable(
        [this]() -> StateValue { return relative_air_velocity_body_.z(); },
        "relativeAirVelocityBodyZ", "m/s");
    addStateVariable(atmosphere_clamped_, "atmosphereClamped");
}

void FlightConditions::addAerodynamics(Aerodynamics& aerodynamics) {
    aerodynamics_.push_back(&aerodynamics);
}

void FlightConditions::addParachute(Parachute& parachute) {
    parachutes_.push_back(&parachute);
}

void FlightConditions::update() {
    altitude_ =
        (vehicle_.kinematics().positionInParent() - altitude_origin_).dot(up_);
    atmosphere_clamped_ = atmosphere_.outsideAltitudeRange(altitude_);

    const Vector3d wind_in_parent = atmosphere_.wind(altitude_);
    const Vector3d relative_air_velocity_parent =
        vehicle_.kinematics().velocityInParent() - wind_in_parent;
    relative_air_velocity_body_ =
        vehicle_.kinematics().rotationIntoParent().conjugate() *
        relative_air_velocity_parent;
    const double speed = relative_air_velocity_body_.norm();
    density_ = atmosphere_.density(altitude_);
    const double sound_speed = atmosphere_.speedOfSound(altitude_);
    const double viscosity = atmosphere_.dynamicViscosity(altitude_);

    // The runtime database currently documents this shared lookup envelope.
    // Preserve the physical Mach number here; the database clamps only its
    // own interpolation coordinate when it evaluates aerodynamic data.
    mach_ = speed / sound_speed;
    reynolds_ = std::clamp(
        density_ * std::max(speed, 1.0) * reference_length_ / viscosity,
        1.0e5, 3.0e7);
    dynamic_pressure_ = 0.5 * density_ * speed * speed;

    const Vector3d velocity_in_parent =
        vehicle_.kinematics().velocityInParent();
    const double ground_speed = velocity_in_parent.norm();
    if (ground_speed > 1.0e-12) {
        const double cosine_from_vertical = std::clamp(
            velocity_in_parent.dot(up_) / ground_speed, -1.0, 1.0);
        velocity_angle_from_vertical_deg_ =
            std::acos(cosine_from_vertical) * 180.0 / std::acos(-1.0);
    } else {
        velocity_angle_from_vertical_deg_ = 0.0;
    }

    if (speed > 1.0e-12) {
        const double radians_to_degrees = 180.0 / std::acos(-1.0);
        alpha_deg_ = std::atan2(
            std::hypot(relative_air_velocity_body_.y(),
                       relative_air_velocity_body_.z()),
            relative_air_velocity_body_.x()) * radians_to_degrees;
        phi_aero_deg_ = std::atan2(
            relative_air_velocity_body_.z(),
            relative_air_velocity_body_.y()) * radians_to_degrees;
    } else {
        alpha_deg_ = 0.0;
        phi_aero_deg_ = 0.0;
    }

    const AerodynamicConditions conditions{
        mach_, alpha_deg_, phi_aero_deg_, reynolds_, dynamic_pressure_,
        relative_air_velocity_body_};
    for (Aerodynamics* aerodynamics : aerodynamics_) {
        aerodynamics->setConditions(conditions);
    }
    for (Parachute* parachute : parachutes_) {
        parachute->setAtmosphere(density_, wind_in_parent);
    }
}

void FlightConditions::initialize() {
    update();
    // Populate setup-time aerodynamic diagnostics before the t=0 truth row.
    for (Aerodynamics* aerodynamics : aerodynamics_) {
        aerodynamics->update();
    }
}

double FlightConditions::altitude() const noexcept {
    return altitude_;
}

double FlightConditions::mach() const noexcept { return mach_; }

}  // namespace rsim
