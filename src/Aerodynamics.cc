#include "Aerodynamics.h"
#include "utility/MatFileReader.h"
#include "Vehicle.h"

#include <algorithm>
#include <cmath>
#include <regex>
#include <stdexcept>
#include <utility>

namespace rsim {

namespace {

double jsonNumber(const std::string& json, const std::string& key) {
    const std::regex pattern(
        "\\\"" + key + "\\\"\\s*:\\s*([-+0-9.eE]+)");
    std::smatch match;
    if (!std::regex_search(json, match, pattern)) {
        throw std::runtime_error("aerodynamic metadata is missing " + key);
    }
    return std::stod(match[1].str());
}

}  // namespace

AerodynamicDatabase::AerodynamicDatabase(const std::string& path) {
    utility::MatFileReader file(path);
    if (file.readUint32("format_version") != 2U) {
        throw std::runtime_error("unsupported aerodynamic database format version");
    }

    Interpolator::Axes axes{
        std::move(file.readDoubleArray("axis__mach").values),
        std::move(file.readDoubleArray("axis__alpha_deg").values),
        std::move(file.readDoubleArray("axis__phi_aero_deg").values),
        std::move(file.readDoubleArray("axis__reynolds").values)};
    interpolator_.emplace(std::move(axes));

    const std::array<std::string, 6> names{
        "CA", "CY", "CZ", "Cmx", "Cmy", "Cmz"};
    const std::vector<std::size_t> expected_dimensions{
        interpolator_->axis(0).size(), interpolator_->axis(1).size(),
        interpolator_->axis(2).size(), interpolator_->axis(3).size()};
    for (std::size_t index = 0U; index < names.size(); ++index) {
        const std::string variable_name = "coefficient__" + names[index];
        utility::MatDoubleArray coefficient = file.readDoubleArray(variable_name);
        if (coefficient.dimensions != expected_dimensions) {
            throw std::runtime_error("coefficient dimensions do not match axes: " +
                                     names[index]);
        }
        if (!std::all_of(coefficient.values.begin(), coefficient.values.end(),
                         [](double value) { return std::isfinite(value); })) {
            throw std::runtime_error("coefficient contains invalid values: " +
                                     names[index]);
        }
        coefficients_[index] = std::move(coefficient.values);
    }
    utility::MatDoubleArray center_pressure =
        file.readDoubleArray("coefficient__xCP");
    if (center_pressure.dimensions != expected_dimensions ||
        !std::all_of(center_pressure.values.begin(), center_pressure.values.end(),
                     [](double value) { return std::isfinite(value); })) {
        throw std::runtime_error("xCP dimensions or values are invalid");
    }
    center_pressure_x_ = std::move(center_pressure.values);

    const std::string metadata = file.readUtf8("metadata_json_utf8");
    reference_area_ = jsonNumber(metadata, "coefficient_reference_area");
    reference_length_ = jsonNumber(metadata, "coefficient_reference_length");
    moment_reference_x_ = jsonNumber(metadata, "moment_reference_x");
    stability_reference_length_ =
        jsonNumber(metadata, "stability_reference_length");
    phi_period_deg_ = jsonNumber(metadata, "phi_aero_period_deg");
    if (reference_area_ <= 0.0 || reference_length_ <= 0.0 ||
        phi_period_deg_ <= 0.0 || stability_reference_length_ <= 0.0) {
        throw std::runtime_error("aerodynamic reference quantities must be positive");
    }
}

double AerodynamicDatabase::evaluateCenterPressureX(
    double mach, double alpha_deg, double phi_aero_deg,
    double reynolds) const {
    return interpolate(center_pressure_x_, mach, alpha_deg, phi_aero_deg,
                       reynolds);
}

double AerodynamicDatabase::interpolate(
    const std::vector<double>& values,
    double mach,
    double alpha_deg,
    double phi_aero_deg,
    double reynolds) const {
    const std::vector<double>& phi_axis = interpolator_->axis(2);
    mach = std::clamp(
        mach, interpolator_->axis(0).front(), interpolator_->axis(0).back());
    alpha_deg = std::clamp(
        alpha_deg, interpolator_->axis(1).front(),
        interpolator_->axis(1).back());
    reynolds = std::clamp(
        reynolds, interpolator_->axis(3).front(),
        interpolator_->axis(3).back());
    if (phi_axis.size() == 1U) {
        phi_aero_deg = phi_axis.front();
    } else {
        phi_aero_deg = phi_axis.front() + std::fmod(
            std::fmod(phi_aero_deg - phi_axis.front(), phi_period_deg_) +
                phi_period_deg_,
            phi_period_deg_);
    }
    return interpolator_->evaluate(
        values, {mach, alpha_deg, phi_aero_deg, reynolds});
}

AerodynamicCoefficients AerodynamicDatabase::evaluate(
    double mach, double alpha_deg, double phi_aero_deg, double reynolds) const {
    const double phi_origin = interpolator_->axis(2).front();
    const double sector = std::floor(
        (phi_aero_deg - phi_origin) / phi_period_deg_);
    const double canonical_phi =
        phi_aero_deg - sector * phi_period_deg_;
    AerodynamicCoefficients result;
    double* outputs[] = {&result.axial, &result.side_y, &result.side_z,
                         &result.roll_moment, &result.pitch_moment,
                         &result.yaw_moment};
    for (std::size_t index = 0U; index < coefficients_.size(); ++index) {
        *outputs[index] = interpolate(coefficients_[index], mach, alpha_deg,
                                      canonical_phi, reynolds);
    }
    // Fin symmetry makes coefficient magnitudes periodic, but transverse
    // body-axis components must still rotate into the requested phi sector.
    const double rotation =
        sector * phi_period_deg_ * std::acos(-1.0) / 180.0;
    const double cosine = std::cos(rotation);
    const double sine = std::sin(rotation);
    const double side_y = result.side_y;
    const double side_z = result.side_z;
    result.side_y = cosine * side_y - sine * side_z;
    result.side_z = sine * side_y + cosine * side_z;
    const double pitch = result.pitch_moment;
    const double yaw = result.yaw_moment;
    result.pitch_moment = cosine * pitch - sine * yaw;
    result.yaw_moment = sine * pitch + cosine * yaw;
    return result;
}

double AerodynamicDatabase::referenceArea() const noexcept { return reference_area_; }
double AerodynamicDatabase::referenceLength() const noexcept { return reference_length_; }
double AerodynamicDatabase::momentReferenceX() const noexcept { return moment_reference_x_; }
double AerodynamicDatabase::minimumAlphaDeg() const noexcept {
    return interpolator_->axis(1).front();
}
double AerodynamicDatabase::maximumAlphaDeg() const noexcept {
    return interpolator_->axis(1).back();
}
double AerodynamicDatabase::stabilityReferenceLength() const noexcept {
    return stability_reference_length_;
}

Aerodynamics::Aerodynamics(std::string name,
                           AerodynamicDatabase database)
    : ForceTorque(std::move(name)),
      database_(std::move(database)) {
    addStateVariable(center_pressure_x_, "centerPressureX", "m");
    addStateVariable(static_stability_margin_, "staticStabilityMargin");
}

void Aerodynamics::setConditions(AerodynamicConditions conditions) {
    if (conditions.mach < 0.0 || conditions.reynolds <= 0.0 ||
        conditions.dynamic_pressure < 0.0) {
        throw std::invalid_argument(
            "Mach and dynamic pressure must be nonnegative and Reynolds positive");
    }
    conditions_ = conditions;
}

const AerodynamicConditions& Aerodynamics::conditions() const noexcept {
    return conditions_;
}

const AerodynamicCoefficients& Aerodynamics::coefficients() const noexcept {
    return coefficients_;
}

const AerodynamicLoad& Aerodynamics::load() const noexcept { return load_; }
double Aerodynamics::centerPressureX() const noexcept {
    return center_pressure_x_;
}
double Aerodynamics::staticStabilityMargin() const noexcept {
    return static_stability_margin_;
}

AerodynamicLoad Aerodynamics::calculateLoad(
    const AerodynamicConditions& conditions) const {
    if (conditions.mach < 0.0 || conditions.reynolds <= 0.0 ||
        conditions.dynamic_pressure < 0.0) {
        throw std::invalid_argument(
            "Mach and dynamic pressure must be nonnegative and Reynolds positive");
    }
    const AerodynamicCoefficients coefficients = database_.evaluate(
        conditions.mach, conditions.alpha_deg,
        conditions.phi_aero_deg, conditions.reynolds);
    const double force_scale =
        conditions.dynamic_pressure * database_.referenceArea();
    const double moment_scale = force_scale * database_.referenceLength();
    return {
        Vector3d{-force_scale * coefficients.axial,
                 -force_scale * coefficients.side_y,
                 -force_scale * coefficients.side_z},
        Vector3d{moment_scale * coefficients.roll_moment,
                 -moment_scale * coefficients.pitch_moment,
                 -moment_scale * coefficients.yaw_moment}};
}

void Aerodynamics::update() {
    coefficients_ = database_.evaluate(
        conditions_.mach, conditions_.alpha_deg,
        conditions_.phi_aero_deg, conditions_.reynolds);
    const double force_scale =
        conditions_.dynamic_pressure * database_.referenceArea();
    const double moment_scale = force_scale * database_.referenceLength();
    const bool outside_alpha_envelope =
        conditions_.alpha_deg < database_.minimumAlphaDeg() ||
        conditions_.alpha_deg > database_.maximumAlphaDeg();
    const double relative_speed =
        conditions_.relative_air_velocity_in_body.norm();
    const double database_cp = database_.evaluateCenterPressureX(
        conditions_.mach, conditions_.alpha_deg,
        conditions_.phi_aero_deg, conditions_.reynolds);
    // Aero x runs nose-to-tail; vehicle x runs engine-to-nose.
    center_pressure_x_ = database_.referenceLength() - database_cp;
    if (isAttachedToVehicle()) {
        static_stability_margin_ =
            (attachedVehicle().rigidBody().cg().x() - center_pressure_x_) /
            database_.stabilityReferenceLength();
    } else {
        static_stability_margin_ = 0.0;
    }
    if (outside_alpha_envelope && relative_speed > 1.0e-12) {
        // Outside the validated coefficient envelope, use dissipative drag at
        // the CG. This cannot add energy or create an unvalidated aero moment.
        const double drag_coefficient =
            std::clamp(std::abs(coefficients_.axial), 0.05, 2.0);
        load_ = {
            -force_scale * drag_coefficient *
                conditions_.relative_air_velocity_in_body / relative_speed,
            Vector3d::Zero()};
        setLoad(load_.force, load_.torque, attachedVehicle().rigidBody().cg());
    } else {
        // CY/CZ are incidence magnitudes; aerodynamic force opposes the
        // corresponding body-relative air velocity. Apply that force at CP,
        // avoiding mixed moment-axis and axial-origin conventions.
        load_ = {
            Vector3d{-force_scale * coefficients_.axial,
                     -force_scale * coefficients_.side_y,
                     -force_scale * coefficients_.side_z},
            Vector3d{moment_scale * coefficients_.roll_moment, 0.0, 0.0}};
        setLoad(load_.force, load_.torque, Vector3d{
            isAttachedToVehicle() ? center_pressure_x_
                                  : database_.momentReferenceX(),
            0.0, 0.0});
    }
}

}  // namespace rsim
