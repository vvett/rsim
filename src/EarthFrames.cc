#include "EarthFrames.h"

#include <Eigen/Geometry>

#include <cmath>

namespace rsim {

namespace {

constexpr double j2000_julian_date = 2451545.0;
constexpr double seconds_per_day = 86400.0;
constexpr double degrees_to_radians =
    3.141592653589793238462643383279502884 / 180.0;
constexpr double sidereal_angle_at_j2000_degrees = 280.46061837;
constexpr double sidereal_degrees_per_day = 360.98564736629;
constexpr double earth_rotation_radians_per_second =
    sidereal_degrees_per_day * degrees_to_radians / seconds_per_day;

double siderealAngle(const JulianDate& date) {
    const double days_since_j2000 =
        (date.first_part - j2000_julian_date) + date.second_part;
    const double angle_degrees =
        sidereal_angle_at_j2000_degrees +
        sidereal_degrees_per_day * days_since_j2000;
    return std::remainder(angle_degrees, 360.0) * degrees_to_radians;
}

}  // namespace

J2000EarthRotation::J2000EarthRotation(const JulianDate& current_date)
    : current_date_(current_date) {
    update();
}

void J2000EarthRotation::update() {
    const Eigen::Quaterniond j2000_from_ecef{
        Eigen::AngleAxisd{siderealAngle(current_date_),
                          Eigen::Vector3d::UnitZ()}
    };

    j2000_from_ecef_ = FrameMotion{
        Transform{j2000_from_ecef, Eigen::Vector3d::Zero()},
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d{0.0, 0.0, earth_rotation_radians_per_second},
        Eigen::Vector3d::Zero()
    };
}

const FrameMotion& J2000EarthRotation::motionIntoParent() const {
    return j2000_from_ecef_;
}

}  // namespace rsim
