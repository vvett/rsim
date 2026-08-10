#include "EarthFrames.h"

#include <Eigen/Core>

#include <cassert>
#include <cmath>
#include <memory>

namespace {

constexpr double degrees_to_radians =
    3.141592653589793238462643383279502884 / 180.0;

bool nearlyEqual(double left, double right, double tolerance) {
    return std::abs(left - right) < tolerance;
}

}  // namespace

int main() {
    rsim::JulianDate current_date{2451545.0, 0.0};

    auto j2000 = std::make_shared<rsim::Frame>(
        "J2000", rsim::InertialStatus::Inertial);
    auto ecef = std::make_shared<rsim::Frame>(
        "ECEF",
        rsim::InertialStatus::NonInertial,
        j2000,
        std::make_unique<rsim::J2000EarthRotation>(current_date));

    rsim::FrameGraph frames;
    frames.addFrame(j2000);
    frames.addFrame(ecef);
    frames.update();

    const Eigen::Matrix3d j2000_from_ecef =
        frames.transform(*j2000, *ecef).rotation().toRotationMatrix();
    const double angle = 280.46061837 * degrees_to_radians;
    assert(nearlyEqual(j2000_from_ecef(0, 0), std::cos(angle), 1.0e-12));
    assert(nearlyEqual(j2000_from_ecef(1, 0), std::sin(angle), 1.0e-12));

    const rsim::MotionState fixed_on_earth{
        Eigen::Vector3d{6378137.0, 0.0, 0.0},
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Zero()
    };
    const rsim::MotionState state_in_j2000 =
        frames.convertState(*j2000, *ecef, fixed_on_earth);

    assert(state_in_j2000.velocity.norm() > 400.0);

    const rsim::MotionState round_trip =
        frames.convertState(*ecef, *j2000, state_in_j2000);
    assert((round_trip.position - fixed_on_earth.position).norm() < 1.0e-6);
    assert(round_trip.velocity.norm() < 1.0e-6);
    assert(round_trip.acceleration.norm() < 1.0e-6);

    current_date.second_part += 1.0;
    frames.update();
    const Eigen::Vector3d next_x_axis =
        frames.transform(*j2000, *ecef).applyVector(Eigen::Vector3d::UnitX());
    assert((next_x_axis - j2000_from_ecef * Eigen::Vector3d::UnitX()).norm() >
           1.0e-3);
}
