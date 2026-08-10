#include "Frames.h"

#include <Eigen/Geometry>

#include <cassert>
#include <cmath>
#include <memory>

namespace {

constexpr double tolerance = 1.0e-12;

bool nearlyEqual(const Eigen::Vector3d& left, const Eigen::Vector3d& right) {
    return (left - right).norm() < tolerance;
}

class MovingFrameDefinition final : public rsim::FrameDefinition {
public:
    void update() override {
        position_ += 1.0;
        motion_ = rsim::FrameMotion::fixed(rsim::Transform{
            Eigen::Quaterniond::Identity(),
            Eigen::Vector3d{position_, 0.0, 0.0}});
    }

    const rsim::FrameMotion& motionIntoParent() const override {
        return motion_;
    }

private:
    double position_ = 0.0;
    rsim::FrameMotion motion_;
};

}  // namespace

int main() {
    auto root = std::make_shared<rsim::Frame>(
        "root", rsim::InertialStatus::Inertial);
    auto ecef = std::make_shared<rsim::Frame>(
        "ecef",
        rsim::InertialStatus::NonInertial,
        root,
        std::make_unique<rsim::FixedFrameDefinition>(rsim::Transform{
            Eigen::Quaterniond::Identity(),
            Eigen::Vector3d{1.0, 0.0, 0.0}
        }));
    auto ned = std::make_shared<rsim::Frame>(
        "ned",
        rsim::InertialStatus::NonInertial,
        ecef,
        std::make_unique<rsim::FixedFrameDefinition>(rsim::Transform{
            Eigen::Quaterniond::Identity(),
            Eigen::Vector3d{0.0, 2.0, 0.0}
        }));
    auto moving = std::make_shared<rsim::Frame>(
        "moving",
        rsim::InertialStatus::NonInertial,
        root,
        std::make_unique<MovingFrameDefinition>());

    assert(root->isInertial());
    assert(!ecef->isInertial());
    assert(ecef->inertialStatus() == rsim::InertialStatus::NonInertial);

    rsim::FrameGraph graph;
    graph.addFrame(ned);
    graph.addFrame(root);
    graph.addFrame(ecef);
    graph.addFrame(moving);

    graph.update();

    const Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    assert(nearlyEqual(graph.transform(*root, *ned).applyPosition(origin),
                       Eigen::Vector3d{1.0, 2.0, 0.0}));
    assert(nearlyEqual(graph.transform(*ned, *root).applyPosition(origin),
                       Eigen::Vector3d{-1.0, -2.0, 0.0}));
    assert(nearlyEqual(graph.transform(*root, *moving).applyPosition(origin),
                       Eigen::Vector3d{1.0, 0.0, 0.0}));

    graph.update();
    assert(nearlyEqual(graph.transform(*root, *moving).applyPosition(origin),
                       Eigen::Vector3d{2.0, 0.0, 0.0}));

    const rsim::Transform quarter_turn{
        Eigen::Quaterniond{Eigen::AngleAxisd{
            std::acos(-1.0) / 2.0, Eigen::Vector3d::UnitZ()}},
        Eigen::Vector3d::Zero()
    };
    assert(nearlyEqual(quarter_turn.applyVector(Eigen::Vector3d::UnitX()),
                       Eigen::Vector3d::UnitY()));

    const rsim::FrameMotion rotating_frame{
        rsim::Transform::identity(),
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d{0.0, 0.0, 2.0},
        Eigen::Vector3d::Zero()
    };
    const rsim::MotionState stationary_in_rotating_frame{
        Eigen::Vector3d{3.0, 0.0, 0.0},
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Zero()
    };
    const rsim::MotionState state_in_parent =
        rotating_frame.convertState(stationary_in_rotating_frame);
    assert(nearlyEqual(state_in_parent.velocity,
                       Eigen::Vector3d{0.0, 6.0, 0.0}));
    assert(nearlyEqual(state_in_parent.acceleration,
                       Eigen::Vector3d{-12.0, 0.0, 0.0}));

    const rsim::MotionState round_trip =
        rotating_frame.inverse().convertState(state_in_parent);
    assert(nearlyEqual(round_trip.position,
                       stationary_in_rotating_frame.position));
    assert(nearlyEqual(round_trip.velocity,
                       stationary_in_rotating_frame.velocity));
    assert(nearlyEqual(round_trip.acceleration,
                       stationary_in_rotating_frame.acceleration));
}
