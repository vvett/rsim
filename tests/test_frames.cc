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

class MovingProvider final : public rsim::TransformProvider {
public:
    void update() override {
        position_ += 1.0;
        transform_ = rsim::Transform{
            Eigen::Quaterniond::Identity(),
            Eigen::Vector3d{position_, 0.0, 0.0}
        };
    }

    const rsim::Transform& parentFromChild() const override {
        return transform_;
    }

private:
    double position_ = 0.0;
    rsim::Transform transform_;
};

}  // namespace

int main() {
    auto root = std::make_shared<rsim::Frame>("root");
    auto ecef = std::make_shared<rsim::Frame>(
        "ecef",
        root,
        std::make_unique<rsim::FixedTransformProvider>(rsim::Transform{
            Eigen::Quaterniond::Identity(),
            Eigen::Vector3d{1.0, 0.0, 0.0}
        }));
    auto ned = std::make_shared<rsim::Frame>(
        "ned",
        ecef,
        std::make_unique<rsim::FixedTransformProvider>(rsim::Transform{
            Eigen::Quaterniond::Identity(),
            Eigen::Vector3d{0.0, 2.0, 0.0}
        }));
    auto moving = std::make_shared<rsim::Frame>(
        "moving", root, std::make_unique<MovingProvider>());

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
}
