#include "Frames.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace rsim {

namespace {

struct RootMotion {
    const Frame* root;
    FrameMotion root_from_frame;
};

RootMotion findRootMotion(const Frame& frame) {
    const Frame* current = &frame;
    FrameMotion current_from_frame;
    std::unordered_set<const Frame*> visited;

    while (!current->isRoot()) {
        if (!visited.insert(current).second) {
            throw std::logic_error("cycle detected in frame graph");
        }

        current_from_frame = current->motionIntoParent() * current_from_frame;
        current = current->parent().get();
    }

    return {current, current_from_frame};
}

void updateFrame(Frame& frame,
                 std::unordered_set<const Frame*>& updating,
                 std::unordered_set<const Frame*>& updated) {
    if (updated.count(&frame) != 0U) {
        return;
    }
    if (!updating.insert(&frame).second) {
        throw std::logic_error("cycle detected in frame graph");
    }

    if (!frame.isRoot()) {
        updateFrame(*frame.parent(), updating, updated);
    }

    frame.update();
    updating.erase(&frame);
    updated.insert(&frame);
}

}  // namespace

Transform::Transform()
    : rotation_(Eigen::Quaterniond::Identity()),
      translation_(Eigen::Vector3d::Zero()) {}

Transform::Transform(Eigen::Quaterniond rotation, Eigen::Vector3d translation)
    : rotation_(std::move(rotation)), translation_(std::move(translation)) {
    if (rotation_.norm() <= 0.0) {
        throw std::invalid_argument("transform rotation must be nonzero");
    }
    rotation_.normalize();
}

Transform Transform::identity() {
    return {};
}

const Eigen::Quaterniond& Transform::rotation() const noexcept {
    return rotation_;
}

const Eigen::Vector3d& Transform::translation() const noexcept {
    return translation_;
}

Eigen::Vector3d Transform::applyPosition(
    const Eigen::Vector3d& position) const {
    return rotation_ * position + translation_;
}

Eigen::Vector3d Transform::applyVector(const Eigen::Vector3d& vector) const {
    return rotation_ * vector;
}

Transform Transform::inverse() const {
    const Eigen::Quaterniond inverse_rotation = rotation_.conjugate();
    return {inverse_rotation, -(inverse_rotation * translation_)};
}

Transform operator*(const Transform& a_from_b, const Transform& b_from_c) {
    return {
        a_from_b.rotation() * b_from_c.rotation(),
        a_from_b.applyPosition(b_from_c.translation())
    };
}

FrameMotion::FrameMotion()
    : parent_from_child_(Transform::identity()),
      child_origin_velocity_(Eigen::Vector3d::Zero()),
      child_origin_acceleration_(Eigen::Vector3d::Zero()),
      child_angular_velocity_(Eigen::Vector3d::Zero()),
      child_angular_acceleration_(Eigen::Vector3d::Zero()) {}

FrameMotion::FrameMotion(
    Transform parent_from_child,
    Eigen::Vector3d child_origin_velocity,
    Eigen::Vector3d child_origin_acceleration,
    Eigen::Vector3d child_angular_velocity,
    Eigen::Vector3d child_angular_acceleration)
    : parent_from_child_(std::move(parent_from_child)),
      child_origin_velocity_(std::move(child_origin_velocity)),
      child_origin_acceleration_(std::move(child_origin_acceleration)),
      child_angular_velocity_(std::move(child_angular_velocity)),
      child_angular_acceleration_(std::move(child_angular_acceleration)) {}

FrameMotion FrameMotion::fixed(Transform parent_from_child) {
    return {std::move(parent_from_child),
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d::Zero()};
}

const Transform& FrameMotion::transform() const noexcept {
    return parent_from_child_;
}

const Eigen::Vector3d& FrameMotion::childOriginVelocity() const noexcept {
    return child_origin_velocity_;
}

const Eigen::Vector3d& FrameMotion::childOriginAcceleration() const noexcept {
    return child_origin_acceleration_;
}

const Eigen::Vector3d& FrameMotion::childAngularVelocity() const noexcept {
    return child_angular_velocity_;
}

const Eigen::Vector3d& FrameMotion::childAngularAcceleration() const noexcept {
    return child_angular_acceleration_;
}

MotionState FrameMotion::convertState(const MotionState& child_state) const {
    const Eigen::Vector3d relative_position =
        parent_from_child_.applyVector(child_state.position);
    const Eigen::Vector3d relative_velocity =
        parent_from_child_.applyVector(child_state.velocity);

    return {
        relative_position + parent_from_child_.translation(),
        relative_velocity +
            child_angular_velocity_.cross(relative_position) +
            child_origin_velocity_,
        parent_from_child_.applyVector(child_state.acceleration) +
            child_origin_acceleration_ +
            child_angular_acceleration_.cross(relative_position) +
            child_angular_velocity_.cross(
                child_angular_velocity_.cross(relative_position)) +
            2.0 * child_angular_velocity_.cross(relative_velocity)
    };
}

FrameMotion FrameMotion::inverse() const {
    const Transform child_from_parent = parent_from_child_.inverse();
    const Eigen::Vector3d& omega = child_angular_velocity_;
    const Eigen::Vector3d& alpha = child_angular_acceleration_;
    const Eigen::Vector3d& translation = parent_from_child_.translation();

    return {
        child_from_parent,
        child_from_parent.applyVector(
            -child_origin_velocity_ + omega.cross(translation)),
        child_from_parent.applyVector(
            -child_origin_acceleration_ + alpha.cross(translation) +
            2.0 * omega.cross(child_origin_velocity_) -
            omega.cross(omega.cross(translation))),
        child_from_parent.applyVector(-omega),
        child_from_parent.applyVector(-alpha)
    };
}

FrameMotion operator*(const FrameMotion& a_from_b,
                      const FrameMotion& b_from_c) {
    const Eigen::Vector3d c_origin_from_b_in_a =
        a_from_b.transform().applyVector(b_from_c.transform().translation());
    const Eigen::Vector3d c_origin_velocity_in_a =
        a_from_b.transform().applyVector(b_from_c.childOriginVelocity());
    const Eigen::Vector3d c_angular_velocity_in_a =
        a_from_b.transform().applyVector(b_from_c.childAngularVelocity());

    return {
        a_from_b.transform() * b_from_c.transform(),
        a_from_b.childOriginVelocity() + c_origin_velocity_in_a +
            a_from_b.childAngularVelocity().cross(c_origin_from_b_in_a),
        a_from_b.childOriginAcceleration() +
            a_from_b.transform().applyVector(
                b_from_c.childOriginAcceleration()) +
            a_from_b.childAngularAcceleration().cross(c_origin_from_b_in_a) +
            a_from_b.childAngularVelocity().cross(
                a_from_b.childAngularVelocity().cross(c_origin_from_b_in_a)) +
            2.0 * a_from_b.childAngularVelocity().cross(
                c_origin_velocity_in_a),
        a_from_b.childAngularVelocity() + c_angular_velocity_in_a,
        a_from_b.childAngularAcceleration() +
            a_from_b.transform().applyVector(
                b_from_c.childAngularAcceleration()) +
            a_from_b.childAngularVelocity().cross(c_angular_velocity_in_a)
    };
}

FixedFrameDefinition::FixedFrameDefinition(Transform parent_from_child)
    : motion_(FrameMotion::fixed(std::move(parent_from_child))) {}

FixedFrameDefinition::FixedFrameDefinition(FrameMotion parent_from_child)
    : motion_(std::move(parent_from_child)) {}

void FixedFrameDefinition::update() {}

const FrameMotion& FixedFrameDefinition::motionIntoParent() const {
    return motion_;
}

Frame::Frame(std::string name, InertialStatus inertial_status)
    : name_(std::move(name)), inertial_status_(inertial_status) {
    if (name_.empty()) {
        throw std::invalid_argument("frame name must not be empty");
    }
}

Frame::Frame(std::string name,
             InertialStatus inertial_status,
             std::shared_ptr<Frame> parent,
             std::unique_ptr<FrameDefinition> definition)
    : name_(std::move(name)),
      inertial_status_(inertial_status),
      parent_(std::move(parent)),
      definition_(std::move(definition)) {
    if (name_.empty()) {
        throw std::invalid_argument("frame name must not be empty");
    }
    if (!parent_) {
        throw std::invalid_argument("child frame must have a parent");
    }
    if (!definition_) {
        throw std::invalid_argument("child frame must have a frame definition");
    }
}

void Frame::update() {
    if (definition_) {
        definition_->update();
    }
}

const std::string& Frame::name() const noexcept {
    return name_;
}

InertialStatus Frame::inertialStatus() const noexcept {
    return inertial_status_;
}

bool Frame::isInertial() const noexcept {
    return inertial_status_ == InertialStatus::Inertial;
}

bool Frame::isRoot() const noexcept {
    return !parent_;
}

const std::shared_ptr<Frame>& Frame::parent() const noexcept {
    return parent_;
}

const FrameMotion& Frame::motionIntoParent() const {
    if (!definition_) {
        throw std::logic_error("root frame has no parent motion");
    }
    return definition_->motionIntoParent();
}

void FrameGraph::addFrame(std::shared_ptr<Frame> frame) {
    if (!frame) {
        throw std::invalid_argument("cannot add a null frame");
    }

    const auto duplicate = std::find_if(
        frames_.begin(), frames_.end(),
        [&frame](const std::shared_ptr<Frame>& existing) {
            return existing.get() == frame.get() ||
                   existing->name() == frame->name();
        });
    if (duplicate != frames_.end()) {
        throw std::invalid_argument("frame already exists: " + frame->name());
    }

    frames_.push_back(std::move(frame));
}

void FrameGraph::update() {
    std::unordered_set<const Frame*> updating;
    std::unordered_set<const Frame*> updated;
    for (const auto& frame : frames_) {
        updateFrame(*frame, updating, updated);
    }
}

FrameMotion FrameGraph::motion(const Frame& destination,
                               const Frame& source) const {
    const RootMotion root_from_destination = findRootMotion(destination);
    const RootMotion root_from_source = findRootMotion(source);

    if (root_from_destination.root != root_from_source.root) {
        throw std::invalid_argument(
            "cannot transform between disconnected frame trees");
    }

    return root_from_destination.root_from_frame.inverse() *
           root_from_source.root_from_frame;
}

Transform FrameGraph::transform(const Frame& destination,
                                const Frame& source) const {
    return motion(destination, source).transform();
}

MotionState FrameGraph::convertState(
    const Frame& destination,
    const Frame& source,
    const MotionState& source_state) const {
    return motion(destination, source).convertState(source_state);
}

}  // namespace rsim
