#include "Frames.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace rsim {

namespace {

struct RootTransform {
    const Frame* root;
    Transform root_from_frame;
};

RootTransform findRootTransform(const Frame& frame) {
    const Frame* current = &frame;
    Transform current_from_frame = Transform::identity();
    std::unordered_set<const Frame*> visited;

    while (!current->isRoot()) {
        if (!visited.insert(current).second) {
            throw std::logic_error("cycle detected in frame graph");
        }

        current_from_frame = current->parentFromThis() * current_from_frame;
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

FixedTransformProvider::FixedTransformProvider(Transform parent_from_child)
    : parent_from_child_(std::move(parent_from_child)) {}

void FixedTransformProvider::update() {}

const Transform& FixedTransformProvider::parentFromChild() const {
    return parent_from_child_;
}

Frame::Frame(std::string name) : name_(std::move(name)) {
    if (name_.empty()) {
        throw std::invalid_argument("frame name must not be empty");
    }
}

Frame::Frame(std::string name,
             std::shared_ptr<Frame> parent,
             std::unique_ptr<TransformProvider> provider)
    : name_(std::move(name)),
      parent_(std::move(parent)),
      provider_(std::move(provider)) {
    if (name_.empty()) {
        throw std::invalid_argument("frame name must not be empty");
    }
    if (!parent_) {
        throw std::invalid_argument("child frame must have a parent");
    }
    if (!provider_) {
        throw std::invalid_argument("child frame must have a transform provider");
    }
}

void Frame::update() {
    if (provider_) {
        provider_->update();
    }
}

const std::string& Frame::name() const noexcept {
    return name_;
}

bool Frame::isRoot() const noexcept {
    return !parent_;
}

const std::shared_ptr<Frame>& Frame::parent() const noexcept {
    return parent_;
}

const Transform& Frame::parentFromThis() const {
    if (!provider_) {
        throw std::logic_error("root frame has no parent transform");
    }
    return provider_->parentFromChild();
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

Transform FrameGraph::transform(const Frame& destination,
                                const Frame& source) const {
    const RootTransform root_from_destination =
        findRootTransform(destination);
    const RootTransform root_from_source = findRootTransform(source);

    if (root_from_destination.root != root_from_source.root) {
        throw std::invalid_argument(
            "cannot transform between disconnected frame trees");
    }

    return root_from_destination.root_from_frame.inverse() *
           root_from_source.root_from_frame;
}

}  // namespace rsim
