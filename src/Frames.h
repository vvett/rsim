#pragma once

#include "Model.h"

#include <Eigen/Geometry>

#include <memory>
#include <string>
#include <vector>

namespace rsim {

// A rigid transform named `a_from_b` converts positions and vectors expressed
// in frame B into frame A.
class Transform {
public:
    // Construct an identity transform.
    Transform();

    // Construct a transform from a rotation and destination-frame translation.
    Transform(Eigen::Quaterniond rotation, Eigen::Vector3d translation);

    // Return a transform with no rotation or translation.
    [[nodiscard]] static Transform identity();

    // Return the rotation from the source axes into the destination axes.
    [[nodiscard]] const Eigen::Quaterniond& rotation() const noexcept;

    // Return the source origin's position expressed in the destination frame.
    [[nodiscard]] const Eigen::Vector3d& translation() const noexcept;

    // Rotate and translate a source-frame position into the destination frame.
    [[nodiscard]] Eigen::Vector3d applyPosition(
        const Eigen::Vector3d& position) const;

    // Rotate a free vector into the destination frame without translating it.
    [[nodiscard]] Eigen::Vector3d applyVector(
        const Eigen::Vector3d& vector) const;

    // Return the transform that performs the opposite coordinate conversion.
    [[nodiscard]] Transform inverse() const;

private:
    Eigen::Quaterniond rotation_;
    Eigen::Vector3d translation_;
};

// Compose `a_from_b * b_from_c` to obtain `a_from_c`.
[[nodiscard]] Transform operator*(const Transform& a_from_b,
                                  const Transform& b_from_c);

// Position, velocity, and acceleration of one object in one frame.
struct MotionState {
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
};

// The pose and relative motion of a child frame with respect to its parent.
// Every vector field is expressed in the parent frame.
class FrameMotion {
public:
    // Construct an identity relationship with all relative rates set to zero.
    FrameMotion();

    // Construct a complete child-to-parent pose and relative-motion relationship.
    FrameMotion(Transform parent_from_child,
                Eigen::Vector3d child_origin_velocity,
                Eigen::Vector3d child_origin_acceleration,
                Eigen::Vector3d child_angular_velocity,
                Eigen::Vector3d child_angular_acceleration);

    // Construct a pose whose relative linear and angular rates are all zero.
    [[nodiscard]] static FrameMotion fixed(Transform parent_from_child);

    // Return the instantaneous transform from the child into the parent.
    [[nodiscard]] const Transform& transform() const noexcept;

    // Return the child origin's velocity expressed in the parent frame.
    [[nodiscard]] const Eigen::Vector3d& childOriginVelocity() const noexcept;

    // Return the child origin's acceleration expressed in the parent frame.
    [[nodiscard]] const Eigen::Vector3d& childOriginAcceleration() const noexcept;

    // Return the child axes' angular velocity expressed in the parent frame.
    [[nodiscard]] const Eigen::Vector3d& childAngularVelocity() const noexcept;

    // Return the child axes' angular acceleration expressed in the parent frame.
    [[nodiscard]] const Eigen::Vector3d& childAngularAcceleration() const noexcept;

    // Convert a full position/velocity/acceleration state into the parent.
    [[nodiscard]] MotionState convertState(const MotionState& child_state) const;

    // Return the same frame relationship in the opposite direction.
    [[nodiscard]] FrameMotion inverse() const;

private:
    Transform parent_from_child_;
    Eigen::Vector3d child_origin_velocity_;
    Eigen::Vector3d child_origin_acceleration_;
    Eigen::Vector3d child_angular_velocity_;
    Eigen::Vector3d child_angular_acceleration_;
};

// Compose `a_from_b * b_from_c` to obtain `a_from_c`, including rates.
[[nodiscard]] FrameMotion operator*(const FrameMotion& a_from_b,
                                    const FrameMotion& b_from_c);

// Defines how a child frame is positioned and moving relative to its parent.
// A definition owns the latest calculated FrameMotion.
class FrameDefinition {
public:
    // Destroy a concrete frame definition correctly through the base interface.
    virtual ~FrameDefinition() = default;

    // Recalculate the frame relationship from current simulation inputs.
    virtual void update() = 0;

    // Return the latest child-to-parent relationship.
    [[nodiscard]] const FrameMotion& motionIntoParent() const noexcept;

protected:
    // Construct a definition with an identity motion snapshot.
    FrameDefinition() = default;

    // Construct a definition with an initial motion snapshot.
    explicit FrameDefinition(FrameMotion motion);

    FrameMotion motion_;
};

// Defines a parent/child relationship that does not change with time.
class FixedFrameDefinition final : public FrameDefinition {
public:
    // Construct an unchanging frame relationship from a fixed transform.
    explicit FixedFrameDefinition(Transform parent_from_child);

    // Construct an unchanging frame relationship from a supplied motion snapshot.
    explicit FixedFrameDefinition(FrameMotion parent_from_child);

    // Leave the fixed frame relationship unchanged.
    void update() override;
};

// Declares whether a frame is suitable for applying inertial equations.
enum class InertialStatus {
    Inertial,
    NonInertial
};

// A named node in a parent-linked tree of Cartesian coordinate frames.
class Frame {
public:
    // Construct a root frame.
    Frame(std::string name, InertialStatus inertial_status);

    // Construct a child frame and its definition relative to the parent.
    Frame(std::string name,
          InertialStatus inertial_status,
          std::shared_ptr<Frame> parent,
          std::unique_ptr<FrameDefinition> definition);

    // Recalculate this frame's definition when it has one.
    void update();

    // Return the frame's identifying name.
    [[nodiscard]] const std::string& name() const noexcept;

    // Return whether the frame is classified as inertial or noninertial.
    [[nodiscard]] InertialStatus inertialStatus() const noexcept;

    // Return true when the frame is classified as inertial.
    [[nodiscard]] bool isInertial() const noexcept;

    // Return true when this frame has no parent.
    [[nodiscard]] bool isRoot() const noexcept;

    // Return the parent frame pointer, or an empty pointer for a root frame.
    [[nodiscard]] const std::shared_ptr<Frame>& parent() const noexcept;

    // Return this child frame's latest motion relationship into its parent.
    [[nodiscard]] const FrameMotion& motionIntoParent() const;

private:
    std::string name_;
    InertialStatus inertial_status_;
    std::shared_ptr<Frame> parent_;
    std::unique_ptr<FrameDefinition> definition_;
};

// Updates a frame tree and converts transforms or motion states within it.
class FrameGraph : public Model {
public:
    // Construct a schedulable frame graph.
    explicit FrameGraph(
        std::string name = "frames",
        UpdateRate update_rate = UpdateRate::everyStep());

    // Register a uniquely named frame with the graph.
    void addFrame(std::shared_ptr<Frame> frame);

    // Return a registered frame by its unique name.
    [[nodiscard]] std::shared_ptr<Frame> frame(const std::string& name) const;

    // Update every registered definition in parent-before-child order.
    void update() override;

    // Return `destination_from_source`, including relative frame rates.
    [[nodiscard]] FrameMotion motion(const Frame& destination,
                                     const Frame& source) const;

    // Return only the instantaneous pose from source into destination.
    [[nodiscard]] Transform transform(const Frame& destination,
                                      const Frame& source) const;

    // Convert an object's full motion state between connected frames.
    [[nodiscard]] MotionState convertState(const Frame& destination,
                                           const Frame& source,
                                           const MotionState& source_state) const;

private:
    std::vector<std::shared_ptr<Frame>> frames_;
};

}  // namespace rsim
