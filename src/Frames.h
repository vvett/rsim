#pragma once

#include <Eigen/Geometry>

#include <memory>
#include <string>
#include <vector>

namespace rsim {

/**
 * A rigid transform named `a_from_b` converts positions and vectors expressed
 * in frame B into frame A.
 */
class Transform {
public:
    Transform();
    Transform(Eigen::Quaterniond rotation, Eigen::Vector3d translation);

    [[nodiscard]] static Transform identity();
    [[nodiscard]] const Eigen::Quaterniond& rotation() const noexcept;
    [[nodiscard]] const Eigen::Vector3d& translation() const noexcept;
    [[nodiscard]] Eigen::Vector3d applyPosition(
        const Eigen::Vector3d& position) const;
    [[nodiscard]] Eigen::Vector3d applyVector(
        const Eigen::Vector3d& vector) const;
    [[nodiscard]] Transform inverse() const;

private:
    Eigen::Quaterniond rotation_;
    Eigen::Vector3d translation_;
};

/** Compose `a_from_b * b_from_c` to obtain `a_from_c`. */
[[nodiscard]] Transform operator*(const Transform& a_from_b,
                                  const Transform& b_from_c);

/** Position, velocity, and acceleration of one object in one frame. */
struct MotionState {
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
};

/**
 * The pose and relative motion of a child frame with respect to its parent.
 * Every vector field is expressed in the parent frame.
 */
class FrameMotion {
public:
    FrameMotion();
    FrameMotion(Transform parent_from_child,
                Eigen::Vector3d child_origin_velocity,
                Eigen::Vector3d child_origin_acceleration,
                Eigen::Vector3d child_angular_velocity,
                Eigen::Vector3d child_angular_acceleration);

    [[nodiscard]] static FrameMotion fixed(Transform parent_from_child);

    [[nodiscard]] const Transform& transform() const noexcept;
    [[nodiscard]] const Eigen::Vector3d& childOriginVelocity() const noexcept;
    [[nodiscard]] const Eigen::Vector3d& childOriginAcceleration() const noexcept;
    [[nodiscard]] const Eigen::Vector3d& childAngularVelocity() const noexcept;
    [[nodiscard]] const Eigen::Vector3d& childAngularAcceleration() const noexcept;

    /** Convert a full position/velocity/acceleration state into the parent. */
    [[nodiscard]] MotionState convertState(const MotionState& child_state) const;

    /** Return the same frame relationship in the opposite direction. */
    [[nodiscard]] FrameMotion inverse() const;

private:
    Transform parent_from_child_;
    Eigen::Vector3d child_origin_velocity_;
    Eigen::Vector3d child_origin_acceleration_;
    Eigen::Vector3d child_angular_velocity_;
    Eigen::Vector3d child_angular_acceleration_;
};

/** Compose `a_from_b * b_from_c` to obtain `a_from_c`, including rates. */
[[nodiscard]] FrameMotion operator*(const FrameMotion& a_from_b,
                                    const FrameMotion& b_from_c);

/**
 * Defines how a child frame is positioned and moving relative to its parent.
 * A definition owns the latest calculated FrameMotion.
 */
class FrameDefinition {
public:
    virtual ~FrameDefinition() = default;

    /** Recalculate the frame relationship from current simulation inputs. */
    virtual void update() = 0;

    /** Return the latest child-to-parent relationship. */
    [[nodiscard]] virtual const FrameMotion& motionIntoParent() const = 0;
};

/** Defines a parent/child relationship that does not change with time. */
class FixedFrameDefinition final : public FrameDefinition {
public:
    explicit FixedFrameDefinition(Transform parent_from_child);
    explicit FixedFrameDefinition(FrameMotion parent_from_child);

    void update() override;
    [[nodiscard]] const FrameMotion& motionIntoParent() const override;

private:
    FrameMotion motion_;
};

/** Declares whether a frame is suitable for applying inertial equations. */
enum class InertialStatus {
    Inertial,
    NonInertial
};

/** A named node in a parent-linked tree of Cartesian coordinate frames. */
class Frame {
public:
    /** Construct a root frame. */
    Frame(std::string name, InertialStatus inertial_status);

    /** Construct a child frame and its definition relative to the parent. */
    Frame(std::string name,
          InertialStatus inertial_status,
          std::shared_ptr<Frame> parent,
          std::unique_ptr<FrameDefinition> definition);

    void update();
    [[nodiscard]] const std::string& name() const noexcept;
    [[nodiscard]] InertialStatus inertialStatus() const noexcept;
    [[nodiscard]] bool isInertial() const noexcept;
    [[nodiscard]] bool isRoot() const noexcept;
    [[nodiscard]] const std::shared_ptr<Frame>& parent() const noexcept;
    [[nodiscard]] const FrameMotion& motionIntoParent() const;

private:
    std::string name_;
    InertialStatus inertial_status_;
    std::shared_ptr<Frame> parent_;
    std::unique_ptr<FrameDefinition> definition_;
};

/** Updates a frame tree and converts transforms or motion states within it. */
class FrameGraph {
public:
    void addFrame(std::shared_ptr<Frame> frame);

    /** Update every registered definition in parent-before-child order. */
    void update();

    /** Return `destination_from_source`, including relative frame rates. */
    [[nodiscard]] FrameMotion motion(const Frame& destination,
                                     const Frame& source) const;

    /** Return only the instantaneous pose from source into destination. */
    [[nodiscard]] Transform transform(const Frame& destination,
                                      const Frame& source) const;

    /** Convert an object's full motion state between connected frames. */
    [[nodiscard]] MotionState convertState(const Frame& destination,
                                           const Frame& source,
                                           const MotionState& source_state) const;

private:
    std::vector<std::shared_ptr<Frame>> frames_;
};

}  // namespace rsim
