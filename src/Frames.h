#pragma once

#include <Eigen/Geometry>

#include <memory>
#include <string>
#include <vector>

namespace rsim {

/**
 * @brief A rigid-body transform between two Cartesian coordinate frames.
 *
 * A transform stores a unit quaternion @f$q@f$ and a translation @f$t@f$.
 * The convention throughout this API is `destination_from_source`: a transform
 * named `a_from_b` maps coordinates expressed in frame B into frame A.
 * For a position @f$p_b@f$,
 *
 * @f[
 *     p_a = R(q_{ab})p_b + t_{ab}.
 * @f]
 *
 * Translation is the location of B's origin expressed in A. Direction vectors
 * do not represent a location, so they are rotated but never translated.
 *
 * This class represents only rigid Cartesian transformations. Latitude,
 * longitude, and altitude require a nonlinear geodetic conversion before they
 * can participate in this frame graph.
 */
class Transform {
public:
    /**
     * @brief Construct the identity transform.
     *
     * The identity has quaternion `(w, x, y, z) = (1, 0, 0, 0)` and zero
     * translation, so applying it leaves every position and vector unchanged.
     */
    Transform();

    /**
     * @brief Construct a transform from a rotation and translation.
     * @param rotation Quaternion that rotates source coordinates into the
     *        destination frame. It is normalized by the constructor.
     * @param translation Position of the source origin expressed in the
     *        destination frame.
     * @throws std::invalid_argument if @p rotation has zero norm.
     */
    Transform(Eigen::Quaterniond rotation, Eigen::Vector3d translation);

    /** @brief Return an identity transform. */
    [[nodiscard]] static Transform identity();

    /**
     * @brief Return the normalized source-to-destination rotation.
     * @return A reference valid for the lifetime of this transform.
     */
    [[nodiscard]] const Eigen::Quaterniond& rotation() const noexcept;

    /**
     * @brief Return the source origin expressed in destination coordinates.
     * @return A reference valid for the lifetime of this transform.
     */
    [[nodiscard]] const Eigen::Vector3d& translation() const noexcept;

    /**
     * @brief Transform a position from source coordinates to destination
     *        coordinates.
     * @param position Source-frame position @f$p_s@f$.
     * @return @f$R_{ds}p_s + t_{ds}@f$.
     */
    [[nodiscard]] Eigen::Vector3d applyPosition(
        const Eigen::Vector3d& position) const;

    /**
     * @brief Transform a free vector from source coordinates to destination
     *        coordinates.
     * @param vector Source-frame vector @f$v_s@f$.
     * @return @f$R_{ds}v_s@f$. Translation is intentionally omitted.
     *
     * Use this for quantities such as a direction, displacement, force, or
     * angular velocity. Use applyPosition() for a point with a location.
     */
    [[nodiscard]] Eigen::Vector3d applyVector(
        const Eigen::Vector3d& vector) const;

    /**
     * @brief Return the transform in the opposite direction.
     *
     * If this transform is @f$T_{ab}=(R_{ab},t_{ab})@f$, its inverse is
     *
     * @f[
     *     T_{ba}=(R_{ab}^{T},-R_{ab}^{T}t_{ab}).
     * @f]
     *
     * For a unit quaternion, conjugation produces @f$R^T=R^{-1}@f$.
     */
    [[nodiscard]] Transform inverse() const;

private:
    Eigen::Quaterniond rotation_;
    Eigen::Vector3d translation_;
};

/**
 * @brief Compose two compatible rigid transforms.
 * @param a_from_b Transform from B coordinates into A coordinates.
 * @param b_from_c Transform from C coordinates into B coordinates.
 * @return The transform from C coordinates directly into A coordinates.
 *
 * Composition follows `a_from_b * b_from_c = a_from_c`. Its components are
 *
 * @f[
 * R_{ac}=R_{ab}R_{bc}, \qquad
 * t_{ac}=R_{ab}t_{bc}+t_{ab}.
 * @f]
 *
 * As with matrix multiplication, the right-hand transform is applied first.
 */
[[nodiscard]] Transform operator*(const Transform& a_from_b,
                                  const Transform& b_from_c);

/**
 * @brief Interface for a possibly time-varying parent-from-child transform.
 *
 * A provider separates frame topology from the model used to compute an edge.
 * Implementations may wrap GeographicLib, SOFA/ERFA, vehicle simulation state,
 * or a fixed calibration. Dependencies such as clocks and state objects should
 * be supplied to an implementation when it is constructed, allowing update()
 * to retain the required no-argument interface.
 */
class TransformProvider {
public:
    /** @brief Destroy a provider through its interface. */
    virtual ~TransformProvider() = default;

    /**
     * @brief Refresh the cached transform from the provider's dependencies.
     *
     * FrameGraph calls this in parent-before-child order. Implementations
     * should update their cached Transform but must not change frame topology.
     */
    virtual void update() = 0;

    /**
     * @brief Return the latest transform from the child into its parent.
     * @return A reference owned by the provider and valid at least until the
     *         next update or provider destruction.
     */
    [[nodiscard]] virtual const Transform& parentFromChild() const = 0;
};

/**
 * @brief Provider for a frame relationship that never changes.
 *
 * Typical uses include a sensor mounting offset, a launch-site tangent frame,
 * or any other rigid calibration that is constant for the simulation.
 */
class FixedTransformProvider final : public TransformProvider {
public:
    /**
     * @brief Store a fixed parent-from-child transform.
     * @param parent_from_child Transform returned for the provider's lifetime.
     */
    explicit FixedTransformProvider(Transform parent_from_child);

    /** @brief Perform no work because the transform is constant. */
    void update() override;

    /** @brief Return the fixed parent-from-child transform. */
    [[nodiscard]] const Transform& parentFromChild() const override;

private:
    Transform parent_from_child_;
};

/**
 * @brief A named node in a tree of Cartesian coordinate frames.
 *
 * Each non-root frame owns one TransformProvider describing the edge from
 * itself to its parent. A frame does not cache a transform to the root;
 * FrameGraph composes the required edges when a conversion is requested.
 *
 * Parent ownership is shared so a child cannot outlive its ancestry. Frames do
 * not own their children, avoiding an ownership cycle.
 */
class Frame {
public:
    /**
     * @brief Construct a root frame.
     * @param name Nonempty name used to identify the frame in a graph.
     * @throws std::invalid_argument if @p name is empty.
     *
     * A root has no parent and no transform provider. Its coordinate system is
     * the reference against which all descendants in that tree are connected.
     */
    explicit Frame(std::string name);

    /**
     * @brief Construct a child frame.
     * @param name Nonempty name used to identify the frame in a graph.
     * @param parent Parent node in the frame tree.
     * @param provider Owner of the transform from this frame into @p parent.
     * @throws std::invalid_argument if the name, parent, or provider is invalid.
     */
    Frame(std::string name,
          std::shared_ptr<Frame> parent,
          std::unique_ptr<TransformProvider> provider);

    /**
     * @brief Refresh this frame's parent transform.
     *
     * Delegates to the provider for a child and does nothing for a root. Prefer
     * FrameGraph::update() when multiple related frames must be synchronized.
     */
    void update();

    /** @brief Return the frame's immutable, nonempty name. */
    [[nodiscard]] const std::string& name() const noexcept;

    /** @brief Return true when this frame has no parent. */
    [[nodiscard]] bool isRoot() const noexcept;

    /**
     * @brief Return this frame's parent pointer.
     * @return A null shared pointer for a root frame.
     */
    [[nodiscard]] const std::shared_ptr<Frame>& parent() const noexcept;

    /**
     * @brief Return the transform that maps this frame into its parent.
     * @throws std::logic_error if called on a root frame.
     */
    [[nodiscard]] const Transform& parentFromThis() const;

private:
    std::string name_;
    std::shared_ptr<Frame> parent_;
    std::unique_ptr<TransformProvider> provider_;
};

/**
 * @brief Owns a collection of frame nodes and evaluates their relationships.
 *
 * The graph may contain multiple disconnected trees, but a transformation can
 * only be computed between frames sharing the same root. Every registered frame
 * is updated at most once per update() call, even if it is also an ancestor of
 * another registered frame.
 */
class FrameGraph {
public:
    /**
     * @brief Register a frame for graph-wide updates.
     * @param frame Non-null frame to register.
     * @throws std::invalid_argument for a null frame, duplicate object, or
     *         duplicate frame name.
     *
     * A frame's ancestors need not be registered separately; update() follows
     * parent links and updates them as needed.
     */
    void addFrame(std::shared_ptr<Frame> frame);

    /**
     * @brief Refresh every registered frame in parent-before-child order.
     * @throws std::logic_error if a cycle is detected.
     *
     * The ordering ensures that a child provider may safely depend on state
     * calculated by an ancestor during the same update pass.
     */
    void update();

    /**
     * @brief Compute the transform between two connected frames.
     * @param destination Frame in which the result will be expressed.
     * @param source Frame in which the input is expressed.
     * @return `destination_from_source`.
     * @throws std::invalid_argument if the frames have different roots.
     * @throws std::logic_error if a cycle is detected while walking a parent
     *         chain.
     *
     * If @f$T_{rd}@f$ maps destination into the common root and @f$T_{rs}@f$
     * maps source into that root, then
     *
     * @f[
     *     T_{ds}=T_{rd}^{-1}T_{rs}.
     * @f]
     */
    [[nodiscard]] Transform transform(const Frame& destination,
                                      const Frame& source) const;

private:
    std::vector<std::shared_ptr<Frame>> frames_;
};

}  // namespace rsim
