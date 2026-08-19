#pragma once

#include "Model.h"
#include "defs.h"

#include <cstddef>
#include <string>
#include <vector>

namespace rsim {

// Separates proper-force contributions from gravity for ideal accelerometers.
enum class ForceClassification { Applied, Gravity };

// Base class for one C++ force model and its most recently calculated load.
class ForceTorque : public Model {
public:
    explicit ForceTorque(std::string name);
    ~ForceTorque() override = default;

    // Return this model's current force expressed in vehicle body coordinates.
    [[nodiscard]] const Vector3d& force() const noexcept;

    // Return this model's torque about an arbitrary body-coordinate point.
    [[nodiscard]] Vector3d torqueAbout(
        const Vector3d& reference_point) const noexcept;

    // Return this model's force followed by torque about reference_point.
    [[nodiscard]] Vector6d forceTorqueAbout(
        const Vector3d& reference_point) const noexcept;

    // Return the point where this model's force acts in body coordinates.
    [[nodiscard]] const Vector3d& applicationPoint() const noexcept;

    // Identify whether this load contributes to measured specific force.
    [[nodiscard]] virtual ForceClassification classification() const noexcept;

protected:
    // Publish the load calculated by a concrete force model's update method.
    void setLoad(Vector3d force,
                 Vector3d torque,
                 Vector3d application_point = Vector3d::Zero()) noexcept;

private:
    Vector3d force_ = Vector3d::Zero();
    Vector3d torque_ = Vector3d::Zero();
    Vector3d application_point_ = Vector3d::Zero();
};

// Supplies a prescribed load that remains fixed until explicitly replaced.
class ConstantForceTorque final : public ForceTorque {
public:
    explicit ConstantForceTorque(std::string name = "constant-load");

    // Replace the prescribed force, intrinsic torque, and application point.
    void setLoad(Vector3d force,
                 Vector3d torque,
                 Vector3d application_point = Vector3d::Zero()) noexcept;

    // Keep the prescribed load unchanged during registry updates.
    void update() override;
};

// Updates an ordered collection of force models and exposes their summed load.
class ForceTorqueRegistry final : public Model {
public:
    explicit ForceTorqueRegistry(
        std::string name = "forces",
        UpdateRate update_rate = UpdateRate::everyStep());

    // Add a non-owning force-model pointer at the end of the update order.
    void add(ForceTorque& force_model);

    // Remove a force model, throwing if it is not currently registered.
    void remove(ForceTorque& force_model);

    // Return whether this exact force-model object has been registered.
    [[nodiscard]] bool contains(const ForceTorque& force_model) const noexcept;

    // Return the number of force models in the ordered registry.
    [[nodiscard]] std::size_t size() const noexcept;

    // Return a registered force model by its insertion-order index.
    [[nodiscard]] ForceTorque& at(std::size_t index);
    [[nodiscard]] const ForceTorque& at(std::size_t index) const;

    // Update every enabled force model in deterministic insertion order.
    void update() override;

    // Return the sum of all enabled models' current forces.
    [[nodiscard]] Vector3d getForce() const noexcept;

    // Return the sum of all enabled models' torques about reference_point.
    [[nodiscard]] Vector3d getTorque(
        const Vector3d& reference_point) const noexcept;

    // Return net force followed by net torque about reference_point.
    [[nodiscard]] Vector6d getForceTorque(
        const Vector3d& reference_point) const noexcept;

private:
    // Pointers are non-owning; Vehicle or the C++ caller controls lifetimes.
    std::vector<ForceTorque*> force_models_;
};

}  // namespace rsim
