#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <variant>
#include <vector>

namespace rsim {

// Describes an exact integer or reciprocal-integer update multiplier.
class UpdateRate {
public:
    // Update once during every global simulator step.
    [[nodiscard]] static UpdateRate everyStep();

    // Update count times during every global simulator step.
    [[nodiscard]] static UpdateRate timesPerStep(std::uint32_t count);

    // Update once after count global simulator steps have accumulated.
    [[nodiscard]] static UpdateRate everyNSteps(std::uint32_t count);

    // Return how many updates occur within one global step.
    [[nodiscard]] std::uint32_t updatesPerStep() const noexcept;

    // Return how many global steps are accumulated before one update.
    [[nodiscard]] std::uint32_t stepsPerUpdate() const noexcept;

    // Compare exact scheduling ratios.
    [[nodiscard]] bool operator==(const UpdateRate& other) const noexcept;
    [[nodiscard]] bool operator!=(const UpdateRate& other) const noexcept;

private:
    UpdateRate(std::uint32_t updates_per_step,
               std::uint32_t steps_per_update);

    std::uint32_t updates_per_step_;
    std::uint32_t steps_per_update_;
};

// Timing information installed by Simulator immediately before an update.
struct UpdateContext {
    double simulation_time = 0.0;
    double time_step = 0.0;
    std::uint64_t global_step = 0;
    std::uint32_t substep = 0;
};

class Simulator;
class ForceTorqueRegistry;
class Vehicle;

enum class ModelPhase {
    Environment, Sensors, FlightControl, Actuators, Plant,
    MassProperties, Forces, Kinematics, Truth
};

using StateValue = std::variant<double, std::int64_t, bool>;

// Describes one stable, read-only value published by a native model.
struct StateVariable {
    std::string alias;
    std::string units;
    std::function<StateValue()> read;
};

// Common interface and scheduling state for an independently updated model.
class Model {
public:
    Model(std::string name,
          UpdateRate update_rate = UpdateRate::everyStep());
    virtual ~Model() = default;

    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) = delete;
    Model& operator=(Model&&) = delete;

    // Recalculate this model using the update context supplied by Simulator.
    virtual void update() = 0;

    // Populate setup-time outputs without advancing dynamics or control state.
    virtual void initialize();

    // Return the model's human-readable identifier.
    [[nodiscard]] const std::string& name() const noexcept;

    // Return the validated machine-readable namespace used by state channels.
    [[nodiscard]] const std::string& stateNamespace() const noexcept;

    // Convert a display label into one valid state-path identifier segment.
    [[nodiscard]] static std::string makeStateIdentifier(
        const std::string& label);

    // Return or replace the model's exact update multiplier.
    [[nodiscard]] const UpdateRate& updateRate() const noexcept;
    void setUpdateRate(UpdateRate update_rate) noexcept;

    // Return or change whether Simulator should update this model.
    [[nodiscard]] bool enabled() const noexcept;
    void setEnabled(bool enabled) noexcept;

    [[nodiscard]] ModelPhase phase() const noexcept;
    void setPhase(ModelPhase phase) noexcept;

    // Return the beginning time of the current or most recent model update.
    [[nodiscard]] double simulationTime() const noexcept;

    // Return the time interval assigned to the current or most recent update.
    [[nodiscard]] double timeStep() const noexcept;

    // Return the global step associated with the current or most recent update.
    [[nodiscard]] std::uint64_t globalStep() const noexcept;

    // Return the substep index for models updated multiple times per global step.
    [[nodiscard]] std::uint32_t substep() const noexcept;

    // Return the native values this model makes available to logging and events.
    [[nodiscard]] const std::vector<StateVariable>& stateVariables() const noexcept;

protected:
    // Publish a stable scalar member without exposing write access to consumers.
    void addStateVariable(const double& value, std::string alias,
                          std::string units = {});
    void addStateVariable(const std::int64_t& value, std::string alias,
                          std::string units = {});
    void addStateVariable(const bool& value, std::string alias,
                          std::string units = {});

    // Publish a value calculated by native code when a consumer samples it.
    void addComputedStateVariable(std::function<StateValue()> read,
                                  std::string alias,
                                  std::string units = {});
    // React when Vehicle connects this model to its native runtime state.
    virtual void onVehicleAttached(Vehicle& vehicle);

    // Release any borrowed vehicle state before the vehicle is destroyed.
    virtual void onVehicleDetached() noexcept;

    // Return the vehicle connected during setup, or throw if none is connected.
    [[nodiscard]] Vehicle& attachedVehicle();
    [[nodiscard]] const Vehicle& attachedVehicle() const;

    // Return whether setup has connected this model to a vehicle.
    [[nodiscard]] bool isAttachedToVehicle() const noexcept;

private:
    friend class Simulator;
    friend class ForceTorqueRegistry;
    friend class Vehicle;

    // Install scheduler-owned timing immediately before calling update().
    void setUpdateContext(UpdateContext context) noexcept;

    // Return a counter changed whenever scheduling configuration changes.
    [[nodiscard]] std::uint64_t scheduleRevision() const noexcept;

    // Connect and disconnect models only through Vehicle's setup API.
    void attachToVehicle(Vehicle& vehicle);
    void detachFromVehicle() noexcept;

    std::string name_;
    std::string state_namespace_;
    UpdateRate update_rate_;
    bool enabled_ = true;
    ModelPhase phase_ = ModelPhase::Plant;
    std::uint64_t schedule_revision_ = 0U;
    UpdateContext update_context_;
    Vehicle* attached_vehicle_ = nullptr;
    std::vector<StateVariable> state_variables_;
};

}  // namespace rsim
