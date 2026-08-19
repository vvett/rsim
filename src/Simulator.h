#pragma once

#include "Model.h"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace rsim {

class Vehicle;
class DataLogger;
class TruthEvent;
enum class VehicleState;

struct RunResult {
    std::string reason;
    std::string event_name;
    double final_time = 0.0;
    std::uint64_t step_count = 0U;
};

// Advances global time and calls registered models at their requested rates.
class Simulator {
public:
    using VehicleStateTimeSteps = std::map<VehicleState, double>;

    explicit Simulator(double global_time_step);

    // Register a non-owning model reference in deterministic update order.
    void addModel(Model& model);

    // Register a vehicle and its owned models in dependency order.
    void addVehicle(Vehicle& vehicle);

    // Replace one vehicle's state-specific timestep policy during setup.
    void setVehicleStateTimeSteps(
        Vehicle& vehicle, VehicleStateTimeSteps time_steps);

    // Remove one vehicle's policy so every state falls back to the base timestep.
    void clearVehicleStateTimeSteps(const Vehicle& vehicle) noexcept;

    // Advance the simulation by one global time step.
    void step();

    // Advance the simulation by a fixed number of global steps.
    void runSteps(std::uint64_t count);

    // Run entirely in native code until an event or a mandatory safety limit stops it.
    RunResult run(double maximum_time, std::uint64_t maximum_steps);

    // Attach one native logger and any number of native truth events before running.
    void addDataLogger(DataLogger& logger);
    void addTruthEvent(TruthEvent& event);

    [[nodiscard]] double time() const noexcept;

    // Return or replace the positive finite fallback timestep.
    [[nodiscard]] double globalTimeStep() const noexcept;
    void setGlobalTimeStep(double global_time_step);

    // Return the timestep selected for the most recently started step.
    [[nodiscard]] double lastTimeStep() const noexcept;

    [[nodiscard]] std::uint64_t stepNumber() const noexcept;

    [[nodiscard]] bool stopRequested() const noexcept;

private:
    // Keeps per-registration counters next to the model they schedule.
    struct ScheduledModel {
        Model* model;
        std::uint64_t observed_schedule_revision;
        std::uint32_t accumulated_steps = 0U;
        double accumulated_time = 0.0;
        double interval_start = 0.0;
    };

    // Associates a borrowed vehicle with a fully native state-to-timestep map.
    struct VehicleTimeStepPolicy {
        Vehicle* vehicle;
        VehicleStateTimeSteps time_steps;
    };

    // Select the smallest active vehicle request, or the base timestep.
    [[nodiscard]] double selectTimeStep() const noexcept;
    void initializeRuntime();
    void evaluateEvents();
    [[nodiscard]] std::vector<Model*> discoverModelGraph() const;

    double time_ = 0.0;
    double global_time_step_;
    double last_time_step_;
    std::uint64_t step_number_ = 0U;
    std::vector<ScheduledModel> models_;
    std::vector<VehicleTimeStepPolicy> vehicle_time_step_policies_;
    std::vector<DataLogger*> data_loggers_;
    std::vector<TruthEvent*> truth_events_;
    bool runtime_initialized_ = false;
    bool stop_requested_ = false;
    std::string stop_event_name_;
};

}  // namespace rsim
