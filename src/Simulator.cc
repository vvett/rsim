#include "Simulator.h"

#include "Vehicle.h"
#include "DataLogger.h"
#include "TruthEvent.h"
#include "FlightComputer.h"
#include "ForceTorque.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <unordered_map>
#include <unordered_set>

namespace rsim {

Simulator::Simulator(double global_time_step)
    : global_time_step_(global_time_step),
      last_time_step_(global_time_step) {
    setGlobalTimeStep(global_time_step);
}

void Simulator::addModel(Model& model) {
    const auto duplicate = std::find_if(
        models_.begin(), models_.end(),
        [&model](const ScheduledModel& scheduled) {
            return scheduled.model == &model;
        });
    if (duplicate != models_.end()) {
        throw std::invalid_argument("model is already registered: " + model.name());
    }
    if (runtime_initialized_) throw std::logic_error("models cannot be added after runtime initialization");
    models_.push_back(
        {&model, model.scheduleRevision(), 0U, 0.0, time_});
}

void Simulator::addVehicle(Vehicle& vehicle) {
    // Registration order is the dependency order used within each global step.
    addModel(vehicle);
    for (std::size_t index = 0; index < vehicle.modelCount(); ++index) {
        addModel(vehicle.modelAt(index));
    }
    addModel(vehicle.rigidBody());
    addModel(vehicle.forces());
    addModel(vehicle.kinematics());
}

void Simulator::setVehicleStateTimeSteps(
    Vehicle& vehicle, VehicleStateTimeSteps time_steps) {
    for (const auto& [state, time_step] : time_steps) {
        static_cast<void>(state);
        if (!std::isfinite(time_step) || time_step <= 0.0) {
            throw std::invalid_argument(
                "vehicle state time steps must be finite and positive");
        }
    }

    if (time_steps.empty()) {
        clearVehicleStateTimeSteps(vehicle);
        return;
    }

    const auto found = std::find_if(
        vehicle_time_step_policies_.begin(),
        vehicle_time_step_policies_.end(),
        [&vehicle](const VehicleTimeStepPolicy& policy) {
            return policy.vehicle == &vehicle;
        });
    if (found != vehicle_time_step_policies_.end()) {
        found->time_steps = std::move(time_steps);
        return;
    }
    vehicle_time_step_policies_.push_back(
        {&vehicle, std::move(time_steps)});
}

void Simulator::clearVehicleStateTimeSteps(
    const Vehicle& vehicle) noexcept {
    vehicle_time_step_policies_.erase(
        std::remove_if(
            vehicle_time_step_policies_.begin(),
            vehicle_time_step_policies_.end(),
            [&vehicle](const VehicleTimeStepPolicy& policy) {
                return policy.vehicle == &vehicle;
            }),
        vehicle_time_step_policies_.end());
}

double Simulator::selectTimeStep() const noexcept {
    if (vehicle_time_step_policies_.empty()) {
        return global_time_step_;
    }

    double selected = std::numeric_limits<double>::infinity();
    for (const VehicleTimeStepPolicy& policy : vehicle_time_step_policies_) {
        const auto requested = policy.time_steps.find(policy.vehicle->state());
        // An unmapped state requests the base timestep for that vehicle.
        const double vehicle_time_step =
            requested == policy.time_steps.end()
                ? global_time_step_
                : requested->second;
        selected = std::min(selected, vehicle_time_step);
    }
    return selected;
}

void Simulator::step() {
    initializeRuntime();
    // State changes during this loop deliberately affect the following step.
    const double current_time_step = selectTimeStep();
    last_time_step_ = current_time_step;

    constexpr ModelPhase phases[] = {
        ModelPhase::Environment, ModelPhase::Sensors, ModelPhase::FlightControl,
        ModelPhase::Actuators, ModelPhase::Plant, ModelPhase::MassProperties,
        ModelPhase::Forces, ModelPhase::Kinematics, ModelPhase::Truth};
    for (ModelPhase phase : phases) for (ScheduledModel& scheduled : models_) {
        Model& model = *scheduled.model;
        if (model.phase() != phase) continue;
        if (scheduled.observed_schedule_revision !=
            model.scheduleRevision()) {
            // Rate or enabled changes start a fresh interval and schedule phase.
            scheduled.observed_schedule_revision = model.scheduleRevision();
            scheduled.accumulated_steps = 0U;
            scheduled.accumulated_time = 0.0;
            scheduled.interval_start = time_;
        }
        if (!model.enabled()) {
            scheduled.accumulated_steps = 0U;
            scheduled.accumulated_time = 0.0;
            scheduled.interval_start = time_ + current_time_step;
            continue;
        }

        const UpdateRate& rate = model.updateRate();
        if (rate.updatesPerStep() > 1U) {
            const double model_dt =
                current_time_step / rate.updatesPerStep();
            for (std::uint32_t substep = 0U;
                 substep < rate.updatesPerStep();
                 ++substep) {
                model.setUpdateContext({
                    time_ + substep * model_dt,
                    model_dt,
                    step_number_,
                    substep
                });
                model.update();
            }
            continue;
        }

        if (scheduled.accumulated_steps == 0U) {
            scheduled.interval_start = time_;
        }
        ++scheduled.accumulated_steps;
        scheduled.accumulated_time += current_time_step;
        if (scheduled.accumulated_steps < rate.stepsPerUpdate()) {
            continue;
        }
        model.setUpdateContext({
            scheduled.interval_start,
            scheduled.accumulated_time,
            step_number_,
            0U
        });
        model.update();
        scheduled.accumulated_steps = 0U;
        scheduled.accumulated_time = 0.0;
        scheduled.interval_start = time_ + current_time_step;
    }

    time_ += current_time_step;
    ++step_number_;
    evaluateEvents();
    for (DataLogger* logger : data_loggers_) logger->record(time_);
}

void Simulator::runSteps(std::uint64_t count) {
    for (std::uint64_t index = 0U; index < count; ++index) {
        step();
    }
}

void Simulator::addDataLogger(DataLogger& logger) {
    if (runtime_initialized_) {
        throw std::logic_error(
            "data loggers cannot be added after runtime initialization");
    }
    if (std::find(data_loggers_.begin(), data_loggers_.end(), &logger) !=
        data_loggers_.end()) {
        throw std::invalid_argument("data logger is already registered");
    }
    data_loggers_.push_back(&logger);
}

void Simulator::addTruthEvent(TruthEvent& event) {
    if (runtime_initialized_) {
        throw std::logic_error(
            "truth events cannot be added after runtime initialization");
    }
    if (std::find(truth_events_.begin(), truth_events_.end(), &event) !=
        truth_events_.end()) {
        throw std::invalid_argument("truth event is already registered");
    }
    truth_events_.push_back(&event);
}

void Simulator::initializeRuntime() {
    if (runtime_initialized_) return;
    const std::vector<Model*> models = discoverModelGraph();
    // Initialize only observational phases; flight control and plant remain idle.
    constexpr ModelPhase initialization_phases[] = {
        ModelPhase::Environment, ModelPhase::Truth, ModelPhase::Sensors};
    for (ModelPhase phase : initialization_phases) {
        for (ScheduledModel& scheduled : models_) {
            if (scheduled.model->phase() != phase) {
                continue;
            }
            scheduled.model->setUpdateContext({time_, 0.0, step_number_, 0U});
            scheduled.model->initialize();
        }
    }
    std::unordered_map<const void*, std::string> actuator_owners;
    for (const Model* model : models) {
        const auto* computer = dynamic_cast<const FlightComputer*>(model);
        if (!computer) continue;
        for (const void* actuator : computer->actuatorIdentities()) {
            const auto [found, inserted] = actuator_owners.emplace(actuator, computer->name());
            if (!inserted) throw std::invalid_argument(
                "actuator is owned by multiple flight computers: " +
                found->second + " and " + computer->name());
        }
    }
    std::unordered_map<std::string, std::size_t> namespace_counts;
    for (const Model* model : models) {
        if (!model->stateVariables().empty()) {
            ++namespace_counts[model->stateNamespace()];
        }
    }
    for (const auto& [model_namespace, count] : namespace_counts) {
        if (count > 1U) {
            throw std::invalid_argument(
                "duplicate model state namespace: " + model_namespace);
        }
    }

    std::unordered_map<std::string, std::size_t> counts;
    for (const Model* model : models) {
        for (const auto& variable : model->stateVariables()) {
            ++counts[variable.alias];
        }
    }
    std::unordered_map<std::string, std::function<StateValue()>> values;
    values.emplace("simulationTime", [this]() -> StateValue { return time_; });
    values.emplace("timeStep", [this]() -> StateValue { return last_time_step_; });
    values.emplace("globalStep", [this]() -> StateValue {
        return static_cast<std::int64_t>(step_number_); });
    for (const Model* model : models) {
        for (const auto& variable : model->stateVariables()) {
            const std::string qualified =
                model->stateNamespace() + "." + variable.alias;
            if (!values.emplace(qualified, variable.read).second) {
                throw std::invalid_argument(
                    "duplicate qualified state-variable alias: " + qualified);
            }
            if (counts[variable.alias] == 1U) values.emplace(variable.alias, variable.read);
        }
    }
    for (TruthEvent* event : truth_events_) event->bind(values);
    for (DataLogger* logger : data_loggers_) {
        logger->bind(models);
        if (logger->includeInitialSample()) logger->record(time_, true);
    }
    runtime_initialized_ = true;
    evaluateEvents();
}

std::vector<Model*> Simulator::discoverModelGraph() const {
    std::vector<Model*> discovered;
    std::unordered_set<Model*> visited;
    std::function<void(Model*)> visit = [&](Model* model) {
        if (!model || !visited.insert(model).second) {
            return;
        }
        discovered.push_back(model);

        if (auto* vehicle = dynamic_cast<Vehicle*>(model)) {
            for (std::size_t index = 0; index < vehicle->modelCount(); ++index) {
                visit(&vehicle->modelAt(index));
            }
            visit(&vehicle->rigidBody());
            visit(&vehicle->forces());
            visit(&vehicle->kinematics());
        }
        if (auto* registry = dynamic_cast<ForceTorqueRegistry*>(model)) {
            for (std::size_t index = 0; index < registry->size(); ++index) {
                visit(&registry->at(index));
            }
        }
    };
    for (const ScheduledModel& scheduled : models_) {
        visit(scheduled.model);
    }
    return discovered;
}

void Simulator::evaluateEvents() {
    for (TruthEvent* event : truth_events_) {
        if (event->evaluate() && event->stop_simulation_) {
            stop_requested_ = true;
            stop_event_name_ = event->name();
        }
    }
}

RunResult Simulator::run(double maximum_time, std::uint64_t maximum_steps) {
    if ((!std::isfinite(maximum_time) || maximum_time <= time_) &&
        maximum_steps == 0U) {
        throw std::invalid_argument(
            "run requires a future maximum time or a positive maximum step count");
    }
    initializeRuntime();
    const std::uint64_t starting_step = step_number_;
    while (!stop_requested_) {
        if (maximum_steps > 0U &&
            step_number_ - starting_step >= maximum_steps) {
            return {"maximum_steps", {}, time_,
                    step_number_ - starting_step};
        }
        if (std::isfinite(maximum_time) && time_ >= maximum_time) {
            return {"maximum_time", {}, time_,
                    step_number_ - starting_step};
        }
        step();
    }
    return {"truth_event", stop_event_name_, time_, step_number_ - starting_step};
}

bool Simulator::stopRequested() const noexcept { return stop_requested_; }

double Simulator::time() const noexcept {
    return time_;
}

double Simulator::globalTimeStep() const noexcept {
    return global_time_step_;
}

void Simulator::setGlobalTimeStep(double global_time_step) {
    if (!std::isfinite(global_time_step) || global_time_step <= 0.0) {
        throw std::invalid_argument(
            "global time step must be finite and positive");
    }
    global_time_step_ = global_time_step;
    if (step_number_ == 0U) {
        last_time_step_ = global_time_step;
    }
}

double Simulator::lastTimeStep() const noexcept {
    return last_time_step_;
}

std::uint64_t Simulator::stepNumber() const noexcept {
    return step_number_;
}

}  // namespace rsim
