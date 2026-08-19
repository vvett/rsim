#pragma once

#include "Model.h"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace rsim {

class Simulator;
class Parachute;
class IdealActuator;

// Evaluates a setup-time parsed expression against registered native values.
class TruthEvent {
public:
    struct Expression;
    TruthEvent(std::string name, std::string trigger,
               bool stop_simulation = false, bool one_shot = true);
    ~TruthEvent();
    TruthEvent(TruthEvent&&) noexcept;
    TruthEvent& operator=(TruthEvent&&) noexcept;
    TruthEvent(const TruthEvent&) = delete;
    TruthEvent& operator=(const TruthEvent&) = delete;

    [[nodiscard]] const std::string& name() const noexcept;
    [[nodiscard]] const std::string& trigger() const noexcept;
    [[nodiscard]] bool fired() const noexcept;
    [[nodiscard]] std::uint64_t fireCount() const noexcept;

    // Queue setup-time native actions; they execute in registration order.
    void addModelEnableCommand(Model& model, bool enabled);
    void addParachuteDeployCommand(Parachute& parachute);

    // Send a scalar command through the same isolated path as flight software.
    void addActuatorCommand(IdealActuator& actuator, double value);

private:
    friend class Simulator;
    void bind(const std::unordered_map<std::string, std::function<StateValue()>>& values);
    bool evaluate();

    std::string name_;
    std::string trigger_;
    bool stop_simulation_;
    bool one_shot_;
    bool fired_ = false;
    std::uint64_t fire_count_ = 0;
    std::unique_ptr<Expression> expression_;
    std::vector<std::function<void()>> commands_;
};

}  // namespace rsim
