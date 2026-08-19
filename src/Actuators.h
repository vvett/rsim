#pragma once

#include "Model.h"

#include <string>

namespace rsim {

class Valve;
class Parachute;

// Common command boundary used by native flight-computer plugins.
class IdealActuator : public Model {
public:
    explicit IdealActuator(std::string name);
    virtual void command(double value) = 0;
    [[nodiscard]] virtual double actualValue() const noexcept = 0;
};

// Applies a held opening command to one native propulsion valve.
class ValveActuator final : public IdealActuator {
public:
    ValveActuator(std::string name, Valve& valve);
    void command(double opening) override;
    void update() override;
    [[nodiscard]] double commanded() const noexcept;
    [[nodiscard]] double actual() const noexcept;
    [[nodiscard]] double actualValue() const noexcept override;
private:
    Valve& valve_; double commanded_ = 0.0; double actual_ = 0.0;
};

// Latches a deployment command and applies it to one parachute.
class ParachuteActuator final : public IdealActuator {
public:
    ParachuteActuator(std::string name, Parachute& parachute);
    void commandDeploy() noexcept;
    void command(double value) override;
    void update() override;
    [[nodiscard]] bool commanded() const noexcept;
    [[nodiscard]] bool actual() const noexcept;
    [[nodiscard]] double actualValue() const noexcept override;
private:
    Parachute& parachute_; bool commanded_ = false; bool actual_ = false;
};

// Applies a Boolean enable command to an arbitrary native model.
class ModelEnableActuator final : public IdealActuator {
public:
    ModelEnableActuator(std::string name, Model& target);
    void command(double value) override;
    void update() override;
    [[nodiscard]] bool commanded() const noexcept;
    [[nodiscard]] bool actual() const noexcept;
    [[nodiscard]] double actualValue() const noexcept override;
private:
    Model& target_;
    bool commanded_ = true;
    bool actual_ = true;
};

}  // namespace rsim
