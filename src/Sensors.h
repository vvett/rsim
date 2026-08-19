#pragma once

#include "Model.h"
#include "defs.h"

#include <string>

namespace rsim {
class Vehicle;

// Samples one published scalar and holds it until the sensor updates again.
class ScalarSensor final : public Model {
public:
    ScalarSensor(std::string name, const Model& source, std::string alias,
                 UpdateRate update_rate = UpdateRate::everyStep());
    void update() override;
    void initialize() override;
    [[nodiscard]] double value() const noexcept;
    [[nodiscard]] bool valid() const noexcept;
    [[nodiscard]] double sampleTime() const noexcept;
    [[nodiscard]] std::uint64_t sequence() const noexcept;
private:
    std::function<StateValue()> read_;
    double value_ = 0.0;
    bool valid_ = false;
    double sample_time_ = 0.0;
    std::int64_t sequence_ = 0;
};

// Reports body-axis non-gravitational force divided by current vehicle mass.
class SpecificForceSensor final : public Model {
public:
    SpecificForceSensor(
        std::string name,
        Vehicle& vehicle,
        UpdateRate update_rate = UpdateRate::everyStep());
    void update() override;
    void initialize() override;
    [[nodiscard]] const Vector3d& specificForceInBody() const noexcept;
    [[nodiscard]] bool valid() const noexcept;
    [[nodiscard]] double sampleTime() const noexcept;
    [[nodiscard]] std::uint64_t sequence() const noexcept;
private:
    Vehicle& vehicle_;
    Vector3d specific_force_in_body_ = Vector3d::Zero();
    bool valid_ = false;
    double sample_time_ = 0.0;
    std::int64_t sequence_ = 0;
};

}  // namespace rsim
