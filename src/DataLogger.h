#pragma once

#include "Model.h"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace rsim {

class Simulator;
using StateSamples = std::variant<std::vector<double>,
                                  std::vector<std::int64_t>,
                                  std::vector<bool>>;
enum class StateScalarType { Float64, Int64, Boolean };
struct StateChannelMetadata {
    std::string name;
    std::string units;
    StateScalarType scalar_type;
};

// Stores registered native state variables without calling Python during a run.
class DataLogger {
public:
    explicit DataLogger(double sample_period = 0.0,
                        bool include_initial_sample = true);

    void reserveSamples(std::size_t count);
    void clear() noexcept;
    [[nodiscard]] std::size_t sampleCount() const noexcept;
    [[nodiscard]] const std::vector<double>& times() const noexcept;
    [[nodiscard]] std::vector<std::string> channelNames() const;
    [[nodiscard]] const std::vector<double>& channel(const std::string& name) const;
    [[nodiscard]] const StateSamples& channelValues(const std::string& name) const;
    [[nodiscard]] std::vector<StateChannelMetadata> metadata() const;
    [[nodiscard]] double samplePeriod() const noexcept;
    [[nodiscard]] bool includeInitialSample() const noexcept;

private:
    friend class Simulator;
    struct Channel {
        std::string name;
        std::string units;
        std::function<StateValue()> read;
        StateSamples samples;
    };

    void bind(const std::vector<Model*>& models);
    void record(double time, bool force = false);

    double sample_period_;
    bool include_initial_sample_;
    double next_sample_time_ = 0.0;
    bool bound_ = false;
    std::vector<double> times_;
    std::vector<Channel> channels_;
    std::unordered_map<std::string, std::size_t> channel_indices_;
};

}  // namespace rsim
