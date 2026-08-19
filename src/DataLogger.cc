#include "DataLogger.h"

#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace rsim {

DataLogger::DataLogger(double sample_period, bool include_initial_sample)
    : sample_period_(sample_period), include_initial_sample_(include_initial_sample) {
    if (!std::isfinite(sample_period) || sample_period < 0.0) {
        throw std::invalid_argument("sample period must be finite and nonnegative");
    }
}

void DataLogger::bind(const std::vector<Model*>& models) {
    channels_.clear();
    channel_indices_.clear();
    for (const Model* model : models) {
        for (const auto& variable : model->stateVariables()) {
            const std::string name =
                model->stateNamespace() + "." + variable.alias;
            if (channel_indices_.count(name)) {
                throw std::invalid_argument("duplicate state-variable name: " + name);
            }
            const StateValue initial_value = variable.read();
            StateSamples samples;
            if (std::holds_alternative<double>(initial_value)) {
                samples = std::vector<double>{};
            } else if (std::holds_alternative<std::int64_t>(initial_value)) {
                samples = std::vector<std::int64_t>{};
            } else {
                samples = std::vector<bool>{};
            }
            channel_indices_[name] = channels_.size();
            channels_.push_back(
                {name, variable.units, variable.read, std::move(samples)});
        }
    }
    bound_ = true;
    next_sample_time_ = 0.0;
}

void DataLogger::record(double time, bool force) {
    if (!bound_) throw std::logic_error("data logger is not bound to a simulator");
    constexpr double tolerance = 1e-12;
    if (!force && sample_period_ > 0.0 && time + tolerance < next_sample_time_) return;
    times_.push_back(time);
    for (Channel& channel : channels_) {
        const StateValue value = channel.read();
        if (value.index() != channel.samples.index()) {
            throw std::runtime_error(
                "state-variable scalar type changed during simulation: " +
                channel.name);
        }
        std::visit([&value](auto& samples) {
            using Sample = typename std::decay_t<decltype(samples)>::value_type;
            samples.push_back(std::get<Sample>(value));
        }, channel.samples);
    }
    if (sample_period_ > 0.0) {
        do { next_sample_time_ += sample_period_; }
        while (next_sample_time_ <= time + tolerance);
    }
}

void DataLogger::reserveSamples(std::size_t count) {
    times_.reserve(count);
    for (auto& channel : channels_) {
        std::visit([count](auto& samples) { samples.reserve(count); },
                   channel.samples);
    }
}
void DataLogger::clear() noexcept {
    times_.clear();
    for (auto& channel : channels_) {
        std::visit([](auto& samples) { samples.clear(); }, channel.samples);
    }
    next_sample_time_ = 0.0;
}
std::size_t DataLogger::sampleCount() const noexcept { return times_.size(); }
const std::vector<double>& DataLogger::times() const noexcept { return times_; }
std::vector<std::string> DataLogger::channelNames() const {
    std::vector<std::string> names;
    for (const auto& channel : channels_) names.push_back(channel.name);
    return names;
}
const std::vector<double>& DataLogger::channel(const std::string& name) const {
    return std::get<std::vector<double>>(channelValues(name));
}
const StateSamples& DataLogger::channelValues(const std::string& name) const {
    const auto found = channel_indices_.find(name);
    if (found == channel_indices_.end()) throw std::out_of_range("unknown logged channel: " + name);
    return channels_[found->second].samples;
}
std::vector<StateChannelMetadata> DataLogger::metadata() const {
    std::vector<StateChannelMetadata> result;
    result.reserve(channels_.size());
    for (const Channel& channel : channels_) {
        StateScalarType scalar_type = StateScalarType::Float64;
        if (std::holds_alternative<std::vector<std::int64_t>>(
                channel.samples)) {
            scalar_type = StateScalarType::Int64;
        } else if (std::holds_alternative<std::vector<bool>>(
                       channel.samples)) {
            scalar_type = StateScalarType::Boolean;
        }
        result.push_back({channel.name, channel.units, scalar_type});
    }
    return result;
}
double DataLogger::samplePeriod() const noexcept { return sample_period_; }
bool DataLogger::includeInitialSample() const noexcept { return include_initial_sample_; }

}  // namespace rsim
