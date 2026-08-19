#include "Model.h"

#include <stdexcept>
#include <utility>
#include <algorithm>
#include <cctype>

namespace rsim {
namespace {

std::string sanitizeIdentifier(const std::string& name) {
    std::string result;
    result.reserve(name.size());
    bool previous_was_separator = false;
    for (const unsigned char character : name) {
        if (std::isalnum(character) || character == '_') {
            result.push_back(static_cast<char>(character));
            previous_was_separator = false;
        } else if (!previous_was_separator) {
            result.push_back('_');
            previous_was_separator = true;
        }
    }
    while (!result.empty() && result.back() == '_') {
        result.pop_back();
    }
    if (result.empty() || std::isdigit(static_cast<unsigned char>(result.front()))) {
        result.insert(result.begin(), '_');
    }
    return result;
}

}  // namespace

UpdateRate::UpdateRate(std::uint32_t updates_per_step,
                       std::uint32_t steps_per_update)
    : updates_per_step_(updates_per_step),
      steps_per_update_(steps_per_update) {
    if (updates_per_step == 0U || steps_per_update == 0U) {
        throw std::invalid_argument("update-rate counts must be positive");
    }
    if (updates_per_step != 1U && steps_per_update != 1U) {
        throw std::invalid_argument(
            "update rate must be an integer or reciprocal integer");
    }
}

UpdateRate UpdateRate::everyStep() {
    return {1U, 1U};
}

UpdateRate UpdateRate::timesPerStep(std::uint32_t count) {
    return {count, 1U};
}

UpdateRate UpdateRate::everyNSteps(std::uint32_t count) {
    return {1U, count};
}

std::uint32_t UpdateRate::updatesPerStep() const noexcept {
    return updates_per_step_;
}

std::uint32_t UpdateRate::stepsPerUpdate() const noexcept {
    return steps_per_update_;
}

bool UpdateRate::operator==(const UpdateRate& other) const noexcept {
    return updates_per_step_ == other.updates_per_step_ &&
           steps_per_update_ == other.steps_per_update_;
}

bool UpdateRate::operator!=(const UpdateRate& other) const noexcept {
    return !(*this == other);
}

Model::Model(std::string name, UpdateRate update_rate)
    : name_(std::move(name)),
      state_namespace_(sanitizeIdentifier(name_)),
      update_rate_(update_rate) {
    if (name_.empty()) {
        throw std::invalid_argument("model name must not be empty");
    }
}

const std::string& Model::name() const noexcept {
    return name_;
}
void Model::initialize() {}

const std::string& Model::stateNamespace() const noexcept {
    return state_namespace_;
}
std::string Model::makeStateIdentifier(const std::string& label) {
    return sanitizeIdentifier(label);
}

const UpdateRate& Model::updateRate() const noexcept {
    return update_rate_;
}

void Model::setUpdateRate(UpdateRate update_rate) noexcept {
    if (update_rate_ != update_rate) {
        ++schedule_revision_;
    }
    update_rate_ = update_rate;
}

bool Model::enabled() const noexcept {
    return enabled_;
}

void Model::setEnabled(bool enabled) noexcept {
    if (enabled_ != enabled) {
        ++schedule_revision_;
    }
    enabled_ = enabled;
}
ModelPhase Model::phase() const noexcept { return phase_; }
void Model::setPhase(ModelPhase phase) noexcept { phase_ = phase; }

double Model::simulationTime() const noexcept {
    return update_context_.simulation_time;
}

double Model::timeStep() const noexcept {
    return update_context_.time_step;
}

std::uint64_t Model::globalStep() const noexcept {
    return update_context_.global_step;
}

std::uint32_t Model::substep() const noexcept {
    return update_context_.substep;
}

const std::vector<StateVariable>& Model::stateVariables() const noexcept {
    return state_variables_;
}

void Model::addStateVariable(const double& value, std::string alias,
                             std::string units) {
    addComputedStateVariable([&value]() -> StateValue { return value; },
                             std::move(alias), std::move(units));
}

void Model::addStateVariable(const std::int64_t& value, std::string alias,
                             std::string units) {
    addComputedStateVariable([&value]() -> StateValue { return value; },
                             std::move(alias), std::move(units));
}

void Model::addStateVariable(const bool& value, std::string alias,
                             std::string units) {
    addComputedStateVariable([&value]() -> StateValue { return value; },
                             std::move(alias), std::move(units));
}

void Model::addComputedStateVariable(std::function<StateValue()> read,
                                     std::string alias,
                                     std::string units) {
    if (alias.empty() || !read) {
        throw std::invalid_argument("state-variable alias and reader must be valid");
    }
    for (const char character : alias) {
        if (!std::isalnum(static_cast<unsigned char>(character)) &&
            character != '_' && character != '.') {
            throw std::invalid_argument(
                "state-variable alias must be machine-readable: " + alias);
        }
    }
    const auto duplicate = std::find_if(
        state_variables_.begin(), state_variables_.end(),
        [&alias](const StateVariable& variable) { return variable.alias == alias; });
    if (duplicate != state_variables_.end()) {
        throw std::invalid_argument("duplicate state-variable alias: " + alias);
    }
    state_variables_.push_back({std::move(alias), std::move(units), std::move(read)});
}

void Model::onVehicleAttached(Vehicle&) {}

void Model::onVehicleDetached() noexcept {}

Vehicle& Model::attachedVehicle() {
    if (!attached_vehicle_) {
        throw std::logic_error(
            "model is not attached to a vehicle: " + name());
    }
    return *attached_vehicle_;
}

const Vehicle& Model::attachedVehicle() const {
    if (!attached_vehicle_) {
        throw std::logic_error(
            "model is not attached to a vehicle: " + name());
    }
    return *attached_vehicle_;
}

bool Model::isAttachedToVehicle() const noexcept {
    return attached_vehicle_ != nullptr;
}

void Model::attachToVehicle(Vehicle& vehicle) {
    if (attached_vehicle_ && attached_vehicle_ != &vehicle) {
        throw std::logic_error(
            "model is already attached to another vehicle: " + name());
    }
    if (attached_vehicle_ == &vehicle) {
        return;
    }
    attached_vehicle_ = &vehicle;
    try {
        onVehicleAttached(vehicle);
    } catch (...) {
        attached_vehicle_ = nullptr;
        throw;
    }
}

void Model::detachFromVehicle() noexcept {
    if (!attached_vehicle_) {
        return;
    }
    onVehicleDetached();
    attached_vehicle_ = nullptr;
}

void Model::setUpdateContext(UpdateContext context) noexcept {
    update_context_ = context;
}

std::uint64_t Model::scheduleRevision() const noexcept {
    return schedule_revision_;
}

}  // namespace rsim
