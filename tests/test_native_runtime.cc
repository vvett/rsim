#include "DataLogger.h"
#include "ForceTorque.h"
#include "Sensors.h"
#include "Simulator.h"
#include "TruthEvent.h"

#include <cassert>
#include <cmath>

namespace {
class TestTruth final : public rsim::Model {
public:
    explicit TestTruth(std::string name = "truth values")
        : Model(std::move(name)) {
        addStateVariable(altitude, "truthAltitude", "m");
        addStateVariable(armed, "armed");
    }
    void update() override { altitude += 10.0 * timeStep(); armed = true; }
    double altitude = 0.0;
    bool armed = false;
};
}

int main() {
    TestTruth truth;
    truth.setPhase(rsim::ModelPhase::Truth);
    rsim::ScalarSensor sensor("altimeter", truth, "truthAltitude");
    rsim::DataLogger logger(0.1, true);
    rsim::TruthEvent stop("altitude limit",
        "truthAltitude greater_than 2.5 and armed equal_to true", true);
    rsim::Simulator simulator(0.1);
    simulator.addModel(truth);
    simulator.addModel(sensor);
    simulator.addTruthEvent(stop);
    simulator.addDataLogger(logger);
    const rsim::RunResult result = simulator.run(5.0, 100U);
    assert(result.reason == "truth_event");
    assert(result.event_name == "altitude limit");
    assert(stop.fired());
    assert(logger.sampleCount() == 4U);
    assert(logger.times().front() == 0.0);
    assert(std::abs(logger.channel("truth_values.truthAltitude").back() - 3.0) < 1e-12);
    assert(std::holds_alternative<std::vector<bool>>(
        logger.channelValues("truth_values.armed")));
    assert(std::holds_alternative<std::vector<std::int64_t>>(
        logger.channelValues("altimeter.sequence")));
    assert(sensor.sequence() == 4U);

    // Registry children are recursively discoverable by canonical alias.
    rsim::ConstantForceTorque force("test force");
    force.setLoad({12.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    rsim::ForceTorqueRegistry registry("force registry");
    registry.add(force);
    rsim::TruthEvent force_stop(
        "force observed", "test_force.forceX greater_than 10", true);
    rsim::DataLogger force_logger;
    rsim::Simulator force_simulator(0.1);
    force_simulator.addModel(registry);
    force_simulator.addTruthEvent(force_stop);
    force_simulator.addDataLogger(force_logger);
    const auto force_result = force_simulator.run(1.0, 10U);
    assert(force_result.event_name == "force observed");
    assert(force_logger.channel("test_force.forceX").front() == 12.0);

    // Display-name spaces are sanitized, and namespace collisions reject.
    TestTruth first_collision("same model");
    TestTruth second_collision("same-model");
    rsim::Simulator collision_simulator(0.1);
    collision_simulator.addModel(first_collision);
    collision_simulator.addModel(second_collision);
    bool collision_rejected = false;
    try {
        collision_simulator.run(0.2, 2U);
    } catch (const std::invalid_argument&) {
        collision_rejected = true;
    }
    assert(collision_rejected);
}
