#include "Integration.h"
#include "PressureController.h"
#include "Propulsion.h"
#include "Simulator.h"

#include <Eigen/Dense>

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
bool near(double a, double b, double tolerance = 1e-8) {
    return std::abs(a - b) <= tolerance;
}

double splitStep(double value, double step) {
    // y' = y - 10y: exact diagonal solve proves the IMEX utility is generic
    // and independent from propulsion or a particular nonlinear solver.
    return rsim::imexAdditiveStep(
        value, 0.0, step,
        [](double, double y) { return y; },
        [](double, double y) { return -10.0 * y; },
        [](double, double rhs, double coefficient, const auto&) {
            return rhs / (1.0 + 10.0 * coefficient);
        });
}
}  // namespace

int main() {
    const rsim::ImexTableau ars = rsim::ImexTableau::ars222();
    ars.validate();
    assert(ars.c.size() == 3U);
    assert(ars.A_implicit[1][1] > 0.0);
    assert(ars.A_implicit[2][2] > 0.0);

    double coarse = 1.0;
    for (int i = 0; i < 10; ++i) coarse = splitStep(coarse, 0.01);
    double fine = 1.0;
    for (int i = 0; i < 20; ++i) fine = splitStep(fine, 0.005);
    const double exact = std::exp(-0.9);
    assert(std::abs(fine - exact) < std::abs(coarse - exact));

    const std::vector<double> pc{1.0e4, 5.0e6};
    const std::vector<double> mr{0.5, 4.0};
    const std::vector<double> temperature(4, 3400.0);
    const std::vector<double> cstar(4, 1500.0);
    const std::vector<double> cf(4, 1.7);
    const std::vector<double> isp(4, 260.0);
    auto table = std::make_shared<rsim::CombustionTableData>(
        pc, mr, temperature, cstar, cf, isp, 5.0, "N2O", "ethanol");
    assert(near(table->cstar(2.0e6, 2.0), 1500.0));
    bool extrapolation_rejected = false;
    try {
        static_cast<void>(table->cstar(6.0e6, 2.0));
    } catch (const std::out_of_range&) {
        extrapolation_rejected = true;
    }
    assert(extrapolation_rejected);

    auto ambient = std::make_shared<rsim::Ambient>("vacuum", 0.0);
    auto oxidizer = std::make_shared<rsim::Tank>(
        "oxidizer", 0.2, 1.0, 0.0, 900.0, 80.0, 3.0e6);
    auto fuel = std::make_shared<rsim::Tank>(
        "fuel", 0.2, 1.0, 1.2, 800.0, 70.0, 3.0e6);
    auto engine = std::make_shared<rsim::Engine>(
        "engine", table, 1.0e-3, 5.0, Eigen::Vector3d::UnitX(),
        Eigen::Vector3d(0.0, 0.1, 0.0), ambient, 5.0e5, 1.0);
    auto network = std::make_shared<rsim::Propulsion>();
    network->connect(oxidizer, engine, std::make_shared<rsim::Valve>(1.0e-7),
                     rsim::FluidPhase::Liquid, rsim::InletRole::Oxidizer);
    network->connect(fuel, engine, std::make_shared<rsim::Orifice>(1.0e-7),
                     rsim::FluidPhase::Liquid, rsim::InletRole::Fuel);
    network->finalize();
    const double mass_before = network->totalStoredMass();
    rsim::Simulator simulator(0.01);
    simulator.addModel(*network);
    simulator.addModel(*engine);
    simulator.step();
    assert(network->totalStoredMass() < mass_before);
    assert(engine->thrust() > 0.0);
    assert(engine->force().x() > 0.0);
    assert(engine->torqueAbout(Eigen::Vector3d::Zero()).z() < 0.0);
    assert(network->lastSolveSummary().find("converged") != std::string::npos);

    // The last partial step consumes only the remaining limiting propellant.
    // On the following step, the missing inlet flames out the engine natively.
    fuel->setPropellantMass(1.0e-8);
    simulator.step();
    assert(near(fuel->propellantMass(), 0.0, 1e-14));
    simulator.step();
    assert(engine->pressure() < 10.0 * 6894.757293168);
    assert(near(engine->oxidizerMassFlow(), 0.0));
    assert(near(engine->fuelMassFlow(), 0.0));
    assert(near(engine->thrust(), 0.0));
    assert(engine->force().isZero());

    rsim::RigidBody body(10.0, Eigen::Matrix3d::Identity());
    body.addVariableMassSource(*oxidizer);
    bool duplicate_rejected = false;
    try {
        body.addVariableMassSource(*oxidizer);
    } catch (const std::invalid_argument&) {
        duplicate_rejected = true;
    }
    assert(duplicate_rejected);
    body.update();
    assert(body.mass() > 10.0);

    rsim::Propulsion skeletons;
    auto junction = std::make_shared<rsim::Junction>("junction");
    skeletons.connect(oxidizer, junction, std::make_shared<rsim::Pump>());
    bool pump_rejected = false;
    try {
        skeletons.finalize();
    } catch (const std::logic_error&) {
        pump_rejected = true;
    }
    assert(pump_rejected);

    // One gas source may conservatively feed any number of tank branches.
    auto nitrogen_source = std::make_shared<rsim::Tank>(
        "shared-N2", 0.2, 1.0, 0.0, 1000.0, 0.0, 5.0e6, 293.15,
        0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.2, "none", "Nitrogen");
    auto target_a = std::make_shared<rsim::Tank>(
        "target-a", 0.2, 1.0, 0.0, 1000.0, 60.0, 1.0e6, 293.15,
        0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0, "water", "Nitrogen");
    auto target_b = std::make_shared<rsim::Tank>(
        "target-b", 0.2, 1.0, 0.0, 1000.0, 60.0, 1.5e6, 293.15,
        0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0, "water", "Nitrogen");
    auto target_closed = std::make_shared<rsim::Tank>(
        "target-closed", 0.2, 1.0, 0.0, 1000.0, 60.0, 2.0e6, 293.15,
        0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0, "water", "Nitrogen");
    auto branches = std::make_shared<rsim::Propulsion>("three-way-N2");
    branches->connect(nitrogen_source, target_a,
                      std::make_shared<rsim::Orifice>(2.0e-7),
                      rsim::FluidPhase::Gas);
    branches->connect(nitrogen_source, target_b,
                      std::make_shared<rsim::Orifice>(7.0e-8),
                      rsim::FluidPhase::Gas);
    branches->connect(nitrogen_source, target_closed,
                      std::make_shared<rsim::Valve>(1.0e-6, 0.0),
                      rsim::FluidPhase::Gas);
    branches->finalize();
    const double branch_mass_before = branches->totalStoredMass();
    const double branch_energy_before = nitrogen_source->pressurantInternalEnergy() +
        target_a->pressurantInternalEnergy() +
        target_b->pressurantInternalEnergy() +
        target_closed->pressurantInternalEnergy();
    const double source_before = nitrogen_source->pressurantMass();
    const double a_before = target_a->pressurantMass();
    const double b_before = target_b->pressurantMass();
    const double closed_before = target_closed->pressurantMass();
    rsim::Simulator branch_simulator(0.01);
    branch_simulator.addModel(*branches);
    branch_simulator.step();
    const double flow_a = branches->connectionMassFlow(0U);
    const double flow_b = branches->connectionMassFlow(1U);
    assert(flow_a > flow_b && flow_b > 0.0);
    assert(near(branches->connectionMassFlow(2U), 0.0));
    assert(near(target_a->pressurantMass() - a_before, 0.01 * flow_a, 1e-12));
    assert(near(target_b->pressurantMass() - b_before, 0.01 * flow_b, 1e-12));
    assert(near(nitrogen_source->pressurantMass() - source_before,
                -0.01 * (flow_a + flow_b), 1e-12));
    assert(near(target_closed->pressurantMass(), closed_before, 1e-12));
    assert(near(branches->totalStoredMass(), branch_mass_before, 1e-12));
    const double branch_energy_after = nitrogen_source->pressurantInternalEnergy() +
        target_a->pressurantInternalEnergy() +
        target_b->pressurantInternalEnergy() +
        target_closed->pressurantInternalEnergy();
    assert(near(branch_energy_after, branch_energy_before, 1e-7));

    // A declared connection also conserves a compatible reverse flow.
    auto reverse_low = std::make_shared<rsim::Tank>(
        "reverse-low", 0.2, 1.0, 0.0, 1000.0, 60.0, 1.0e6, 293.15,
        0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0, "water", "Nitrogen");
    auto reverse_high = std::make_shared<rsim::Tank>(
        "reverse-high", 0.2, 1.0, 0.0, 1000.0, 60.0, 4.0e6, 293.15,
        0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0, "water", "Nitrogen");
    rsim::Propulsion reverse("reverse-N2");
    reverse.connect(reverse_low, reverse_high,
                    std::make_shared<rsim::Orifice>(1.0e-7),
                    rsim::FluidPhase::Gas);
    reverse.finalize();
    const double reverse_total = reverse.totalStoredMass();
    const double low_before = reverse_low->pressurantMass();
    rsim::Simulator reverse_simulator(0.01);
    reverse_simulator.addModel(reverse);
    reverse_simulator.step();
    assert(reverse.connectionMassFlow(0U) < 0.0);
    assert(reverse_low->pressurantMass() > low_before);
    assert(near(reverse.totalStoredMass(), reverse_total, 1e-12));

    bool species_rejected = false;
    try {
        auto helium = std::make_shared<rsim::Tank>(
            "helium", 0.2, 1.0, 0.0, 1000.0, 60.0, 1.0e6, 293.15,
            0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0,
            "water", "Helium");
        rsim::Propulsion incompatible;
        incompatible.connect(nitrogen_source, helium,
                             std::make_shared<rsim::Orifice>(1.0e-7),
                             rsim::FluidPhase::Gas);
    } catch (const std::invalid_argument&) {
        species_rejected = true;
    }
    assert(species_rejected);

    const auto pressure_tank = [](const std::string& name, double pressure) {
        return std::make_shared<rsim::Tank>(
            name, 0.2, 1.0, 0.0, 1000.0, 60.0, pressure, 293.15,
            0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0,
            "water", "Nitrogen");
    };
    auto low_pressure_tank = pressure_tank("controller-low", 0.8e6);
    auto low_valve = std::make_shared<rsim::Valve>(1.0e-7, 0.0);
    rsim::BangBangPressureController open_controller(
        "open-low", low_pressure_tank, low_valve, 1.0e6, 5.0e4, false);
    open_controller.update();
    assert(open_controller.valveOpen());
    assert(near(low_valve->opening(), 1.0));

    auto high_pressure_tank = pressure_tank("controller-high", 1.2e6);
    auto high_valve = std::make_shared<rsim::Valve>(1.0e-7, 1.0);
    rsim::BangBangPressureController close_controller(
        "close-high", high_pressure_tank, high_valve, 1.0e6, 5.0e4, true);
    close_controller.update();
    assert(!close_controller.valveOpen());
    assert(near(high_valve->opening(), 0.0));

    auto deadband_tank = pressure_tank("controller-deadband", 1.0e6);
    auto memory_closed_valve = std::make_shared<rsim::Valve>(1.0e-7, 1.0);
    auto memory_open_valve = std::make_shared<rsim::Valve>(1.0e-7, 0.0);
    rsim::BangBangPressureController memory_closed(
        "memory-closed", deadband_tank, memory_closed_valve,
        1.0e6, 5.0e4, false);
    rsim::BangBangPressureController memory_open(
        "memory-open", deadband_tank, memory_open_valve,
        1.0e6, 5.0e4, true);
    memory_closed.update();
    memory_open.update();
    assert(!memory_closed.valveOpen());
    assert(memory_open.valveOpen());

    auto controlled_source = std::make_shared<rsim::Tank>(
        "controlled-shared-N2", 0.2, 1.0, 0.0, 1000.0, 0.0, 5.0e6,
        293.15, 0.0, 1.0, rsim::BlowdownMode::Polytropic, 1.0,
        "none", "Nitrogen");
    auto controlled_low = pressure_tank("controlled-low", 0.8e6);
    auto controlled_high = pressure_tank("controlled-high", 1.2e6);
    auto controlled_low_valve = std::make_shared<rsim::Valve>(1.0e-7, 0.0);
    auto controlled_high_valve = std::make_shared<rsim::Valve>(1.0e-7, 1.0);
    auto low_controller = std::make_shared<rsim::BangBangPressureController>(
        "low-branch-controller", controlled_low, controlled_low_valve,
        1.0e6, 5.0e4, false);
    auto high_controller = std::make_shared<rsim::BangBangPressureController>(
        "high-branch-controller", controlled_high, controlled_high_valve,
        1.0e6, 5.0e4, true);
    auto controlled_network = std::make_shared<rsim::Propulsion>(
        "controlled-branches");
    controlled_network->connect(controlled_source, controlled_low,
                                controlled_low_valve, rsim::FluidPhase::Gas);
    controlled_network->connect(controlled_source, controlled_high,
                                controlled_high_valve, rsim::FluidPhase::Gas);
    controlled_network->finalize();
    const double controlled_source_before = controlled_source->pressurantMass();
    rsim::Simulator controlled_simulator(0.01);
    controlled_simulator.addModel(*low_controller);
    controlled_simulator.addModel(*high_controller);
    controlled_simulator.addModel(*controlled_network);
    controlled_simulator.step();
    assert(low_controller->valveOpen());
    assert(!high_controller->valveOpen());
    assert(controlled_network->connectionMassFlow(0U) > 0.0);
    assert(near(controlled_network->connectionMassFlow(1U), 0.0));
    assert(near(controlled_source_before - controlled_source->pressurantMass(),
                0.01 * controlled_network->connectionMassFlow(0U), 1e-12));

    return 0;
}
