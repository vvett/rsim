#include "MassProperties.h"
#include "Simulator.h"

#include <cassert>
#include <cmath>
#include <stdexcept>

int main() {
    rsim::Matrix3d inertia;
    inertia << 1.0, 0.0, 0.0,
               0.0, 2.0, 0.0,
               0.0, 0.0, 3.0;
    rsim::RigidBody body(10.0, inertia);
    assert(body.mass() == 10.0);
    assert(body.momentOfInertia().isApprox(inertia));

    body.setMass(12.0);
    assert(body.mass() == 12.0);

    bool rejected_invalid_mass = false;
    try {
        body.setMass(0.0);
    } catch (const std::invalid_argument&) {
        rejected_invalid_mass = true;
    }
    assert(rejected_invalid_mass);

    constexpr double pi = 3.14159265358979323846;
    rsim::CylindricalPropellantTank tank(
        1.0, 2.0, 0.0, 1.0, 2.0 * pi, pi);
    body.addPropellantTank(tank);

    // At time zero the full tank extends from x=0 to x=2.
    body.update();
    assert(std::abs(body.mass() - (12.0 + 2.0 * pi)) < 1e-12);
    assert(std::abs(body.cg().x() - 2.0 * pi / body.mass()) < 1e-12);

    // Absolute simulator time removes half the propellant after one second.
    rsim::Simulator simulator(1.0);
    simulator.addModel(body);
    simulator.step();
    simulator.step();
    assert(std::abs(body.propellantTank(0).propellantMass() - pi) < 1e-12);
    assert(std::abs(body.propellantTank(0).centerOfMass().x() - 0.5) < 1e-12);
    assert(std::abs(body.mass() - (12.0 + pi)) < 1e-12);
}
