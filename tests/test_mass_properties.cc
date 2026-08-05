#include "MassProperties.h"

#include <cassert>
#include <stdexcept>

int main() {
    const rsim::Body::Matrix3 inertia{{{1.0, 0.0, 0.0},
                                         {0.0, 2.0, 0.0},
                                         {0.0, 0.0, 3.0}}};
    rsim::Body body(10.0, inertia);
    assert(body.mass() == 10.0);
    assert(body.momentOfInertia() == inertia);

    body.setMass(12.0);
    assert(body.mass() == 12.0);

    bool rejected_invalid_mass = false;
    try {
        body.setMass(0.0);
    } catch (const std::invalid_argument&) {
        rejected_invalid_mass = true;
    }
    assert(rejected_invalid_mass);
}
