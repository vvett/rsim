#include "MassProperties.h"
#include <iostream>
#include <array>

int main() {

    std::array<std::array<double, 3>, 3> moi = {{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}};
    double m = 10;

    rsim::Body testBody(m, moi);

    testBody.printMass();

    return 0;
}
