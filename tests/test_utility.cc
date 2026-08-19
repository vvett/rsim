#include "utility/RegularGridInterpolator.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

bool near(double left, double right, double tolerance = 1e-12) {
    return std::abs(left - right) <= tolerance;
}

}  // namespace

int main() {
    // The one-dimensional specialization performs ordinary linear interpolation.
    const rsim::utility::RegularGridInterpolator<1> line({{{0.0, 1.0}}});
    assert(near(line.evaluate({0.0, 2.0}, {0.25}), 0.5));

    // The two-dimensional specialization reads the first axis as the fast axis.
    const rsim::utility::RegularGridInterpolator<2> surface(
        {{{0.0, 1.0}, {0.0, 2.0}}});
    const std::vector<double> surface_values{0.0, 1.0, 4.0, 5.0};
    assert(near(surface.evaluate(surface_values, {0.25, 0.5}), 1.25));

    // Singleton axes are constant dimensions and do not require duplicate points.
    const rsim::utility::RegularGridInterpolator<2> singleton(
        {{{3.0}, {0.0, 2.0}}});
    assert(near(singleton.evaluate({10.0, 14.0}, {3.0, 0.5}), 11.0));

    // Three and higher dimensions use mlinterp through compile-time expansion.
    const rsim::utility::RegularGridInterpolator<3> volume(
        {{{0.0, 1.0}, {0.0, 2.0}, {0.0, 4.0}}});
    std::vector<double> volume_values;
    for (double z : {0.0, 4.0}) {
        for (double y : {0.0, 2.0}) {
            for (double x : {0.0, 1.0}) {
                volume_values.push_back(x + 2.0 * y + 3.0 * z);
            }
        }
    }
    assert(near(volume.evaluate(volume_values, {0.25, 0.5, 1.0}), 4.25));

    bool outside_rejected = false;
    try {
        const double unused = line.evaluate({0.0, 2.0}, {2.0});
        (void)unused;
    } catch (const std::out_of_range&) {
        outside_rejected = true;
    }
    assert(outside_rejected);

    return 0;
}
