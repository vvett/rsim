#pragma once

#include <mlinterp/mlinterp.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace rsim::utility {

// Interpolates scalar values on a regular N-dimensional rectilinear grid.
template <std::size_t Dimensions>
class RegularGridInterpolator {
public:
    static_assert(Dimensions > 0U, "an interpolator needs at least one dimension");

    using Axes = std::array<std::vector<double>, Dimensions>;
    using Point = std::array<double, Dimensions>;

    explicit RegularGridInterpolator(Axes axes)
        : axes_(std::move(axes)) {
        for (const auto& axis : axes_) {
            const bool finite = std::all_of(
                axis.begin(), axis.end(),
                [](double value) { return std::isfinite(value); });
            const bool increasing = std::adjacent_find(
                axis.begin(), axis.end(), std::greater_equal<>()) == axis.end();
            if (axis.empty() || !finite || !increasing) {
                throw std::invalid_argument(
                    "interpolation axes must be finite and strictly increasing");
            }
        }
    }

    // Interpolate one column-major value table at the supplied point.
    [[nodiscard]] double evaluate(
        const std::vector<double>& values,
        const Point& point) const {
        validate(values, point);
        if constexpr (Dimensions == 1U) {
            return evaluate1D(values, point[0]);
        } else if constexpr (Dimensions == 2U) {
            return evaluate2D(values, point);
        } else {
            return evaluateND(values, point,
                              std::make_index_sequence<Dimensions>{});
        }
    }

    // Return one axis for bounds, periodic wrapping, or diagnostics.
    [[nodiscard]] const std::vector<double>& axis(
        std::size_t dimension) const {
        return axes_.at(dimension);
    }

private:
    struct Interval {
        std::size_t lower;
        double fraction;
    };

    // Locate the lower point and fraction used by the fast 1-D/2-D paths.
    [[nodiscard]] static Interval interval(
        const std::vector<double>& axis, double value) {
        if (axis.size() == 1U) {
            return {0U, 0.0};
        }
        if (value == axis.back()) {
            return {axis.size() - 2U, 1.0};
        }
        const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
        const std::size_t lower =
            static_cast<std::size_t>(upper - axis.begin() - 1);
        return {lower,
                (value - axis[lower]) / (axis[lower + 1U] - axis[lower])};
    }

    // Dedicated linear interpolation avoids the general corner traversal.
    [[nodiscard]] double evaluate1D(
        const std::vector<double>& values, double point) const {
        const Interval location = interval(axes_[0], point);
        if (axes_[0].size() == 1U) {
            return values[0];
        }
        return values[location.lower] * (1.0 - location.fraction) +
               values[location.lower + 1U] * location.fraction;
    }

    // Dedicated bilinear interpolation performs only the required four blends.
    [[nodiscard]] double evaluate2D(
        const std::vector<double>& values, const Point& point) const {
        const Interval first = interval(axes_[0], point[0]);
        const Interval second = interval(axes_[1], point[1]);
        const std::size_t first_upper =
            axes_[0].size() == 1U ? first.lower : first.lower + 1U;
        const std::size_t second_upper =
            axes_[1].size() == 1U ? second.lower : second.lower + 1U;
        const auto value = [&](std::size_t i, std::size_t j) {
            return values[i + axes_[0].size() * j];
        };
        const double lower = value(first.lower, second.lower) *
                                 (1.0 - first.fraction) +
                             value(first_upper, second.lower) * first.fraction;
        const double upper = value(first.lower, second_upper) *
                                 (1.0 - first.fraction) +
                             value(first_upper, second_upper) * first.fraction;
        return lower * (1.0 - second.fraction) + upper * second.fraction;
    }

    // Expand axes and point addresses into mlinterp's compile-time N-D API.
    template <std::size_t... Indices>
    [[nodiscard]] double evaluateND(
        const std::vector<double>& values,
        const Point& point,
        std::index_sequence<Indices...>) const {
        const std::array<int, Dimensions> sizes{
            static_cast<int>(axes_[Indices].size())...};
        const auto interpolation_arguments = std::tuple_cat(
            std::make_tuple(axes_[Indices].data(), &point[Indices])...);
        double result = 0.0;
        // rnatord matches MATLAB column-major storage with the first axis fastest.
        std::apply(
            [&](const auto*... arguments) {
                mlinterp::interp<mlinterp::rnatord>(
                    sizes.data(), 1, values.data(), &result, arguments...);
            },
            interpolation_arguments);
        return result;
    }

    void validate(const std::vector<double>& values, const Point& point) const {
        const std::size_t expected = std::accumulate(
            axes_.begin(), axes_.end(), std::size_t{1},
            [](std::size_t product, const std::vector<double>& axis) {
                return product * axis.size();
            });
        if (values.size() != expected) {
            throw std::invalid_argument(
                "interpolation value count does not match the grid");
        }
        for (std::size_t dimension = 0U; dimension < Dimensions; ++dimension) {
            if (!std::isfinite(point[dimension]) ||
                point[dimension] < axes_[dimension].front() ||
                point[dimension] > axes_[dimension].back()) {
                throw std::out_of_range(
                    "interpolation point is outside dimension " +
                    std::to_string(dimension));
            }
        }
    }

    Axes axes_;
};

}  // namespace rsim::utility
