#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rsim {

// Integrate one state step using RK4
template <typename State, typename DerivativeFunction>
State rungeKutta4(const State& initial_state,
                  double initial_time,
                  double time_step,
                  DerivativeFunction&& derivative) {
    // std::forward preserves whether the supplied callable was temporary or named.
    auto&& evaluate = std::forward<DerivativeFunction>(derivative);
    const State k1 = evaluate(initial_time, initial_state);
    const State k2 = evaluate(
        initial_time + 0.5 * time_step,
        initial_state + 0.5 * time_step * k1);
    const State k3 = evaluate(
        initial_time + 0.5 * time_step,
        initial_state + 0.5 * time_step * k2);
    const State k4 = evaluate(
        initial_time + time_step,
        initial_state + time_step * k3);
    return initial_state +
           (time_step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// Coefficients for an additive Runge--Kutta method y'=f_explicit+f_implicit.
// Matrices are stored row-major.  The utility intentionally has no dependency
// on a particular nonlinear solver: callers provide the diagonal-stage solve.
struct ImexTableau {
    std::vector<std::vector<double>> A_explicit;
    std::vector<std::vector<double>> A_implicit;
    std::vector<double> b_explicit;
    std::vector<double> b_implicit;
    std::vector<double> c;

    void validate() const {
        const std::size_t stages = c.size();
        if (stages == 0U || A_explicit.size() != stages ||
            A_implicit.size() != stages || b_explicit.size() != stages ||
            b_implicit.size() != stages) {
            throw std::invalid_argument("inconsistent IMEX tableau dimensions");
        }
        for (std::size_t i = 0; i < stages; ++i) {
            if (A_explicit[i].size() != stages ||
                A_implicit[i].size() != stages) {
                throw std::invalid_argument("IMEX matrices must be square");
            }
            for (std::size_t j = i; j < stages; ++j) {
                if (std::abs(A_explicit[i][j]) > 1e-15) {
                    throw std::invalid_argument(
                        "IMEX explicit matrix must be strictly lower triangular");
                }
            }
            for (std::size_t j = i + 1U; j < stages; ++j) {
                if (std::abs(A_implicit[i][j]) > 1e-15) {
                    throw std::invalid_argument(
                        "IMEX implicit matrix must be lower triangular");
                }
            }
        }
    }

    // Three-row representation of canonical ARS(2,2,2).  The zero first row
    // is the shared initial stage; rows one and two are the two DIRK solves.
    static ImexTableau ars222() {
        const double gamma = 1.0 - 1.0 / std::sqrt(2.0);
        const double delta = 1.0 - 1.0 / (2.0 * gamma);
        return {
            {{0.0, 0.0, 0.0},
             {gamma, 0.0, 0.0},
             {delta, 1.0 - delta, 0.0}},
            {{0.0, 0.0, 0.0},
             {0.0, gamma, 0.0},
             {0.0, 1.0 - gamma, gamma}},
            {delta, 1.0 - delta, 0.0},
            {0.0, 1.0 - gamma, gamma},
            {0.0, gamma, 1.0}};
    }
};

// Advance one additive RK step.  State must support addition and scalar
// multiplication. solve_implicit(time, rhs, h_aii, implicit_derivative)
// returns the stage state satisfying Y=rhs+h_aii*f_implicit(time,Y).
template <typename State, typename ExplicitFunction,
          typename ImplicitFunction, typename ImplicitStageSolver>
State imexAdditiveStep(const State& initial_state,
                       double initial_time,
                       double time_step,
                       ExplicitFunction&& explicit_derivative,
                       ImplicitFunction&& implicit_derivative,
                       ImplicitStageSolver&& solve_implicit,
                       const ImexTableau& tableau = ImexTableau::ars222()) {
    tableau.validate();
    const std::size_t count = tableau.c.size();
    std::vector<State> explicit_slopes;
    std::vector<State> implicit_slopes;
    explicit_slopes.reserve(count);
    implicit_slopes.reserve(count);

    for (std::size_t i = 0; i < count; ++i) {
        State rhs = initial_state;
        for (std::size_t j = 0; j < i; ++j) {
            rhs = rhs + time_step *
                (tableau.A_explicit[i][j] * explicit_slopes[j] +
                 tableau.A_implicit[i][j] * implicit_slopes[j]);
        }
        const double stage_time = initial_time + tableau.c[i] * time_step;
        State stage = rhs;
        const double diagonal = tableau.A_implicit[i][i];
        if (std::abs(diagonal) > 1e-15) {
            stage = solve_implicit(stage_time, rhs, time_step * diagonal,
                                   implicit_derivative);
        }
        explicit_slopes.push_back(explicit_derivative(stage_time, stage));
        implicit_slopes.push_back(implicit_derivative(stage_time, stage));
    }

    State result = initial_state;
    for (std::size_t i = 0; i < count; ++i) {
        result = result + time_step *
            (tableau.b_explicit[i] * explicit_slopes[i] +
             tableau.b_implicit[i] * implicit_slopes[i]);
    }
    return result;
}

}  // namespace rsim
