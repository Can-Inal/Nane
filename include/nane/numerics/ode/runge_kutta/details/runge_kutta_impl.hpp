#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"
#include "nane/numerics/nonlinear/newton.hpp"
#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/order.hpp"
#include "nane/numerics/ode/solution.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace nane::ode
{
    namespace detail
    {
        inline constexpr double adaptive_tolerance = 1e-6;

        /**
         * @brief Common Runge-Kutta implementation for vector-valued states.
         *
         * Scalar problems are converted to one-dimensional vector problems
         * by the public scalar overload.
         *
         * A grid with a prescribed point count uses fixed step sizes.
         * Otherwise, the step size is chosen adaptively using an embedded
         * Runge-Kutta method.
         *
         * @tparam Derivative Type of the right-hand-side function.
         * @tparam Dimension State dimension.
         *
         * @param derivative Right-hand-side function.
         * @param initial_value Initial state.
         * @param time_grid Time grid specification.
         * @param table Butcher table defining the Runge-Kutta method.
         *
         * @return Time points and numerical states.
         */
        template <typename Derivative, std::size_t Dimension>
        [[nodiscard]] nane::ode::solution<nane::matrix<Dimension, nane::dynamic>>
        runge_kutta_impl(const Derivative& derivative, const nane::vector<Dimension>& initial_value, const nane::geometry::grid<1>& time_grid,
                         const nane::ode::butcher_table& table)
        {
            const auto& alpha = table.alpha();
            const auto& beta = table.beta();
            const auto& gamma = table.gamma();

            const std::size_t stage_count = table.stages();
            const std::size_t dimension = initial_value.size();

            // Create a zero vector with the same dimension as the state.
            //
            // Fixed-size vectors already know their dimension, while dynamic
            // vectors need it at runtime.
            const auto make_state = [&]()
            {
                nane::vector<Dimension> state;

                if constexpr (Dimension == nane::dynamic)
                    state = nane::vector<>(dimension);

                state.fill(0.0);

                return state;
            };

            // Compute all Runge-Kutta stages for one time step.
            //
            // Each column of the stage matrix contains one stage vector.
            const auto compute_stages = [&](double current_time, const nane::vector<Dimension>& current_state, double tau)
            {
                nane::matrix<Dimension, nane::dynamic> stage(dimension, stage_count);

                stage.fill(0.0);

                // Explicit methods compute the stages one after another.
                if (table.is_explicit())
                {
                    for (std::size_t j = 0; j < stage_count; ++j)
                    {
                        // Weighted sum of the previously computed stages:
                        //
                        // stage_sum = sum_l beta_jl * k_l
                        auto stage_sum = make_state();

                        for (std::size_t l = 0; l < j; ++l)
                        {
                            for (std::size_t i = 0; i < dimension; ++i)
                                stage_sum[i] += beta(j, l) * stage(i, l);
                        }

                        // k_j = f(
                        //     t_n + alpha_j * tau,
                        //     x_n + tau * sum_l beta_jl * k_l
                        // )
                        const auto stage_state = current_state + tau * stage_sum;

                        const auto stage_value = derivative(current_time + alpha[j] * tau, stage_state);

                        NANE_REQUIRE(stage_value.size() == dimension, "ODE derivative dimension must match state dimension");

                        // Store k_j in column j.
                        for (std::size_t i = 0; i < dimension; ++i)
                            stage(i, j) = stage_value[i];
                    }
                }
                else
                {
                    // Implicit stages depend on each other, so all stages are
                    // solved simultaneously.
                    //
                    // The stage vectors are flattened into one vector:
                    //
                    // K = [k_0, k_1, ..., k_s].
                    const auto mapping = [&](const nane::vector<>& current_stage)
                    {
                        nane::vector<> next_stage(dimension * stage_count);

                        for (std::size_t j = 0; j < stage_count; ++j)
                        {
                            // Weighted sum of all stage vectors:
                            //
                            // stage_sum = sum_l beta_jl * k_l
                            auto stage_sum = make_state();

                            for (std::size_t l = 0; l < stage_count; ++l)
                            {
                                for (std::size_t i = 0; i < dimension; ++i)
                                {
                                    stage_sum[i] += beta(j, l) * current_stage[l * dimension + i];
                                }
                            }

                            // k_j = f(
                            //     t_n + alpha_j * tau,
                            //     x_n + tau * sum_l beta_jl * k_l
                            // )
                            const auto stage_state = current_state + tau * stage_sum;

                            const auto stage_value = derivative(current_time + alpha[j] * tau, stage_state);

                            NANE_REQUIRE(stage_value.size() == dimension, "ODE derivative dimension must match state dimension");

                            for (std::size_t i = 0; i < dimension; ++i)
                            {
                                next_stage[j * dimension + i] = stage_value[i];
                            }
                        }

                        return next_stage;
                    };

                    // The implicit stage equations are
                    //
                    // K = mapping(K),
                    //
                    // so Newton's method is applied to
                    //
                    // residual(K) = K - mapping(K) = 0.
                    const auto residual = [&](const nane::vector<>& current_stage) -> nane::vector<>
                    { return current_stage - mapping(current_stage); };

                    // Start Newton's method with all stages equal to zero.
                    const nane::vector<> initial_stage(dimension * stage_count, 0.0);

                    const auto solved_stage = nane::nonlinear::newton(residual, initial_stage);

                    // Unpack the flattened stages back into matrix columns.
                    for (std::size_t j = 0; j < stage_count; ++j)
                    {
                        for (std::size_t i = 0; i < dimension; ++i)
                        {
                            stage(i, j) = solved_stage[j * dimension + i];
                        }
                    }
                }

                return stage;
            };

            // A prescribed point count means a fixed, uniformly spaced grid.
            if (time_grid.has_count(0))
            {
                const std::size_t count = *time_grid.count(0);

                const double tau = *time_grid.spacing(0);

                std::vector<double> time(count);

                for (std::size_t i = 0; i < count; ++i)
                    time[i] = time_grid.axis(0)[i];

                // Each column contains the state at one grid point.
                nane::matrix<Dimension, nane::dynamic> values(dimension, count);

                // The first column contains the initial state.
                for (std::size_t i = 0; i < dimension; ++i)
                    values(i, 0) = initial_value[i];

                auto current_state = initial_value;

                for (std::size_t index = 0; index + 1 < count; ++index)
                {
                    const double current_time = time[index];

                    const auto stage = compute_stages(current_time, current_state, tau);

                    // Runge-Kutta update:
                    //
                    // x_{n+1}
                    // =
                    // x_n + tau * sum_j gamma_j * k_j.
                    //
                    // stage * gamma is the weighted sum of all stages.
                    current_state = current_state + tau * (stage * gamma);

                    for (std::size_t i = 0; i < dimension; ++i)
                        values(i, index + 1) = current_state[i];
                }

                return {
                    time,
                    values,
                };
            }

            // Adaptive integration needs a second set of weights. Both
            // approximations share the same stages, and their difference
            // estimates the local error.
            NANE_REQUIRE(table.has_embedded_gamma(), "adaptive Runge-Kutta integration requires embedded gamma coefficients");

            const auto& embedded_gamma = table.embedded_gamma();

            // The method order determines how the next step size is chosen.
            const int method_order = nane::ode::order(table);

            NANE_REQUIRE(method_order > 0, "adaptive Runge-Kutta integration requires a method with positive order");

            const double start = time_grid.start(0);

            const double stop = time_grid.stop(0);

            // Adaptive integration does not know the number of accepted
            // points beforehand, so they are collected dynamically.
            std::vector<double> time_values;
            std::vector<nane::vector<Dimension>> state_values;

            time_values.push_back(start);
            state_values.push_back(initial_value);

            // Nothing to integrate if both interval bounds are equal.
            if (start == stop)
            {
                nane::matrix<Dimension, nane::dynamic> values(dimension, 1);

                for (std::size_t i = 0; i < dimension; ++i)
                    values(i, 0) = initial_value[i];

                return {
                    {start},
                    values,
                };
            }

            // Positive for forward integration, negative for backward.
            const double direction = stop > start ? 1.0 : -1.0;

            double current_time = start;

            auto current_state = initial_value;

            // Use one hundred equal parts of the interval as the initial guess.
            double tau = (stop - start) / 100.0;

            // Keep the next step slightly below the estimated error boundary.
            constexpr double safety = 0.9;

            // Limit how much the step size can change after one step.
            constexpr double minimum_factor = 0.2;
            constexpr double maximum_factor = 5.0;

            while (direction * (stop - current_time) > 0.0)
            {
                // Stop if floating-point arithmetic can no longer advance time.
                const double minimum_step = std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(current_time));

                NANE_FAIL(std::abs(tau) <= minimum_step, "adaptive Runge-Kutta step size became too small");

                // Shorten the final step so the last point is exactly stop.
                if (direction * (current_time + tau - stop) > 0.0)
                    tau = stop - current_time;

                // Both embedded approximations use the same stages.
                const auto stage = compute_stages(current_time, current_state, tau);

                // Primary approximation:
                //
                // x_{n+1}
                // =
                // x_n + tau * sum_j gamma_j * k_j.
                const auto next_state = current_state + tau * (stage * gamma);

                // Embedded approximation:
                //
                // x_hat_{n+1}
                // =
                // x_n + tau * sum_j gamma_hat_j * k_j.
                const auto embedded_state = current_state + tau * (stage * embedded_gamma);

                // Difference between both approximations estimates local error.
                const double error = (next_state - embedded_state).norm();

                if (error <= adaptive_tolerance)
                {
                    current_time += tau;
                    current_state = next_state;

                    time_values.push_back(current_time);
                    state_values.push_back(current_state);
                }

                // Start with the maximum increase. If the error is nonzero,
                // replace it with the usual adaptive step-size factor.
                double factor = maximum_factor;

                if (error > 0.0)
                {
                    // factor
                    // =
                    // safety * (tolerance / error)^(1 / p)
                    //
                    // where p is the order of the primary method.
                    factor = safety * std::pow(adaptive_tolerance / error, 1.0 / static_cast<double>(method_order));

                    factor = std::clamp(factor, minimum_factor, maximum_factor);
                }

                tau *= factor;
            }

            // Convert the collected states into the matrix used by solution.
            nane::matrix<Dimension, nane::dynamic> values(dimension, state_values.size());

            for (std::size_t column = 0; column < state_values.size(); ++column)
            {
                for (std::size_t row = 0; row < dimension; ++row)
                    values(row, column) = state_values[column][row];
            }

            return {
                time_values,
                values,
            };
        }

    } // namespace detail
} // namespace nane::ode
