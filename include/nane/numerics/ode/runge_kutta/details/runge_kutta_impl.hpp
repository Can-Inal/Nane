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
         * @param tolerance Tolerance for adaptive ODE control.
         *
         * @return Time points and numerical states.
         */
        template <typename Derivative, std::size_t Dimension>
        [[nodiscard]] nane::ode::solution<nane::matrix<Dimension, nane::dynamic>>
        runge_kutta_impl(const Derivative& derivative, const nane::vector<Dimension>& initial_value, const nane::grid<1>& time_grid,
                         const nane::ode::butcher_table& table, double tolerance = 1e-6)
        {
            NANE_REQUIRE(tolerance > 0.0, "Runge-Kutta tolerance must be positive");

            const auto& alpha = table.alpha();
            const auto& beta = table.beta();
            const auto& gamma = table.gamma();

            const std::size_t stage_count = table.stages();
            const std::size_t dimension = initial_value.size();

            // Compute all Runge-Kutta stages for one time step.
            //
            // Each column of the stage matrix contains one stage vector.
            const auto compute_stages = [&](double current_time, const nane::vector<Dimension>& current_state, double tau)
            {
                nane::matrix<Dimension, nane::dynamic> stage(dimension, stage_count, 0.0);

                // Explicit methods compute the stages one after another.
                if (table.is_explicit())
                {
                    // k_j = f(t + alpha_j * tau, x[n] + tau * sum(beta_jl * k_l))
                    // For an explicit RK method, stage j depends only on previously computed stages l<j.
                    for (std::size_t j = 0; j < stage_count; ++j)
                    {
                        nane::vector<Dimension> stage_sum(dimension, 0.0);

                        for (std::size_t l = 0; l < j; ++l)
                        {
                            for (std::size_t i = 0; i < dimension; ++i)
                                stage_sum[i] += beta(j, l) * stage(i, l);
                        }

                        const auto stage_state = current_state + tau * stage_sum;
                        const auto stage_value = derivative(current_time + alpha[j] * tau, stage_state);

                        NANE_REQUIRE(stage_value.size() == dimension, "ODE derivative dimension must match state dimension");

                        // Store k_j in column j.
                        stage.set_column(j, stage_value);
                    }
                }
                else
                {
                    // Implicit stages depend on each other, so all stages are solved simultaneously.
                    const auto mapping = [&](const nane::vector<>& current_stage)
                    {
                        nane::matrix<Dimension, nane::dynamic> next_stage(dimension, stage_count);

                        for (std::size_t j = 0; j < stage_count; ++j)
                        {
                            nane::vector<Dimension> stage_sum(dimension, 0.0);

                            for (std::size_t l = 0; l < stage_count; ++l)
                            {
                                // sum(beta_jl * k_l)
                                for (std::size_t i = 0; i < dimension; ++i)
                                    stage_sum[i] += beta(j, l) * current_stage[l * dimension + i];
                            }

                            const auto stage_state = current_state + tau * stage_sum;
                            const auto stage_value = derivative(current_time + alpha[j] * tau, stage_state);

                            NANE_REQUIRE(stage_value.size() == dimension, "ODE derivative dimension must match state dimension");

                            next_stage.set_column(j, stage_value);
                        }

                        return next_stage.flatten_columns();
                    };

                    // clang-format off
                    const auto residual = [&](const nane::vector<>& current_stage) 
		    { 
			return current_stage - mapping(current_stage); 
		    };
                    // clang-format on

                    const nane::vector<> initial_stage(dimension * stage_count, 0.0);
                    const auto solved_stage = nane::nonlinear::newton(residual, initial_stage);

                    // Unpack the flattened stages back into matrix columns.
                    stage = nane::matrix<Dimension, nane::dynamic>::from_columns(solved_stage, stage_count);
                }

                return stage;
            };

            // Fixed step size section.
            // A prescribed point count means a fixed, uniformly spaced grid.
            if (time_grid.has_count(0))
            {
                const std::size_t count = *time_grid.count(0);
                const double tau = *time_grid.spacing(0);
                std::vector<double> time(time_grid.axis(0).begin(), time_grid.axis(0).end());

                // Each column contains the state at one grid point.
                nane::matrix<Dimension, nane::dynamic> values(dimension, count);
                values.set_column(0, initial_value);

                auto current_state = initial_value;

                for (std::size_t index = 0; index + 1 < count; ++index)
                {
                    const auto stage = compute_stages(time[index], current_state, tau);

                    // x[n+1] = x[n] + tau * sum_j gamma_j * k_j
                    current_state = current_state + tau * (stage * gamma);
                    values.set_column(index + 1, current_state);
                }

                return {time, values};
            }

            // Adaptive step size section.

            NANE_REQUIRE(table.has_embedded_gamma(), "adaptive Runge-Kutta integration requires embedded gamma coefficients");
            const auto& embedded_gamma = table.embedded_gamma();

            // The method order determines how the next step size is chosen.
            const int method_order = nane::ode::order(table);

            const double start = time_grid.start(0);
            const double stop = time_grid.stop(0);

            // Adaptive integration does not know the number of accepted
            // points beforehand, so they are collected dynamically.
            std::vector<double> time_values;
            std::vector<nane::vector<Dimension>> state_values;

            time_values.push_back(start);
            state_values.push_back(initial_value);

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

                // Compute states for both methods.
                const auto next_state = current_state + tau * (stage * gamma);
                const auto embedded_state = current_state + tau * (stage * embedded_gamma);

                // Difference between both approximations estimates local error.
                const double error = (next_state - embedded_state).norm();

                if (error <= tolerance)
                {
                    current_time += tau;
                    current_state = next_state;

                    time_values.push_back(current_time);
                    state_values.push_back(current_state);
                }

                // Start with the maximum increase.
                double factor = maximum_factor;

                // If the error is nonzero, replace it with the usual adaptive step-size factor.
                if (error > 0.0)
                {
                    factor = safety * std::pow(tolerance / error, 1.0 / static_cast<double>(method_order));
                    factor = std::clamp(factor, minimum_factor, maximum_factor);
                }

                tau *= factor;
            }

            // Convert the collected states into the matrix used by solution.
            const auto values = nane::matrix<Dimension, nane::dynamic>::from_columns(state_values);

            return {time_values, values};
        }

    } // namespace detail
} // namespace nane::ode
