#pragma once

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
#include <stdexcept>
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

            // Creates a zero vector with the same dimension as the state.
            //
            // Fixed-size vectors already know their dimension at compile time,
            // while dynamic vectors require it at runtime.
            const auto make_state = [&]()
            {
                nane::vector<Dimension> state;

                if constexpr (Dimension == nane::dynamic)
                    state = nane::vector<>(dimension);

                state.fill(0.0);

                return state;
            };

            // Computes all Runge-Kutta stage vectors for one time step.
            //
            // The stage matrix contains one stage vector in each column.
            // For example, with dimension 3 and 4 stages we get a 3x4 matrix.
            const auto compute_stages = [&](double current_time, const nane::vector<Dimension>& current_state, double tau)
            {
                nane::matrix<Dimension, nane::dynamic> stage(dimension, stage_count);

                stage.fill(0.0);

                // Explicit Runge-Kutta methods allow the stages to be computed
                // one after another.
                if (table.is_explicit())
                {
                    // Going over the Runge-Kutta stages.
                    for (std::size_t j = 0; j < stage_count; ++j)
                    {
                        // Weighted sum of the previous stage vectors.
                        //
                        // stage_sum = sum_l beta_jl * k_l
                        auto stage_sum = make_state();

                        // The first stage has no previous stages.
                        //
                        // After that we use all already computed stages.
                        for (std::size_t l = 0; l < j; ++l)
                        {
                            for (std::size_t i = 0; i < dimension; ++i)
                                stage_sum[i] += beta(j, l) * stage(i, l);
                        }

                        // Computes
                        //
                        // k_j = f(
                        //     t_n + alpha_j * tau,
                        //     x_n + tau * sum_l beta_jl * k_l
                        // )
                        const auto stage_state = current_state + tau * stage_sum;

                        const auto stage_value = derivative(current_time + alpha[j] * tau, stage_state);

                        if (stage_value.size() != dimension)
                        {
                            throw std::invalid_argument("ODE derivative dimension must match state dimension");
                        }

                        // Store k_j in column j of the stage matrix.
                        for (std::size_t i = 0; i < dimension; ++i)
                            stage(i, j) = stage_value[i];
                    }
                }
                else
                {
                    // Implicit Runge-Kutta stages depend on each other.
                    //
                    // Therefore all stages must be solved simultaneously.
                    //
                    // The stage vectors are flattened into one vector:
                    //
                    // K = [k_0, k_1, ..., k_s].
                    const auto mapping = [&](const nane::vector<>& current_stage)
                    {
                        nane::vector<> next_stage(dimension * stage_count);

                        // Going over all Runge-Kutta stages.
                        for (std::size_t j = 0; j < stage_count; ++j)
                        {
                            // Weighted sum of all stage vectors.
                            //
                            // stage_sum = sum_l beta_jl * k_l
                            auto stage_sum = make_state();

                            // Implicit methods can depend on previous, current,
                            // and later stages, so every stage is included.
                            for (std::size_t l = 0; l < stage_count; ++l)
                            {
                                for (std::size_t i = 0; i < dimension; ++i)
                                {
                                    stage_sum[i] += beta(j, l) * current_stage[l * dimension + i];
                                }
                            }

                            // Computes
                            //
                            // k_j = f(
                            //     t_n + alpha_j * tau,
                            //     x_n + tau * sum_l beta_jl * k_l
                            // )
                            const auto stage_state = current_state + tau * stage_sum;

                            const auto stage_value = derivative(current_time + alpha[j] * tau, stage_state);

                            if (stage_value.size() != dimension)
                            {
                                throw std::invalid_argument("ODE derivative dimension must match state dimension");
                            }

                            // Write k_j back into the corresponding part of
                            // the flattened stage vector.
                            for (std::size_t i = 0; i < dimension; ++i)
                            {
                                next_stage[j * dimension + i] = stage_value[i];
                            }
                        }

                        return next_stage;
                    };

                    // The implicit stage equations have the form
                    //
                    // K = mapping(K).
                    //
                    // Newton's method solves equations of the form
                    //
                    // residual(K) = 0.
                    //
                    // Therefore
                    //
                    // residual(K) = K - mapping(K).
                    const auto residual = [&](const nane::vector<>& current_stage) -> nane::vector<>
                    { return current_stage - mapping(current_stage); };

                    // Use a zero vector as the initial approximation for all
                    // flattened stage vectors.
                    const nane::vector<> initial_stage(dimension * stage_count, 0.0);

                    const auto solved_stage = nane::nonlinear::newton(residual, initial_stage);

                    // Unpack the flattened stage vector back into the stage
                    // matrix. Every section becomes one column.
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

            // A prescribed point count means that the time grid is fixed and
            // uniformly spaced.
            if (time_grid.has_count(0))
            {
                const std::size_t count = *time_grid.count(0);

                const double tau = *time_grid.spacing(0);

                // Contains all prescribed time points.
                std::vector<double> time(count);

                for (std::size_t i = 0; i < count; ++i)
                    time[i] = time_grid.axis(0)[i];

                // Solution matrix.
                //
                // Each column contains the numerical state at one grid point.
                //
                // For example, with dimension 3 and 100 grid points we get
                // a 3x100 matrix.
                nane::matrix<Dimension, nane::dynamic> values(dimension, count);

                // The first column contains the initial state.
                for (std::size_t i = 0; i < dimension; ++i)
                    values(i, 0) = initial_value[i];

                auto current_state = initial_value;

                // Going over the fixed time steps.
                for (std::size_t index = 0; index + 1 < count; ++index)
                {
                    const double current_time = time[index];

                    // Compute all Runge-Kutta stages for the current step.
                    const auto stage = compute_stages(current_time, current_state, tau);

                    // Computes the Runge-Kutta update
                    //
                    // x_{n+1}
                    // =
                    // x_n + tau * sum_j gamma_j * k_j.
                    //
                    // stage * gamma computes the weighted sum of all
                    // stage vectors.
                    current_state = current_state + tau * (stage * gamma);

                    // Store the new state in the next solution column.
                    for (std::size_t i = 0; i < dimension; ++i)
                        values(i, index + 1) = current_state[i];
                }

                return {
                    time,
                    values,
                };
            }

            // Adaptive integration requires a second set of solution weights.
            //
            // Both approximations use the same Runge-Kutta stage vectors and
            // their difference provides an estimate of the local error.
            if (!table.has_embedded_gamma())
            {
                throw std::invalid_argument("adaptive Runge-Kutta integration requires embedded gamma coefficients");
            }

            const auto& embedded_gamma = table.embedded_gamma();

            // The convergence order is used when adjusting the next step size.
            const int method_order = nane::ode::order(table);

            if (method_order <= 0)
            {
                throw std::invalid_argument("adaptive Runge-Kutta integration requires a method with positive order");
            }

            const double start = time_grid.start(0);

            const double stop = time_grid.stop(0);

            // Adaptive integration does not know the number of accepted
            // time points beforehand, so values are collected dynamically.
            std::vector<double> time_values;
            std::vector<nane::vector<Dimension>> state_values;

            time_values.push_back(start);
            state_values.push_back(initial_value);

            // There is nothing to integrate if both interval bounds are equal.
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

            // Direction is positive for forward integration and negative for
            // backward integration.
            const double direction = stop > start ? 1.0 : -1.0;

            double current_time = start;

            auto current_state = initial_value;

            // Use one hundred equal parts of the interval as the initial
            // step-size guess.
            double tau = (stop - start) / 100.0;

            // Safety factor prevents the new step size from being chosen
            // directly at the estimated error boundary.
            constexpr double safety = 0.9;

            // Limits how strongly the step size may change after one step.
            constexpr double minimum_factor = 0.2;
            constexpr double maximum_factor = 5.0;

            // Continue until the end of the integration interval is reached.
            while (direction * (stop - current_time) > 0.0)
            {
                // Prevent the adaptive procedure from reducing the step size
                // until floating-point arithmetic can no longer advance time.
                const double minimum_step = std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(current_time));

                if (std::abs(tau) <= minimum_step)
                {
                    throw std::runtime_error("adaptive Runge-Kutta step size became too small");
                }

                // Shorten the final step if necessary so that the final
                // accepted time point is exactly the interval bound.
                if (direction * (current_time + tau - stop) > 0.0)
                    tau = stop - current_time;

                // Compute all Runge-Kutta stages once.
                //
                // Both embedded approximations use these same stages.
                const auto stage = compute_stages(current_time, current_state, tau);

                // Primary Runge-Kutta approximation
                //
                // x_{n+1}
                // =
                // x_n + tau * sum_j gamma_j * k_j.
                const auto next_state = current_state + tau * (stage * gamma);

                // Embedded Runge-Kutta approximation
                //
                // x_hat_{n+1}
                // =
                // x_n + tau * sum_j gamma_hat_j * k_j.
                const auto embedded_state = current_state + tau * (stage * embedded_gamma);

                // The difference between both approximations estimates the
                // local error of the current step.
                const double error = (next_state - embedded_state).norm();

                // Accept the step if the estimated error is small enough.
                if (error <= adaptive_tolerance)
                {
                    current_time += tau;
                    current_state = next_state;

                    time_values.push_back(current_time);
                    state_values.push_back(current_state);
                }

                // Choose how much the next step size should change.
                //
                // If the estimated error is zero we can use the maximum
                // allowed increase directly.
                double factor = maximum_factor;

                if (error > 0.0)
                {
                    // Computes
                    //
                    // factor
                    // =
                    // safety * (tolerance / error)^(1 / p),
                    //
                    // where p is the order of the primary Runge-Kutta method.
                    factor = safety * std::pow(adaptive_tolerance / error, 1.0 / static_cast<double>(method_order));

                    // Prevent excessively large changes in the step size.
                    factor = std::clamp(factor, minimum_factor, maximum_factor);
                }

                tau *= factor;
            }

            // Convert the dynamically collected state vectors into the matrix
            // used by the solution representation.
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
