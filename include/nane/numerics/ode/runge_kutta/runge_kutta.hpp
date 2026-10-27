#pragma once

#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/details/runge_kutta_impl.hpp"
#include "nane/numerics/ode/solution.hpp"

#include <cstddef>

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Solves a scalar initial value problem using a Runge-Kutta method.
     *
     * @tparam Derivative Type of the right-hand-side function.
     *
     * @param problem Scalar initial value problem.
     * @param table Butcher table defining the Runge-Kutta method.
     * @param time_grid Time grid specification.
     * @param tolerance Tolerance for adaptive ODE control.
     *
     * @return Time points and numerical solution values.
     */
    template <typename Derivative>
    [[nodiscard]] nane::ode::solution<nane::vector<>> runge_kutta(const nane::ode::ivp<Derivative, double>& problem,
                                                                  const nane::ode::butcher_table& table, const nane::geometry::grid<1>& time_grid,
                                                                  double tolerance = 1e-6)
    {
        const nane::vector<1> initial_value{
            problem.initial_value(),
        };

        const auto derivative = [&](double time, const nane::vector<1>& state)
        {
            return nane::vector<1>{
                problem.derivative()(time, state[0]),
            };
        };

        const auto result = nane::ode::detail::runge_kutta_impl(derivative, initial_value, time_grid, table, tolerance);

        nane::vector<> values(result.values.columns());

        for (std::size_t i = 0; i < result.values.columns(); ++i)
        {
            values[i] = result.values(0, i);
        }

        return {
            result.time,
            values,
        };
    }

    /**
     * @ingroup ode
     *
     * @brief Solves a vector-valued initial value problem using a
     * Runge-Kutta method.
     *
     * @tparam Derivative Type of the right-hand-side function.
     * @tparam Dimension State dimension.
     *
     * @param problem Vector-valued initial value problem.
     * @param table Butcher table defining the Runge-Kutta method.
     * @param time_grid Time grid specification.
     * @param tolerance Tolerance for adaptive ODE control.
     *
     * @return Time points and numerical solution states.
     */
    template <typename Derivative, std::size_t Dimension>
    [[nodiscard]] nane::ode::solution<nane::matrix<Dimension, nane::dynamic>>
    runge_kutta(const nane::ode::ivp<Derivative, nane::vector<Dimension>>& problem, const nane::ode::butcher_table& table,
                const nane::geometry::grid<1>& time_grid, double tolerance = 1e-6)
    {
        return nane::ode::detail::runge_kutta_impl(problem.derivative(), problem.initial_value(), time_grid, table, tolerance);
    }

} // namespace nane::ode
