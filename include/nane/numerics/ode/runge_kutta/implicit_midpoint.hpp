#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/runge_kutta.hpp"

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Solves an initial value problem using the implicit midpoint method.
     *
     * Approximates
     *
     * @f[
     * \dot{x}(t) = f(t, x(t))
     * @f]
     *
     * using
     *
     * @f[
     * x_{n+1}
     * =
     * x_n
     * +
     * \tau_n
     * f\left(
     * t_n + \frac{\tau_n}{2},
     * \frac{x_n + x_{n+1}}{2}
     * \right).
     * @f]
     *
     * The resulting implicit stage equation is solved by the
     * Runge-Kutta implementation using Newton's method.
     *
     * The method is second-order accurate in time. For a sufficiently smooth
     * solution, the global discretization error satisfies
     *
     * @f[
     * O(\tau^2).
     * @f]
     *
     * @tparam Derivative Type of the right-hand-side function.
     * @tparam InitialValue Type of the initial value.
     *
     * @param problem Initial value problem containing the derivative and
     * initial value.
     * @param time_grid Time discretization.
     *
     * @return Numerical solution at all time-grid points.
     */
    template <typename Derivative, typename InitialValue>
    [[nodiscard]] auto implicit_midpoint(const nane::ode::ivp<Derivative, InitialValue>& problem, const nane::geometry::uniform_grid<1>& time_grid)
    {
        // Butcher coefficients for implicit midpoint:
        //
        // alpha = [1/2]
        // beta  = [1/2]
        // gamma = [1]

        static const nane::ode::butcher_table table{
            {0.5},
            {
                {0.5},
            },
            {1.0},
        };

        return nane::ode::runge_kutta(problem, time_grid, table);
    }

} // namespace nane::ode
