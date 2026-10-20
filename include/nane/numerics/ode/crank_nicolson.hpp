#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/butcher_table.hpp"
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/runge_kutta.hpp"

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Solves an initial value problem using the Crank-Nicolson method.
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
     * \frac{\tau_n}{2}
     * \left(
     * f(t_n, x_n)
     * +
     * f(t_{n+1}, x_{n+1})
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
    [[nodiscard]] auto crank_nicolson(const nane::ode::ivp<Derivative, InitialValue>& problem, const nane::geometry::uniform_grid<1>& time_grid)
    {
        // Butcher coefficients for Crank-Nicolson:
        //
        // alpha = [0, 1]^T
        //
        //        [0   0  ]
        // beta = [1/2 1/2]
        //
        // gamma = [1/2, 1/2]^T

        static const nane::ode::butcher_table table{
            {0.0, 1.0},
            {
                {0.0, 0.0},
                {0.5, 0.5},
            },
            {0.5, 0.5},
        };

        return nane::ode::runge_kutta(problem, time_grid, table);
    }

} // namespace nane::ode
