#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/butcher_table.hpp"
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/runge_kutta.hpp"

namespace nane
{
    /**
     * @ingroup ode
     *
     * @brief Solves an initial value problem using modified Euler.
     *
     * Approximates
     *
     * @f[
     * \dot{x}(t) = f(t, x(t))
     * @f]
     *
     * using the midpoint predictor
     *
     * @f[
     * \tilde{x}_{n+\frac{1}{2}}
     * =
     * x_n
     * +
     * \frac{\tau_n}{2} f(t_n, x_n)
     * @f]
     *
     * and the update
     *
     * @f[
     * x_{n+1}
     * =
     * x_n
     * +
     * \tau_n
     * f\left(
     * t_n + \frac{\tau_n}{2},
     * \tilde{x}_{n+\frac{1}{2}}
     * \right).
     * @f]
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
    [[nodiscard]] auto modified_euler(const nane::ivp<Derivative, InitialValue>& problem, const nane::uniform_grid<1>& time_grid)
    {
        // Butcher coefficients for modified Euler:
        //
        // alpha = [0, 1/2]^T
        //
        //        [0   0]
        // beta = [1/2 0]
        //
        // gamma = [0, 1]^T

        static const nane::butcher_table table{
            {0.0, 0.5},
            {
                {0.0, 0.0},
                {0.5, 0.0},
            },
            {0.0, 1.0},
        };

        return nane::runge_kutta(problem, time_grid, table);
    }

} // namespace nane
