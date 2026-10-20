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
     * @brief Solves an initial value problem using the classical fourth-order
     * Runge-Kutta method.
     *
     * Approximates
     *
     * @f[
     * \dot{x}(t) = f(t, x(t))
     * @f]
     *
     * using the four stages
     *
     * @f[
     * k_1 = f(t_n, x_n),
     * @f]
     *
     * @f[
     * k_2 =
     * f\left(
     * t_n + \frac{\tau_n}{2},
     * x_n + \frac{\tau_n}{2} k_1
     * \right),
     * @f]
     *
     * @f[
     * k_3 =
     * f\left(
     * t_n + \frac{\tau_n}{2},
     * x_n + \frac{\tau_n}{2} k_2
     * \right),
     * @f]
     *
     * @f[
     * k_4 =
     * f\left(
     * t_n + \tau_n,
     * x_n + \tau_n k_3
     * \right),
     * @f]
     *
     * followed by the update
     *
     * @f[
     * x_{n+1}
     * =
     * x_n
     * +
     * \frac{\tau_n}{6}
     * \left(
     * k_1 + 2k_2 + 2k_3 + k_4
     * \right).
     * @f]
     *
     * The method is fourth-order accurate in time. For a sufficiently smooth
     * solution, the global discretization error satisfies
     *
     * @f[
     * O(\tau^4).
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
    [[nodiscard]] auto rk4(const nane::ode::ivp<Derivative, InitialValue>& problem, const nane::geometry::uniform_grid<1>& time_grid)
    {
        // Butcher coefficients for classical RK4:
        //
        // alpha = [0, 1/2, 1/2, 1]^T
        //
        //        [0   0   0  0]
        //        [1/2 0   0  0]
        // beta = [0   1/2 0  0]
        //        [0   0   1  0]
        //
        // gamma = [1/6, 1/3, 1/3, 1/6]^T

        static const nane::ode::butcher_table table{
            {0.0, 0.5, 0.5, 1.0},
            {
                {0.0, 0.0, 0.0, 0.0},
                {0.5, 0.0, 0.0, 0.0},
                {0.0, 0.5, 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0},
            },
            {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0},
        };

        return nane::ode::runge_kutta(problem, time_grid, table);
    }

} // namespace nane::ode
