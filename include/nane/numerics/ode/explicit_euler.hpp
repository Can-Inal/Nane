#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/runge_kutta.hpp"

#include <Eigen/Core>

namespace nane
{
    /**
     * @ingroup ode
     *
     * @brief Solves an initial value problem using explicit Euler.
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
     * \tau_n f(t_n, x_n).
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
    [[nodiscard]] auto explicit_euler(const nane::ivp<Derivative, InitialValue>& problem, const nane::uniform_grid<1>& time_grid)
    {
        // Butcher coefficients for explicit Euler:
        //
        // alpha = [0]
        // beta  = [0]
        // gamma = [1]

        Eigen::VectorXd alpha(1);
        alpha << 0.0;

        Eigen::MatrixXd beta(1, 1);
        beta << 0.0;

        Eigen::VectorXd gamma(1);
        gamma << 1.0;

        return nane::runge_kutta(problem, time_grid, alpha, beta, gamma);
    }

} // namespace nane
