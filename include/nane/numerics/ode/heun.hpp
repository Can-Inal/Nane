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
     * @brief Solves an initial value problem using Heun's method.
     *
     * Approximates
     *
     * @f[
     * \dot{x}(t) = f(t, x(t))
     * @f]
     *
     * using the predictor
     *
     * @f[
     * \tilde{x}_{n+1}
     * =
     * x_n
     * +
     * \tau_n f(t_n, x_n)
     * @f]
     *
     * and the corrected update
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
     * f(t_{n+1}, \tilde{x}_{n+1})
     * \right).
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
    [[nodiscard]] auto heun(const nane::ivp<Derivative, InitialValue>& problem, const nane::uniform_grid<1>& time_grid)
    {
        // Butcher coefficients for Heun's method:
        //
        // alpha = [0, 1]^T
        //
        //        [0  0]
        // beta = [1  0]
        //
        // gamma = [1/2, 1/2]^T

        Eigen::VectorXd alpha(2);
        alpha << 0.0, 1.0;

        // clang-format off
        Eigen::MatrixXd beta(2, 2);
        beta << 0.0, 0.0,
                1.0, 0.0;
        // clang-format on

        Eigen::VectorXd gamma(2);
        gamma << 0.5, 0.5;

        return nane::runge_kutta(problem, time_grid, alpha, beta, gamma);
    }

} // namespace nane
