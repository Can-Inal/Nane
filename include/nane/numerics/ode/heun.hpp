#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/runge_kutta.hpp"

#include <Eigen/Core>
#include <utility>

namespace nane
{
    /**
     * @ingroup ode
     *
     * @brief Solves a scalar initial value problem using Heun's method.
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
     * @tparam Function Type of the right-hand-side function.
     *
     * @param function Right-hand-side function @f$f(t,x)@f$.
     * @param initial_value Initial value.
     * @param time_grid Time discretization.
     *
     * @return Numerical solution at all time-grid points.
     */
    template <typename Function>
    [[nodiscard]] Eigen::VectorXd heun(Function&& function, double initial_value, const nane::uniform_grid<1>& time_grid)
    {
        // Butcher table for heun's method.
        // alpha(0) = 0, alpha(1) = 1.
        // Beta(2, 1) = 1, others 0.
        // gamma(1) = gamma(2) = 1/2.

        Eigen::VectorXd alpha(2);
        alpha << 0.0, 1.0;

        // clang-format off
        Eigen::MatrixXd beta(2, 2);
        beta << 0.0, 0.0,
                1.0, 0.0;
        // clang-format on

        Eigen::VectorXd gamma(2);
        gamma << 0.5, 0.5;

        return nane::runge_kutta(std::forward<Function>(function), initial_value, time_grid, alpha, beta, gamma);
    }

    /**
     * @ingroup ode
     *
     * @brief Solves a vector-valued initial value problem using Heun's method.
     *
     * Approximates
     *
     * @f[
     * \dot{\mathbf{x}}(t)
     * =
     * \mathbf{f}(t, \mathbf{x}(t))
     * @f]
     *
     * using the predictor
     *
     * @f[
     * \tilde{\mathbf{x}}_{n+1}
     * =
     * \mathbf{x}_n
     * +
     * \tau_n
     * \mathbf{f}(t_n, \mathbf{x}_n)
     * @f]
     *
     * and the corrected update
     *
     * @f[
     * \mathbf{x}_{n+1}
     * =
     * \mathbf{x}_n
     * +
     * \frac{\tau_n}{2}
     * \left(
     * \mathbf{f}(t_n, \mathbf{x}_n)
     * +
     * \mathbf{f}(t_{n+1}, \tilde{\mathbf{x}}_{n+1})
     * \right).
     * @f]
     *
     * Each column of the returned matrix contains the numerical state at
     * one time-grid point.
     *
     * @tparam Function Type of the right-hand-side function.
     *
     * @param function Right-hand-side function.
     * @param initial_value Initial state vector.
     * @param time_grid Time discretization.
     *
     * @return Matrix whose columns contain the numerical states.
     */
    template <typename Function>
    [[nodiscard]] Eigen::MatrixXd heun(Function&& function, const Eigen::VectorXd& initial_value, const nane::uniform_grid<1>& time_grid)
    {
        // Butcher table for heun's method.
        // alpha(0) = 0, alpha(1) = 1.
        // Beta(2, 1) = 1, others 0.
        // gamma(1) = gamma(2) = 1/2.

        Eigen::VectorXd alpha(2);
        alpha << 0.0, 1.0;

        // clang-format off
        Eigen::MatrixXd beta(2, 2);
        beta << 0.0, 0.0,
                1.0, 0.0;
        // clang-format on

        Eigen::VectorXd gamma(2);
        gamma << 0.5, 0.5;

        return nane::runge_kutta(std::forward<Function>(function), initial_value, time_grid, alpha, beta, gamma);
    }

} // namespace nane
