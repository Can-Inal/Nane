#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/runge_kutta.hpp"

#include <Eigen/Core>
#include <utility>

namespace nane
{
    /**
     * @brief Solves a scalar initial value problem using implicit Euler.
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
     * f(t_{n+1}, x_{n+1}).
     * @f]
     *
     * The resulting implicit equation is solved through the Runge-Kutta
     * formulation using fixed-point iteration for the stage equation.
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
    [[nodiscard]] Eigen::VectorXd implicit_euler(Function&& function, double initial_value, const nane::uniform_grid<1>& time_grid)
    {
        // x(n+1) = x(n) + tau*f(t(n+1), x(n+1))
        // x(n+1) = x(n) + tau*k
        // k = f(t(n+1), x(n+1))
        // k = f(t(n+1), x(n) + tau*k)

        // Butcher table for implicit euler.
        // alpha = [1], Beta = [1], gamma = [1]

        Eigen::VectorXd alpha(1);
        alpha << 1.0;

        Eigen::MatrixXd beta(1, 1);
        beta << 1.0;

        Eigen::VectorXd gamma(1);
        gamma << 1.0;

        return nane::runge_kutta(std::forward<Function>(function), initial_value, time_grid, alpha, beta, gamma);
    }

    /**
     * @brief Solves a vector-valued initial value problem using implicit Euler.
     *
     * Approximates
     *
     * @f[
     * \dot{\mathbf{x}}(t)
     * =
     * \mathbf{f}(t, \mathbf{x}(t))
     * @f]
     *
     * using
     *
     * @f[
     * \mathbf{x}_{n+1}
     * =
     * \mathbf{x}_n
     * +
     * \tau_n
     * \mathbf{f}(t_{n+1}, \mathbf{x}_{n+1}).
     * @f]
     *
     * The resulting implicit system is solved through the Runge-Kutta
     * formulation using fixed-point iteration for the stage equation.
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
    [[nodiscard]] Eigen::MatrixXd implicit_euler(Function&& function, const Eigen::VectorXd& initial_value, const nane::uniform_grid<1>& time_grid)
    {
        // x(n+1) = x(n) + tau*f(t(n+1), x(n+1))
        // x(n+1) = x(n) + tau*k
        // k = f(t(n+1), x(n+1))
        // k = f(t(n+1), x(n) + tau*k)

        // Butcher table for implicit euler.
        // alpha = [1], Beta = [1], gamma = [1]

        Eigen::VectorXd alpha(1);
        alpha << 1.0;

        Eigen::MatrixXd beta(1, 1);
        beta << 1.0;

        Eigen::VectorXd gamma(1);
        gamma << 1.0;

        return nane::runge_kutta(std::forward<Function>(function), initial_value, time_grid, alpha, beta, gamma);
    }

} // namespace nane
