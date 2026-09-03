#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/linalg/vector.hpp"

#include <Eigen/Core>
#include <cstddef>

namespace nane
{
    /**
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
    [[nodiscard]] Eigen::VectorXd heun(Function&& function, double initial_value, const uniform_grid<1>& time_grid)
    {
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        Eigen::VectorXd solution = nane::vector(count);
        solution[0] = initial_value;

        for (auto i = 0; i < (int)count - 1; ++i)
        {
            const double current_time = time_grid.axis(0)[i];
            const double next_time = time_grid.axis(0)[i + 1];

            const auto f_n = function(current_time, solution[i]);

            const auto approximate_next = solution[i] + (tau * f_n);
            const auto f_next = function(next_time, approximate_next);

            solution[i + 1] = solution[i] + 0.5 * tau * (f_n + f_next);
        }

        return solution;
    }

    /**
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
    [[nodiscard]] Eigen::MatrixXd heun(Function&& function, const Eigen::VectorXd& initial_value, const uniform_grid<1>& time_grid)
    {
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        Eigen::MatrixXd solution(initial_value.size(), count);
        solution.col(0) = initial_value;

        for (auto i = 0; i < (int)count - 1; ++i)
        {
            const double current_time = time_grid.axis(0)[i];
            const double next_time = time_grid.axis(0)[i + 1];

            const Eigen::VectorXd current_state = solution.col(i);
            const Eigen::VectorXd f_n = function(current_time, current_state);

            const Eigen::VectorXd approximate_next = current_state + tau * f_n;
            const Eigen::VectorXd f_next = function(next_time, approximate_next);

            solution.col(i + 1) = current_state + 0.5 * tau * (f_n + f_next);
        }

        return solution;
    }

} // namespace nane
