#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/linalg/vector.hpp"

#include <Eigen/Core>
#include <cstddef>

namespace nane
{
    /**
     * @brief Solves a scalar initial value problem using explicit Euler.
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
     * @tparam Function Type of the right-hand-side function.
     *
     * @param function Right-hand-side function @f$f(t,x)@f$.
     * @param initial_value Initial value.
     * @param time_grid Time discretization.
     *
     * @return Numerical solution at all time-grid points.
     */
    template <typename Function>
    [[nodiscard]] Eigen::VectorXd explicit_euler(Function&& function, double initial_value, const uniform_grid<1>& time_grid)
    {
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        Eigen::VectorXd solution = nane::vector(count);
        solution[0] = initial_value;

        for (auto i = 0; i < (int)count - 1; ++i)
        {
            const double t = time_grid.axis(0)[i];
            solution[i + 1] = solution[i] + tau * function(t, solution[i]);
        }

        return solution;
    }

    /**
     * @brief Solves a vector-valued initial value problem using explicit Euler.
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
     * \mathbf{f}(t_n, \mathbf{x}_n).
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
    [[nodiscard]] Eigen::MatrixXd explicit_euler(Function&& function, const Eigen::VectorXd& initial_value, const uniform_grid<1>& time_grid)
    {
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        const auto dimension = initial_value.size();

        Eigen::MatrixXd solution(dimension, count);
        solution.col(0) = initial_value;

        for (auto i = 0; i < (int)count - 1; ++i)
        {
            const double t = time_grid.axis(0)[i];
            solution.col(i + 1) = solution.col(i) + tau * function(t, solution.col(i));
        }

        return solution;
    }
} // namespace nane
