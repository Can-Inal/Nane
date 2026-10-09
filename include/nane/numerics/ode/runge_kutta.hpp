#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/butcher_table.hpp"
#include "nane/numerics/ode/details/runge_kutta_impl.hpp"
#include "nane/numerics/ode/ivp.hpp"

#include <Eigen/Core>

namespace nane
{
    /**
     * @ingroup ode
     *
     * @brief Solves a scalar initial value problem using a Runge-Kutta method.
     *
     * @tparam Derivative Type of the right-hand-side function.
     *
     * @param problem Scalar initial value problem.
     * @param time_grid Time discretization.
     * @param table Butcher table defining the Runge-Kutta method.
     *
     * @return Numerical solution at all time-grid points.
     */
    template <typename Derivative>
    [[nodiscard]] Eigen::VectorXd runge_kutta(const nane::ivp<Derivative, double>& problem, const nane::uniform_grid<1>& time_grid,
                                              const nane::butcher_table& table)
    {
        Eigen::VectorXd initial_value(1);
        initial_value[0] = problem.initial_value();

        // Adapt the scalar derivative to the vector interface used internally.
        const auto derivative = [&](double time, const Eigen::VectorXd& state)
        {
            Eigen::VectorXd value(1);
            value[0] = problem.derivative()(time, state[0]);

            return value;
        };

        const Eigen::MatrixXd solution = nane::detail::runge_kutta_impl(derivative, initial_value, time_grid, table);

        return solution.row(0).transpose();
    }

    /**
     * @ingroup ode
     *
     * @brief Solves a vector-valued initial value problem using a
     * Runge-Kutta method.
     *
     * @tparam Derivative Type of the right-hand-side system.
     *
     * @param problem Vector-valued initial value problem.
     * @param time_grid Time discretization.
     * @param table Butcher table defining the Runge-Kutta method.
     *
     * @return Matrix whose columns contain the numerical states.
     */
    template <typename Derivative>
    [[nodiscard]] Eigen::MatrixXd runge_kutta(const nane::ivp<Derivative, Eigen::VectorXd>& problem, const nane::uniform_grid<1>& time_grid,
                                              const nane::butcher_table& table)
    {
        return nane::detail::runge_kutta_impl(problem.derivative(), problem.initial_value(), time_grid, table);
    }

} // namespace nane
