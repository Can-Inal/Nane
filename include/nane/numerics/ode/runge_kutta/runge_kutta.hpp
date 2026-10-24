#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/details/runge_kutta_impl.hpp"

#include <Eigen/Core>

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Solves a scalar initial value problem using a Runge-Kutta method.
     *
     * The Runge-Kutta method is defined by the supplied Butcher table.
     * The order of accuracy depends on the coefficients of the Butcher table.
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
    [[nodiscard]] Eigen::VectorXd runge_kutta(const nane::ode::ivp<Derivative, double>& problem, const nane::geometry::uniform_grid<1>& time_grid,
                                              const nane::ode::butcher_table& table)
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

        const Eigen::MatrixXd solution = nane::ode::detail::runge_kutta_impl(derivative, initial_value, time_grid, table);

        return solution.row(0).transpose();
    }

    /**
     * @ingroup ode
     *
     * @brief Solves a vector-valued initial value problem using a
     * Runge-Kutta method.
     *
     * The Runge-Kutta method is defined by the supplied Butcher table.
     * The order of accuracy depends on the coefficients of the Butcher table.
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
    [[nodiscard]] Eigen::MatrixXd runge_kutta(const nane::ode::ivp<Derivative, Eigen::VectorXd>& problem,
                                              const nane::geometry::uniform_grid<1>& time_grid, const nane::ode::butcher_table& table)
    {
        return nane::ode::detail::runge_kutta_impl(problem.derivative(), problem.initial_value(), time_grid, table);
    }

} // namespace nane::ode
