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
     * @brief Solves an initial value problem using implicit Euler.
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
     * The resulting implicit stage equation is solved by the
     * Runge-Kutta implementation using fixed-point iteration.
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
    [[nodiscard]] auto implicit_euler(const nane::ivp<Derivative, InitialValue>& problem, const nane::uniform_grid<1>& time_grid)
    {
        // Implicit Euler satisfies
        //
        // x_(n+1) = x_n + tau * f(t_(n+1), x_(n+1)).
        //
        // Introducing one Runge-Kutta stage k gives
        //
        // k = f(t_n + tau, x_n + tau * k),
        //
        // followed by
        //
        // x_(n+1) = x_n + tau * k.
        //
        // Therefore its Butcher coefficients are:
        //
        // alpha = [1]
        // beta  = [1]
        // gamma = [1]

        Eigen::VectorXd alpha(1);
        alpha << 1.0;

        Eigen::MatrixXd beta(1, 1);
        beta << 1.0;

        Eigen::VectorXd gamma(1);
        gamma << 1.0;

        return nane::runge_kutta(problem, time_grid, alpha, beta, gamma);
    }

} // namespace nane
