#pragma once

#include "nane/details/comparison.hpp"
#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"

#include <Eigen/Core>

namespace nane::ode
{

    /**
     * @ingroup ode
     *
     * @brief Determines the convergence order of a Runge-Kutta method.
     *
     * Determines the order of a Runge-Kutta method from its Butcher table
     * by checking the classical Runge-Kutta order conditions.
     *
     * Before checking the order conditions, the method is required to be
     * invariant against autonomization:
     *
     * @f[
     * \alpha_i = \sum_j \beta_{ij}.
     * @f]
     *
     * If this condition is not satisfied, the function returns zero.
     *
     * The currently implemented order conditions determine convergence
     * orders up to and including order four.
     *
     * The order conditions are checked successively. Therefore, a method
     * can only be classified as order @f$p@f$ if all conditions for the
     * lower orders are also satisfied.
     *
     * Floating-point comparisons are performed approximately.
     *
     * @param table Butcher table defining the Runge-Kutta method.
     *
     * @return Convergence order of the Runge-Kutta method.
     *
     * @note Orders higher than four are currently not detected.
     */
    [[nodiscard]] inline int order(const butcher_table& table)
    {
        int order = 0;

        const auto& alpha = table.alpha();
        const auto& gamma = table.gamma();
        const auto& beta = table.beta();

        // Invariance against autonomization
        {
            const Eigen::VectorXd row_sums = beta.rowwise().sum();

            if (!nane::detail::is_equal(alpha, row_sums))
                return 0;
        }

        // Convergence order of 1.
        {
            if (nane::detail::is_equal(gamma.sum(), 1.0))
                order = 1;
        }

        // Convergence order of 2.
        {
            // sum_i (gamma_i * alpha_i) = 1/2
            const auto condition1 = nane::detail::is_equal(gamma.dot(alpha), 0.5);

            if (order == 1 && condition1)
                order = 2;
        }

        // Convergence order of 3.
        {
            // sum_i (gamma_i * alpha_i^2) = 1/3
            const auto condition1 = nane::detail::is_equal(gamma.dot(alpha.array().square().matrix()), 1.0 / 3.0);

            // sum_i sum_j (gamma_i * beta_ij * alpha_j) = 1/6
            const auto condition2 = nane::detail::is_equal(gamma.dot(beta * alpha), 1.0 / 6.0);

            if (order == 2 && condition1 && condition2)
                order = 3;
        }

        // Convergence order of 4.
        {
            // sum_i (gamma_i * alpha_i^3) = 1/4
            const auto condition1 = nane::detail::is_equal(gamma.dot(alpha.array().cube().matrix()), 1.0 / 4.0);

            // sum_j sum_i (gamma_i * alpha_i * beta_ij * alpha_j) = 1/8
            const auto condition2 = nane::detail::is_equal((gamma.cwiseProduct(alpha)).dot(beta * alpha), 1.0 / 8.0);

            // sum_j sum_i (gamma_i * beta_ij * alpha_j^2) = 1/12
            const auto condition3 = nane::detail::is_equal(gamma.dot(beta * alpha.array().square().matrix()), 1.0 / 12.0);

            // sum_k sum_j sum_i (gamma_i * beta_ij * beta_jk * alpha_k) = 1/24
            const auto condition4 = nane::detail::is_equal(gamma.dot(beta * beta * alpha), 1.0 / 24.0);

            if (order == 3 && condition1 && condition2 && condition3 && condition4)
                order = 4;
        }

        return order;
    }
} // namespace nane::ode
