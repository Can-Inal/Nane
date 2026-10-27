#pragma once

#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"

#include <Eigen/Core>
#include <Eigen/LU>
#include <complex>
#include <cstddef>

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Evaluates the Runge-Kutta stability function.
     *
     * The stability function is
     *
     * @f[
     * R(z)
     * =
     * 1
     * +
     * z
     * \gamma^T
     * (I - z\beta)^{-1}
     * \mathbf{1}.
     * @f]
     *
     * @param table Butcher table defining the Runge-Kutta method.
     * @param z Complex stability parameter.
     *
     * @return Value of the stability function.
     */
    [[nodiscard]] inline std::complex<double> stability_function(const butcher_table& table, std::complex<double> z)
    {
        const auto& beta = table.beta();
        const auto& gamma = table.gamma();

        const std::size_t stage_count = table.stages();

        const auto stages = static_cast<Eigen::Index>(stage_count);

        Eigen::MatrixXcd beta_complex(stages, stages);

        Eigen::VectorXcd gamma_complex(stages);

        for (std::size_t i = 0; i < stage_count; ++i)
        {
            gamma_complex[static_cast<Eigen::Index>(i)] = gamma[i];

            for (std::size_t j = 0; j < stage_count; ++j)
            {
                beta_complex(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = beta(i, j);
            }
        }

        const Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(stages, stages);

        const Eigen::VectorXcd ones = Eigen::VectorXcd::Ones(stages);

        const Eigen::VectorXcd solution = (identity - z * beta_complex).partialPivLu().solve(ones);

        return 1.0 + z * gamma_complex.dot(solution);
    }

    /**
     * @ingroup ode
     *
     * @brief Determines whether a point lies inside the absolute
     * stability region.
     *
     * @param table Butcher table defining the Runge-Kutta method.
     * @param z Complex stability parameter.
     *
     * @return True if @f$|R(z)| \leq 1@f$.
     */
    [[nodiscard]] inline bool is_stable(const butcher_table& table, std::complex<double> z)
    {
        return std::abs(stability_function(table, z)) <= 1.0;
    }

} // namespace nane::ode
