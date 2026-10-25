#pragma once

#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"

#include <Eigen/Core>
#include <Eigen/LU>
#include <complex>

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Evaluates the stability function of a Runge-Kutta method.
     *
     * For a Runge-Kutta method defined by the Butcher coefficients
     * @f$\beta@f$ and @f$\gamma@f$, the stability function is given by
     *
     * @f[
     * R(z)
     * =
     * 1
     * +
     * z \gamma^T
     * (I - z\beta)^{-1}
     * \mathbf{1}.
     * @f]
     *
     * Here, @f$I@f$ denotes the identity matrix and @f$\mathbf{1}@f$
     * denotes the vector whose entries are all equal to one.
     *
     * The complex argument @f$z@f$ typically represents
     * @f$z = \tau \lambda@f$ when analyzing the linear test equation
     * @f$\dot{x} = \lambda x@f$.
     *
     * @param table Butcher table defining the Runge-Kutta method.
     * @param z Complex argument at which the stability function is evaluated.
     *
     * @return Value of the stability function @f$R(z)@f$.
     */
    [[nodiscard]] inline std::complex<double> stability_function(const butcher_table& table, std::complex<double> z)
    {
        const auto& gamma = table.gamma();
        const auto& beta = table.beta();

        const Eigen::Index stages = beta.rows();

        const Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(stages, stages);
        const Eigen::VectorXcd ones = Eigen::VectorXcd::Ones(stages);

        const Eigen::VectorXcd gamma_complex = gamma.cast<std::complex<double>>();
        const Eigen::MatrixXcd beta_complex = beta.cast<std::complex<double>>();

        const auto value = gamma_complex.transpose() * (identity - z * beta_complex).inverse() * ones;

        return 1.0 + z * value(0, 0);
    }

    /**
     * @ingroup ode
     *
     * @brief Checks whether a Runge-Kutta method is absolutely stable at a given point.
     *
     * A Runge-Kutta method is absolutely stable at @f$z@f$ if
     *
     * @f[
     * |R(z)| \leq 1,
     * @f]
     *
     * where @f$R(z)@f$ is the stability function of the method.
     *
     * @param table Butcher table defining the Runge-Kutta method.
     * @param z Complex point in the stability plane.
     *
     * @return True if @f$z@f$ lies inside the absolute stability region.
     */
    [[nodiscard]] inline bool is_stable(const butcher_table& table, std::complex<double> z)
    {
        return std::abs(stability_function(table, z)) <= 1.0;
    }

} // namespace nane::ode
