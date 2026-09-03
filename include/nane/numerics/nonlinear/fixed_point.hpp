#pragma once

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>

namespace nane
{
    /**
     * @ingroup nonlinear
     *
     * @brief Computes a scalar fixed point using fixed-point iteration.
     *
     * Solves the nonlinear equation
     *
     * @f[
     * x = F(x)
     * @f]
     *
     * by the iteration
     *
     * @f[
     * x^{(m+1)}
     * =
     * F\left(x^{(m)}\right).
     * @f]
     *
     * The iteration is considered converged when
     *
     * @f[
     * \left|
     * x^{(m+1)} - x^{(m)}
     * \right|
     * <
     * \varepsilon,
     * @f]
     *
     * where @f$\varepsilon@f$ is the specified tolerance.
     *
     * @tparam Function Type of the fixed-point mapping.
     *
     * @param function Fixed-point mapping @f$F(x)@f$.
     * @param initial_guess Initial approximation of the fixed point.
     * @param tolerance Convergence tolerance.
     * @param max_iteration Maximum number of iterations.
     *
     * @return Approximation of the fixed point.
     *
     * @throws std::runtime_error If the iteration does not converge within
     * the maximum number of iterations.
     */
    template <typename Function>
    [[nodiscard]] double fixed_point(Function&& function, double initial_guess, double tolerance = 1e-10, int max_iteration = 100)
    {
        auto current = initial_guess;

        for (auto i = 0; i < max_iteration; i++)
        {
            auto next = function(current);

            if (std::abs(next - current) < tolerance)
                return next;

            current = next;
        }

        throw std::runtime_error("fixed-point iteration did not converge.");
    }

    /**
     * @ingroup nonlinear
     *
     * @brief Computes a vector-valued fixed point using fixed-point iteration.
     *
     * Solves the nonlinear system
     *
     * @f[
     * \mathbf{x}
     * =
     * \mathbf{F}(\mathbf{x})
     * @f]
     *
     * by the iteration
     *
     * @f[
     * \mathbf{x}^{(m+1)}
     * =
     * \mathbf{F}\left(\mathbf{x}^{(m)}\right).
     * @f]
     *
     * The iteration is considered converged when
     *
     * @f[
     * \left\|
     * \mathbf{x}^{(m+1)}
     * -
     * \mathbf{x}^{(m)}
     * \right\|_2
     * <
     * \varepsilon,
     * @f]
     *
     * where @f$\varepsilon@f$ is the specified tolerance.
     *
     * @tparam Function Type of the fixed-point mapping.
     *
     * @param function Fixed-point mapping @f$\mathbf{F}(\mathbf{x})@f$.
     * @param initial_guess Initial approximation of the fixed point.
     * @param tolerance Convergence tolerance.
     * @param max_iterations Maximum number of iterations.
     *
     * @return Approximation of the fixed point.
     *
     * @throws std::runtime_error If the iteration does not converge within
     * the maximum number of iterations.
     */
    template <typename Function>
    [[nodiscard]] Eigen::VectorXd fixed_point(Function&& function, const Eigen::VectorXd& initial_guess, double tolerance = 1e-10,
                                              int max_iterations = 100)
    {
        Eigen::VectorXd current = initial_guess;

        for (auto i = 0; i < max_iterations; i++)
        {
            Eigen::VectorXd next = function(current);

            if ((next - current).norm() < tolerance)
                return next;

            current = next;
        }

        throw std::runtime_error("fixed-point iteration did not converge.");
    }

} // namespace nane
