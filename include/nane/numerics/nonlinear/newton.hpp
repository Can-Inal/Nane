#pragma once

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <stdexcept>

namespace nane
{
    /**
     * @ingroup nonlinear
     *
     * @brief Computes a scalar root using Newton's method.
     *
     * Solves the nonlinear equation
     *
     * @f[
     * f(x) = 0
     * @f]
     *
     * by the Newton iteration
     *
     * @f[
     * x^{(m+1)}
     * =
     * x^{(m)}
     * -
     * \frac{
     * f\left(x^{(m)}\right)
     * }{
     * f'\left(x^{(m)}\right)
     * }.
     * @f]
     *
     * The derivative is approximated by the central difference
     *
     * @f[
     * f'(x)
     * \approx
     * \frac{
     * f(x+h) - f(x-h)
     * }{
     * 2h
     * },
     * @f]
     *
     * where @f$h@f$ is a small finite-difference step.
     *
     * The iteration is considered converged when either
     *
     * @f[
     * \left|
     * f\left(x^{(m)}\right)
     * \right|
     * <
     * \varepsilon
     * @f]
     *
     * or
     *
     * @f[
     * \left|
     * x^{(m+1)}
     * -
     * x^{(m)}
     * \right|
     * <
     * \varepsilon,
     * @f]
     *
     * where @f$\varepsilon@f$ is the specified tolerance.
     *
     * @tparam Function Type of the nonlinear function.
     *
     * @param function Nonlinear function @f$f(x)@f$.
     * @param initial_guess Initial approximation of the root.
     * @param tolerance Convergence tolerance.
     * @param max_iteration Maximum number of iterations.
     *
     * @return Approximation of the root.
     *
     * @throws std::runtime_error If a near-zero derivative is encountered.
     * @throws std::runtime_error If the iteration does not converge within
     * the maximum number of iterations.
     */
    template <typename Function>
    [[nodiscard]] double newton(Function&& function, double initial_guess, double tolerance = 1e-10, int max_iteration = 100)
    {
        constexpr double step_size = 1e-6;
        auto current = initial_guess;

        for (int i = 0; i < max_iteration; i++)
        {
            // x[n]
            const auto value = function(current);

            if (std::abs(value) < tolerance)
                return current;

            // central difference: f'(n) = (f(n + h) - f(n - h)) / 2 * h
            const auto derivative = (function(current + step_size) - function(current - step_size)) / (2 * step_size);

            if (std::abs(derivative) < tolerance)
                throw std::runtime_error("Newton iteration encountered a near-zero derivative.");

            // newton step: x[n + 1] = x[n] - (f(n)/f'(n))
            const auto next = current - (value / derivative);

            if (std::abs(next - current) < tolerance)
                return next;

            current = next;
        }

        throw std::runtime_error("Newton iteration did not converge.");
    }

    /**
     * @ingroup nonlinear
     *
     * @brief Computes a vector-valued root using Newton's method.
     *
     * Solves the nonlinear system
     *
     * @f[
     * \mathbf{F}(\mathbf{x})
     * =
     * \mathbf{0}
     * @f]
     *
     * by the Newton iteration
     *
     * @f[
     * J_{\mathbf{F}}
     * \left(
     * \mathbf{x}^{(m)}
     * \right)
     * \Delta \mathbf{x}^{(m)}
     * =
     * -
     * \mathbf{F}
     * \left(
     * \mathbf{x}^{(m)}
     * \right),
     * @f]
     *
     * followed by
     *
     * @f[
     * \mathbf{x}^{(m+1)}
     * =
     * \mathbf{x}^{(m)}
     * +
     * \Delta \mathbf{x}^{(m)}.
     * @f]
     *
     * The Jacobian is approximated column-wise using central differences,
     *
     * @f[
     * \frac{
     * \partial \mathbf{F}
     * }{
     * \partial x_j
     * }
     * \approx
     * \frac{
     * \mathbf{F}
     * \left(
     * \mathbf{x} + h\mathbf{e}_j
     * \right)
     * -
     * \mathbf{F}
     * \left(
     * \mathbf{x} - h\mathbf{e}_j
     * \right)
     * }{
     * 2h
     * },
     * @f]
     *
     * where @f$\mathbf{e}_j@f$ is the @f$j@f$-th Cartesian basis vector
     * and @f$h@f$ is a small finite-difference step.
     *
     * The iteration is considered converged when either
     *
     * @f[
     * \left\|
     * \mathbf{F}
     * \left(
     * \mathbf{x}^{(m)}
     * \right)
     * \right\|_2
     * <
     * \varepsilon
     * @f]
     *
     * or
     *
     * @f[
     * \left\|
     * \Delta \mathbf{x}^{(m)}
     * \right\|_2
     * <
     * \varepsilon,
     * @f]
     *
     * where @f$\varepsilon@f$ is the specified tolerance.
     *
     * @tparam Function Type of the nonlinear vector-valued function.
     *
     * @param function Nonlinear function @f$\mathbf{F}(\mathbf{x})@f$.
     * @param initial_guess Initial approximation of the root.
     * @param tolerance Convergence tolerance.
     * @param max_iteration Maximum number of iterations.
     *
     * @return Approximation of the root.
     *
     * @throws std::runtime_error If the iteration does not converge within
     * the maximum number of iterations.
     */
    template <typename Function>
    [[nodiscard]] Eigen::VectorXd newton(Function&& function, const Eigen::VectorXd& initial_guess, double tolerance = 1e-10, int max_iteration = 100)
    {
        constexpr double step_size = 1e-6;
        Eigen::VectorXd current = initial_guess;

        for (int i = 0; i < max_iteration; i++)
        {
            // x[n]
            const Eigen::VectorXd value = function(current);

            if (value.norm() < tolerance)
                return current;

            // Building the jacobi matrix.
            Eigen::MatrixXd jacobian(value.size(), current.size());

            for (Eigen::Index j = 0; j < current.size(); ++j)
            {
                Eigen::VectorXd plus = current;
                Eigen::VectorXd minus = current;

                plus[j] += step_size;
                minus[j] -= step_size;

                // central difference: f'(n) = (f(n + h) - f(n - h)) / 2 * h
                jacobian.col(j) = (function(plus) - function(minus)) / (2.0 * step_size);
            }

            // newton step: x[n + 1] = x[n] - (f(n)/f'(n))
            const Eigen::VectorXd delta = jacobian.partialPivLu().solve(-value);
            const Eigen::VectorXd next = current + delta;

            if (delta.norm() < tolerance)
                return next;

            current = next;
        }

        throw std::runtime_error("Newton iteration did not converge.");
    }

} // namespace nane
