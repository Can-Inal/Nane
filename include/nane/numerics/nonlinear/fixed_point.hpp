#pragma once

#include "nane/error.hpp"
#include "nane/math/vector.hpp"

#include <cmath>
#include <cstddef>

namespace nane::nonlinear
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
    [[nodiscard]] inline double fixed_point(auto&& function, double initial_guess, double tolerance = 1e-10, int max_iteration = 100)
    {
        auto current = initial_guess;

        for (auto i = 0; i < max_iteration; i++)
        {
            auto next = function(current);

            if (std::abs(next - current) < tolerance)
                return next;

            current = next;
        }

        NANE_FAIL(true, "fixed-point iteration did not converge.");
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
     * @tparam Dimension Vector dimension.
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
    template <std::size_t Dimension>
    [[nodiscard]] nane::vector<Dimension> fixed_point(auto&& function, const nane::vector<Dimension>& initial_guess, double tolerance = 1e-10,
                                                      int max_iterations = 100)
    {
        auto current = initial_guess;

        for (auto i = 0; i < max_iterations; i++)
        {
            auto next = function(current);

            if ((next - current).norm() < tolerance)
                return next;

            current = next;
        }

        NANE_FAIL(true, "fixed-point iteration did not converge.");
    }

} // namespace nane::nonlinear
