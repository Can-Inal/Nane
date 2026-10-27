#pragma once

#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"
#include "nane/numerics/calculus/central_difference.hpp"

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace nane::nonlinear
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
     * using Newton iteration.
     *
     * @tparam Function Type of the nonlinear function.
     *
     * @param function Nonlinear function.
     * @param initial_guess Initial approximation of the root.
     * @param tolerance Convergence tolerance.
     * @param max_iteration Maximum number of iterations.
     *
     * @return Approximation of the root.
     */
    template <typename Function>
    [[nodiscard]] double newton(Function&& function, double initial_guess, double tolerance = 1e-10, int max_iteration = 100)
    {
        constexpr double step_size = 1e-6;

        auto current = initial_guess;

        for (int i = 0; i < max_iteration; ++i)
        {
            const auto value = function(current);

            if (std::abs(value) < tolerance)
                return current;

            const auto derivative = nane::calculus::central_difference(function, current, step_size);

            if (std::abs(derivative) < tolerance)
            {
                throw std::runtime_error("Newton iteration encountered a near-zero derivative.");
            }

            const auto next = current - value / derivative;

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
     * Solves
     *
     * @f[
     * \mathbf{F}(\mathbf{x}) = \mathbf{0}
     * @f]
     *
     * using Newton iteration. The Jacobian is approximated using
     * central differences.
     *
     * @tparam Function Type of the nonlinear function.
     * @tparam Dimension Vector dimension.
     *
     * @param function Nonlinear function.
     * @param initial_guess Initial approximation of the root.
     * @param tolerance Convergence tolerance.
     * @param max_iteration Maximum number of iterations.
     *
     * @return Approximation of the root.
     */
    template <typename Function, std::size_t Dimension>
    [[nodiscard]] nane::vector<Dimension> newton(Function&& function, const nane::vector<Dimension>& initial_guess, double tolerance = 1e-10,
                                                 int max_iteration = 100)
    {
        constexpr double step_size = 1e-6;

        const std::size_t dimension = initial_guess.size();

        nane::vector<Dimension> current = initial_guess;

        for (int iteration = 0; iteration < max_iteration; ++iteration)
        {
            const auto value = function(current);

            if (value.size() != dimension)
            {
                throw std::invalid_argument("Newton function dimension must match initial guess dimension");
            }

            if (value.norm() < tolerance)
                return current;

            nane::matrix<> jacobian(dimension, dimension);

            // J_ij =
            // (F_i(x + h e_j) - F_i(x - h e_j)) / (2h)
            for (std::size_t j = 0; j < dimension; ++j)
            {
                auto positive = current;
                auto negative = current;

                positive[j] += step_size;
                negative[j] -= step_size;

                const auto positive_value = function(positive);

                const auto negative_value = function(negative);

                if (positive_value.size() != dimension || negative_value.size() != dimension)
                {
                    throw std::invalid_argument("Newton function dimension must match initial guess dimension");
                }

                for (std::size_t i = 0; i < dimension; ++i)
                {
                    jacobian(i, j) = (positive_value[i] - negative_value[i]) / (2.0 * step_size);
                }
            }

            nane::vector<Dimension> delta;

            if constexpr (Dimension == nane::dynamic)
                delta = nane::vector<>(dimension);

            using eigen_matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

            const auto eigen_dimension = static_cast<Eigen::Index>(dimension);

            Eigen::Map<const eigen_matrix> jacobian_map(jacobian.data(), eigen_dimension, eigen_dimension);

            Eigen::Map<const Eigen::VectorXd> value_map(value.data(), eigen_dimension);

            Eigen::Map<Eigen::VectorXd> delta_map(delta.data(), eigen_dimension);

            // J(x[n]) * delta = -F(x[n])
            //
            // Eigen operates directly on Nane-owned memory here.
            delta_map = jacobian_map.partialPivLu().solve(-value_map);

            const auto next = current + delta;

            if (delta.norm() < tolerance)
                return next;

            current = next;
        }

        throw std::runtime_error("Newton iteration did not converge.");
    }

} // namespace nane::nonlinear
