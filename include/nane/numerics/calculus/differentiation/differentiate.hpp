#pragma once

#include "nane/error.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"
#include "nane/numerics/calculus/differentiation/stencil.hpp"

#include <cmath>
#include <cstddef>
#include <type_traits>

namespace nane::calculus
{
    /**
     * @ingroup calculus
     *
     * @brief Approximates a scalar derivative using a finite-difference stencil.
     *
     * For a stencil of derivative order @f$m@f$,
     *
     * @f[
     * f^{(m)}(x)
     * \approx
     * \frac{1}{h^m}
     * \sum_{i=0}^{s-1}
     * w_i f(x + c_i h),
     * @f]
     *
     * where @f$c_i@f$ are the stencil offsets, @f$w_i@f$ are the
     * corresponding weights, and @f$h@f$ is the step size.
     *
     * @tparam Function Type of the function.
     *
     * @param function Function to differentiate.
     * @param point Point at which the derivative is approximated.
     * @param stencil Finite-difference stencil.
     *
     * @return Approximation of the derivative.
     */
    template <typename Function>
    [[nodiscard]] double differentiate(Function&& function, double point, const nane::calculus::stencil& stencil)
    {
        const auto& offsets = stencil.offsets();
        const auto& weights = stencil.weights();

        const double step_size = stencil.step_size();

        double sum = 0.0;

        for (std::size_t i = 0; i < stencil.size(); ++i)
        {
            const double x = point + static_cast<double>(offsets[i]) * step_size;
            sum += weights[i] * function(x);
        }

        const double scale = std::pow(step_size, static_cast<double>(stencil.derivative_order()));
        return sum / scale;
    }

    /**
     * @ingroup calculus
     *
     * @brief Approximates the derivative of a vector-valued function.
     *
     * For a function
     *
     * @f[
     * \mathbf{F}:\mathbb{R}^n \rightarrow \mathbb{R}^m,
     * @f]
     *
     * the derivative is the Jacobian matrix
     *
     * @f[
     * J_{ij}
     * =
     * \frac{\partial F_i}{\partial x_j}.
     * @f]
     *
     * Each column is approximated by applying the finite-difference
     * stencil in the corresponding coordinate direction.
     *
     * @tparam Function Type of the vector-valued function.
     * @tparam Dimension Input dimension.
     *
     * @param function Vector-valued function.
     * @param point Point at which the derivative is approximated.
     * @param stencil First-derivative finite-difference stencil.
     *
     * @return Approximation of the Jacobian matrix.
     */
    template <typename Function, std::size_t Dimension>
    [[nodiscard]] auto differentiate(Function&& function, const nane::vector<Dimension>& point, const nane::calculus::stencil& stencil)
    {
        NANE_REQUIRE(stencil.derivative_order() == 1, "vector differentiation requires a first-derivative stencil");

        const auto value = function(point);
        constexpr std::size_t output_dimension = std::remove_cvref_t<decltype(value)>::dimension();
        nane::matrix<output_dimension, Dimension> result(value.size(), point.size());

        result.fill(0.0);

        const auto& offsets = stencil.offsets();
        const auto& weights = stencil.weights();

        const double step_size = stencil.step_size();

        // Differentiate with respect to each input component.
        for (std::size_t j = 0; j < point.size(); ++j)
        {
            for (std::size_t k = 0; k < stencil.size(); ++k)
            {
                auto shifted = point;

                shifted[j] += static_cast<double>(offsets[k]) * step_size;

                const auto shifted_value = function(shifted);

                NANE_REQUIRE(shifted_value.size() == value.size(), "function output dimension must remain constant");

                for (std::size_t i = 0; i < value.size(); ++i)
                {
                    result(i, j) += weights[k] * shifted_value[i];
                }
            }

            for (std::size_t i = 0; i < value.size(); ++i)
                result(i, j) /= step_size;
        }

        return result;
    }

} // namespace nane::calculus
