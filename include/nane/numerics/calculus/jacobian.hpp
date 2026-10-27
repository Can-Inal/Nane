#pragma once

#include "nane/error.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

#include <cstddef>
#include <type_traits>

namespace nane::calculus
{
    /**
     * @ingroup calculus
     *
     * @brief Approximates the Jacobian matrix using central differences.
     *
     * Each column is approximated using a central difference,
     *
     * @f[
     * \frac{\partial \mathbf{F}}{\partial x_j}
     * \approx
     * \frac{
     * \mathbf{F}(\mathbf{x} + h\mathbf{e}_j)
     * -
     * \mathbf{F}(\mathbf{x} - h\mathbf{e}_j)
     * }{2h}.
     * @f]
     *
     * The approximation is second-order accurate with respect to
     * the finite-difference step,
     *
     * @f[
     * O(h^2).
     * @f]
     *
     * @tparam Function Type of the vector-valued function.
     * @tparam Dimension Input dimension.
     *
     * @param function Vector-valued function.
     * @param current Point at which the Jacobian is approximated.
     * @param step_size Finite-difference step.
     *
     * @return Approximation of the Jacobian matrix.
     */
    template <typename Function, std::size_t Dimension>
    [[nodiscard]] auto jacobian(Function&& function, const nane::vector<Dimension>& current, double step_size)
    {
        const auto value = function(current);
        constexpr std::size_t output_dimension = std::remove_cvref_t<decltype(value)>::dimension();

        nane::matrix<output_dimension, Dimension> result(value.size(), current.size());

        for (std::size_t j = 0; j < current.size(); ++j)
        {
            auto plus = current;
            auto minus = current;

            plus[j] += step_size;
            minus[j] -= step_size;

            const auto plus_value = function(plus);
            const auto minus_value = function(minus);

            NANE_REQUIRE(plus_value.size() == value.size() && minus_value.size() == value.size(), "function output dimension must remain constant");

            for (std::size_t i = 0; i < value.size(); ++i)
            {
                result(i, j) = (plus_value[i] - minus_value[i]) / (2.0 * step_size);
            }
        }

        return result;
    }

} // namespace nane::calculus
