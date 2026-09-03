#pragma once

#include <nane/core/details/factory.hpp>
#include <utility>

namespace nane
{
    /**
     * @ingroup symbolic
     *
     * @brief Creates a scalar-valued function from a symbolic expression.
     *
     * The returned expression is callable and can be passed directly to
     * numerical algorithms accepting a scalar-valued function object.
     *
     * For example,
     *
     * @code{.cpp}
     * const nane::symbol t{0};
     * const nane::symbol x{1};
     *
     * const auto derivative = nane::function(
     *     -t * t * x
     * );
     * @endcode
     *
     * represents
     *
     * @f[
     * f(t,x) = -t^2x.
     * @f]
     *
     * @tparam Expression Type of the symbolic expression.
     *
     * @param expression Scalar symbolic expression.
     *
     * @return Callable scalar-valued expression.
     */
    template <details::expression_operand Expression>
    [[nodiscard]] auto function(Expression&& expression)
    {
        return details::make_expression(std::forward<Expression>(expression));
    }

} // namespace nane
