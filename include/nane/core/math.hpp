#pragma once

#include <nane/core/details/factory.hpp>
#include <utility>

namespace nane
{
    /**
     * @ingroup symbolic
     *
     * @brief Creates the sine of a symbolic expression.
     *
     * @param expression Symbolic expression.
     *
     * @return Expression representing @f$\sin(x)@f$.
     */
    template <details::symbolic_expression Expression>
    [[nodiscard]] auto sin(Expression&& expression)
    {
        return details::make_unary_expression<details::sine>(std::forward<Expression>(expression));
    }

    /**
     * @ingroup symbolic
     *
     * @brief Creates the cosine of a symbolic expression.
     *
     * @param expression Symbolic expression.
     *
     * @return Expression representing @f$\cos(x)@f$.
     */
    template <details::symbolic_expression Expression>
    [[nodiscard]] auto cos(Expression&& expression)
    {
        return details::make_unary_expression<details::cosine>(std::forward<Expression>(expression));
    }

    /**
     * @ingroup symbolic
     *
     * @brief Creates the tangent of a symbolic expression.
     *
     * @param expression Symbolic expression.
     *
     * @return Expression representing @f$\tan(x)@f$.
     */
    template <details::symbolic_expression Expression>
    [[nodiscard]] auto tan(Expression&& expression)
    {
        return details::make_unary_expression<details::tangent>(std::forward<Expression>(expression));
    }

    /**
     * @ingroup symbolic
     *
     * @brief Creates the exponential of a symbolic expression.
     *
     * @param expression Symbolic expression.
     *
     * @return Expression representing @f$e^x@f$.
     */
    template <details::symbolic_expression Expression>
    [[nodiscard]] auto exp(Expression&& expression)
    {
        return details::make_unary_expression<details::exponential>(std::forward<Expression>(expression));
    }

    /**
     * @ingroup symbolic
     *
     * @brief Creates the natural logarithm of a symbolic expression.
     *
     * @param expression Symbolic expression.
     *
     * @return Expression representing @f$\log(x)@f$.
     */
    template <details::symbolic_expression Expression>
    [[nodiscard]] auto log(Expression&& expression)
    {
        return details::make_unary_expression<details::logarithm>(std::forward<Expression>(expression));
    }

    /**
     * @ingroup symbolic
     *
     * @brief Creates the square root of a symbolic expression.
     *
     * @param expression Symbolic expression.
     *
     * @return Expression representing @f$\sqrt{x}@f$.
     */
    template <details::symbolic_expression Expression>
    [[nodiscard]] auto sqrt(Expression&& expression)
    {
        return details::make_unary_expression<details::square_root>(std::forward<Expression>(expression));
    }

    /**
     * @ingroup symbolic
     *
     * @brief Creates the absolute value of a symbolic expression.
     *
     * @param expression Symbolic expression.
     *
     * @return Expression representing @f$|x|@f$.
     */
    template <details::symbolic_expression Expression>
    [[nodiscard]] auto abs(Expression&& expression)
    {
        return details::make_unary_expression<details::absolute_value>(std::forward<Expression>(expression));
    }

    /**
     * @ingroup symbolic
     *
     * @brief Raises a symbolic expression to a power.
     *
     * At least one of the supplied operands must be symbolic.
     *
     * For example,
     *
     * @code{.cpp}
     * const nane::symbol x{0};
     *
     * const auto square = nane::pow(x, 2.0);
     * @endcode
     *
     * represents
     *
     * @f[
     * x^2.
     * @f]
     *
     * @tparam Base Type of the base expression.
     * @tparam Exponent Type of the exponent expression.
     *
     * @param base Base expression.
     * @param exponent Exponent expression.
     *
     * @return Expression representing the power.
     */
    template <details::expression_operand Base, details::expression_operand Exponent>
        requires(details::symbolic_expression<Base> || details::symbolic_expression<Exponent>)
    [[nodiscard]] auto pow(Base&& base, Exponent&& exponent)
    {
        return details::make_binary_expression<details::power>(std::forward<Base>(base), std::forward<Exponent>(exponent));
    }

} // namespace nane
