#pragma once

#include <cstddef>
#include <nane/core/details/factory.hpp>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace nane
{
    /**
     * @ingroup symbolic
     *
     * @brief Represents an argument of a symbolic expression.
     *
     * A symbol refers to one of the arguments supplied when evaluating a
     * symbolic function. Argument indices are zero-based.
     *
     * Symbols are normally created with nane::symbols():
     *
     * @code{.cpp}
     * const auto [t, x] = nane::symbols<2>();
     *
     * const auto derivative = nane::function(
     *     -t * t * x
     * );
     * @endcode
     *
     * This represents
     *
     * @f[
     * f(t,x) = -t^2x.
     * @f]
     *
     * A symbol referring to a vector-valued argument can be indexed:
     *
     * @code{.cpp}
     * const auto [t, x] = nane::symbols<2>();
     *
     * const auto expression = x[0] + x[1];
     * @endcode
     *
     * @note Symbol and component indices are zero-based.
     */
    class symbol : public details::expression_marker
    {
    public:
        /**
         * @brief Constructs a symbol referring to an argument.
         *
         * @param index Zero-based argument index.
         */
        explicit constexpr symbol(std::size_t index) noexcept : index_(index)
        {
        }

        /**
         * @brief Returns the argument index represented by the symbol.
         *
         * @return Zero-based argument index.
         */
        [[nodiscard]] constexpr std::size_t index() const noexcept
        {
            return index_;
        }

        /**
         * @brief Creates an expression referring to a component of a
         * vector-valued argument.
         *
         * @param component Zero-based component index.
         *
         * @return Symbolic expression representing the selected component.
         */
        [[nodiscard]] auto operator[](std::size_t component) const
        {
            return expression<details::component_node>(details::component_node(index_, component));
        }

        /**
         * @brief Evaluates the symbol for the supplied arguments.
         *
         * @tparam Arguments Types of the supplied arguments.
         *
         * @param arguments Arguments used to evaluate the expression.
         *
         * @return Scalar value represented by the symbol.
         *
         * @throws std::out_of_range If the argument index does not exist.
         * @throws std::invalid_argument If the selected argument is not
         * scalar-valued.
         */
        template <typename... Arguments>
        [[nodiscard]] double operator()(Arguments&&... arguments) const
        {
            auto values = std::forward_as_tuple(arguments...);

            return details::evaluate_scalar_argument(index_, values);
        }

    private:
        std::size_t index_;
    };

    /**
     * @ingroup symbolic
     *
     * @brief Creates an ordered collection of symbolic arguments.
     *
     * The returned tuple contains symbols with consecutive zero-based
     * argument indices.
     *
     * Structured bindings provide a convenient way to name the symbols:
     *
     * @code{.cpp}
     * const auto [t, x] = nane::symbols<2>();
     * @endcode
     *
     * Here, @c t represents argument zero and @c x represents argument one.
     * The names themselves have no special meaning; their position determines
     * the corresponding argument index.
     *
     * For example,
     *
     * @code{.cpp}
     * const auto [t, x] = nane::symbols<2>();
     *
     * const auto derivative = nane::function(
     *     -t * t * x
     * );
     * @endcode
     *
     * @tparam Count Number of symbols to create.
     *
     * @return Tuple containing symbols with indices from zero to
     * @p Count minus one.
     */
    template <std::size_t Count>
    [[nodiscard]] constexpr auto symbols()
    {
        static_assert(Count > 0, "at least one symbol must be requested.");

        return []<std::size_t... Indices>(std::index_sequence<Indices...>)
        {
            return std::tuple{
                symbol{Indices}...,
            };
        }(std::make_index_sequence<Count>{});
    }

    template <details::symbolic_expression Expression>
    [[nodiscard]] auto operator-(Expression&& expression)
    {
        return details::make_unary_expression<details::negation>(std::forward<Expression>(expression));
    }

    template <details::expression_operand Left, details::expression_operand Right>
        requires(details::symbolic_expression<Left> || details::symbolic_expression<Right>)
    [[nodiscard]] auto operator+(Left&& left, Right&& right)
    {
        return details::make_binary_expression<details::addition>(std::forward<Left>(left), std::forward<Right>(right));
    }

    template <details::expression_operand Left, details::expression_operand Right>
        requires(details::symbolic_expression<Left> || details::symbolic_expression<Right>)
    [[nodiscard]] auto operator-(Left&& left, Right&& right)
    {
        return details::make_binary_expression<details::subtraction>(std::forward<Left>(left), std::forward<Right>(right));
    }

    template <details::expression_operand Left, details::expression_operand Right>
        requires(details::symbolic_expression<Left> || details::symbolic_expression<Right>)
    [[nodiscard]] auto operator*(Left&& left, Right&& right)
    {
        return details::make_binary_expression<details::multiplication>(std::forward<Left>(left), std::forward<Right>(right));
    }

    template <details::expression_operand Left, details::expression_operand Right>
        requires(details::symbolic_expression<Left> || details::symbolic_expression<Right>)
    [[nodiscard]] auto operator/(Left&& left, Right&& right)
    {
        return details::make_binary_expression<details::division>(std::forward<Left>(left), std::forward<Right>(right));
    }

} // namespace nane
