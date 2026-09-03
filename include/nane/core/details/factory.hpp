#pragma once

#include <nane/core/expression.hpp>
#include <type_traits>
#include <utility>

namespace nane::details
{
    template <typename Type>
    [[nodiscard]] auto make_expression(Type&& value)
    {
        if constexpr (symbolic_expression<Type>)
        {
            return std::remove_cvref_t<Type>(std::forward<Type>(value));
        }
        else
        {
            return nane::expression<constant_node>(constant_node(static_cast<double>(value)));
        }
    }

    template <typename Operation, typename Operand>
    [[nodiscard]] auto make_unary_expression(Operand&& operand)
    {
        auto expression = make_expression(std::forward<Operand>(operand));

        using node_type = unary_node<Operation, decltype(expression)>;

        return nane::expression<node_type>(node_type(std::move(expression)));
    }

    template <typename Operation, typename Left, typename Right>
    [[nodiscard]] auto make_binary_expression(Left&& left, Right&& right)
    {
        auto left_expression = make_expression(std::forward<Left>(left));
        auto right_expression = make_expression(std::forward<Right>(right));

        using node_type = binary_node<Operation, decltype(left_expression), decltype(right_expression)>;

        return nane::expression<node_type>(node_type(std::move(left_expression), std::move(right_expression)));
    }

} // namespace nane::details
