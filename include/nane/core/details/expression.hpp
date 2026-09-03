#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

namespace nane::details
{
    struct expression_marker
    {
    };

    template <typename Type>
    concept symbolic_expression = std::derived_from<std::remove_cvref_t<Type>, expression_marker>;

    template <typename Type>
    concept expression_operand = symbolic_expression<Type> || std::is_arithmetic_v<std::remove_cvref_t<Type>>;

    template <std::size_t Position = 0, typename Tuple>
    [[nodiscard]] double evaluate_scalar_argument(std::size_t index, Tuple& arguments)
    {
        if constexpr (Position == std::tuple_size_v<std::remove_reference_t<Tuple>>)
        {
            throw std::out_of_range("symbol index exceeds the number of supplied arguments.");
        }
        else
        {
            if (index == Position)
            {
                auto&& value = std::get<Position>(arguments);

                if constexpr (std::convertible_to<decltype(value), double>)
                {
                    return static_cast<double>(value);
                }
                else
                {
                    throw std::invalid_argument("symbol does not refer to a scalar argument.");
                }
            }

            return evaluate_scalar_argument<Position + 1>(index, arguments);
        }
    }

    template <std::size_t Position = 0, typename Tuple>
    [[nodiscard]] double evaluate_component_argument(std::size_t index, std::size_t component, Tuple& arguments)
    {
        if constexpr (Position == std::tuple_size_v<std::remove_reference_t<Tuple>>)
        {
            throw std::out_of_range("symbol index exceeds the number of supplied arguments.");
        }
        else
        {
            if (index == Position)
            {
                auto&& value = std::get<Position>(arguments);

                if constexpr (requires {
                                  value.size();
                                  value[0];
                              })
                {
                    const auto size = value.size();
                    using index_type = std::remove_cvref_t<decltype(size)>;

                    if constexpr (std::is_signed_v<index_type>)
                    {
                        if (size < 0)
                            throw std::out_of_range("symbolic vector argument has an invalid size.");
                    }

                    if (component >= static_cast<std::size_t>(size))
                        throw std::out_of_range("symbol component index exceeds the argument size.");

                    auto&& result = value[static_cast<index_type>(component)];

                    if constexpr (std::convertible_to<decltype(result), double>)
                    {
                        return static_cast<double>(result);
                    }
                    else
                    {
                        throw std::invalid_argument("symbol component is not scalar-valued.");
                    }
                }
                else
                {
                    throw std::invalid_argument("indexed symbol does not refer to a vector-valued argument.");
                }
            }

            return evaluate_component_argument<Position + 1>(index, component, arguments);
        }
    }

    class constant_node
    {
    public:
        explicit constexpr constant_node(double value) noexcept : value_(value)
        {
        }

        template <typename... Arguments>
        [[nodiscard]] constexpr double operator()(Arguments&&...) const noexcept
        {
            return value_;
        }

    private:
        double value_;
    };

    class component_node
    {
    public:
        constexpr component_node(std::size_t symbol_index, std::size_t component_index) noexcept
            : symbol_index_(symbol_index), component_index_(component_index)
        {
        }

        template <typename... Arguments>
        [[nodiscard]] double operator()(Arguments&&... arguments) const
        {
            auto values = std::forward_as_tuple(arguments...);

            return evaluate_component_argument(symbol_index_, component_index_, values);
        }

    private:
        std::size_t symbol_index_;
        std::size_t component_index_;
    };

    struct addition
    {
        [[nodiscard]] constexpr double operator()(double left, double right) const noexcept
        {
            return left + right;
        }
    };

    struct subtraction
    {
        [[nodiscard]] constexpr double operator()(double left, double right) const noexcept
        {
            return left - right;
        }
    };

    struct multiplication
    {
        [[nodiscard]] constexpr double operator()(double left, double right) const noexcept
        {
            return left * right;
        }
    };

    struct division
    {
        [[nodiscard]] constexpr double operator()(double left, double right) const noexcept
        {
            return left / right;
        }
    };

    struct negation
    {
        [[nodiscard]] constexpr double operator()(double value) const noexcept
        {
            return -value;
        }
    };

    struct sine
    {
        [[nodiscard]] double operator()(double value) const
        {
            return std::sin(value);
        }
    };

    struct cosine
    {
        [[nodiscard]] double operator()(double value) const
        {
            return std::cos(value);
        }
    };

    struct tangent
    {
        [[nodiscard]] double operator()(double value) const
        {
            return std::tan(value);
        }
    };

    struct exponential
    {
        [[nodiscard]] double operator()(double value) const
        {
            return std::exp(value);
        }
    };

    struct logarithm
    {
        [[nodiscard]] double operator()(double value) const
        {
            return std::log(value);
        }
    };

    struct square_root
    {
        [[nodiscard]] double operator()(double value) const
        {
            return std::sqrt(value);
        }
    };

    struct absolute_value
    {
        [[nodiscard]] double operator()(double value) const
        {
            return std::abs(value);
        }
    };

    struct power
    {
        [[nodiscard]] double operator()(double base, double exponent) const
        {
            return std::pow(base, exponent);
        }
    };

    template <typename Operation, typename Operand>
    class unary_node
    {
    public:
        explicit unary_node(Operand operand) : operand_(std::move(operand))
        {
        }

        template <typename... Arguments>
        [[nodiscard]] double operator()(Arguments&&... arguments) const
        {
            return Operation{}(operand_(arguments...));
        }

    private:
        Operand operand_;
    };

    template <typename Operation, typename Left, typename Right>
    class binary_node
    {
    public:
        binary_node(Left left, Right right) : left_(std::move(left)), right_(std::move(right))
        {
        }

        template <typename... Arguments>
        [[nodiscard]] double operator()(Arguments&&... arguments) const
        {
            return Operation{}(left_(arguments...), right_(arguments...));
        }

    private:
        Left left_;
        Right right_;
    };

} // namespace nane::details
