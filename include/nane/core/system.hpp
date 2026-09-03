#pragma once

#include <nane/core/details/factory.hpp>
#include <nane/core/details/system.hpp>
#include <utility>

namespace nane
{
    /**
     * @ingroup symbolic
     *
     * @brief Creates a vector-valued function from symbolic expressions.
     *
     * Each supplied expression defines one component of the returned vector.
     * Components may depend on any component of the input state, allowing
     * fully coupled systems to be represented.
     *
     * For example,
     *
     * @f[
     * \begin{aligned}
     * x_0' &= x_1, \\
     * x_1' &= -x_0
     * \end{aligned}
     * @f]
     *
     * can be written as
     *
     * @code{.cpp}
     * const nane::symbol x{1};
     *
     * const auto derivative = nane::system(
     *     x[1],
     *     -x[0]
     * );
     * @endcode
     *
     * @tparam Expressions Types of the component expressions.
     *
     * @param expressions Scalar symbolic expressions defining the components.
     *
     * @return Callable vector-valued function.
     */
    template <details::expression_operand... Expressions>
        requires(sizeof...(Expressions) > 0)
    [[nodiscard]] auto system(Expressions&&... expressions)
    {
        return details::system_function<decltype(details::make_expression(std::forward<Expressions>(expressions)))...>(
            details::make_expression(std::forward<Expressions>(expressions))...);
    }

} // namespace nane
