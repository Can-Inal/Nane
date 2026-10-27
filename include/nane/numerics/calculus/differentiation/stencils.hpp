#pragma once

#include "nane/numerics/calculus/differentiation/stencil.hpp"

namespace nane::calculus::stencils
{
    /**
     * @ingroup calculus
     *
     * @brief Creates a centered first-derivative stencil.
     *
     * Uses
     *
     * @f[
     * f'(x)
     * \approx
     * \frac{f(x+h)-f(x-h)}{2h}.
     * @f]
     *
     * @param step_size Finite-difference step size.
     *
     * @return Centered first-derivative stencil.
     */
    [[nodiscard]] inline nane::calculus::stencil central_first(double step_size)
    {
        return {
            1,
            step_size,
            {
                {-1, -0.5},
                {1, 0.5},
            },
        };
    }

    /**
     * @ingroup calculus
     *
     * @brief Creates a centered second-derivative stencil.
     *
     * Uses
     *
     * @f[
     * f''(x)
     * \approx
     * \frac{
     * f(x-h)-2f(x)+f(x+h)
     * }{h^2}.
     * @f]
     *
     * @param step_size Finite-difference step size.
     *
     * @return Centered second-derivative stencil.
     */
    [[nodiscard]] inline nane::calculus::stencil central_second(double step_size)
    {
        return {
            2,
            step_size,
            {
                {-1, 1.0},
                {0, -2.0},
                {1, 1.0},
            },
        };
    }

    /**
     * @ingroup calculus
     *
     * @brief Creates a forward first-derivative stencil.
     *
     * Uses
     *
     * @f[
     * f'(x)
     * \approx
     * \frac{f(x+h)-f(x)}{h}.
     * @f]
     *
     * @param step_size Finite-difference step size.
     *
     * @return Forward first-derivative stencil.
     */
    [[nodiscard]] inline nane::calculus::stencil forward_first(double step_size)
    {
        return {
            1,
            step_size,
            {
                {0, -1.0},
                {1, 1.0},
            },
        };
    }

    /**
     * @ingroup calculus
     *
     * @brief Creates a backward first-derivative stencil.
     *
     * Uses
     *
     * @f[
     * f'(x)
     * \approx
     * \frac{f(x)-f(x-h)}{h}.
     * @f]
     *
     * @param step_size Finite-difference step size.
     *
     * @return Backward first-derivative stencil.
     */
    [[nodiscard]] inline nane::calculus::stencil backward_first(double step_size)
    {
        return {
            1,
            step_size,
            {
                {-1, -1.0},
                {0, 1.0},
            },
        };
    }

} // namespace nane::calculus::stencils
