#pragma once

namespace nane
{
    /**
     * @ingroup calculus
     *
     * @brief Approximates the first derivative using a central difference.
     *
     * Approximates the derivative of a scalar function using
     *
     * @f[
     * f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}.
     * @f]
     *
     * The approximation has second-order accuracy,
     *
     * @f[
     * f'(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2).
     * @f]
     *
     * @tparam Function Type of the scalar function.
     *
     * @param function Function @f$f(x)@f$ whose derivative is approximated.
     * @param current Point @f$x@f$ at which the derivative is approximated.
     * @param step_size Finite-difference step @f$h@f$.
     *
     * @return Approximation of @f$f'(x)@f$.
     */
    template <typename Function>
    [[nodiscard]] double central_difference(Function&& function, double current, double step_size)
    {
        return (function(current + step_size) - function(current - step_size)) / (2.0 * step_size);
    }

} // namespace nane
