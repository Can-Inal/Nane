#pragma once

namespace nane
{
    /**
     * @ingroup calculus
     *
     * @brief Approximates the first derivative using a forward difference.
     *
     * Approximates the derivative of a scalar function
     *
     * @f[
     * f'(x)
     * @f]
     *
     * using the forward finite-difference formula
     *
     * @f[
     * f'(x)
     * \approx
     * \frac{
     * f(x+h) - f(x)
     * }{
     * h
     * }.
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
    [[nodiscard]] double forward_difference(Function&& function, double current, double step_size)
    {
        return (function(current + step_size) - function(current)) / step_size;
    }

} // namespace nane
