#pragma once

#include "nane/numerics/calculus/integration/quadrature_rule.hpp"

#include <cmath>

namespace nane::calculus::rules
{
    /**
     * @ingroup calculus
     *
     * @brief Trapezoidal quadrature rule.
     *
     * On the reference interval @f$[0,1]@f$,
     *
     * @f[
     * \int_0^1 f(x)\,dx
     * \approx
     * \frac{1}{2}f(0)
     * +
     * \frac{1}{2}f(1).
     * @f]
     */
    inline const nane::calculus::quadrature_rule trapezoidal{
        {0.0, 0.5},
        {1.0, 0.5},
    };

    /**
     * @ingroup calculus
     *
     * @brief Simpson quadrature rule.
     *
     * On the reference interval @f$[0,1]@f$,
     *
     * @f[
     * \int_0^1 f(x)\,dx
     * \approx
     * \frac{1}{6}f(0)
     * +
     * \frac{4}{6}f\left(\frac{1}{2}\right)
     * +
     * \frac{1}{6}f(1).
     * @f]
     */
    inline const nane::calculus::quadrature_rule simpson{
        {0.0, 1.0 / 6.0},
        {0.5, 4.0 / 6.0},
        {1.0, 1.0 / 6.0},
    };

    /**
     * @ingroup calculus
     *
     * @brief Two-point Gauss-Legendre quadrature rule.
     *
     * The rule has degree of exactness three.
     */
    inline const nane::calculus::quadrature_rule gauss_legendre_2{
        {0.5 * (1.0 - 1.0 / std::sqrt(3.0)), 0.5},
        {0.5 * (1.0 + 1.0 / std::sqrt(3.0)), 0.5},
    };

    /**
     * @ingroup calculus
     *
     * @brief Three-point Gauss-Legendre quadrature rule.
     *
     * The rule has degree of exactness five.
     */
    inline const nane::calculus::quadrature_rule gauss_legendre_3{
        {0.5 * (1.0 - std::sqrt(3.0 / 5.0)), 5.0 / 18.0},
        {0.5, 4.0 / 9.0},
        {0.5 * (1.0 + std::sqrt(3.0 / 5.0)), 5.0 / 18.0},
    };

} // namespace nane::calculus::rules
