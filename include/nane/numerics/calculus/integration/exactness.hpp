#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/numerics/calculus/integration/integrate.hpp"
#include "nane/numerics/calculus/integration/quadrature_rule.hpp"

#include <cmath>

namespace nane::calculus
{

    /**
     * @brief Determines the degree of exactness of a quadrature rule.
     *
     * For an interval @f$[0, 1]@f$, the quadrature rule is applied to the test function
     *
     * @f[
     * x^k
     * @f]
     *
     * whose analytical integral is
     *
     * @f[
     * \int_0^1 x^k\,dx = \frac{1}{k + 1}.
     * @f]
     *
     * If the numerical result differs by more than the numerical tolerance for degree @f$k@f$,
     * the quadrature rule has degree of exactness @f$k-1@f$.
     *
     * @param rule Quadrature rule.
     * @param check_limit Highest polynomial degree to test.
     *
     * @return Degree of exactness for the provided quadrature rule.
     */
    [[nodiscard]] inline int degree_of_exactness(const nane::calculus::quadrature_rule& rule, int check_limit = 6)
    {
        NANE_REQUIRE(check_limit >= 0, "Check limit must be non-negative.");

        const auto tolerance = 1e-12;

        for (auto i = 0; i <= check_limit; i++)
        {
            const auto test_function = [i](double x) { return std::pow(x, i); };
            const auto analytic_result = 1.0 / (static_cast<double>(i) + 1.0);

            const auto numeric_result = nane::calculus::integrate(test_function, rule, {0.0, 1.0});

            if (std::abs(numeric_result - analytic_result) >= tolerance)
                return i - 1;
        }

        NANE_FAIL(true, "Exactness could not be determined within the provided limit");
    }

} // namespace nane::calculus
