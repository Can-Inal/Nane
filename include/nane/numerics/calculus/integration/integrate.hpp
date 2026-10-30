#pragma once

#include "nane/geometry/grid.hpp"
#include "nane/numerics/calculus/integration/quadrature_rule.hpp"

#include <cstddef>

namespace nane::calculus
{
    /**
     * @ingroup calculus
     *
     * @brief Approximates a definite integral using a quadrature rule.
     *
     * For an interval @f$[a,b]@f$, the quadrature rule is applied as
     *
     * @f[
     * \int_a^b f(x)\,dx
     * \approx
     * (b-a)
     * \sum_{i=0}^{s-1}
     * w_i
     * f\left(a + (b-a)c_i\right),
     * @f]
     *
     * where @f$c_i@f$ are the quadrature nodes and @f$w_i@f$ are
     * the corresponding weights.
     *
     * If the grid has a prescribed point count, the rule is applied
     * separately on each grid subinterval and the results are summed.
     *
     * @param function Function to integrate.
     * @param rule Quadrature rule defined on the reference interval @f$[0,1]@f$.
     * @param interval Integration interval or prescribed integration grid.
     *
     * @return Approximation of the definite integral.
     */
    [[nodiscard]] inline double integrate(auto&& function, const nane::calculus::quadrature_rule& rule, const nane::geometry::grid<1>& interval)
    {
        const auto& nodes = rule.nodes();
        const auto& weights = rule.weights();

        // Apply the quadrature rule once on a given interval.
        const auto apply_rule = [&](double start, double stop)
        {
            const double width = stop - start;

            double sum = 0.0;

            for (std::size_t i = 0; i < rule.size(); ++i)
            {
                const double x = start + width * nodes[i];
                sum += weights[i] * function(x);
            }

            return width * sum;
        };

        // Apply the rule on each grid subinterval.
        if (interval.has_count(0))
        {
            const std::size_t count = *interval.count(0);
            const auto& axis = interval.axis(0);

            double result = 0.0;

            for (std::size_t i = 0; i + 1 < count; ++i)
            {
                result += apply_rule(axis[i], axis[i + 1]);
            }

            return result;
        }

        // Apply the rule once on the whole interval.
        return apply_rule(interval.start(0), interval.stop(0));
    }

} // namespace nane::calculus
