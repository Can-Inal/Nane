#pragma once

#include "nane/error.hpp"
#include "nane/math/vector.hpp"

#include <cstddef>
#include <initializer_list>
#include <utility>

namespace nane::calculus
{
    /**
     * @ingroup calculus
     *
     * @brief Represents a numerical quadrature rule.
     *
     * A quadrature rule approximates an integral on the reference interval
     * @f$[0,1]@f$ using nodes @f$c_i@f$ and weights @f$w_i@f$,
     *
     * @f[
     * \int_0^1 f(x)\,dx
     * \approx
     * \sum_{i=0}^{s-1}
     * w_i f(c_i).
     * @f]
     *
     * Each coefficient pair is given as
     *
     * @f[
     * (c_i, w_i).
     * @f]
     */
    class quadrature_rule
    {
    public:
        /**
         * @brief Constructs a quadrature rule from node-weight pairs.
         *
         * Each pair contains a quadrature node followed by its corresponding
         * weight.
         *
         * @param coefficients Node-weight pairs @f$(c_i, w_i)@f$.
         */
        quadrature_rule(std::initializer_list<std::pair<double, double>> coefficients) : m_Nodes(coefficients.size()), m_Weights(coefficients.size())
        {
            NANE_REQUIRE(coefficients.size() > 0, "quadrature rule must contain at least one node");

            std::size_t i = 0;

            for (const auto& [node, weight] : coefficients)
            {
                m_Nodes[i] = node;
                m_Weights[i] = weight;

                ++i;
            }
        }

        /**
         * @brief Returns the quadrature nodes.
         */
        [[nodiscard]] const nane::vector<>& nodes() const noexcept
        {
            return m_Nodes;
        }

        /**
         * @brief Returns the quadrature weights.
         */
        [[nodiscard]] const nane::vector<>& weights() const noexcept
        {
            return m_Weights;
        }

        /**
         * @brief Returns the number of quadrature nodes.
         */
        [[nodiscard]] std::size_t size() const noexcept
        {
            return m_Nodes.size();
        }

    private:
        nane::vector<> m_Nodes;
        nane::vector<> m_Weights;
    };

} // namespace nane::calculus
