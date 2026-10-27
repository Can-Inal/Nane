#pragma once

#include "nane/error.hpp"

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <vector>

namespace nane::calculus
{
    /**
     * @ingroup calculus
     *
     * @brief Represents a finite-difference stencil.
     *
     * A stencil approximates a derivative of order @f$m@f$ as
     *
     * @f[
     * f^{(m)}(x)
     * \approx
     * \frac{1}{h^m}
     * \sum_{i=0}^{s-1}
     * w_i f(x + c_i h),
     * @f]
     *
     * where @f$c_i@f$ are integer offsets, @f$w_i@f$ are the
     * corresponding weights, and @f$h@f$ is the step size.
     *
     * Each coefficient pair is given as
     *
     * @f[
     * (c_i, w_i).
     * @f]
     */
    class stencil
    {
    public:
        /**
         * @brief Constructs a finite-difference stencil.
         *
         * @param derivative_order Order of the approximated derivative.
         * @param step_size Finite-difference step size.
         * @param coefficients Offset-weight pairs @f$(c_i, w_i)@f$.
         */
        stencil(std::size_t derivative_order, double step_size, std::initializer_list<std::pair<std::ptrdiff_t, double>> coefficients)
            : m_DerivativeOrder(derivative_order), m_StepSize(step_size)
        {
            NANE_REQUIRE(derivative_order > 0, "stencil derivative order must be positive");
            NANE_REQUIRE(step_size > 0.0, "stencil step size must be positive");
            NANE_REQUIRE(coefficients.size() > 0, "stencil must contain at least one coefficient");

            m_Offsets.reserve(coefficients.size());
            m_Weights.reserve(coefficients.size());

            for (const auto& [offset, weight] : coefficients)
            {
                m_Offsets.push_back(offset);
                m_Weights.push_back(weight);
            }
        }

        /**
         * @brief Returns the derivative order.
         */
        [[nodiscard]] std::size_t derivative_order() const noexcept
        {
            return m_DerivativeOrder;
        }

        /**
         * @brief Returns the finite-difference step size.
         */
        [[nodiscard]] double step_size() const noexcept
        {
            return m_StepSize;
        }

        /**
         * @brief Returns the stencil offsets.
         */
        [[nodiscard]] const std::vector<std::ptrdiff_t>& offsets() const noexcept
        {
            return m_Offsets;
        }

        /**
         * @brief Returns the stencil weights.
         */
        [[nodiscard]] const std::vector<double>& weights() const noexcept
        {
            return m_Weights;
        }

        /**
         * @brief Returns the number of stencil coefficients.
         */
        [[nodiscard]] std::size_t size() const noexcept
        {
            return m_Offsets.size();
        }

    private:
        std::size_t m_DerivativeOrder;
        double m_StepSize;

        std::vector<std::ptrdiff_t> m_Offsets;
        std::vector<double> m_Weights;
    };

} // namespace nane::calculus
