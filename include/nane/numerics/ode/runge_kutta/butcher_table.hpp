#pragma once

#include "nane/error.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

#include <cstddef>
#include <initializer_list>
#include <optional>

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Represents the Butcher coefficients of a Runge-Kutta method.
     *
     * An optional second weight vector defines an embedded Runge-Kutta
     * approximation for error estimation.
     */
    class butcher_table
    {
    public:
        /**
         * @brief Constructs a Butcher table from its coefficients.
         *
         * @param alpha Stage-time coefficients.
         * @param beta Stage coupling coefficients.
         * @param gamma Solution weights.
         * @param embedded_gamma Optional embedded solution weights.
         */
        butcher_table(std::initializer_list<double> alpha, std::initializer_list<std::initializer_list<double>> beta,
                      std::initializer_list<double> gamma, std::initializer_list<double> embedded_gamma = {})
            : m_Alpha(alpha.size()), m_Beta(alpha.size(), alpha.size()), m_Gamma(gamma.size())
        {
            const std::size_t stage_count = alpha.size();

            NANE_REQUIRE(stage_count > 0, "Butcher table must contain at least one stage");
            NANE_REQUIRE(beta.size() == stage_count, "Butcher table beta row count must match stage count");
            NANE_REQUIRE(gamma.size() == stage_count, "Butcher table gamma size must match stage count");
            NANE_REQUIRE(embedded_gamma.size() == 0 || embedded_gamma.size() == stage_count,
                         "Butcher table embedded gamma size must match stage count");

            std::size_t i = 0;

            for (const double value : alpha)
                m_Alpha[i++] = value;

            i = 0;

            for (const auto& row : beta)
            {
                NANE_REQUIRE(row.size() == stage_count, "Butcher table beta must be square");

                std::size_t j = 0;

                for (const double value : row)
                    m_Beta(i, j++) = value;

                ++i;
            }

            i = 0;

            for (const double value : gamma)
                m_Gamma[i++] = value;

            if (embedded_gamma.size() != 0)
            {
                m_EmbeddedGamma.emplace(stage_count);

                i = 0;

                for (const double value : embedded_gamma)
                    (*m_EmbeddedGamma)[i++] = value;
            }
        }

        [[nodiscard]] const nane::vector<>& alpha() const noexcept
        {
            return m_Alpha;
        }

        [[nodiscard]] const nane::matrix<>& beta() const noexcept
        {
            return m_Beta;
        }

        [[nodiscard]] const nane::vector<>& gamma() const noexcept
        {
            return m_Gamma;
        }

        [[nodiscard]] bool has_embedded_gamma() const noexcept
        {
            return m_EmbeddedGamma.has_value();
        }

        [[nodiscard]] const nane::vector<>& embedded_gamma() const
        {
            NANE_ENSURE(m_EmbeddedGamma.has_value(), "Butcher table has no embedded gamma coefficients");

            return *m_EmbeddedGamma;
        }

        [[nodiscard]] std::size_t stages() const noexcept
        {
            return m_Alpha.size();
        }

        [[nodiscard]] bool is_explicit() const noexcept
        {
            for (std::size_t i = 0; i < m_Beta.rows(); ++i)
            {
                for (std::size_t j = i; j < m_Beta.columns(); ++j)
                {
                    if (m_Beta(i, j) != 0.0)
                        return false;
                }
            }

            return true;
        }

    private:
        nane::vector<> m_Alpha;
        nane::matrix<> m_Beta;
        nane::vector<> m_Gamma;
        std::optional<nane::vector<>> m_EmbeddedGamma;
    };

} // namespace nane::ode
