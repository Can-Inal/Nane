#pragma once

#include "nane/error.hpp"

#include <array>
#include <cstddef>
#include <initializer_list>
#include <variant>

namespace nane::pde
{
    /**
     * @ingroup pde
     *
     * @brief Represents a Dirichlet boundary condition.
     *
     * A Dirichlet boundary condition prescribes the value of the solution
     * at the boundary.
     */
    struct dirichlet
    {
        double value;
    };

    /**
     * @ingroup pde
     *
     * @brief Represents a Neumann boundary condition.
     *
     * A Neumann boundary condition prescribes the normal derivative of the
     * solution at the boundary.
     */
    struct neumann
    {
        double value;
    };

    /**
     * @ingroup pde
     *
     * @brief Boundary condition type supporting Dirichlet and Neumann conditions.
     */
    using boundary_condition = std::variant<dirichlet, neumann>;

    /**
     * @ingroup pde
     *
     * @brief Stores lower and upper boundary conditions for each spatial dimension.
     *
     * For each spatial dimension, two boundary conditions are stored:
     * the lower boundary and the upper boundary.
     *
     * @tparam Dimension Number of spatial dimensions.
     */
    template <std::size_t Dimension>
    class boundaries
    {
    public:
        /**
         * @brief Constructs boundary conditions for each spatial dimension.
         *
         * Each inner initializer list must contain exactly two conditions,
         * corresponding to the lower and upper boundaries of one axis.
         *
         * @param conditions Boundary condition pairs for each spatial dimension.
         */
        explicit boundaries(std::initializer_list<std::initializer_list<boundary_condition>> conditions)
        {
            NANE_REQUIRE(conditions.size() == Dimension, "number of boundary pairs must match dimension");

            std::size_t axis = 0;

            for (const auto& pair : conditions)
            {
                NANE_REQUIRE(pair.size() == 2, "each dimension must have a lower and upper boundary condition");

                std::size_t side = 0;

                for (const auto& condition : pair)
                    m_Conditions[axis][side++] = condition;

                ++axis;
            }
        }

        /**
         * @brief Returns the lower boundary condition of an axis.
         *
         * @param axis Spatial axis.
         *
         * @return Lower boundary condition.
         */
        [[nodiscard]] const boundary_condition& lower(std::size_t axis) const
        {
            NANE_REQUIRE(axis < Dimension, "boundary axis is out of range");
            return m_Conditions[axis][0];
        }

        /**
         * @brief Returns the upper boundary condition of an axis.
         *
         * @param axis Spatial axis.
         *
         * @return Upper boundary condition.
         */
        [[nodiscard]] const boundary_condition& upper(std::size_t axis) const
        {
            NANE_REQUIRE(axis < Dimension, "boundary axis is out of range");
            return m_Conditions[axis][1];
        }

    private:
        std::array<std::array<boundary_condition, 2>, Dimension> m_Conditions;
    };

} // namespace nane::pde
