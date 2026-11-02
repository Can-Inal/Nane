#pragma once

#include "nane/numerics/pde/boundary_conditions.hpp"

#include <cstddef>
#include <utility>

namespace nane::pde
{
    /**
     * @ingroup pde
     *
     * @brief Represents an initial-boundary value problem for a partial
     * differential equation.
     *
     * An initial-boundary value problem consists of an evolution equation
     *
     * @f[
     * \frac{\partial u}{\partial t} = F(t, u),
     * @f]
     *
     * together with an initial state and boundary conditions.
     *
     * The initial time is determined by the first point of the time grid
     * supplied to the numerical method. The spatial grid is supplied
     * separately to the PDE solver.
     *
     * @tparam Equation Type of the right-hand-side function.
     * @tparam InitialState Type of the initial spatial state.
     * @tparam Dimension Number of spatial dimensions.
     */
    template <typename Equation, typename InitialState, std::size_t Dimension>
    class ibvp
    {
    public:
        /**
         * @brief Constructs an initial-boundary value problem.
         *
         * @param equation Right-hand-side function of the partial differential equation.
         * @param initial_state Initial spatial state of the solution.
         * @param boundary_conditions Boundary conditions of the problem.
         */
        ibvp(Equation equation, InitialState initial_state, boundaries<Dimension> boundary_conditions)
            : m_Equation(std::move(equation)), m_InitialState(std::move(initial_state)), m_BoundaryConditions(std::move(boundary_conditions))
        {
        }

        /**
         * @brief Returns the right-hand-side function.
         *
         * @return Constant reference to the equation.
         */
        [[nodiscard]] const Equation& equation() const noexcept
        {
            return m_Equation;
        }

        /**
         * @brief Returns the initial spatial state.
         *
         * @return Constant reference to the initial state.
         */
        [[nodiscard]] const InitialState& initial_state() const noexcept
        {
            return m_InitialState;
        }

        /**
         * @brief Returns the boundary conditions.
         *
         * @return Constant reference to the boundary conditions.
         */
        [[nodiscard]] const boundaries<Dimension>& boundary_conditions() const noexcept
        {
            return m_BoundaryConditions;
        }

    private:
        Equation m_Equation;
        InitialState m_InitialState;
        boundaries<Dimension> m_BoundaryConditions;
    };

} // namespace nane::pde
