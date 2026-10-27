#pragma once

#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"

namespace nane::ode::methods
{
    /**
     * @brief Explicit Euler method.
     */
    inline const nane::ode::butcher_table explicit_euler{
        {0.0},
        {
            {0.0},
        },
        {1.0},
    };

    /**
     * @brief Implicit Euler method.
     */
    inline const nane::ode::butcher_table implicit_euler{
        {1.0},
        {
            {1.0},
        },
        {1.0},
    };

    /**
     * @brief Heun's second-order Runge-Kutta method.
     */
    inline const nane::ode::butcher_table heun{
        {0.0, 1.0},
        {
            {0.0, 0.0},
            {1.0, 0.0},
        },
        {0.5, 0.5},
    };

    /**
     * @brief Modified Euler method.
     */
    inline const nane::ode::butcher_table modified_euler{
        {0.0, 0.5},
        {
            {0.0, 0.0},
            {0.5, 0.0},
        },
        {0.0, 1.0},
    };

    /**
     * @brief Implicit midpoint method.
     */
    inline const nane::ode::butcher_table implicit_midpoint{
        {0.5},
        {
            {0.5},
        },
        {1.0},
    };

    /**
     * @brief Crank-Nicolson method.
     */
    inline const nane::ode::butcher_table crank_nicolson{
        {0.0, 1.0},
        {
            {0.0, 0.0},
            {0.5, 0.5},
        },
        {0.5, 0.5},
    };

    /**
     * @brief Classical fourth-order Runge-Kutta method.
     */
    inline const nane::ode::butcher_table rk4{
        {0.0, 0.5, 0.5, 1.0},
        {
            {0.0, 0.0, 0.0, 0.0},
            {0.5, 0.0, 0.0, 0.0},
            {0.0, 0.5, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0},
        },
        {
            1.0 / 6.0,
            1.0 / 3.0,
            1.0 / 3.0,
            1.0 / 6.0,
        },
    };

} // namespace nane::ode::methods
