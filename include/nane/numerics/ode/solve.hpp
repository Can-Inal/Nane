#pragma once

#include "nane/geometry/grid.hpp"
#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/runge_kutta.hpp"

#include <utility>

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Solves an initial value problem using a Runge-Kutta method.
     *
     * The supplied Butcher table defines the Runge-Kutta method used
     * to advance the solution in time.
     *
     * A prescribed point count in the time grid produces a fixed time
     * grid. Without a point count, an embedded Butcher table may be used
     * for adaptive time stepping.
     *
     * @tparam Problem Type of the initial value problem.
     *
     * @param problem Initial value problem to solve.
     * @param method Butcher table defining the Runge-Kutta method.
     * @param time_grid Time interval or prescribed time grid.
     * @param tolerance Error tolerance used for adaptive time stepping.
     *
     * @return Numerical solution of the initial value problem.
     */
    template <typename Problem>
    [[nodiscard]] auto solve(Problem&& problem, const nane::ode::butcher_table& method, const nane::geometry::grid<1>& time_grid,
                             double tolerance = 1e-6)
    {
        return nane::ode::runge_kutta(std::forward<Problem>(problem), method, time_grid, tolerance);
    }

} // namespace nane::ode
