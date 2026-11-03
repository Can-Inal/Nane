#pragma once

#include "nane/geometry/grid.hpp"
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/solve.hpp"
#include "nane/numerics/pde/assemble.hpp"
#include "nane/numerics/pde/boundary_conditions.hpp"
#include "nane/numerics/pde/ibvp.hpp"
#include "nane/numerics/pde/solution.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace nane::pde
{
    template <typename Equation, typename InitialState>
    [[nodiscard]] auto solve(const nane::pde::ibvp<Equation, InitialState, 2>& problem, const nane::geometry::grid<2>& spatial_grid,
                             const nane::ode::butcher_table& time_method, const nane::geometry::grid<1>& time_grid, double tolerance = 1e-6)
    {
        const auto equation = problem.equation()(spatial_grid);
        const auto system = nane::pde::assemble(equation, spatial_grid, problem.boundary_conditions());

        auto initial = problem.initial_state();
        const auto initial_vector = initial.flatten_columns();

        nane::pde::apply_dirichlet(initial, problem.boundary_conditions());

        // clang-format off
        const auto system_function = [&](double, const nane::vector<>& state)
	{
	    return system.matrix * state + system.rhs;
	};
        // clang-format on

        const auto ode_problem = nane::ode::ivp(system_function, initial_vector);
        const auto ode_solution = nane::ode::solve(ode_problem, time_method, time_grid, tolerance);

        // each ODE solution column is being reshaped back into a PDE matrix.
        std::vector<nane::matrix<>> values;
        values.reserve(ode_solution.values.columns());

        for (std::size_t n = 0; n < ode_solution.values.columns(); ++n)
            values.push_back(nane::matrix<>::from_columns(ode_solution.values.column(n), initial.columns()));

        return nane::pde::solution<2, nane::matrix<>>{spatial_grid, ode_solution.time, std::move(values)};
    }

} // namespace nane::pde
