#pragma once

#include "nane/geometry/grid.hpp"
#include "nane/numerics/pde/boundary_conditions.hpp"
#include "nane/numerics/pde/finite_difference.hpp"

#include <cstddef>

namespace nane::pde
{
    struct discrete_system
    {
        nane::matrix<> matrix;
        nane::vector<> rhs;
    };

    template <std::size_t Rows, std::size_t Columns>
    [[nodiscard]] discrete_system assemble(const nane::pde::finite_difference<Rows, Columns>& equation, const nane::geometry::grid<2>& grid,
                                           const nane::pde::boundaries<2>& boundary_conditions)
    {
        const auto& kernel = equation.kernel();

        const int center_row = static_cast<int>(kernel.rows() / 2);
        const int center_column = static_cast<int>(kernel.columns() / 2);

        const auto& x_axis = grid.axis(0);
        const auto& y_axis = grid.axis(1);

        const auto Nx = x_axis.size();
        const auto Ny = y_axis.size();

        const auto hx = *grid.spacing(0);
        const auto hy = *grid.spacing(1);

        const auto N = Nx * Ny;

        // clang-format off
	const auto is_dirichlet_boundary = [&](std::size_t i, std::size_t j)
        {
            return (i == 0 && is_dirichlet(boundary_conditions.lower(0))) || (i == x_axis.size() - 1 && is_dirichlet(boundary_conditions.upper(0))) ||
                   (j == 0 && is_dirichlet(boundary_conditions.lower(1))) || (j == y_axis.size() - 1 && is_dirichlet(boundary_conditions.upper(1)));
        };

        const auto index = [&](std::size_t i, std::size_t j)
	{ 
	    return j * x_axis.size() + i; 
	};
        // clang-format on

        nane::matrix A(N, N, 0.0);
        nane::vector b(N, 0.0);

        for (std::size_t j = 0; j < Ny; ++j)
        {
            for (std::size_t i = 0; i < Nx; ++i)
            {
                const auto row = index(i, j);

                // Dirichlet boundary nodes have zero derivative.
                if (is_dirichlet_boundary(i, j))
                    continue;

                for (std::size_t kj = 0; kj < kernel.columns(); ++kj)
                {
                    for (std::size_t ki = 0; ki < kernel.rows(); ++ki)
                    {
                        if (kernel(ki, kj) == 0.0)
                            continue;

                        const int di = static_cast<int>(ki) - center_row;
                        const int dj = static_cast<int>(kj) - center_column;

                        const int ti = static_cast<int>(i) + di;
                        const int tj = static_cast<int>(j) + dj;

                        const bool inside = ti >= 0 && ti < static_cast<int>(Nx) && tj >= 0 && tj < static_cast<int>(Ny);

                        if (inside)
                        {
                            A(row, index(static_cast<std::size_t>(ti), static_cast<std::size_t>(tj))) += kernel(ki, kj);
                            continue;
                        }

                        // Neumann X-axis
                        {
                            // Lower boundary.
                            if (ti < 0)
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.lower(0)).value;

                                A(row, index(1, j)) += kernel(ki, kj);
                                b[row] += kernel(ki, kj) * 2 * hx * q;
                            }

                            // Upper boundary.
                            else if (ti >= static_cast<int>(Nx))
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.upper(0)).value;

                                A(row, index(Nx - 2, j)) += kernel(ki, kj);
                                b[row] += kernel(ki, kj) * 2 * hx * q;
                            }
                        }

                        // Neumann Y-axis
                        {
                            // Lower boundary.
                            if (tj < 0)
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.lower(1)).value;

                                A(row, index(i, 1)) += kernel(ki, kj);
                                b[row] += kernel(ki, kj) * 2 * hy * q;
                            }

                            // Upper boundary.
                            else if (tj >= static_cast<int>(Ny))
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.upper(1)).value;

                                A(row, index(i, Ny - 2)) += kernel(ki, kj);
                                b[row] += kernel(ki, kj) * 2 * hy * q;
                            }
                        }
                    }
                }
            }
        }

        return discrete_system{A, b};
    }

} // namespace nane::pde
