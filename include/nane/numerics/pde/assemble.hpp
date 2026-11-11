#pragma once

#include "nane/geometry/grid.hpp"
#include "nane/math/sparse_matrix.hpp"
#include "nane/numerics/pde/boundary_conditions.hpp"
#include "nane/numerics/pde/finite_difference.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace nane::pde
{
    struct discrete_system
    {
        nane::sparse_matrix<> matrix;
        nane::vector<> rhs;
    };

    // NOTE: Nane currently supports only 1 ghost node. (3x3 kernel)
    template <std::size_t Rows, std::size_t Columns>
    [[nodiscard]] discrete_system assemble(const nane::pde::finite_difference<Rows, Columns>& equation, const nane::grid<2>& grid,
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
            return (i == 0 && is_dirichlet(boundary_conditions.lower(0))) || (i == Nx - 1 && is_dirichlet(boundary_conditions.upper(0))) ||
                   (j == 0 && is_dirichlet(boundary_conditions.lower(1))) || (j == Ny - 1 && is_dirichlet(boundary_conditions.upper(1)));
        };

        const auto index = [&](std::size_t i, std::size_t j)
        {
            return j * Nx + i;
        };
        // clang-format on

        std::vector<nane::sparse_matrix<>::entry> entries;
        entries.reserve(N * kernel.size());

        nane::vector b(N, 0.0);

        for (std::size_t j = 0; j < Ny; ++j)
        {
            for (std::size_t i = 0; i < Nx; ++i)
            {
                const auto row = index(i, j);

                // Dirichlet boundary values are applied to the initial state in solve().
                // Zero derivative keeps them fixed during time integration.
                if (is_dirichlet_boundary(i, j))
                    continue;

                for (std::size_t kj = 0; kj < kernel.columns(); ++kj)
                {
                    for (std::size_t ki = 0; ki < kernel.rows(); ++ki)
                    {
                        const double coefficient = kernel(ki, kj);

                        if (coefficient == 0.0)
                            continue;

                        const int di = static_cast<int>(ki) - center_row;
                        const int dj = static_cast<int>(kj) - center_column;

                        const int ti = static_cast<int>(i) + di;
                        const int tj = static_cast<int>(j) + dj;

                        const bool inside = ti >= 0 && ti < static_cast<int>(Nx) && tj >= 0 && tj < static_cast<int>(Ny);

                        if (inside)
                        {
                            entries.push_back({row, index(static_cast<std::size_t>(ti), static_cast<std::size_t>(tj)), coefficient});
                            continue;
                        }

                        // TODO: support more ghost nodes for bigger kernels. (Can-Inal)
                        // TODO: support corners. (Can-Inal)

                        // Neumann X-axis.
                        {
                            if (ti < 0)
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.lower(0)).value;

                                entries.push_back({row, index(1, j), coefficient});
                                b[row] += coefficient * 2.0 * hx * q;
                            }
                            else if (ti >= static_cast<int>(Nx))
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.upper(0)).value;

                                entries.push_back({row, index(Nx - 2, j), coefficient});
                                b[row] += coefficient * 2.0 * hx * q;
                            }
                        }

                        // Neumann Y-axis.
                        {
                            if (tj < 0)
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.lower(1)).value;

                                entries.push_back({row, index(i, 1), coefficient});
                                b[row] += coefficient * 2.0 * hy * q;
                            }
                            else if (tj >= static_cast<int>(Ny))
                            {
                                const auto q = std::get<nane::pde::neumann>(boundary_conditions.upper(1)).value;

                                entries.push_back({row, index(i, Ny - 2), coefficient});
                                b[row] += coefficient * 2.0 * hy * q;
                            }
                        }
                    }
                }
            }
        }

        nane::sparse_matrix<> A(N, N);
        A.set_entries(entries);

        return discrete_system{
            std::move(A),
            std::move(b),
        };
    }

} // namespace nane::pde
