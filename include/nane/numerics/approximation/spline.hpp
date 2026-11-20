#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/sparse_matrix.hpp"
#include "nane/math/vector.hpp"

#include <cstddef>

namespace nane::approximation
{
    /**
     * @brief Evaluates a natural cubic spline on a uniform one-dimensional grid.
     *
     * The spline has zero second derivative at both endpoints. Outside the
     * grid, it extrapolates using the first or last cubic segment.
     *
     * @param values Sample values in grid order, one per grid point.
     * @param position Coordinate at which to evaluate the spline.
     * @param grid Uniform grid with a prescribed count of at least two points.
     * @return The interpolated or extrapolated value at position.
     * @throws std::invalid_argument If the grid has no prescribed count, has
     *         fewer than two points, or has a different number of points from values.
     */
    [[nodiscard]] inline double spline(const nane::vector<>& values, double position, const nane::grid<1>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "grid must have a prescribed point count");

        const auto& axis = grid.axis(0);
        const auto n = axis.size();

        NANE_REQUIRE(n >= 2, "grid must have at least 2 points");
        NANE_REQUIRE(values.size() == n, "number of values must match grid points");

        // c[i] is the second derivative scaled by h² / 6.
        // The two endpoint entries are zero by the natural boundary conditions.
        nane::vector<> c(n, 0.0);
        const auto interior_count = n - 2;

        if (interior_count > 0)
        {
            const auto difference_matrix = nane::sparse_matrix<>::tridiagonal(n, {1.0, -2.0, 1.0});
            const auto differences = difference_matrix * values;
            const auto system = nane::sparse_matrix<>::tridiagonal(interior_count, {1.0, 4.0, 1.0});

            nane::vector<> diagonal(interior_count);
            nane::vector<> rhs(interior_count);

            // System row r corresponds to grid point r + 1.
            for (std::size_t row = 0; row < interior_count; ++row)
            {
                diagonal[row] = system(row, row);
                rhs[row] = differences[row + 1];
            }

            // Eliminate the lower diagonal.
            for (std::size_t row = 1; row < interior_count; ++row)
            {
                const double factor = system(row, row - 1) / diagonal[row - 1];
                diagonal[row] -= factor * system(row - 1, row);
                rhs[row] -= factor * rhs[row - 1];
            }

            // Recover the interior curvatures from right to left.
            for (std::size_t remaining = interior_count; remaining > 0; --remaining)
            {
                const std::size_t row = remaining - 1;
                double residual = rhs[row];

                if (row + 1 < interior_count)
                {
                    residual -= system(row, row + 1) * c[row + 2];
                }

                c[row + 1] = residual / diagonal[row];
            }
        }

        // Select the containing interval or an end interval for extrapolation.
        std::size_t index = 0;
        while (index + 2 < n && axis[index + 1] < position)
        {
            ++index;
        }

        const double t = (position - axis[index]) / (axis[index + 1] - axis[index]);
        const double a = 1.0 - t;

        return a * values[index] + t * values[index + 1] + (a * a * a - a) * c[index] + (t * t * t - t) * c[index + 1];
    }

    /**
     * @brief Evaluates a natural cubic spline on a uniform two-dimensional grid.
     *
     * Interpolates each row along the x-axis, then interpolates the resulting
     * values along the y-axis. Each pass applies natural endpoint conditions.
     * Outside either axis, that pass extrapolates using its first or last
     * cubic segment.
     *
     * @param values Sample matrix. The entry at row j and column i corresponds
     *               to the grid point (x[i], y[j]).
     * @param position Query coordinates in (x, y) order.
     * @param grid Uniform grid with a prescribed count of at least two points
     *             on each axis.
     * @return The interpolated or extrapolated value at position.
     * @throws std::invalid_argument If either axis has no prescribed point
     *         count or fewer than two points, or if the matrix dimensions
     *         do not match the grid.
     */
    [[nodiscard]] inline double spline(const nane::matrix<>& values, const nane::vector<2>& position, const nane::grid<2>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "x axis must have a prescribed point count");
        NANE_REQUIRE(grid.has_count(1), "y axis must have a prescribed point count");

        const auto& x_axis = grid.axis(0);
        const auto& y_axis = grid.axis(1);

        NANE_REQUIRE(x_axis.size() >= 2, "x axis must have at least 2 points");
        NANE_REQUIRE(y_axis.size() >= 2, "y axis must have at least 2 points");
        NANE_REQUIRE(values.columns() == x_axis.size(), "number of columns must match x axis points");
        NANE_REQUIRE(values.rows() == y_axis.size(), "number of rows must match y axis points");

        // Give the 1D overload a grid for each axis.
        const nane::grid<1> x_grid(grid.start(0), grid.stop(0), x_axis.size());
        const nane::grid<1> y_grid(grid.start(1), grid.stop(1), y_axis.size());

        nane::vector<> row_values(x_axis.size());
        nane::vector<> intermediate(y_axis.size());

        const auto interpolate_row = [&](std::size_t row)
        {
            // Reuse the same buffer to pass each row to the 1D overload.
            for (std::size_t column = 0; column < x_axis.size(); ++column)
            {
                row_values[column] = values(row, column);
            }

            return spline(row_values, position[0], x_grid);
        };

        // Entry j is the spline value at (position[0], y_axis[j]).
        for (std::size_t row = 0; row < y_axis.size(); ++row)
        {
            intermediate[row] = interpolate_row(row);
        }

        // Evaluate the intermediate values along y at the requested coordinate.
        return spline(intermediate, position[1], y_grid);
    }

} // namespace nane::approximation
