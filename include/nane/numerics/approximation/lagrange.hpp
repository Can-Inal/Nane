
#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

#include <algorithm>
#include <cstddef>

namespace nane::approximation
{
    /**
     * @brief Interpolates or extrapolates 1D samples with a local Lagrange polynomial.
     *
     * Selects consecutive samples around the interval bracketing @p position.
     * The window shifts near either endpoint to stay within the grid. For
     * extrapolation, the polynomial through the nearest end window is evaluated
     * at the original query coordinate.
     *
     * @tparam Points Number of samples in the window. Must be at least two.
     *                The polynomial has degree at most @c Points - 1.
     * @param values Sample values in the order of the grid coordinates.
     * @param position Coordinate at which to evaluate the polynomial.
     * @param grid One-dimensional grid with at least @p Points coordinates.
     * @return The interpolated or extrapolated value.
     * @pre Grid coordinates are finite and strictly increasing.
     * @throws std::invalid_argument If the grid lacks a point count, has too few
     *         points, or has a different number of points than @p values.
     */
    template <std::size_t Points = 4>
        requires(Points >= 2)
    [[nodiscard]] double lagrange(const nane::vector<>& values, double position, const nane::grid<1>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "grid must have a prescribed point count");
        NANE_REQUIRE(*grid.count(0) >= Points, "grid has too few points for the requested interpolation");

        const auto& axis = grid.axis(0);
        NANE_REQUIRE(values.size() == axis.size(), "number of values must match grid points");

        // Find the bracketing interval, stopping at the last interval.
        std::size_t index = 0;
        while (index + 2 < axis.size() && axis[index + 1] < position)
        {
            ++index;
        }

        // Place the window around that interval and keep it within the grid.
        constexpr std::size_t left_offset = (Points - 2) / 2;
        const std::size_t start = std::min(index >= left_offset ? index - left_offset : std::size_t{0}, axis.size() - Points);

        double result = 0.0;
        for (std::size_t i = 0; i < Points; ++i)
        {
            const double x_i = axis[start + i];
            double weight = 1.0;

            // Build the basis weight for sample i using the other nodes.
            for (std::size_t j = 0; j < Points; ++j)
            {
                if (j != i)
                {
                    const double x_j = axis[start + j];
                    weight *= (position - x_j) / (x_i - x_j);
                }
            }

            result += values[start + i] * weight;
        }

        return result;
    }

    /**
     * @brief Evaluates a local Lagrange interpolant on a two-dimensional grid.
     *
     * Applies the one-dimensional Lagrange interpolant along the x-axis for
     * each row, then interpolates those results along the y-axis. A query
     * outside either axis uses one-dimensional extrapolation along that axis.
     *
     * @tparam Points Number of interpolation points per axis. Must be at least two.
     * @param values Sample matrix. The entry at row j and column i corresponds
     *               to the grid point (x[i], y[j]).
     * @param position Query coordinates in (x, y) order.
     * @param grid Sampling grid with a prescribed point count on each axis.
     * @return The interpolated or extrapolated value at position.
     * @throws std::invalid_argument If either axis has no prescribed point count
     *         or fewer than Points coordinates, or if the matrix dimensions do
     *         not match the grid.
     */
    template <std::size_t Points = 4>
        requires(Points >= 2)
    [[nodiscard]] double lagrange(const nane::matrix<>& values, const nane::vector<2>& position, const nane::grid<2>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "x axis must have a prescribed point count");
        NANE_REQUIRE(grid.has_count(1), "y axis must have a prescribed point count");

        const auto& x_axis = grid.axis(0);
        const auto& y_axis = grid.axis(1);

        NANE_REQUIRE(x_axis.size() >= Points, "x axis has too few points for Lagrange interpolation");
        NANE_REQUIRE(y_axis.size() >= Points, "y axis has too few points for Lagrange interpolation");
        NANE_REQUIRE(values.columns() == x_axis.size(), "number of columns must match x axis points");
        NANE_REQUIRE(values.rows() == y_axis.size(), "number of rows must match y axis points");

        // Recreate each axis as a grid for the 1D overload.
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

            return lagrange<Points>(row_values, position[0], x_grid);
        };

        // Entry j is the interpolated value at (position[0], y_axis[j]).
        for (std::size_t row = 0; row < y_axis.size(); ++row)
        {
            intermediate[row] = interpolate_row(row);
        }

        // The second pass selects Points of these values along the y-axis.
        return lagrange<Points>(intermediate, position[1], y_grid);
    }

} // namespace nane::approximation
