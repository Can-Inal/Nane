#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

#include <cstddef>

namespace nane::approximation
{
    /**
     * @brief Interpolates or extrapolates sampled values on a 1D grid.
     *
     * Uses the line through the two samples that bracket the query. For a
     * query outside the grid, extends the line through the nearest end interval.
     *
     * @param values Sampled values in grid order. Entry @c values[i] corresponds to
     *               @f$x_i = \mathrm{grid.axis}(0)[i]@f$.
     * @param position Coordinate at which to evaluate the line.
     * @param grid One-dimensional grid with at least two points.
     * @return The interpolated or extrapolated value.
     * @pre The grid coordinates are finite and strictly increasing.
     * @throws std::invalid_argument If the grid lacks a point count, has fewer than
     *         two points, or has a different number of points than @p values.
     */
    [[nodiscard]] inline double linear(const nane::vector<>& values, double position, const nane::grid<1>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "grid must have a prescribed point count");
        NANE_REQUIRE(*grid.count(0) >= 2, "grid must have at least 2 points");

        const auto& axis = grid.axis(0);
        NANE_REQUIRE(values.size() == axis.size(), "number of values must match grid points");

        std::size_t index = 0;
        while (index + 2 < axis.size() && axis[index + 1] < position)
        {
            ++index;
        }

        const double x_0 = axis[index];
        const double x_1 = axis[index + 1];
        const double y_0 = values[index];
        const double y_1 = values[index + 1];

        return y_0 + ((position - x_0) / (x_1 - x_0)) * (y_1 - y_0);
    }

    /**
     * @brief Interpolates or extrapolates sampled values on a 2D grid.
     *
     * Uses bilinear interpolation through the four corners of the selected
     * grid cell. The interval on each axis is selected independently. If a
     * coordinate lies outside its axis, uses the nearest end interval and
     * evaluates the same bilinear polynomial at the original coordinates.
     *
     * @param values Sampled values with X on rows and Y on columns, so
     *               @c values(i,j) corresponds to @f$f(x_i,y_j)@f$.
     * @param position Query coordinates in X, Y order.
     * @param grid Two-dimensional grid with at least two points per axis.
     * @return The interpolated or extrapolated value.
     * @pre The coordinates on each axis are finite and strictly increasing.
     * @throws std::invalid_argument If either axis lacks a point count, has fewer
     *         than two points, or has a size that disagrees with @p values.
     */
    [[nodiscard]] inline double linear(const nane::matrix<>& values, const nane::vector<2>& position, const nane::grid<2>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "x axis must have a prescribed point count");
        NANE_REQUIRE(*grid.count(0) >= 2, "x axis must have at least 2 points");

        NANE_REQUIRE(grid.has_count(1), "y axis must have a prescribed point count");
        NANE_REQUIRE(*grid.count(1) >= 2, "y axis must have at least 2 points");

        const auto& x_axis = grid.axis(0);
        const auto& y_axis = grid.axis(1);

        NANE_REQUIRE(values.rows() == x_axis.size(), "number of rows must match x axis points");
        NANE_REQUIRE(values.columns() == y_axis.size(), "number of columns must match y axis points");

        const double x = position[0];
        const double y = position[1];

        std::size_t i = 0;
        std::size_t j = 0;

        while (i + 2 < x_axis.size() && x_axis[i + 1] < x)
        {
            ++i;
        }

        while (j + 2 < y_axis.size() && y_axis[j + 1] < y)
        {
            ++j;
        }

        const double x_0 = x_axis[i];
        const double x_1 = x_axis[i + 1];
        const double y_0 = y_axis[j];
        const double y_1 = y_axis[j + 1];

        const double t = (x - x_0) / (x_1 - x_0);
        const double u = (y - y_0) / (y_1 - y_0);

        const double f_00 = values(i, j);
        const double f_01 = values(i, j + 1);
        const double f_10 = values(i + 1, j);
        const double f_11 = values(i + 1, j + 1);

        return (1.0 - t) * (1.0 - u) * f_00 + t * (1.0 - u) * f_10 + (1.0 - t) * u * f_01 + t * u * f_11;
    }
} // namespace nane::approximation
