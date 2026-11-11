#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

namespace nane::approximation
{
    [[nodiscard]] inline double linear(const nane::vector<>& values, double position, const nane::grid<1>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "grid must have a prescribed point count");
        NANE_REQUIRE(*grid.count(0) >= 2, "grid must have at least 2 points");

        const auto& axis = grid.axis(0);

        NANE_REQUIRE(values.size() == axis.size(), "number of values must match grid points");

        const auto start = axis[0];
        const auto end = axis[axis.size() - 1];

        const auto formula = [&](std::size_t index)
        {
            const double x_0 = axis[index];
            const double x_1 = axis[index + 1];

            const double y_0 = values[index];
            const double y_1 = values[index + 1];

            return y_0 + ((position - x_0) / (x_1 - x_0)) * (y_1 - y_0);
        };

        std::size_t index = 0;

        while (index + 2 < axis.size() && axis[index + 1] < position)
        {
            ++index;
        }

        return formula(index);
    }

    [[nodiscard]] inline double linear(const nane::matrix<>& values, const nane::vector<2>& position, const nane::grid<2>& grid)
    {
        NANE_REQUIRE(grid.has_count(0), "x axis must have a prescribed point count");
        NANE_REQUIRE(*grid.count(0) >= 2, "x axis must have at least 2 points");

        NANE_REQUIRE(grid.has_count(1), "y axis must have a prescribed point count");
        NANE_REQUIRE(*grid.count(1) >= 2, "y axis must have at least 2 points");

        const auto& x_axis = grid.axis(0);
        const auto& y_axis = grid.axis(1);

        NANE_REQUIRE(values.columns() == x_axis.size(), "number of columns must match x axis points");
        NANE_REQUIRE(values.rows() == y_axis.size(), "number of rows must match y axis points");

        const auto x_start = x_axis[0];
        const auto x_end = x_axis[x_axis.size() - 1];

        const auto y_start = y_axis[0];
        const auto y_end = y_axis[y_axis.size() - 1];

        const auto x = position[0];
        const auto y = position[1];

        const auto formula = [&](std::size_t i, std::size_t j)
        {
            const double x_0 = x_axis[i];
            const double x_1 = x_axis[i + 1];

            const double y_0 = y_axis[j];
            const double y_1 = y_axis[j + 1];

            const auto t = (x - x_0) / (x_1 - x_0);
            const auto u = (y - y_0) / (y_1 - y_0);

            const auto f_00 = values(j, i);
            const auto f_01 = values(j + 1, i);
            const auto f_10 = values(j, i + 1);
            const auto f_11 = values(j + 1, i + 1);

            return (1 - t) * (1 - u) * f_00 + t * (1 - u) * f_10 + (1 - t) * u * f_01 + t * u * f_11;
        };

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

        return formula(i, j);
    }

} // namespace nane::approximation
