#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

#include <cstddef>

namespace nane::geometry
{
    /**
     * @brief Samples a scalar function on a one-dimensional grid.
     *
     * For grid points @f$x_i@f$, the sampled values are
     *
     * @f[
     * u_i = f(x_i).
     * @f]
     *
     * @param grid One-dimensional grid on which the function is sampled.
     * @param function Scalar function to sample.
     *
     * @return Vector containing the sampled function values.
     */
    [[nodiscard]] inline nane::vector<> sample(const nane::geometry::grid<1>& grid, auto&& function)
    {
        NANE_REQUIRE(grid.has_count(0), "grid axis must have prescribed count");

        const auto dimension = grid.axis(0).size();

        nane::vector<> result(dimension);

        for (std::size_t i = 0; i < dimension; ++i)
        {
            const double x = grid.axis(0)[i];
            result[i] = function(x);
        }

        return result;
    }

    /**
     * @brief Samples a scalar function on a two-dimensional grid.
     *
     * For grid points @f$(x_i, y_j)@f$, the sampled values are
     *
     * @f[
     * u_{ij} = f(x_i, y_j).
     * @f]
     *
     * @param grid Two-dimensional grid on which the function is sampled.
     * @param function Scalar function to sample.
     *
     * @return Matrix containing the sampled function values.
     */
    [[nodiscard]] inline nane::matrix<> sample(const nane::geometry::grid<2>& grid, auto&& function)
    {
        NANE_REQUIRE((grid.has_count(0) && grid.has_count(1)), "both grid axes must have prescribed counts");

        const auto rows = grid.axis(0).size();
        const auto columns = grid.axis(1).size();

        nane::matrix<> result(rows, columns);

        for (std::size_t i = 0; i < rows; ++i)
        {
            const double x = grid.axis(0)[i];

            for (std::size_t j = 0; j < columns; ++j)
            {
                const double y = grid.axis(1)[j];
                result(i, j) = function(x, y);
            }
        }

        return result;
    }

} // namespace nane::geometry
