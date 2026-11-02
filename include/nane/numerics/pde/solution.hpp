#pragma once

#include "nane/geometry/grid.hpp"
#include "nane/math/vector.hpp"

#include <cstddef>
#include <vector>

namespace nane::pde
{
    /**
     * @ingroup pde
     *
     * @brief Represents the numerical solution of a partial differential equation.
     *
     * Stores the spatial grid, time points, and corresponding spatial states.
     *
     * @tparam Dimension Number of spatial dimensions.
     * @tparam State Type representing the spatial state.
     */
    template <std::size_t Dimension, typename State>
    struct solution
    {
        nane::geometry::grid<Dimension> spatial_grid;
        nane::vector<> time;
        std::vector<State> values;
    };

} // namespace nane::pde
