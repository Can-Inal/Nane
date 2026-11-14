
#pragma once

#include "nane/error.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <utility>

namespace nane
{
    /**
     * @brief Generates a fixed-size vector of uniformly distributed values.
     *
     * @tparam Dimension Number of vector components.
     * @param range Lower (inclusive) and upper (exclusive) sampling bounds.
     * @return A vector with one sample per component.
     * @throws std::invalid_argument If the interval has no positive finite width.
     */
    template <std::size_t Dimension>
        requires(Dimension != dynamic)
    [[nodiscard]] vector<Dimension> random(std::pair<double, double> range)
    {
        const auto [lower, upper] = range;
        NANE_REQUIRE(lower < upper && std::isfinite(upper - lower), "random range must have positive finite width");

        static thread_local std::mt19937 engine{std::random_device{}()};
        std::uniform_real_distribution<double> distribution(lower, upper);

        vector<Dimension> result;
        for (auto& value : result)
        {
            value = distribution(engine);
        }

        return result;
    }

    /**
     * @brief Generates a fixed-size matrix of uniformly distributed values.
     *
     * @tparam Rows Number of matrix rows.
     * @tparam Columns Number of matrix columns.
     * @param range Lower (inclusive) and upper (exclusive) sampling bounds.
     * @return A matrix with one sample per entry.
     * @throws std::invalid_argument If the interval has no positive finite width.
     */
    template <std::size_t Rows, std::size_t Columns>
        requires(Rows != dynamic && Columns != dynamic)
    [[nodiscard]] matrix<Rows, Columns> random(std::pair<double, double> range)
    {
        const auto [lower, upper] = range;
        NANE_REQUIRE(lower < upper && std::isfinite(upper - lower), "random range must have positive finite width");

        static thread_local std::mt19937 engine{std::random_device{}()};
        std::uniform_real_distribution<double> distribution(lower, upper);

        matrix<Rows, Columns> result;
        for (std::size_t row = 0; row < Rows; ++row)
        {
            for (std::size_t column = 0; column < Columns; ++column)
            {
                result(row, column) = distribution(engine);
            }
        }

        return result;
    }
} // namespace nane
