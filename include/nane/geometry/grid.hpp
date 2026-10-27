#pragma once

#include "nane/error.hpp"
#include "nane/math/vector.hpp"

#include <array>
#include <cstddef>
#include <initializer_list>
#include <optional>

namespace nane::geometry
{
    /**
     * @ingroup geometry
     *
     * @brief Describes one grid axis.
     *
     * If a point count is specified, the interval is discretized uniformly.
     * Otherwise, only the interval bounds are prescribed.
     */
    struct grid_axis
    {
        double start;
        double stop;
        std::optional<std::size_t> count = std::nullopt;
    };

    /**
     * @ingroup geometry
     *
     * @brief Represents a Cartesian grid specification.
     *
     * @tparam Dimensions Number of dimensions.
     */
    template <std::size_t Dimensions>
    class grid
    {
    public:
        static_assert(Dimensions > 0, "grid must have at least one dimension");

        /**
         * @brief Constructs a one-dimensional grid without a prescribed count.
         */
        grid(double start, double stop)
            requires(Dimensions == 1)
        {
            m_Axes[0] = {
                start,
                stop,
                std::nullopt,
            };
        }

        /**
         * @brief Constructs a one-dimensional grid with a prescribed count.
         */
        grid(double start, double stop, std::size_t count)
            requires(Dimensions == 1)
        {
            NANE_REQUIRE(count >= 2, "grid count must be at least 2");

            m_Axes[0] = {
                start,
                stop,
                count,
            };

            m_Coordinates[0] = make_coordinates(start, stop, count);
        }

        /**
         * @brief Constructs a grid from its axis descriptions.
         */
        explicit grid(std::initializer_list<grid_axis> axes)
        {
            NANE_REQUIRE(axes.size() == Dimensions, "grid axis count must match its dimensions");

            std::size_t dimension = 0;

            for (const auto& axis : axes)
            {
                NANE_REQUIRE(!axis.count || *axis.count >= 2, "grid count must be at least 2");

                m_Axes[dimension] = axis;

                if (axis.count)
                {
                    m_Coordinates[dimension] = make_coordinates(axis.start, axis.stop, *axis.count);
                }

                ++dimension;
            }
        }

        [[nodiscard]] static constexpr std::size_t dimensions() noexcept
        {
            return Dimensions;
        }

        [[nodiscard]] double start(std::size_t dimension) const
        {
            return m_Axes.at(dimension).start;
        }

        [[nodiscard]] double stop(std::size_t dimension) const
        {
            return m_Axes.at(dimension).stop;
        }

        [[nodiscard]] bool has_count(std::size_t dimension) const
        {
            return m_Axes.at(dimension).count.has_value();
        }

        [[nodiscard]] std::optional<std::size_t> count(std::size_t dimension) const
        {
            return m_Axes.at(dimension).count;
        }

        [[nodiscard]] std::optional<double> spacing(std::size_t dimension) const
        {
            const auto& axis = m_Axes.at(dimension);

            if (!axis.count)
                return std::nullopt;

            return (axis.stop - axis.start) / static_cast<double>(*axis.count - 1);
        }

        [[nodiscard]] const nane::vector<>& axis(std::size_t dimension) const
        {
            NANE_ENSURE(m_Axes.at(dimension).count.has_value(), "grid axis has no prescribed point count");

            return m_Coordinates.at(dimension);
        }

    private:
        [[nodiscard]] static nane::vector<> make_coordinates(double start, double stop, std::size_t count)
        {
            nane::vector<> coordinates(count);

            const double spacing = (stop - start) / static_cast<double>(count - 1);

            for (std::size_t i = 0; i < count; ++i)
            {
                coordinates[i] = start + static_cast<double>(i) * spacing;
            }

            // Ensure that both prescribed interval bounds are represented
            // exactly despite floating-point rounding.
            coordinates[0] = start;
            coordinates[count - 1] = stop;

            return coordinates;
        }

        std::array<grid_axis, Dimensions> m_Axes;
        std::array<nane::vector<>, Dimensions> m_Coordinates;
    };

} // namespace nane::geometry
