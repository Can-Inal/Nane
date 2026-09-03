#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>

namespace nane
{
    /**
     * @brief Describes one uniformly spaced grid axis.
     *
     * A uniform axis is defined by a closed interval and the number of
     * points distributed evenly across that interval.
     *
     * The spacing is:
     *
     * @f[
     * h = \frac{\text{stop} - \text{start}}{\text{count} - 1}
     * @f]
     */
    struct uniform_axis
    {
        /**
         * @brief Lower bound of the axis.
         */
        double start;

        /**
         * @brief Upper bound of the axis.
         */
        double stop;

        /**
         * @brief Number of points along the axis.
         *
         * Must be at least 2.
         */
        std::size_t count;
    };

    /**
     * @brief Represents a uniform Cartesian grid.
     *
     * A uniform grid consists of one uniformly spaced axis for each spatial
     * dimension.
     *
     * Only the coordinate values of the individual axes are stored.
     * Cartesian grid points are implied by combinations of those axis values.
     *
     * Example:
     *
     * ```cpp
     * nane::uniform_grid<2> grid({
     *     { 0.0,  1.0, 6 },
     *     { -1.0, 1.0, 5 },
     * });
     * ```
     *
     * The first axis contains 6 points over `[0, 1]` and the second axis
     * contains 5 points over `[-1, 1]`.
     *
     * @tparam Dimensions Number of spatial dimensions.
     */
    template <std::size_t Dimensions>
    class uniform_grid
    {
        static_assert(Dimensions > 0, "uniform_grid must have at least one dimension");

    public:
        /**
         * @brief Constructs a uniform Cartesian grid.
         *
         * Each element of @p axes describes one spatial dimension.
         *
         * The number of supplied axes must match @p Dimensions.
         *
         * @param axes Axis descriptions.
         *
         * @throws std::invalid_argument If the number of supplied axes does
         * not match @p Dimensions.
         *
         * @throws std::invalid_argument If any axis contains fewer than
         * 2 points.
         */
        explicit uniform_grid(std::initializer_list<uniform_axis> axes)
        {
            if (axes.size() != Dimensions)
            {
                throw std::invalid_argument("uniform_grid axis count must match its dimensions");
            }

            std::size_t dimension = 0;

            for (const auto& axis : axes)
            {
                if (axis.count < 2)
                {
                    throw std::invalid_argument("uniform_grid count must be at least 2");
                }

                axes_[dimension] = axis;
                coordinates_[dimension] = Eigen::ArrayXd::LinSpaced(static_cast<Eigen::Index>(axis.count), axis.start, axis.stop);
                ++dimension;
            }
        }

        /**
         * @brief Returns the number of spatial dimensions.
         *
         * @return Number of dimensions represented by the grid.
         */
        [[nodiscard]] static constexpr std::size_t dimensions() noexcept
        {
            return Dimensions;
        }

        /**
         * @brief Returns the lower bound of an axis.
         *
         * @param dimension Zero-based axis index.
         *
         * @return Lower bound of the selected axis.
         *
         * @throws std::out_of_range If @p dimension is outside the valid
         * range.
         */
        [[nodiscard]] double start(std::size_t dimension) const
        {
            return axes_.at(dimension).start;
        }

        /**
         * @brief Returns the upper bound of an axis.
         *
         * @param dimension Zero-based axis index.
         *
         * @return Upper bound of the selected axis.
         *
         * @throws std::out_of_range If @p dimension is outside the valid
         * range.
         */
        [[nodiscard]] double stop(std::size_t dimension) const
        {
            return axes_.at(dimension).stop;
        }

        /**
         * @brief Returns the number of points along an axis.
         *
         * @param dimension Zero-based axis index.
         *
         * @return Number of points along the selected axis.
         *
         * @throws std::out_of_range If @p dimension is outside the valid
         * range.
         */
        [[nodiscard]] std::size_t count(std::size_t dimension) const
        {
            return axes_.at(dimension).count;
        }

        /**
         * @brief Returns the spacing between neighboring points on an axis.
         *
         * The spacing is computed as:
         *
         * @f[
         * h = \frac{b-a}{N-1}
         * @f]
         *
         * where `a` is the axis start value, `b` is the axis stop value,
         * and `N` is the number of points.
         *
         * @param dimension Zero-based axis index.
         *
         * @return Uniform spacing along the selected axis.
         *
         * @throws std::out_of_range If @p dimension is outside the valid
         * range.
         */
        [[nodiscard]] double spacing(std::size_t dimension) const
        {
            const auto& axis = axes_.at(dimension);

            return (axis.stop - axis.start) / static_cast<double>(axis.count - 1);
        }

        /**
         * @brief Returns the coordinates along an axis.
         *
         * The returned array contains all evenly spaced coordinate values
         * from the axis start value to the axis stop value, including both
         * endpoints.
         *
         * @param dimension Zero-based axis index.
         *
         * @return Constant reference to the coordinate array.
         *
         * @throws std::out_of_range If @p dimension is outside the valid
         * range.
         */
        [[nodiscard]] const Eigen::ArrayXd& axis(std::size_t dimension) const
        {
            return coordinates_.at(dimension);
        }

    private:
        std::array<uniform_axis, Dimensions> axes_;
        std::array<Eigen::ArrayXd, Dimensions> coordinates_;
    };
} // namespace nane
