#pragma once

#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace nane::viz
{
    /**
     * @brief Represents a two-dimensional plot or collection of subplots.
     *
     * Multiple data series can be added to each subplot.
     *
     * The visualization backend is implemented using VTK.
     */
    class plot
    {
    public:
        /**
         * @brief Constructs an empty plot.
         *
         * The plot initially contains a single subplot.
         */
        plot();

        /**
         * @brief Destroys the plot.
         */
        ~plot();

        plot(const plot&) = delete;
        plot& operator=(const plot&) = delete;
        plot(plot&&) noexcept;
        plot& operator=(plot&&) noexcept;

        /**
         * @brief Selects a subplot in a subplot grid.
         *
         * The subplot index is zero-based and follows row-major order,
         * starting from the upper-left subplot.
         *
         * For example, a 2-by-2 grid is indexed as
         *
         * @code
         * 0  1
         * 2  3
         * @endcode
         *
         * Subsequent calls to add_line(), add_points(), title(),
         * x_label(), and y_label() operate on the selected subplot.
         *
         * The subplot grid must be configured before data series are added.
         *
         * @param rows Number of subplot rows.
         * @param columns Number of subplot columns.
         * @param index Zero-based subplot index.
         *
         * @throws std::invalid_argument If the dimensions are zero or the
         * subplot index is outside the grid.
         *
         * @throws std::logic_error If an existing subplot layout is changed
         * after data series have already been added.
         */
        void subplot(std::size_t rows, std::size_t columns, std::size_t index);

        /**
         * @brief Adds a point series to the active subplot.
         *
         * @param x Horizontal coordinate values.
         * @param y Vertical coordinate values.
         * @param label Optional legend label.
         *
         * @throws std::invalid_argument If @p x and @p y have different sizes.
         */
        void add_points(std::span<const double> x, std::span<const double> y, std::string_view label = {});

        /**
         * @brief Adds a line series to the active subplot.
         *
         * @param x Horizontal coordinate values.
         * @param y Vertical coordinate values.
         * @param label Optional legend label.
         *
         * @throws std::invalid_argument If @p x and @p y have different sizes.
         */
        void add_line(std::span<const double> x, std::span<const double> y, std::string_view label = {});

        /**
         * @brief Adds a point series from indexable containers or expressions.
         *
         * The inputs must provide size() and operator[].
         *
         * This allows containers and views such as std::vector, Eigen vectors,
         * and Eigen row or column expressions to be plotted directly.
         *
         * @tparam XValues Horizontal coordinate type.
         * @tparam YValues Vertical coordinate type.
         *
         * @param x Horizontal coordinate values.
         * @param y Vertical coordinate values.
         * @param label Optional legend label.
         *
         * @throws std::invalid_argument If @p x and @p y have different sizes.
         */
        template <typename XValues, typename YValues>
        void add_points(const XValues& x, const YValues& y, std::string_view label = {})
        {
            if (x.size() != y.size())
                throw std::invalid_argument("plot coordinate sizes must match");

            auto x_values = copy_values(x);
            auto y_values = copy_values(y);

            add_points(std::span<const double>(x_values), std::span<const double>(y_values), label);
        }

        /**
         * @brief Adds a line series from indexable containers or expressions.
         *
         * The inputs must provide size() and operator[].
         *
         * This allows containers and views such as std::vector, Eigen vectors,
         * and Eigen row or column expressions to be plotted directly.
         *
         * @tparam XValues Horizontal coordinate type.
         * @tparam YValues Vertical coordinate type.
         *
         * @param x Horizontal coordinate values.
         * @param y Vertical coordinate values.
         * @param label Optional legend label.
         *
         * @throws std::invalid_argument If @p x and @p y have different sizes.
         */
        template <typename XValues, typename YValues>
        void add_line(const XValues& x, const YValues& y, std::string_view label = {})
        {
            if (x.size() != y.size())
                throw std::invalid_argument("plot coordinate sizes must match");

            auto x_values = copy_values(x);
            auto y_values = copy_values(y);

            add_line(std::span<const double>(x_values), std::span<const double>(y_values), label);
        }

        /**
         * @brief Sets the title of the active subplot.
         *
         * @param value Plot title.
         */
        void title(std::string_view value);

        /**
         * @brief Sets the horizontal axis label of the active subplot.
         *
         * @param value Horizontal axis label.
         */
        void x_label(std::string_view value);

        /**
         * @brief Sets the vertical axis label of the active subplot.
         *
         * @param value Vertical axis label.
         */
        void y_label(std::string_view value);

        /**
         * @brief Saves the complete plot window as a PNG image.
         *
         * All subplots are included in the saved image.
         *
         * If no file extension is specified, ".png" is appended
         * automatically.
         *
         * @param filename Output file name.
         *
         * @throws std::invalid_argument If the file name is empty or uses an
         * unsupported file extension.
         */
        void save(std::string_view filename);

        /**
         * @brief Displays the complete plot window.
         *
         * Opens the visualization window and starts the VTK interaction loop.
         *
         * This function blocks until the visualization window is closed.
         */
        void show();

    private:
        template <typename Values>
        [[nodiscard]] static std::vector<double> copy_values(const Values& values)
        {
            std::vector<double> result;
            result.reserve(static_cast<std::size_t>(values.size()));

            for (decltype(values.size()) i = 0; i < values.size(); ++i)
                result.push_back(static_cast<double>(values[i]));

            return result;
        }

    private:
        struct implementation;

        std::unique_ptr<implementation> implementation_;
    };
} // namespace nane::viz
