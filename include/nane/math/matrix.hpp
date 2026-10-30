#pragma once

#include "nane/error.hpp"
#include "nane/math/vector.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <vector>

namespace nane
{
    /**
     * @brief Numerical matrix.
     *
     * One template argument defines a square matrix:
     *
     * @code
     * nane::matrix<3>
     * @endcode
     *
     * Two template arguments define a rectangular matrix:
     *
     * @code
     * nane::matrix<2, 3>
     * @endcode
     *
     * @tparam Rows Compile-time row count or nane::dynamic.
     * @tparam Columns Compile-time column count or nane::dynamic.
     */
    template <std::size_t Rows = dynamic, std::size_t Columns = Rows>
    class matrix
    {
    public:
        matrix() = default;

        /**
         * @brief Constructs a matrix with the given dimensions.
         *
         * For fixed dimensions, the given shape must match the
         * compile-time dimensions.
         *
         * @param rows Number of rows.
         * @param columns Number of columns.
         */
        matrix(std::size_t rows, std::size_t columns)
        {
            set_shape(rows, columns);
        }

        /**
         * @brief Constructs a matrix with the given dimensions
         * and initial value.
         *
         * @param rows Number of rows.
         * @param columns Number of columns.
         * @param value Initial value of every matrix entry.
         */
        matrix(std::size_t rows, std::size_t columns, double value) : matrix(rows, columns)
        {
            fill(value);
        }

        /**
         * @brief Constructs a matrix from nested initializer lists.
         */
        matrix(std::initializer_list<std::initializer_list<double>> values)
        {
            const std::size_t rows = values.size();
            const std::size_t columns = rows == 0 ? 0 : values.begin()->size();

            for (const auto& row : values)
                NANE_REQUIRE(row.size() == columns, "matrix rows must have equal length");

            set_shape(rows, columns);

            std::size_t index = 0;

            for (const auto& row : values)
            {
                for (double value : row)
                    m_Data[index++] = value;
            }
        }

        /**
         * @brief Returns the number of rows.
         */
        [[nodiscard]] constexpr std::size_t rows() const noexcept
        {
            return m_Rows;
        }

        /**
         * @brief Returns the number of columns.
         */
        [[nodiscard]] constexpr std::size_t columns() const noexcept
        {
            return m_Columns;
        }

        /**
         * @brief Returns the number of matrix entries.
         */
        [[nodiscard]] constexpr std::size_t size() const noexcept
        {
            return m_Rows * m_Columns;
        }

        /**
         * @brief Returns the compile-time row count.
         */
        [[nodiscard]] static constexpr std::size_t static_rows() noexcept
        {
            return Rows;
        }

        /**
         * @brief Returns the compile-time column count.
         */
        [[nodiscard]] static constexpr std::size_t static_columns() noexcept
        {
            return Columns;
        }

        /**
         * @brief Returns whether the matrix has a dynamic dimension.
         */
        [[nodiscard]] static constexpr bool is_dynamic() noexcept
        {
            return s_Dynamic;
        }

        [[nodiscard]] double& operator()(std::size_t row, std::size_t column) noexcept
        {
            return m_Data[row * m_Columns + column];
        }

        [[nodiscard]] const double& operator()(std::size_t row, std::size_t column) const noexcept
        {
            return m_Data[row * m_Columns + column];
        }

        /**
         * @brief Adds two matrices.
         */
        [[nodiscard]] matrix operator+(const matrix& other) const
        {
            NANE_REQUIRE(rows() == other.rows() && columns() == other.columns(), "matrix dimensions must match");

            matrix result(rows(), columns());

            for (std::size_t i = 0; i < rows(); ++i)
            {
                for (std::size_t j = 0; j < columns(); ++j)
                    result(i, j) = (*this)(i, j) + other(i, j);
            }

            return result;
        }

        /**
         * @brief Subtracts two matrices.
         */
        [[nodiscard]] matrix operator-(const matrix& other) const
        {
            NANE_REQUIRE(rows() == other.rows() && columns() == other.columns(), "matrix dimensions must match");

            matrix result(rows(), columns());

            for (std::size_t i = 0; i < rows(); ++i)
            {
                for (std::size_t j = 0; j < columns(); ++j)
                    result(i, j) = (*this)(i, j) - other(i, j);
            }

            return result;
        }

        /**
         * @brief Multiplies the matrix by a scalar.
         */
        [[nodiscard]] matrix operator*(double scalar) const
        {
            matrix result(rows(), columns());

            for (std::size_t i = 0; i < rows(); ++i)
            {
                for (std::size_t j = 0; j < columns(); ++j)
                    result(i, j) = (*this)(i, j) * scalar;
            }

            return result;
        }

        /**
         * @brief Divides the matrix by a scalar.
         */
        [[nodiscard]] matrix operator/(double scalar) const
        {
            matrix result(rows(), columns());

            for (std::size_t i = 0; i < rows(); ++i)
            {
                for (std::size_t j = 0; j < columns(); ++j)
                    result(i, j) = (*this)(i, j) / scalar;
            }

            return result;
        }

        /**
         * @brief Multiplies the matrix by a vector.
         */
        template <std::size_t Dimension>
            requires(Columns == dynamic || Dimension == dynamic || Columns == Dimension)
        [[nodiscard]] vector<Rows> operator*(const vector<Dimension>& value) const
        {
            NANE_REQUIRE(columns() == value.size(), "matrix and vector dimensions must match");

            vector<Rows> result(rows(), 0.0);

            for (std::size_t i = 0; i < rows(); ++i)
            {
                for (std::size_t j = 0; j < columns(); ++j)
                    result[i] += (*this)(i, j) * value[j];
            }

            return result;
        }

        /**
         * @brief Sets every matrix entry to the same value.
         */
        void fill(double value)
        {
            std::fill(m_Data.begin(), m_Data.end(), value);
        }

        [[nodiscard]] double* data() noexcept
        {
            return m_Data.data();
        }

        [[nodiscard]] const double* data() const noexcept
        {
            return m_Data.data();
        }

        friend matrix operator*(double scalar, const matrix& value)
        {
            return value * scalar;
        }

    private:
        void set_shape(std::size_t rows, std::size_t columns)
        {
            if constexpr (Rows != dynamic)
                NANE_REQUIRE(rows == Rows, "matrix row count does not match");

            if constexpr (Columns != dynamic)
                NANE_REQUIRE(columns == Columns, "matrix column count does not match");

            m_Rows = rows;
            m_Columns = columns;

            if constexpr (s_Dynamic)
                m_Data.assign(m_Rows * m_Columns, 0.0);
        }

    private:
        static constexpr bool s_Dynamic = Rows == dynamic || Columns == dynamic;
        static constexpr std::size_t s_StaticSize = s_Dynamic ? 0 : Rows * Columns;

        using storage_type = std::conditional_t<s_Dynamic, std::vector<double>, std::array<double, s_StaticSize>>;

        std::size_t m_Rows = Rows == dynamic ? 0 : Rows;
        std::size_t m_Columns = Columns == dynamic ? 0 : Columns;

        storage_type m_Data{};
    };

    matrix(std::initializer_list<std::initializer_list<double>>) -> matrix<dynamic, dynamic>;
    matrix(std::size_t, std::size_t) -> matrix<dynamic, dynamic>;
    matrix(std::size_t, std::size_t, double) -> matrix<dynamic, dynamic>;

} // namespace nane
