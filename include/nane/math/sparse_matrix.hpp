#pragma once

#include "nane/error.hpp"
#include "nane/math/vector.hpp"

#include <Eigen/Sparse>
#include <cstddef>
#include <vector>

namespace nane
{
    /**
     * @brief Numerical sparse matrix.
     *
     * The matrix uses sparse storage and stores only its nonzero entries.
     *
     * @tparam Rows Compile-time row count or nane::dynamic.
     * @tparam Columns Compile-time column count or nane::dynamic.
     */
    template <std::size_t Rows = dynamic, std::size_t Columns = Rows>
    class sparse_matrix
    {
    public:
        /**
         * @brief Sparse matrix entry.
         */
        struct entry
        {
            std::size_t row;
            std::size_t column;
            double value;
        };

        sparse_matrix()
        {
            m_Data.resize(static_cast<Eigen::Index>(m_Rows), static_cast<Eigen::Index>(m_Columns));
        }

        /**
         * @brief Constructs a sparse matrix with the given dimensions.
         *
         * For fixed dimensions, the given shape must match the
         * compile-time dimensions.
         *
         * @param rows Number of rows.
         * @param columns Number of columns.
         */
        sparse_matrix(std::size_t rows, std::size_t columns)
        {
            set_shape(rows, columns);
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
         * @brief Returns the number of stored nonzero entries.
         */
        [[nodiscard]] std::size_t nonzeros() const noexcept
        {
            return static_cast<std::size_t>(m_Data.nonZeros());
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

        /**
         * @brief Sets the sparse matrix entries.
         *
         * Duplicate entries are summed.
         *
         * @param entries Matrix entries.
         */
        void set_entries(const std::vector<entry>& entries)
        {
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(entries.size());

            for (const auto& entry : entries)
            {
                NANE_REQUIRE(entry.row < rows(), "sparse matrix row index out of range");
                NANE_REQUIRE(entry.column < columns(), "sparse matrix column index out of range");

                if (entry.value == 0.0)
                    continue;

                triplets.emplace_back(static_cast<Eigen::Index>(entry.row), static_cast<Eigen::Index>(entry.column), entry.value);
            }

            m_Data.setFromTriplets(triplets.begin(), triplets.end());
            m_Data.makeCompressed();
        }

        /**
         * @brief Returns a matrix entry.
         *
         * Entries not stored in the sparse matrix are zero.
         */
        [[nodiscard]] double operator()(std::size_t row, std::size_t column) const
        {
            return m_Data.coeff(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column));
        }

        /**
         * @brief Multiplies the sparse matrix by a vector.
         */
        template <std::size_t Dimension>
            requires(Columns == dynamic || Dimension == dynamic || Columns == Dimension)
        [[nodiscard]] vector<Rows> operator*(const vector<Dimension>& value) const
        {
            NANE_REQUIRE(columns() == value.size(), "matrix and vector dimensions must match");

            Eigen::VectorXd input(static_cast<Eigen::Index>(value.size()));

            for (std::size_t i = 0; i < value.size(); ++i)
                input[static_cast<Eigen::Index>(i)] = value[i];

            const Eigen::VectorXd output = m_Data * input;

            vector<Rows> result(rows());

            for (std::size_t i = 0; i < rows(); ++i)
                result[i] = output[static_cast<Eigen::Index>(i)];

            return result;
        }

    private:
        void set_shape(std::size_t rows, std::size_t columns)
        {
            if constexpr (Rows != dynamic)
                NANE_REQUIRE(rows == Rows, "sparse matrix row count does not match");

            if constexpr (Columns != dynamic)
                NANE_REQUIRE(columns == Columns, "sparse matrix column count does not match");

            m_Rows = rows;
            m_Columns = columns;

            m_Data.resize(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(columns));
        }

    private:
        static constexpr bool s_Dynamic = Rows == dynamic || Columns == dynamic;

        std::size_t m_Rows = Rows == dynamic ? 0 : Rows;
        std::size_t m_Columns = Columns == dynamic ? 0 : Columns;

        Eigen::SparseMatrix<double, Eigen::RowMajor> m_Data;
    };

} // namespace nane
