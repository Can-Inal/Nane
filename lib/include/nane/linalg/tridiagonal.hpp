#pragma once

#include <Eigen/Dense>
#include <array>
#include <cstddef>

namespace nane
{
    /**
     * @brief Creates a dense tridiagonal matrix.
     *
     * The three supplied values define the lower, main, and upper
     * diagonals respectively.
     *
     * For example:
     *
     * ```cpp
     * auto matrix = nane::tridiagonal(
     *     { -1.0, 2.0, -1.0 },
     *     4
     * );
     * ```
     *
     * produces:
     *
     * ```text
     *  2 -1  0  0
     * -1  2 -1  0
     *  0 -1  2 -1
     *  0  0 -1  2
     * ```
     *
     * @param diagonals Values of the lower, main, and upper diagonals.
     * @param size Number of rows and columns.
     *
     * @return Dense square tridiagonal matrix.
     */
    [[nodiscard]] inline Eigen::MatrixXd tridiagonal(std::array<double, 3> diagonals, std::size_t size)
    {
        const auto eigen_size = static_cast<Eigen::Index>(size);

        Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(eigen_size, eigen_size);

        if (size == 0)
            return matrix;

        matrix.diagonal().setConstant(diagonals[1]);

        if (size > 1)
        {
            matrix.diagonal(-1).setConstant(diagonals[0]);
            matrix.diagonal(1).setConstant(diagonals[2]);
        }

        return matrix;
    }
} // namespace nane
