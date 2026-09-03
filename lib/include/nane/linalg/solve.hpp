#pragma once

#include <Eigen/Dense>
#include <stdexcept>

namespace nane
{
    /**
     * @brief Solves a dense linear system.
     *
     * Solves the linear system
     *
     * @f[
     * Ax = b
     * @f]
     *
     * using an LDLT decomposition.
     *
     * @param matrix Square system matrix.
     * @param rhs Right-hand-side vector.
     *
     * @return Solution vector.
     *
     * @throws std::invalid_argument If @p matrix is not square.
     * @throws std::invalid_argument If the size of @p rhs does not match
     * the number of rows of @p matrix.
     * @throws std::runtime_error If the decomposition or solve fails.
     */
    [[nodiscard]] inline Eigen::VectorXd solve(const Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs)
    {
        if (matrix.rows() != matrix.cols())
            throw std::invalid_argument("solve requires a square matrix");

        if (matrix.rows() != rhs.size())
            throw std::invalid_argument("solve matrix and vector dimensions must match");

        const auto decomposition = matrix.ldlt();

        if (decomposition.info() != Eigen::Success)
            throw std::runtime_error("solve decomposition failed");

        Eigen::VectorXd solution = decomposition.solve(rhs);

        if (decomposition.info() != Eigen::Success)
            throw std::runtime_error("solve failed");

        return solution;
    }
} // namespace nane
