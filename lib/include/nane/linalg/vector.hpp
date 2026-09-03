#pragma once

#include <Eigen/Dense>
#include <cstddef>

namespace nane
{
    /**
     * @brief Creates a zero-initialized dynamic vector.
     *
     * @param size Number of vector elements.
     *
     * @return Dynamic column vector containing @p size zeros.
     */
    [[nodiscard]] inline Eigen::VectorXd vector(std::size_t size)
    {
        return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(size));
    }
} // namespace nane
