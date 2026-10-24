#pragma once

#include <Eigen/Core>
#include <cmath>

namespace nane::detail
{
    /**
     * @brief Checks whether two scalar values are approximately equal.
     *
     * The comparison uses Eigen's default dummy precision for double
     * precision floating-point values.
     *
     * @param left First scalar value.
     * @param right Second scalar value.
     *
     * @return True if the values are approximately equal.
     */
    [[nodiscard]] inline bool is_equal(double left, double right)
    {
        return std::abs(left - right) < Eigen::NumTraits<double>::dummy_precision();
    }

    /**
     * @brief Checks whether two Eigen matrices or vectors are approximately equal.
     *
     * The comparison uses Eigen's default approximation tolerance.
     *
     * @tparam Left Type of the left Eigen expression.
     * @tparam Right Type of the right Eigen expression.
     *
     * @param left First matrix or vector.
     * @param right Second matrix or vector.
     *
     * @return True if the values are approximately equal.
     */
    template <typename Left, typename Right>
    [[nodiscard]] bool is_equal(const Eigen::MatrixBase<Left>& left, const Eigen::MatrixBase<Right>& right)
    {
        return left.isApprox(right);
    }

} // namespace nane::detail
