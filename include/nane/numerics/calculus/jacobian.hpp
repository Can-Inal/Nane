#pragma once

#include <Eigen/Core>

namespace nane::calculus
{
    /**
     * Each column is approximated using a central difference,
     *
     * @f[
     * \frac{\partial \mathbf{F}}{\partial x_j}
     * \approx
     * \frac{
     * \mathbf{F}(\mathbf{x} + h\mathbf{e}_j)
     * -
     * \mathbf{F}(\mathbf{x} - h\mathbf{e}_j)
     * }{2h}.
     * @f]
     *
     * The approximation is second-order accurate with respect to
     * the finite-difference step,
     *
     * @f[
     * O(h^2).
     * @f]
     */
    template <typename Function>
    [[nodiscard]] Eigen::MatrixXd jacobian(Function&& function, const Eigen::VectorXd& current, double step_size)
    {
        const Eigen::VectorXd value = function(current);
        Eigen::MatrixXd result(value.size(), current.size());

        for (Eigen::Index j = 0; j < current.size(); ++j)
        {
            Eigen::VectorXd plus = current;
            Eigen::VectorXd minus = current;

            plus[j] += step_size;
            minus[j] -= step_size;

            result.col(j) = (function(plus) - function(minus)) / (2.0 * step_size);
        }

        return result;
    }

} // namespace nane::calculus
