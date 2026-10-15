#pragma once

#include <Eigen/Core>

namespace nane
{
    /**
     * @ingroup calculus
     *
     * @brief Approximates the Jacobian matrix using central differences.
     *
     * Approximates the Jacobian of a vector-valued function
     *
     * @f[
     * \mathbf{F}
     * :
     * \mathbb{R}^n
     * \rightarrow
     * \mathbb{R}^m
     * @f]
     *
     * at the point @f$\mathbf{x}@f$.
     *
     * The Jacobian is defined as
     *
     * @f[
     * J_{\mathbf{F}}(\mathbf{x})
     * =
     * \begin{pmatrix}
     * \frac{\partial F_1}{\partial x_1}
     * &
     * \cdots
     * &
     * \frac{\partial F_1}{\partial x_n}
     * \\
     * \vdots
     * &
     * \ddots
     * &
     * \vdots
     * \\
     * \frac{\partial F_m}{\partial x_1}
     * &
     * \cdots
     * &
     * \frac{\partial F_m}{\partial x_n}
     * \end{pmatrix}.
     * @f]
     *
     * Each column is approximated using the central difference
     *
     * @f[
     * \frac{
     * \partial \mathbf{F}
     * }{
     * \partial x_j
     * }
     * \approx
     * \frac{
     * \mathbf{F}
     * \left(
     * \mathbf{x} + h\mathbf{e}_j
     * \right)
     * -
     * \mathbf{F}
     * \left(
     * \mathbf{x} - h\mathbf{e}_j
     * \right)
     * }{
     * 2h
     * },
     * @f]
     *
     * where @f$\mathbf{e}_j@f$ is the @f$j@f$-th Cartesian basis vector.
     *
     * @tparam Function Type of the vector-valued function.
     *
     * @param function Vector-valued function @f$\mathbf{F}(\mathbf{x})@f$.
     * @param current Point @f$\mathbf{x}@f$ at which the Jacobian is approximated.
     * @param step_size Finite-difference step @f$h@f$.
     *
     * @return Approximation of the Jacobian matrix.
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

} // namespace nane
