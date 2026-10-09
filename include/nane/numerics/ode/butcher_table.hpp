#pragma once

#include <Eigen/Core>
#include <initializer_list>
#include <stdexcept>

namespace nane
{

    /**
     * @ingroup ode
     *
     * @brief Represents the Butcher coefficients of a Runge-Kutta method.
     *
     * An @f$s@f$-stage Runge-Kutta method is defined by the coefficients
     * @f$\alpha@f$, @f$\beta@f$, and @f$\gamma@f$.
     *
     * The stage values are given by
     *
     * @f[
     * k_i
     * =
     * f\left(
     * t_n + \alpha_i \tau,
     * x_n + \tau \sum_{j=1}^{s} \beta_{ij} k_j
     * \right),
     * @f]
     *
     * followed by the update
     *
     * @f[
     * x_{n+1}
     * =
     * x_n
     * +
     * \tau
     * \sum_{i=1}^{s}
     * \gamma_i k_i.
     * @f]
     *
     * The coefficients are commonly written as the Butcher tableau
     *
     * @f[
     * \begin{array}{c|c}
     * \alpha & \beta \\
     * \hline
     *        & \gamma^{T}
     * \end{array}.
     * @f]
     *
     * The table validates that the coefficient dimensions are consistent
     * when it is constructed.
     */
    class butcher_table
    {
    public:
        /**
         * @brief Constructs a Butcher table from its coefficients.
         *
         * The number of stages is determined by the size of @p alpha.
         * Consequently, @p beta must be a square matrix with the same
         * number of rows and columns, and @p gamma must contain the same
         * number of entries.
         *
         * For example, Heun's method can be represented as
         *
         * @code
         * const nane::butcher_table table{
         *     {0.0, 1.0},
         *     {
         *         {0.0, 0.0},
         *         {1.0, 0.0},
         *     },
         *     {0.5, 0.5},
         * };
         * @endcode
         *
         * @param alpha Stage-time coefficients.
         * @param beta Stage coupling coefficients.
         * @param gamma Solution weights.
         *
         * @throws std::invalid_argument If the table contains no stages.
         * @throws std::invalid_argument If the coefficient dimensions are
         * inconsistent.
         */
        butcher_table(std::initializer_list<double> alpha, std::initializer_list<std::initializer_list<double>> beta,
                      std::initializer_list<double> gamma)
        {
            const auto stage_count = static_cast<Eigen::Index>(alpha.size());

            if (stage_count == 0)
                throw std::invalid_argument("Butcher table must contain at least one stage.");

            if (static_cast<Eigen::Index>(beta.size()) != stage_count)
                throw std::invalid_argument("Butcher table beta row count must match stage count.");

            if (static_cast<Eigen::Index>(gamma.size()) != stage_count)
                throw std::invalid_argument("Butcher table gamma size must match stage count.");

            alpha_.resize(stage_count);
            beta_.resize(stage_count, stage_count);
            gamma_.resize(stage_count);

            Eigen::Index i = 0;

            for (const auto value : alpha)
                alpha_[i++] = value;

            i = 0;

            for (const auto& row : beta)
            {
                if (static_cast<Eigen::Index>(row.size()) != stage_count)
                    throw std::invalid_argument("Butcher table beta must be square.");

                Eigen::Index j = 0;

                for (const auto value : row)
                    beta_(i, j++) = value;

                ++i;
            }

            i = 0;

            for (const auto value : gamma)
                gamma_[i++] = value;
        }

        /**
         * @brief Returns the stage-time coefficients.
         *
         * @return Constant reference to @f$\alpha@f$.
         */
        [[nodiscard]] const Eigen::VectorXd& alpha() const noexcept
        {
            return alpha_;
        }

        /**
         * @brief Returns the stage coupling coefficients.
         *
         * @return Constant reference to @f$\beta@f$.
         */
        [[nodiscard]] const Eigen::MatrixXd& beta() const noexcept
        {
            return beta_;
        }

        /**
         * @brief Returns the solution weights.
         *
         * @return Constant reference to @f$\gamma@f$.
         */
        [[nodiscard]] const Eigen::VectorXd& gamma() const noexcept
        {
            return gamma_;
        }

        /**
         * @brief Returns the number of Runge-Kutta stages.
         *
         * @return Number of stages in the method.
         */
        [[nodiscard]] Eigen::Index stages() const noexcept
        {
            return alpha_.size();
        }

        /**
         * @brief Determines whether the Runge-Kutta method is explicit.
         *
         * A Runge-Kutta method is explicit when @f$\beta@f$ is strictly
         * lower triangular, that is,
         *
         * @f[
         * \beta_{ij} = 0
         * \qquad
         * \text{for } j \geq i.
         * @f]
         *
         * @return @c true if the method is explicit, otherwise @c false.
         */
        [[nodiscard]] bool is_explicit() const noexcept
        {
            for (Eigen::Index i = 0; i < beta_.rows(); ++i)
            {
                for (Eigen::Index j = i; j < beta_.cols(); ++j)
                {
                    if (beta_(i, j) != 0.0)
                        return false;
                }
            }

            return true;
        }

    private:
        Eigen::VectorXd alpha_;
        Eigen::MatrixXd beta_;
        Eigen::VectorXd gamma_;
    };

} // namespace nane
