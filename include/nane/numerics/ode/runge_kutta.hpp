#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/nonlinear/fixed_point.hpp"

#include <Eigen/Core>
#include <cstddef>
#include <stdexcept>

namespace nane
{
    /**
     * @ingroup ode
     *
     * @brief Solves a scalar initial value problem using a Runge-Kutta method.
     *
     * Approximates
     *
     * @f[
     * \dot{x}(t) = f(t, x(t))
     * @f]
     *
     * using an @f$s@f$-stage Runge-Kutta method defined by the
     * Butcher coefficients @f$\alpha@f$, @f$\beta@f$, and @f$\gamma@f$.
     *
     * The stage values are defined by
     *
     * @f[
     * k_i
     * =
     * f\left(
     * t_n + \alpha_i \tau_n,
     * x_n + \tau_n \sum_{j=1}^{s} \beta_{ij} k_j
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
     * \tau_n
     * \sum_{i=1}^{s} \gamma_i k_i.
     * @f]
     *
     * If @f$\beta@f$ is strictly lower triangular, the stages are computed
     * explicitly. Otherwise, the coupled stage equations are solved using
     * fixed-point iteration.
     *
     * @tparam Function Type of the right-hand-side function.
     *
     * @param function Right-hand-side function @f$f(t,x)@f$.
     * @param initial_value Initial value.
     * @param time_grid Time discretization.
     * @param alpha Runge-Kutta stage-time coefficients.
     * @param beta Runge-Kutta stage coefficients.
     * @param gamma Runge-Kutta update coefficients.
     *
     * @return Numerical solution at all time-grid points.
     */
    template <typename Function>
    [[nodiscard]] Eigen::VectorXd runge_kutta(Function&& function, double initial_value, const nane::uniform_grid<1>& time_grid,
                                              const Eigen::VectorXd& alpha, const Eigen::MatrixXd& beta, const Eigen::VectorXd& gamma)
    {
        const double tolerance = 1e-12;

        // we determine the stage count.
        const int stage_count = (int)gamma.size();

        // checking the dimension match.
        if (stage_count == 0)
            throw std::invalid_argument("Runge-Kutta method must have at least one stage.");

        if (alpha.size() != stage_count)
            throw std::invalid_argument("Runge-Kutta alpha size must match stage count.");

        if (beta.rows() != stage_count || beta.cols() != stage_count)
            throw std::invalid_argument("Runge-Kutta beta matrix must be square and match stage count.");

        // we determine if the given butcher table would construct an explicit or implicit Runge-Kutta method.
        const bool is_explicit = beta.isLowerTriangular(tolerance) && beta.diagonal().isZero(tolerance);

        // create the stage vector.
        Eigen::VectorXd stage = Eigen::VectorXd::Zero(stage_count);

        // now the count of time steps and distance between time points.
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        // we construct the solution vector and give the initial value.
        Eigen::VectorXd solution = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(count));
        solution[0] = initial_value;

        for (auto i = 0; i < (int)count - 1; i++)
        {
            const auto current_time = time_grid.axis(0)[i];

            if (is_explicit)
            {
                // for an explicit Runge-Kutta method, every stage only depends on previously computed stages.
                for (auto j = 0; j < stage_count; j++)
                {
                    const auto stage_sum = beta.row(j).head(j).dot(stage.head(j));

                    stage[j] = function(current_time + alpha[j] * tau, solution[i] + tau * stage_sum);
                }
            }
            else
            {
                // for an implicit Runge-Kutta method, all stages have to be solved simultaneously.
                const auto mapping = [&](const Eigen::VectorXd& current_stage)
                {
                    Eigen::VectorXd next_stage(stage_count);

                    for (auto j = 0; j < stage_count; j++)
                    {
                        const auto stage_sum = beta.row(j).dot(current_stage);

                        next_stage[j] = function(current_time + alpha[j] * tau, solution[i] + tau * stage_sum);
                    }

                    return next_stage;
                };

                // we solve the coupled stage system by fixed-point iteration.
                stage = nane::fixed_point(mapping, Eigen::VectorXd::Zero(stage_count));
            }

            // then the update step.
            solution[i + 1] = solution[i] + tau * gamma.dot(stage);
        }

        return solution;
    }

    /**
     * @ingroup ode
     *
     * @brief Solves a vector-valued initial value problem using a Runge-Kutta method.
     *
     * Approximates
     *
     * @f[
     * \dot{\mathbf{x}}(t)
     * =
     * \mathbf{f}(t, \mathbf{x}(t))
     * @f]
     *
     * using an @f$s@f$-stage Runge-Kutta method defined by the
     * Butcher coefficients @f$\alpha@f$, @f$\beta@f$, and @f$\gamma@f$.
     *
     * The stage vectors are defined by
     *
     * @f[
     * \mathbf{k}_i
     * =
     * \mathbf{f}\left(
     * t_n + \alpha_i \tau_n,
     * \mathbf{x}_n
     * +
     * \tau_n
     * \sum_{j=1}^{s}
     * \beta_{ij}
     * \mathbf{k}_j
     * \right),
     * @f]
     *
     * followed by the update
     *
     * @f[
     * \mathbf{x}_{n+1}
     * =
     * \mathbf{x}_n
     * +
     * \tau_n
     * \sum_{i=1}^{s}
     * \gamma_i
     * \mathbf{k}_i.
     * @f]
     *
     * If @f$\beta@f$ is strictly lower triangular, the stages are computed
     * explicitly. Otherwise, the coupled stage equations are solved using
     * fixed-point iteration.
     *
     * Each column of the returned matrix contains the numerical state at
     * one time-grid point.
     *
     * @tparam Function Type of the right-hand-side function.
     *
     * @param function Right-hand-side function.
     * @param initial_value Initial state vector.
     * @param time_grid Time discretization.
     * @param alpha Runge-Kutta stage-time coefficients.
     * @param beta Runge-Kutta stage coefficients.
     * @param gamma Runge-Kutta update coefficients.
     *
     * @return Matrix whose columns contain the numerical states.
     */
    template <typename Function>
    [[nodiscard]] Eigen::MatrixXd runge_kutta(Function&& function, const Eigen::VectorXd& initial_value, const nane::uniform_grid<1>& time_grid,
                                              const Eigen::VectorXd& alpha, const Eigen::MatrixXd& beta, const Eigen::VectorXd& gamma)
    {
        const double tolerance = 1e-12;

        // we determine the stage count.
        const int stage_count = (int)gamma.size();

        // checking the dimension match.
        if (stage_count == 0)
            throw std::invalid_argument("Runge-Kutta method must have at least one stage.");

        if (alpha.size() != stage_count)
            throw std::invalid_argument("Runge-Kutta alpha size must match stage count.");

        if (beta.rows() != stage_count || beta.cols() != stage_count)
            throw std::invalid_argument("Runge-Kutta beta matrix must be square and match stage count.");

        // we determine if the given butcher table would construct an explicit or implicit Runge-Kutta method.
        const bool is_explicit = beta.isLowerTriangular(tolerance) && beta.diagonal().isZero(tolerance);

        // determine the dimension of the state vector.
        const auto dimension = initial_value.size();

        // create the stage matrix, where every column represents one stage vector.
        Eigen::MatrixXd stage(dimension, stage_count);

        // now the count of time steps and distance between time points.
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        // we construct the solution matrix and give the initial value.
        Eigen::MatrixXd solution(dimension, count);
        solution.col(0) = initial_value;

        for (auto i = 0; i < (int)count - 1; i++)
        {
            const auto current_time = time_grid.axis(0)[i];
            const Eigen::VectorXd current_state = solution.col(i);

            if (is_explicit)
            {
                // for an explicit Runge-Kutta method, every stage only depends on previously computed stages.
                for (auto j = 0; j < stage_count; j++)
                {
                    Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                    if (j > 0)
                        stage_sum = stage.leftCols(j) * beta.row(j).head(j).transpose();

                    stage.col(j) = function(current_time + alpha[j] * tau, current_state + tau * stage_sum);
                }
            }
            else
            {
                // for a vector-valued implicit Runge-Kutta method, all stage vectors are coupled.
                //
                // the fixed-point solver operates on a single vector, therefore the stage
                // vectors are stored consecutively:
                //
                // [ k_0 ]
                // [ k_1 ]
                // [ ... ]
                // [ k_s ]
                //
                // giving a vector of dimension * stage_count entries.

                const auto mapping = [&](const Eigen::VectorXd& current_stage)
                {
                    Eigen::VectorXd next_stage(dimension * stage_count);

                    for (auto j = 0; j < stage_count; j++)
                    {
                        Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                        for (auto l = 0; l < stage_count; l++)
                        {
                            stage_sum += beta(j, l) * current_stage.segment(l * dimension, dimension);
                        }

                        next_stage.segment(j * dimension, dimension) = function(current_time + alpha[j] * tau, current_state + tau * stage_sum);
                    }

                    return next_stage;
                };

                // we solve the coupled stage system by fixed-point iteration.
                const Eigen::VectorXd initial_stage = Eigen::VectorXd::Zero(dimension * stage_count);
                const Eigen::VectorXd fixed_stage = nane::fixed_point(mapping, initial_stage);

                // copy the flattened fixed-point result back into the stage matrix.
                for (auto j = 0; j < stage_count; j++)
                {
                    stage.col(j) = fixed_stage.segment(j * dimension, dimension);
                }
            }

            // then the update step.
            solution.col(i + 1) = current_state + tau * stage * gamma;
        }

        return solution;
    }

} // namespace nane
