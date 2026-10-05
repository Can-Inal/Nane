#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/nonlinear/fixed_point.hpp"
#include "nane/numerics/ode/ivp.hpp"

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
     * @tparam Derivative Type of the right-hand-side function.
     *
     * @param problem Scalar initial value problem.
     * @param time_grid Time discretization.
     * @param alpha Runge-Kutta stage-time coefficients.
     * @param beta Runge-Kutta stage coefficients.
     * @param gamma Runge-Kutta update coefficients.
     *
     * @return Numerical solution at all time-grid points.
     */
    template <typename Derivative>
    [[nodiscard]] Eigen::VectorXd runge_kutta(const nane::ivp<Derivative, double>& problem, const nane::uniform_grid<1>& time_grid,
                                              const Eigen::VectorXd& alpha, const Eigen::MatrixXd& beta, const Eigen::VectorXd& gamma)
    {
        const double tolerance = 1e-12;

        // Extract the mathematical data of the IVP. The initial value
        // and derivative are owned by the IVP and remain valid for the
        // complete integration.
        const auto& derivative = problem.derivative();
        const double initial_value = problem.initial_value();

        // The number of Runge-Kutta stages is determined by gamma.
        const int stage_count = static_cast<int>(gamma.size());

        // Check that the supplied Butcher coefficients have compatible
        // dimensions.
        if (stage_count == 0)
            throw std::invalid_argument("Runge-Kutta method must have at least one stage.");

        if (alpha.size() != stage_count)
            throw std::invalid_argument("Runge-Kutta alpha size must match stage count.");

        if (beta.rows() != stage_count || beta.cols() != stage_count)
            throw std::invalid_argument("Runge-Kutta beta matrix must be square and match stage count.");

        // A Runge-Kutta method is explicit when beta is strictly lower
        // triangular. Otherwise the stage equations are coupled and must
        // be solved implicitly.
        const bool is_explicit = beta.isLowerTriangular(tolerance) && beta.diagonal().isZero(tolerance);

        // Stores the scalar stage values k_i.
        Eigen::VectorXd stage = Eigen::VectorXd::Zero(stage_count);

        // The time grid is uniform, so every integration step uses the
        // same step size tau.
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        // One scalar solution value is stored for every point of the
        // time grid.
        Eigen::VectorXd solution = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(count));

        solution[0] = initial_value;

        for (auto i = 0; i < static_cast<int>(count) - 1; ++i)
        {
            const auto current_time = time_grid.axis(0)[i];

            if (is_explicit)
            {
                // For an explicit Runge-Kutta method, stage j depends
                // only on stages 0, ..., j - 1, which have already been
                // computed.
                for (auto j = 0; j < stage_count; ++j)
                {
                    const auto stage_sum = beta.row(j).head(j).dot(stage.head(j));

                    stage[j] = derivative(current_time + alpha[j] * tau, solution[i] + tau * stage_sum);
                }
            }
            else
            {
                // For an implicit Runge-Kutta method, the stages depend
                // on each other. Therefore all stage equations are solved
                // simultaneously as a fixed-point problem.
                const auto mapping = [&](const Eigen::VectorXd& current_stage)
                {
                    Eigen::VectorXd next_stage(stage_count);

                    for (auto j = 0; j < stage_count; ++j)
                    {
                        const auto stage_sum = beta.row(j).dot(current_stage);

                        next_stage[j] = derivative(current_time + alpha[j] * tau, solution[i] + tau * stage_sum);
                    }

                    return next_stage;
                };

                stage = nane::fixed_point(mapping, Eigen::VectorXd::Zero(stage_count));
            }

            // Combine all stage values according to gamma to advance
            // the numerical solution by one time step.
            solution[i + 1] = solution[i] + tau * gamma.dot(stage);
        }

        return solution;
    }

    /**
     * @ingroup ode
     *
     * @brief Solves a vector-valued initial value problem using a
     * Runge-Kutta method.
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
     * @tparam Derivative Type of the right-hand-side system.
     *
     * @param problem Vector-valued initial value problem.
     * @param time_grid Time discretization.
     * @param alpha Runge-Kutta stage-time coefficients.
     * @param beta Runge-Kutta stage coefficients.
     * @param gamma Runge-Kutta update coefficients.
     *
     * @return Matrix whose columns contain the numerical states.
     */
    template <typename Derivative>
    [[nodiscard]] Eigen::MatrixXd runge_kutta(const nane::ivp<Derivative, Eigen::VectorXd>& problem, const nane::uniform_grid<1>& time_grid,
                                              const Eigen::VectorXd& alpha, const Eigen::MatrixXd& beta, const Eigen::VectorXd& gamma)
    {
        const double tolerance = 1e-12;

        // Extract the mathematical data of the IVP.
        const auto& derivative = problem.derivative();
        const auto& initial_value = problem.initial_value();

        // The number of Runge-Kutta stages is determined by gamma.
        const int stage_count = static_cast<int>(gamma.size());

        // Check that the supplied Butcher coefficients have compatible
        // dimensions.
        if (stage_count == 0)
            throw std::invalid_argument("Runge-Kutta method must have at least one stage.");

        if (alpha.size() != stage_count)
            throw std::invalid_argument("Runge-Kutta alpha size must match stage count.");

        if (beta.rows() != stage_count || beta.cols() != stage_count)
            throw std::invalid_argument("Runge-Kutta beta matrix must be square and match stage count.");

        // Determine whether the stage equations are explicit or coupled.
        const bool is_explicit = beta.isLowerTriangular(tolerance) && beta.diagonal().isZero(tolerance);

        // Dimension of the vector-valued state.
        const Eigen::Index dimension = initial_value.size();

        // Every column contains one Runge-Kutta stage vector.
        Eigen::MatrixXd stage(dimension, stage_count);

        // The grid is uniform, so the same step size is used everywhere.
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        // Every column of the solution matrix represents the state at one
        // point of the time grid.
        Eigen::MatrixXd solution(dimension, static_cast<Eigen::Index>(count));

        solution.col(0) = initial_value;

        for (auto i = 0; i < static_cast<int>(count) - 1; ++i)
        {
            const auto current_time = time_grid.axis(0)[i];

            const Eigen::VectorXd current_state = solution.col(i);

            if (is_explicit)
            {
                // In an explicit method, stage j depends only on the
                // previously computed stages.
                for (auto j = 0; j < stage_count; ++j)
                {
                    Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                    if (j > 0)
                    {
                        stage_sum = stage.leftCols(j) * beta.row(j).head(j).transpose();
                    }

                    stage.col(j) = derivative(current_time + alpha[j] * tau, current_state + tau * stage_sum);
                }
            }
            else
            {
                // For an implicit method, all stage vectors are coupled.
                //
                // fixed_point operates on one Eigen::VectorXd, so the
                // complete collection of stage vectors is flattened as
                //
                // [ k_0 ]
                // [ k_1 ]
                // [ ... ]
                // [ k_(s-1) ]
                //
                // with dimension * stage_count total entries.
                const auto mapping = [&](const Eigen::VectorXd& current_stage)
                {
                    Eigen::VectorXd next_stage(dimension * stage_count);

                    for (auto j = 0; j < stage_count; ++j)
                    {
                        Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                        for (auto l = 0; l < stage_count; ++l)
                        {
                            stage_sum += beta(j, l) * current_stage.segment(l * dimension, dimension);
                        }

                        next_stage.segment(j * dimension, dimension) = derivative(current_time + alpha[j] * tau, current_state + tau * stage_sum);
                    }

                    return next_stage;
                };

                // Start fixed-point iteration from zero stage vectors.
                const Eigen::VectorXd initial_stage = Eigen::VectorXd::Zero(dimension * stage_count);

                const Eigen::VectorXd fixed_stage = nane::fixed_point(mapping, initial_stage);

                // Restore the flattened fixed-point result to the stage
                // matrix used by the Runge-Kutta update.
                for (auto j = 0; j < stage_count; ++j)
                {
                    stage.col(j) = fixed_stage.segment(j * dimension, dimension);
                }
            }

            // Combine the stage vectors according to gamma to advance
            // the numerical state by one time step.
            solution.col(i + 1) = current_state + tau * stage * gamma;
        }

        return solution;
    }

} // namespace nane
