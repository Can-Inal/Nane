#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/nonlinear/newton.hpp"
#include "nane/numerics/ode/butcher_table.hpp"

#include <Eigen/Core>
#include <cstddef>

namespace nane::ode
{
    namespace detail
    {

        /**
         * @brief Common Runge-Kutta implementation for vector-valued states.
         *
         * Scalar problems are converted to one-dimensional vector problems
         * by the public scalar overload.
         *
         * @tparam Derivative Type of the right-hand-side function.
         *
         * @param derivative Right-hand-side function.
         * @param initial_value Initial state.
         * @param time_grid Time discretization.
         * @param table Butcher table defining the Runge-Kutta method.
         *
         * @return Matrix whose columns contain the numerical states.
         */
        template <typename Derivative>
        [[nodiscard]] Eigen::MatrixXd runge_kutta_impl(const Derivative& derivative, const Eigen::VectorXd& initial_value,
                                                       const nane::geometry::uniform_grid<1>& time_grid, const nane::ode::butcher_table& table)
        {
            const auto& alpha = table.alpha();
            const auto& beta = table.beta();
            const auto& gamma = table.gamma();

            const Eigen::Index stage_count = table.stages();
            const Eigen::Index dimension = initial_value.size();

            const std::size_t count = time_grid.count(0);
            const double tau = time_grid.spacing(0);

            // stage matrix.
            // each column contains one Runge-Kutta stage vector.
            // for example, with dimension 3 and 4 stages we get a 3x4 matrix.
            Eigen::MatrixXd stage = Eigen::MatrixXd::Zero(dimension, stage_count);

            // solution matrix.
            // each column contains the numerical state at one grid point.
            // for example, with dimension 3 and 100 grid points we get a 3x100 matrix.
            Eigen::MatrixXd solution(dimension, static_cast<Eigen::Index>(count));

            // the first column of the solution is the initial state.
            // this corresponds to the first grid point.
            solution.col(0) = initial_value;

            // going over the time steps.
            for (std::size_t i = 0; i + 1 < count; ++i)
            {
                // index of the current grid point.
                const auto index = static_cast<Eigen::Index>(i);

                // time value of the current grid point.
                // for example t = 2s.
                const double current_time = time_grid.axis(0)[index];

                // state at the current grid point.
                // for example [x1(t), x2(t), x3(t)] for a 3 dimensional system.
                const Eigen::VectorXd current_state = solution.col(index);

                // explicit Runge-Kutta methods allow the stages to be computed one after another.
                if (table.is_explicit())
                {
                    // going over the Runge-Kutta stages.
                    for (std::size_t j = 0; j < static_cast<std::size_t>(stage_count); ++j)
                    {
                        // index of the current stage.
                        const auto stage_index = static_cast<Eigen::Index>(j);

                        // weighted sum of the previous stage vectors.
                        // the vector has the same dimension as the state.
                        Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                        // the first stage has no previous stages.
                        // after that we use all already computed stages.
                        if (j > 0)
                        {
                            // stage.leftCols(...) gives the already computed stage vectors.
                            // beta.row(...).head(...) gives the corresponding beta coefficients.
                            // their product computes sum_l beta(j,l) * k_l.
                            stage_sum = stage.leftCols(stage_index) * beta.row(stage_index).head(stage_index).transpose();
                        }

                        // computes
                        // k_j = f(t_n + alpha_j * tau,
                        //         x_n + tau * sum_l beta(j,l) * k_l)
                        stage.col(stage_index) = derivative(current_time + alpha[stage_index] * tau, current_state + tau * stage_sum);
                    }
                }
                else
                {
                    // implicit Runge-Kutta stages depend on each other.
                    // therefore all stages must be solved simultaneously.
                    // all stage vectors are flattened into one large vector.
                    const auto mapping = [&](const Eigen::VectorXd& current_stage)
                    {
                        // contains the new values of all flattened stage vectors.
                        Eigen::VectorXd next_stage(dimension * stage_count);

                        // going over all stages.
                        for (std::size_t j = 0; j < static_cast<std::size_t>(stage_count); ++j)
                        {
                            // index of the current stage.
                            const auto stage_index = static_cast<Eigen::Index>(j);

                            // weighted sum of all stage vectors.
                            Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                            // going over all stages because implicit methods can depend
                            // on previous, current and later stages.
                            for (std::size_t l = 0; l < static_cast<std::size_t>(stage_count); ++l)
                            {
                                const auto other_stage_index = static_cast<Eigen::Index>(l);

                                // current_stage stores all stage vectors in one flattened vector.
                                // segment(...) extracts the stage vector k_l.
                                stage_sum += beta(stage_index, other_stage_index) * current_stage.segment(other_stage_index * dimension, dimension);
                            }

                            // computes
                            // k_j = f(t_n + alpha_j * tau, x_n + tau * sum_l beta(j,l) * k_l)
                            //
                            // the new stage vector is written back into the corresponding
                            // segment of the flattened vector.
                            next_stage.segment(stage_index * dimension, dimension) =
                                derivative(current_time + alpha[stage_index] * tau, current_state + tau * stage_sum);
                        }

                        return next_stage;
                    };

                    // the implicit stage equations have the form
                    //
                    // K = mapping(K)
                    //
                    // Newton's method solves equations of the form
                    //
                    // residual(K) = 0
                    //
                    // therefore we use
                    //
                    // residual(K) = K - mapping(K).

                    // clang-format off
                    const auto residual = [&](const Eigen::VectorXd& current_stage) -> Eigen::VectorXd
                    { 
			return current_stage - mapping(current_stage); 
		    };
                    // clang-format on

                    // uses Newton's method to solve all implicit stage vectors simultaneously.
                    // the initial guess is a zero vector containing all flattened stages.
                    const Eigen::VectorXd solved_stage = nane::nonlinear::newton(residual, Eigen::VectorXd::Zero(dimension * stage_count));

                    // unpack the flattened stage vector back into the stage matrix.
                    // every segment becomes one column of the stage matrix.
                    for (std::size_t j = 0; j < static_cast<std::size_t>(stage_count); ++j)
                    {
                        const auto stage_index = static_cast<Eigen::Index>(j);

                        stage.col(stage_index) = solved_stage.segment(stage_index * dimension, dimension);
                    }
                }

                // computes the Runge-Kutta update
                //
                // x_{n+1} = x_n + tau * sum_j gamma_j * k_j
                //
                // stage * gamma computes the weighted sum of all stage vectors.
                solution.col(index + 1) = current_state + tau * stage * gamma;
            }

            return solution;
        }

    } // namespace detail
} // namespace nane::ode
