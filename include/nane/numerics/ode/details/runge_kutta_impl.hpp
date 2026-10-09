#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/numerics/nonlinear/fixed_point.hpp"
#include "nane/numerics/ode/butcher_table.hpp"

#include <Eigen/Core>
#include <cstddef>

namespace nane
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
                                                       const nane::uniform_grid<1>& time_grid, const nane::butcher_table& table)
        {
            const auto& alpha = table.alpha();
            const auto& beta = table.beta();
            const auto& gamma = table.gamma();

            const Eigen::Index stage_count = table.stages();
            const Eigen::Index dimension = initial_value.size();

            const std::size_t count = time_grid.count(0);
            const double tau = time_grid.spacing(0);

            Eigen::MatrixXd stage = Eigen::MatrixXd::Zero(dimension, stage_count);

            Eigen::MatrixXd solution(dimension, static_cast<Eigen::Index>(count));

            solution.col(0) = initial_value;

            for (std::size_t i = 0; i + 1 < count; ++i)
            {
                const auto index = static_cast<Eigen::Index>(i);

                const double current_time = time_grid.axis(0)[index];
                const Eigen::VectorXd current_state = solution.col(index);

                if (table.is_explicit())
                {
                    for (std::size_t j = 0; j < static_cast<std::size_t>(stage_count); ++j)
                    {
                        const auto stage_index = static_cast<Eigen::Index>(j);

                        Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                        if (j > 0)
                        {
                            stage_sum = stage.leftCols(stage_index) * beta.row(stage_index).head(stage_index).transpose();
                        }

                        stage.col(stage_index) = derivative(current_time + alpha[stage_index] * tau, current_state + tau * stage_sum);
                    }
                }
                else
                {
                    // Implicit stages are solved simultaneously by flattening
                    // all stage vectors into one fixed-point vector.
                    const auto mapping = [&](const Eigen::VectorXd& current_stage)
                    {
                        Eigen::VectorXd next_stage(dimension * stage_count);

                        for (std::size_t j = 0; j < static_cast<std::size_t>(stage_count); ++j)
                        {
                            const auto stage_index = static_cast<Eigen::Index>(j);

                            Eigen::VectorXd stage_sum = Eigen::VectorXd::Zero(dimension);

                            for (std::size_t l = 0; l < static_cast<std::size_t>(stage_count); ++l)
                            {
                                const auto other_stage_index = static_cast<Eigen::Index>(l);

                                stage_sum += beta(stage_index, other_stage_index) * current_stage.segment(other_stage_index * dimension, dimension);
                            }

                            next_stage.segment(stage_index * dimension, dimension) =
                                derivative(current_time + alpha[stage_index] * tau, current_state + tau * stage_sum);
                        }

                        return next_stage;
                    };

                    const Eigen::VectorXd fixed_stage = nane::fixed_point(mapping, Eigen::VectorXd::Zero(dimension * stage_count));

                    for (std::size_t j = 0; j < static_cast<std::size_t>(stage_count); ++j)
                    {
                        const auto stage_index = static_cast<Eigen::Index>(j);

                        stage.col(stage_index) = fixed_stage.segment(stage_index * dimension, dimension);
                    }
                }

                solution.col(index + 1) = current_state + tau * stage * gamma;
            }

            return solution;
        }

    } // namespace detail
} // namespace nane
