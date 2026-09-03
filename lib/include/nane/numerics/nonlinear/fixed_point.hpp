#pragma once

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>

namespace nane
{

    template <typename Function>
    [[nodiscard]] double fixed_point(Function&& function, double initial_guess, double tolerance = 1e-10, int max_iteration = 100)
    {
        auto current = initial_guess;

        for (auto i = 0; i < max_iteration; i++)
        {
            auto next = function(current);

            if (std::abs(next - current) < tolerance)
                return next;

            current = next;
        }

        throw std::runtime_error("fixed-point iteration did not converge.");
    }

    template <typename Function>
    [[nodiscard]] Eigen::VectorXd fixed_point(Function&& function, const Eigen::VectorXd& initial_guess, double tolerance = 1e-10,
                                              int max_iterations = 100)
    {
        Eigen::VectorXd current = initial_guess;

        for (auto i = 0; i < max_iterations; i++)
        {
            Eigen::VectorXd next = function(current);

            if ((next - current).norm() < tolerance)
                return next;

            current = next;
        }

        throw std::runtime_error("fixed-point iteration did not converge.");
    }

} // namespace nane
