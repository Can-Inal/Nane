#pragma once

#include "nane/geometry/uniform_grid.hpp"
#include "nane/linalg/vector.hpp"
#include "nane/numerics/nonlinear/fixed_point.hpp"

#include <Eigen/Core>
#include <cstddef>

namespace nane
{
    template <typename Function>
    [[nodiscard]] Eigen::VectorXd implicit_euler(Function&& function, double initial_value, const uniform_grid<1>& time_grid)
    {
        // x(n+1) = x(n) + tau*f(t(n+1), x(n+1))
        // x(n+1) = x(n) + tau*k
        // k = f(t(n+1), x(n+1))
        // k = f(t(n+1), x(n) + tau*k)

        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        Eigen::VectorXd solution = nane::vector(count);
        solution[0] = initial_value;

        for (auto i = 0; i < (int)count - 1; i++)
        {
            const auto next_time = time_grid.axis(0)[i + 1];                                                   // t(n+1)
            const auto mapping = [&](double stage) { return function(next_time, solution[i] + tau * stage); }; // f(t(n+1), x(n) + tau*k)

            const auto k = nane::fixed_point(mapping, 0.0);

            solution[i + 1] = solution[i] + tau * k;
        }

        return solution;
    }

    template <typename Function>
    [[nodiscard]] Eigen::MatrixXd implicit_euler(Function&& function, const Eigen::VectorXd& initial_value, const uniform_grid<1>& time_grid)
    {
        const std::size_t count = time_grid.count(0);
        const double tau = time_grid.spacing(0);

        const auto dimension = initial_value.size();

        Eigen::MatrixXd solution(dimension, count);
        solution.col(0) = initial_value;

        for (auto i = 0; i < (int)count - 1; i++)
        {
            const auto next_time = time_grid.axis(0)[i + 1];
            const auto mapping = [&](const Eigen::VectorXd& stage) { return function(next_time, solution.col(i) + tau * stage); };

            const auto k = nane::fixed_point(mapping, nane::vector(static_cast<std::size_t>(dimension)));

            solution.col(i + 1) = solution.col(i) + tau * k;
        }

        return solution;
    }
} // namespace nane
