#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <nane.hpp>
#include <numbers>

namespace
{
    void run_heat_equation(std::size_t spatial_points, std::size_t time_points)
    {
        const double kappa = 0.1;

        const nane::geometry::grid<2> spatial_grid{
            {0.0, 1.0, spatial_points},
            {0.0, 1.0, spatial_points},
        };

        const auto initial =
            nane::geometry::sample(spatial_grid, [](double x, double y) { return std::sin(std::numbers::pi * x) * std::sin(std::numbers::pi * y); });

        const auto equation = [&](const nane::geometry::grid<2>& grid)
        {
            const double h = *grid.spacing(0);

            return nane::pde::finite_difference{
                (kappa / (h * h)) *
                    nane::matrix<3>{
                        {0.0, 1.0, 0.0},
                        {1.0, -4.0, 1.0},
                        {0.0, 1.0, 0.0},
                    },
            };
        };

        const nane::pde::boundaries<2> boundary_conditions{
            {nane::pde::dirichlet{0.0}, nane::pde::dirichlet{0.0}},
            {nane::pde::dirichlet{0.0}, nane::pde::dirichlet{0.0}},
        };

        const auto problem = nane::pde::ibvp(equation, initial, boundary_conditions);

        const auto start = std::chrono::steady_clock::now();

        const auto solution = nane::pde::solve(problem, spatial_grid, nane::ode::methods::heun, nane::geometry::grid<1>{{0.0, 0.01, time_points}});

        const auto end = std::chrono::steady_clock::now();

        const auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << spatial_points << "x" << spatial_points << ", " << time_points << " time points: " << milliseconds << " ms\n";

        REQUIRE(solution.values.size() == time_points);
    }
} // namespace

TEST_CASE("Two-dimensional PDE solve performance", "[pde][performance]")
{
    SECTION("11x11")
    {
        run_heat_equation(11, 101);
    }

    SECTION("21x21")
    {
        run_heat_equation(21, 101);
    }

    SECTION("31x31")
    {
        run_heat_equation(31, 101);
    }

    SECTION("31x31")
    {
        run_heat_equation(31, 101);
    }
}
