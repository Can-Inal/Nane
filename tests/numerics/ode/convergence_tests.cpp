#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <nane/geometry/uniform_grid.hpp>
#include <nane/numerics/ode/crank_nicolson.hpp>
#include <nane/numerics/ode/explicit_euler.hpp>
#include <nane/numerics/ode/heun.hpp>
#include <nane/numerics/ode/implicit_euler.hpp>
#include <nane/numerics/ode/implicit_midpoint.hpp>
#include <nane/numerics/ode/ivp.hpp>
#include <nane/numerics/ode/modified_euler.hpp>
#include <nane/numerics/ode/rk4.hpp>

namespace
{
    template <typename Solver>
    double observed_convergence_order(Solver solver, std::size_t coarse_count, std::size_t fine_count)
    {
        // x' = -x
        const auto derivative = [](double time, double value) { return -value; };
        const auto problem = nane::ivp(derivative, 1.0);

        nane::uniform_grid<1> coarse_grid({
            {0.0, 1.0, coarse_count},
        });

        nane::uniform_grid<1> fine_grid({
            {0.0, 1.0, fine_count},
        });

        const Eigen::VectorXd coarse_solution = solver(problem, coarse_grid);
        const Eigen::VectorXd fine_solution = solver(problem, fine_grid);

        const double exact = std::exp(-1.0);

        const double coarse_value = coarse_solution[coarse_solution.size() - 1];
        const double fine_value = fine_solution[fine_solution.size() - 1];

        const double coarse_error = std::abs(coarse_value - exact);
        const double fine_error = std::abs(fine_value - exact);

        const double step_ratio = coarse_grid.spacing(0) / fine_grid.spacing(0);
        const double order = std::log(coarse_error / fine_error) / std::log(step_ratio);

        return order;
    }

} // namespace

TEST_CASE("Explicit Euler has first-order convergence", "[ode][explicit_euler][convergence]")
{
    const double order =
        observed_convergence_order([](const auto& problem, const auto& time_grid) { return nane::explicit_euler(problem, time_grid); }, 11, 21);

    REQUIRE(order == Catch::Approx(1.0).margin(0.1));
}

TEST_CASE("Implicit Euler has first-order convergence", "[ode][implicit_euler][convergence]")
{
    // clang-format off
    const double order = observed_convergence_order(
	    [](const auto& problem, const auto& time_grid) 
	    { 
		return nane::implicit_euler(problem, time_grid); 
	    }, 
	    11, 21
	);
    // clang-format on

    REQUIRE(order == Catch::Approx(1.0).margin(0.1));
}

TEST_CASE("Heun has second-order convergence", "[ode][heun][convergence]")
{
    // clang-format off
    const double order = observed_convergence_order(
	    [](const auto& problem, const auto& time_grid) 
	    { 
		return nane::heun(problem, time_grid); 
	    }, 
	    11, 21
	);
    // clang-format on

    REQUIRE(order == Catch::Approx(2.0).margin(0.1));
}

TEST_CASE("Modified Euler has second-order convergence", "[ode][modified_euler][convergence]")
{
    // clang-format off
    const double order = observed_convergence_order(
	    [](const auto& problem, const auto& time_grid) 
	    { 
		return nane::modified_euler(problem, time_grid); 
	    }, 
	    11, 21
	);
    // clang-format on

    REQUIRE(order == Catch::Approx(2.0).margin(0.1));
}

TEST_CASE("Implicit midpoint has second-order convergence", "[ode][implicit_midpoint][convergence]")
{
    // clang-format off
    const double order = observed_convergence_order(
	    [](const auto& problem, const auto& time_grid) 
	    { 
		return nane::implicit_midpoint(problem, time_grid); 
	    }, 
	    11, 21
	);
    // clang-format on 

    REQUIRE(order == Catch::Approx(2.0).margin(0.1));
}

TEST_CASE("Crank-Nicolson has second-order convergence", "[ode][crank_nicolson][convergence]")
{
    // clang-format off 
    const double order = observed_convergence_order(
	    [](const auto& problem, const auto& time_grid)
	    { 
		return nane::crank_nicolson(problem, time_grid); 
	    }, 
	    11, 21
	);
    // clang-format on 

    REQUIRE(order == Catch::Approx(2.0).margin(0.1));
}

TEST_CASE("RK4 has fourth-order convergence", "[ode][rk4][convergence]")
{
    // clang-format off 
    const double order = observed_convergence_order(
	    [](const auto& problem, const auto& time_grid) 
	    { 
		return nane::rk4(problem, time_grid); 
	    }, 
	    11, 21
	);
    // clang-format on 

    REQUIRE(order == Catch::Approx(4.0).margin(0.15));
}
