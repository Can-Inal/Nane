#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <nane/core/symbolic.hpp>
#include <nane/geometry/uniform_grid.hpp>
#include <nane/numerics/ode/explicit_euler.hpp>
#include <nane/numerics/ode/heun.hpp>
#include <nane/numerics/ode/implicit_euler.hpp>

TEST_CASE("Explicit Euler preserves a constant scalar solution", "[ode][explicit_euler]")
{
    nane::uniform_grid<1> time_grid({
        {0.0, 1.0, 11},
    });

    const auto derivative = nane::function(0.0);

    const auto solution = nane::explicit_euler(derivative, 5.3, time_grid);

    for (Eigen::Index i = 0; i < solution.size(); i++)
        REQUIRE(solution[i] == Catch::Approx(5.3).margin(1e-12));
}

TEST_CASE("Implicit Euler preserves a constant scalar solution", "[ode][implicit_euler]")
{
    nane::uniform_grid<1> time_grid({
        {0.0, 1.0, 11},
    });

    const auto derivative = nane::function(0.0);

    const auto solution = nane::implicit_euler(derivative, 5.3, time_grid);

    for (Eigen::Index i = 0; i < solution.size(); i++)
        REQUIRE(solution[i] == Catch::Approx(5.3).margin(1e-12));
}

TEST_CASE("Heun preserves a constant scalar solution", "[ode][heun]")
{
    nane::uniform_grid<1> time_grid({
        {0.0, 1.0, 11},
    });

    const auto derivative = nane::function(0.0);

    const auto solution = nane::heun(derivative, 5.3, time_grid);

    for (Eigen::Index i = 0; i < solution.size(); i++)
        REQUIRE(solution[i] == Catch::Approx(5.3).margin(1e-12));
}

TEST_CASE("Implicit Euler supports vector-valued systems", "[ode][implicit_euler][vector]")
{
    nane::uniform_grid<1> time_grid({
        {0.0, 1.0, 11},
    });

    const auto derivative = nane::system(0.0, 0.0);

    Eigen::VectorXd initial_value(2);
    initial_value << 2.0, 1.0;

    const auto solution = nane::implicit_euler(derivative, initial_value, time_grid);

    for (Eigen::Index i = 0; i < solution.cols(); i++)
        REQUIRE((solution.col(i) - initial_value).norm() < 1e-12);
}

TEST_CASE("Heun accepts a symbolic scalar derivative", "[ode][heun][symbolic]")
{
    nane::uniform_grid<1> time_grid({
        {0.0, 1.0, 101},
    });

    const auto [t, x] = nane::symbols<2>();

    const auto derivative = nane::function(-t * t * x);

    const auto solution = nane::heun(derivative, 1.0, time_grid);

    REQUIRE(solution[0] == Catch::Approx(1.0).margin(1e-12));

    REQUIRE(solution[solution.size() - 1] == Catch::Approx(std::exp(-1.0 / 3.0)).margin(1e-4));
}

TEST_CASE("Heun accepts a symbolic coupled system", "[ode][heun][vector][symbolic]")
{
    nane::uniform_grid<1> time_grid({
        {0.0, 1.0, 101},
    });

    const auto [t, x] = nane::symbols<2>();

    const auto derivative = nane::system(x[1], -x[0]);

    Eigen::VectorXd initial_value(2);
    initial_value << 1.0, 0.0;

    const auto solution = nane::heun(derivative, initial_value, time_grid);

    REQUIRE(solution(0, 0) == Catch::Approx(1.0).margin(1e-12));
    REQUIRE(solution(1, 0) == Catch::Approx(0.0).margin(1e-12));

    REQUIRE(solution(0, solution.cols() - 1) == Catch::Approx(std::cos(1.0)).margin(1e-4));

    REQUIRE(solution(1, solution.cols() - 1) == Catch::Approx(-std::sin(1.0)).margin(1e-4));
}
