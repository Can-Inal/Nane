#include "nane/geometry/sample.hpp"
#include "nane/numerics/calculus/differentiation/differentiate.hpp"
#include "nane/numerics/calculus/differentiation/stencils.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Scalar first derivative with different stencils", "[calculus][differentiation][scalar first derivative]")
{
    const auto test_function = [](double x) { return x * x * x; };
    const double point = 2.0;

    const double step_size = 1e-6;

    const auto backward_first = nane::calculus::differentiate(test_function, point, nane::calculus::stencils::backward_first(step_size));
    const auto central_first = nane::calculus::differentiate(test_function, point, nane::calculus::stencils::central_first(step_size));
    const auto forward_first = nane::calculus::differentiate(test_function, point, nane::calculus::stencils::forward_first(step_size));

    REQUIRE(backward_first == Catch::Approx(12.0).margin(1e-5));
    REQUIRE(central_first == Catch::Approx(12.0).margin(1e-5));
    REQUIRE(forward_first == Catch::Approx(12.0).margin(1e-5));
}

TEST_CASE("Scalar second derivative with central difference", "[calculus][differentiation][scalar second derivative]")
{
    const auto test_function = [](double x) { return x * x * x; };
    const double point = 2.0;

    const double step_size = 1e-4;

    const auto central_second = nane::calculus::differentiate(test_function, point, nane::calculus::stencils::central_second(step_size));

    REQUIRE(central_second == Catch::Approx(12.0).margin(1e-5));
}

TEST_CASE("Jacobian of vector-valued function", "[calculus][differentiation][jacobian]")
{
    const auto system = [](const nane::vector<2>& point)
    {
        nane::vector<2> f;

        f[0] = (point[0] * point[0]) + point[1];
        f[1] = point[0] * point[1];

        return f;
    };

    const nane::vector<2> point{2.0, 3.0};

    const double step_size = 1e-4;

    const auto result = nane::calculus::differentiate(system, point, nane::calculus::stencils::central_first(step_size));

    REQUIRE(result(0, 0) == Catch::Approx(4.0).margin(1e-5));
    REQUIRE(result(0, 1) == Catch::Approx(1.0).margin(1e-5));
    REQUIRE(result(1, 0) == Catch::Approx(3.0).margin(1e-5));
    REQUIRE(result(1, 1) == Catch::Approx(2.0).margin(1e-5));
}

TEST_CASE("Derivative of sampled one-dimensional data", "[calculus][differentiation][sampled][1d]")
{
    constexpr std::size_t N = 5;

    const nane::geometry::grid<1> grid{{0.0, 1.0, N}};

    // u(x) = x^2
    const auto values = nane::geometry::sample(grid, [](double x) { return x * x; });

    // u_xx = 2
    const auto result = nane::calculus::differentiate(values, grid, nane::calculus::stencils::central_second(1.0));

    REQUIRE(result[0] == Catch::Approx(0.0));

    for (std::size_t i = 1; i < N - 1; ++i)
        REQUIRE(result[i] == Catch::Approx(2.0).margin(1e-12));

    REQUIRE(result[N - 1] == Catch::Approx(0.0));
}

TEST_CASE("Partial derivative of sampled two-dimensional data along axis 0", "[calculus][differentiation][sampled][2d][axis 0]")
{
    constexpr std::size_t Nx = 5;
    constexpr std::size_t Ny = 4;

    const nane::geometry::grid<2> grid{
        {0.0, 1.0, Nx},
        {0.0, 1.0, Ny},
    };

    // u(x, y) = x^2 + 3y^2
    const auto values = nane::geometry::sample(grid, [](double x, double y) { return x * x + 3.0 * y * y; });

    // u_xx = 2
    const auto result = nane::calculus::differentiate(values, grid, nane::calculus::stencils::central_second(1.0), 0);

    for (std::size_t j = 0; j < Ny; ++j)
    {
        REQUIRE(result(0, j) == Catch::Approx(0.0));
        REQUIRE(result(Nx - 1, j) == Catch::Approx(0.0));
    }

    for (std::size_t i = 1; i < Nx - 1; ++i)
    {
        for (std::size_t j = 0; j < Ny; ++j)
            REQUIRE(result(i, j) == Catch::Approx(2.0).margin(1e-12));
    }
}

TEST_CASE("Partial derivative of sampled two-dimensional data along axis 1", "[calculus][differentiation][sampled][2d][axis 1]")
{
    constexpr std::size_t Nx = 5;
    constexpr std::size_t Ny = 4;

    const nane::geometry::grid<2> grid{
        {0.0, 1.0, Nx},
        {0.0, 1.0, Ny},
    };

    // u(x, y) = x^2 + 3y^2
    const auto values = nane::geometry::sample(grid, [](double x, double y) { return x * x + 3.0 * y * y; });

    // u_yy = 6
    const auto result = nane::calculus::differentiate(values, grid, nane::calculus::stencils::central_second(1.0), 1);

    for (std::size_t i = 0; i < Nx; ++i)
    {
        REQUIRE(result(i, 0) == Catch::Approx(0.0));
        REQUIRE(result(i, Ny - 1) == Catch::Approx(0.0));
    }

    for (std::size_t i = 0; i < Nx; ++i)
    {
        for (std::size_t j = 1; j < Ny - 1; ++j)
            REQUIRE(result(i, j) == Catch::Approx(6.0).margin(1e-12));
    }
}

TEST_CASE("Mixed derivative of sampled two-dimensional data", "[calculus][differentiation][sampled][2d][mixed]")
{
    constexpr std::size_t Nx = 5;
    constexpr std::size_t Ny = 5;

    const nane::geometry::grid<2> grid{
        {0.0, 1.0, Nx},
        {0.0, 1.0, Ny},
    };

    // u(x, y) = xy
    const auto values = nane::geometry::sample(grid, [](double x, double y) { return x * y; });

    // u_xy = 1
    const auto u_x = nane::calculus::differentiate(values, grid, nane::calculus::stencils::central_first(1.0), 0);
    const auto u_xy = nane::calculus::differentiate(u_x, grid, nane::calculus::stencils::central_first(1.0), 1);

    for (std::size_t i = 1; i < Nx - 1; ++i)
    {
        for (std::size_t j = 1; j < Ny - 1; ++j)
            REQUIRE(u_xy(i, j) == Catch::Approx(1.0).margin(1e-12));
    }
}
