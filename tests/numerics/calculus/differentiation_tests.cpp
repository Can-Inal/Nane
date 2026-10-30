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
