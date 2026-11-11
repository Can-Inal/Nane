#include "nane/geometry/sample.hpp"
#include "nane/numerics/approximation/linear.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Linear interpolation on a known linear function", "[approximation][interpolation][linear]")
{
    const nane::grid<1> grid{0.0, 10.0, 1001};
    const nane::vector<> values = nane::sample(grid, [](double x) { return 3.0 * x + 2.0; });

    const double position = 5.234;
    const double expected = 3.0 * position + 2.0;

    const double result = nane::approximation::linear(values, position, grid);

    REQUIRE(result == Catch::Approx(expected).margin(1e-12));
}

TEST_CASE("Linear extrapolation on a known linear function", "[approximation][extrapolation][linear]")
{
    const nane::grid<1> grid{0.0, 10.0, 1001};
    const nane::vector<> values = nane::sample(grid, [](double x) { return 3.0 * x + 2.0; });

    SECTION("Left extrapolation")
    {
        const double position = -2.5;
        const double expected = 3.0 * position + 2.0;

        const double result = nane::approximation::linear(values, position, grid);

        REQUIRE(result == Catch::Approx(expected).margin(1e-12));
    }

    SECTION("Right extrapolation")
    {
        const double position = 12.5;
        const double expected = 3.0 * position + 2.0;

        const double result = nane::approximation::linear(values, position, grid);

        REQUIRE(result == Catch::Approx(expected).margin(1e-12));
    }
}
