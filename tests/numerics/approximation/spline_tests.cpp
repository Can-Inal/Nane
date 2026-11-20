#include "nane/geometry/sample.hpp"
#include "nane/numerics/approximation/spline.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

TEST_CASE("Natural cubic spline matches a known four-point result", "[approximation][spline][interpolation]")
{
    const nane::grid<1> grid{0.0, 3.0, 4};
    const nane::vector<> values{0.0, 1.0, 1.0, 0.0};

    REQUIRE(nane::approximation::spline(values, 1.5, grid) == Catch::Approx(1.15).margin(1e-12));

    // Every given sample must be reproduced at its own grid coordinate.
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        REQUIRE(nane::approximation::spline(values, grid.axis(0)[i], grid) == Catch::Approx(values[i]).margin(1e-12));
    }
}

TEST_CASE("Natural cubic spline reduces to linear interpolation with two points", "[approximation][spline][extrapolation]")
{
    const nane::grid<1> grid{0.0, 2.0, 2};
    const nane::vector<> values{1.0, 5.0};

    REQUIRE(nane::approximation::spline(values, 1.0, grid) == Catch::Approx(3.0));
    REQUIRE(nane::approximation::spline(values, -1.0, grid) == Catch::Approx(-1.0));
    REQUIRE(nane::approximation::spline(values, 3.0, grid) == Catch::Approx(7.0));
}

TEST_CASE("Two-dimensional natural cubic spline on a rectangular grid", "[approximation][spline][2d]")
{
    const nane::grid<2> grid{{0.0, 3.0, 4}, {0.0, 2.0, 3}};

    // At the grid nodes, u(x) = {0, 1, 1, 0} and v(y) = {0, 1, 0}.
    // Sampling a separable function makes the expected 2D result easy to verify.
    const auto values = nane::sample(grid, [](double x, double y) { return 0.5 * x * (3.0 - x) * y * (2.0 - y); });

    SECTION("Interpolation between nodes")
    {
        // The 1D natural spline factors are 1.15 along X and 0.6875 along Y.
        REQUIRE(nane::approximation::spline(values, nane::vector<2>{1.5, 0.5}, grid) == Catch::Approx(0.790625).margin(1e-12));
    }

    SECTION("Values are reproduced at grid points")
    {
        REQUIRE(nane::approximation::spline(values, nane::vector<2>{1.0, 1.0}, grid) == Catch::Approx(1.0).margin(1e-12));
        REQUIRE(nane::approximation::spline(values, nane::vector<2>{3.0, 2.0}, grid) == Catch::Approx(0.0).margin(1e-12));
    }
}

TEST_CASE("Two-dimensional spline preserves an affine function when extrapolating", "[approximation][spline][2d][extrapolation]")
{
    const nane::grid<2> grid{{0.0, 2.0, 3}, {0.0, 3.0, 4}};
    const auto values = nane::sample(grid, [](double x, double y) { return 1.0 + 2.0 * x - 3.0 * y; });

    REQUIRE(nane::approximation::spline(values, nane::vector<2>{-1.0, 4.0}, grid) == Catch::Approx(-13.0).margin(1e-12));
}

TEST_CASE("Two-dimensional spline rejects mismatched matrix dimensions", "[approximation][spline][2d][validation]")
{
    const nane::grid<2> grid{{0.0, 3.0, 4}, {0.0, 2.0, 3}};
    const nane::matrix<> values(3, 4); // Transposed dimensions: expected 4 rows and 3 columns.

    REQUIRE_THROWS_AS(nane::approximation::spline(values, nane::vector<2>{1.5, 0.5}, grid), std::invalid_argument);
}

TEST_CASE("Spline rejects a mismatched number of values", "[approximation][spline][validation]")
{
    const nane::grid<1> grid{0.0, 3.0, 4};
    const nane::vector<> values{0.0, 1.0, 0.0};

    REQUIRE_THROWS_AS(nane::approximation::spline(values, 1.5, grid), std::invalid_argument);
}
