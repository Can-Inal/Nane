#include "nane/geometry/sample.hpp"
#include "nane/numerics/approximation/lagrange.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

TEST_CASE("Lagrange interpolation reproduces a cubic polynomial in 1D", "[approximation][lagrange][interpolation]")
{
    const nane::grid<1> grid{0.0, 5.0, 6};
    const auto values = nane::sample(grid, [](double x) { return x * x * x; });

    REQUIRE(nane::approximation::lagrange(values, 2.5, grid) == Catch::Approx(15.625).margin(1e-12));
}

TEST_CASE("Lagrange extrapolation reproduces a cubic polynomial in 1D", "[approximation][lagrange][extrapolation]")
{
    const nane::grid<1> grid{0.0, 5.0, 6};
    const auto values = nane::sample(grid, [](double x) { return x * x * x; });

    SECTION("Left extrapolation")
    {
        REQUIRE(nane::approximation::lagrange(values, -1.0, grid) == Catch::Approx(-1.0).margin(1e-12));
    }

    SECTION("Right extrapolation")
    {
        REQUIRE(nane::approximation::lagrange(values, 6.0, grid) == Catch::Approx(216.0).margin(1e-12));
    }
}

TEST_CASE("Two-point Lagrange interpolation reproduces a linear function", "[approximation][lagrange][interpolation]")
{
    const nane::grid<1> grid{0.0, 3.0, 4};
    const nane::vector<> values{1.0, 3.0, 5.0, 7.0};

    REQUIRE(nane::approximation::lagrange<2>(values, 1.5, grid) == Catch::Approx(4.0));
}

TEST_CASE("Two-dimensional Lagrange interpolation and extrapolation", "[approximation][lagrange][2d]")
{
    // A rectangular grid distinguishes X rows from Y columns.
    const nane::grid<2> grid{{0.0, 5.0, 6}, {0.0, 4.0, 5}};

    // nane::sample stores f(x_i, y_j) in values(i, j).
    const auto values = nane::sample(grid, [](double x, double y) { return x * x * x + 2.0 * y * y; });

    SECTION("Interpolation")
    {
        REQUIRE(nane::approximation::lagrange(values, nane::vector<2>{2.5, 1.5}, grid) == Catch::Approx(20.125).margin(1e-11));
    }

    SECTION("Extrapolation along both axes")
    {
        REQUIRE(nane::approximation::lagrange(values, nane::vector<2>{-1.0, 5.0}, grid) == Catch::Approx(49.0).margin(1e-10));
    }

    SECTION("Values are reproduced at grid points")
    {
        REQUIRE(nane::approximation::lagrange(values, nane::vector<2>{2.0, 3.0}, grid) == Catch::Approx(26.0).margin(1e-11));
        REQUIRE(nane::approximation::lagrange(values, nane::vector<2>{5.0, 4.0}, grid) == Catch::Approx(157.0).margin(1e-11));
    }
}

TEST_CASE("Two-dimensional Lagrange rejects mismatched matrix dimensions", "[approximation][lagrange][2d][validation]")
{
    const nane::grid<2> grid{{0.0, 5.0, 6}, {0.0, 4.0, 5}};
    const nane::matrix<> values(5, 6); // Transposed dimensions: expected 6 rows and 5 columns.

    REQUIRE_THROWS_AS(nane::approximation::lagrange(values, nane::vector<2>{2.5, 1.5}, grid), std::invalid_argument);
}

TEST_CASE("Lagrange rejects too few grid points", "[approximation][lagrange][validation]")
{
    const nane::grid<1> grid{0.0, 2.0, 3};
    const nane::vector<> values{0.0, 1.0, 2.0};

    REQUIRE_THROWS_AS(nane::approximation::lagrange<4>(values, 1.0, grid), std::invalid_argument);
}
