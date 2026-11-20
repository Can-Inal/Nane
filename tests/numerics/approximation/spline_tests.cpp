#include "nane/numerics/approximation/spline.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

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
    const nane::vector<> u{0.0, 1.0, 1.0, 0.0};
    const nane::vector<> v{0.0, 1.0, 0.0};
    nane::matrix<> values(3, 4); // row = Y, column = X (current approximation overload convention)

    for (std::size_t j = 0; j < values.rows(); ++j)
    {
        for (std::size_t i = 0; i < values.columns(); ++i)
        {
            values(j, i) = u[i] * v[j];
        }
    }

    REQUIRE(nane::approximation::spline(values, nane::vector<2>{1.5, 0.5}, grid) == Catch::Approx(0.790625).margin(1e-12));
}

TEST_CASE("Spline rejects a mismatched number of values", "[approximation][spline][validation]")
{
    const nane::grid<1> grid{0.0, 3.0, 4};
    const nane::vector<> values{0.0, 1.0, 0.0};

    REQUIRE_THROWS_AS(nane::approximation::spline(values, 1.5, grid), std::invalid_argument);
}
