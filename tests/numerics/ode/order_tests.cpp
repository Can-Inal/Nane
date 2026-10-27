#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/order.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cmath>

TEST_CASE("Runge-Kutta order conditions", "[ode][runge-kutta][order]")
{
    SECTION("order zero")
    {
        const nane::ode::butcher_table table{
            {0.5},
            {
                {0.0},
            },
            {1.0},
        };

        REQUIRE(nane::ode::order(table) == 0);
    }

    SECTION("order one")
    {
        const nane::ode::butcher_table table{
            {0.0, 1.0},
            {
                {0.0, 0.0},
                {1.0, 0.0},
            },
            {1.0, 0.0},
        };

        REQUIRE(nane::ode::order(table) == 1);
    }

    SECTION("order two")
    {
        const nane::ode::butcher_table table{
            {0.0, 2.0 / 3.0},
            {
                {0.0, 0.0},
                {2.0 / 3.0, 0.0},
            },
            {1.0 / 4.0, 3.0 / 4.0},
        };

        REQUIRE(nane::ode::order(table) == 2);
    }

    SECTION("order three")
    {
        const nane::ode::butcher_table table{
            {0.0, 0.5, 1.0},
            {
                {0.0, 0.0, 0.0},
                {0.5, 0.0, 0.0},
                {-1.0, 2.0, 0.0},
            },
            {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
        };

        REQUIRE(nane::ode::order(table) == 3);
    }

    SECTION("order four")
    {
        const nane::ode::butcher_table table{
            {0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0},
            {
                {0.0, 0.0, 0.0, 0.0},
                {1.0 / 3.0, 0.0, 0.0, 0.0},
                {-1.0 / 3.0, 1.0, 0.0, 0.0},
                {1.0, -1.0, 1.0, 0.0},
            },
            {1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0},
        };

        REQUIRE(nane::ode::order(table) == 4);
    }

    SECTION("implicit order four")
    {
        const double r = std::sqrt(3.0) / 6.0;

        const nane::ode::butcher_table table{
            {0.5 - r, 0.5 + r},
            {
                {0.25, 0.25 - r},
                {0.25 + r, 0.25},
            },
            {0.5, 0.5},
        };

        REQUIRE(nane::ode::order(table) == 4);
    }
}
