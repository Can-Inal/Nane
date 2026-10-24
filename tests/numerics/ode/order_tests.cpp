#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/order.hpp"

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

TEST_CASE("Ralston's second order method butcher table order test", "[ode][custom][order]")
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

TEST_CASE("Kutta's third order method butcher table order test", "[ode][custom][order]")
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

TEST_CASE("Runge-Kutta 3/8 rule butcher table order test", "[ode][custom][order]")
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

TEST_CASE("Two-stage Gauss-Legendre butcher table order test", "[ode][custom][order]")
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

TEST_CASE("Butcher table violating invariance against autonomization has order zero", "[ode][custom][order]")
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

TEST_CASE("Butcher table satisfying only first order conditions has order one", "[ode][custom][order]")
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
