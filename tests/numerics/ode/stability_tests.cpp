#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/stability.hpp"

#include <catch2/catch_test_macros.hpp>
#include <complex>

TEST_CASE("Explicit Euler's method stability test", "[ode][explicit_euler][stability]")
{
    // Butcher table for explicit Euler.
    const nane::ode::butcher_table table{
        {0.0},
        {
            {0.0},
        },
        {1.0},
    };

    // Explicit Euler should be stable inside its stability region.
    {
        const std::complex<double> z{-1.0, -0.5};
        REQUIRE(nane::ode::is_stable(table, z) == true);
    }

    // Explicit Euler should be stable on the boundary of its stability region.
    {
        const std::complex<double> z{-1.0, -1.0};
        REQUIRE(nane::ode::is_stable(table, z) == true);
    }

    // Explicit Euler should not be stable outside its stability region.
    {
        const std::complex<double> z{5.0, -3.0};
        REQUIRE(nane::ode::is_stable(table, z) == false);
    }
}

TEST_CASE("Implicit Midpoint method stability test", "[ode][implicit_midpoint][stability]")
{
    // Butcher table for Implicit Midpoint.
    const nane::ode::butcher_table table{
        {0.5},
        {
            {0.5},
        },
        {1.0},
    };

    // Implicit Midpoint should be stable inside its stability region.
    {
        const std::complex<double> z{-2.0, 3.0};
        REQUIRE(nane::ode::is_stable(table, z) == true);
    }

    // Implicit Midpoint should be stable on the boundary of its stability region.
    {
        const std::complex<double> z{0.0, 2.0};
        REQUIRE(nane::ode::is_stable(table, z) == true);
    }

    // Implicit Midpoint should not be stable outside its stability region.
    {
        const std::complex<double> z{1.0, -1.0};
        REQUIRE(nane::ode::is_stable(table, z) == false);
    }
}

TEST_CASE("Explicit Euler stability function test", "[ode][explicit_euler][stability]")
{
    const nane::ode::butcher_table table{
        {0.0},
        {
            {0.0},
        },
        {1.0},
    };

    const std::complex<double> z{-1.0, 0.5};

    const auto value = nane::ode::stability_function(table, z);
    const auto expected = 1.0 + z;

    REQUIRE(std::abs(value - expected) < 1e-12);
}
