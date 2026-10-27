#include "nane/numerics/ode/runge_kutta/methods.hpp"
#include "nane/numerics/ode/runge_kutta/stability.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>

TEST_CASE("Runge-Kutta stability", "[ode][runge-kutta][stability]")
{
    using namespace nane::ode;

    REQUIRE(is_stable(methods::explicit_euler, {-1.0, -0.5}));
    REQUIRE(is_stable(methods::explicit_euler, {-1.0, -1.0}));
    REQUIRE_FALSE(is_stable(methods::explicit_euler, {5.0, -3.0}));

    REQUIRE(is_stable(methods::implicit_midpoint, {-2.0, 3.0}));
    REQUIRE(is_stable(methods::implicit_midpoint, {0.0, 2.0}));
    REQUIRE_FALSE(is_stable(methods::implicit_midpoint, {1.0, -1.0}));

    const std::complex<double> z{-1.0, 0.5};

    REQUIRE(std::abs(stability_function(methods::explicit_euler, z) - (1.0 + z)) < 1e-12);
}
