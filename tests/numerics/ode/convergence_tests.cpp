#include "nane/numerics/ode/runge_kutta/methods.hpp"
#include "nane/numerics/ode/runge_kutta/runge_kutta.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>

namespace
{
    double observed_order(const nane::ode::butcher_table& method, std::size_t coarse_count, std::size_t fine_count)
    {
        const auto derivative = [](double, double x) { return -x; };
        const auto problem = nane::ode::ivp(derivative, 1.0);

        const auto coarse = nane::ode::runge_kutta(problem, method, {0.0, 1.0, coarse_count});
        const auto fine = nane::ode::runge_kutta(problem, method, {0.0, 1.0, fine_count});

        const double exact = std::exp(-1.0);

        const double coarse_error = std::abs(coarse.values[coarse.values.size() - 1] - exact);
        const double fine_error = std::abs(fine.values[fine.values.size() - 1] - exact);

        const double coarse_step = 1.0 / static_cast<double>(coarse_count - 1);
        const double fine_step = 1.0 / static_cast<double>(fine_count - 1);

        return std::log(coarse_error / fine_error) / std::log(coarse_step / fine_step);
    }
} // namespace

TEST_CASE("Runge-Kutta methods converge with their expected order", "[ode][runge-kutta][convergence]")
{
    using namespace nane::ode;

    REQUIRE(observed_order(methods::explicit_euler, 11, 21) == Catch::Approx(1.0).margin(0.1));
    REQUIRE(observed_order(methods::implicit_euler, 11, 21) == Catch::Approx(1.0).margin(0.1));
    REQUIRE(observed_order(methods::heun, 11, 21) == Catch::Approx(2.0).margin(0.1));
    REQUIRE(observed_order(methods::modified_euler, 11, 21) == Catch::Approx(2.0).margin(0.1));
    REQUIRE(observed_order(methods::implicit_midpoint, 11, 21) == Catch::Approx(2.0).margin(0.1));
    REQUIRE(observed_order(methods::crank_nicolson, 11, 21) == Catch::Approx(2.0).margin(0.1));
    REQUIRE(observed_order(methods::rk4, 11, 21) == Catch::Approx(4.0).margin(0.15));
}
