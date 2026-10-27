#include <nane/numerics/ode/runge_kutta/methods.hpp>
#include <nane/numerics/ode/runge_kutta/runge_kutta.hpp>

int main(int argc, char* argv[])
{
    const auto derivative = [](double t, double x) { return -t * t * x; };

    const auto problem = nane::ode::ivp(derivative, 5.3);

    [[maybe_unused]] const auto solution = nane::ode::runge_kutta(problem, nane::ode::methods::heun, {0.0, 2.0, 101});
}
