#include <cmath>
#include <nane/math/vector.hpp>
#include <nane/numerics/ode/runge_kutta/methods.hpp>
#include <nane/numerics/ode/runge_kutta/runge_kutta.hpp>

int main(int argc, char* argv[])
{
    const auto derivative = [](double t, const nane::vector<2>& x)
    {
        nane::vector<2> dx;

        dx[0] = x[1];
        dx[1] = -std::sin(x[0]) - 0.1 * x[1];

        return dx;
    };

    const nane::vector<2> initial_value{
        1.0,
        0.0,
    };

    const auto problem = nane::ode::ivp(derivative, initial_value);

    [[maybe_unused]] const auto solution = nane::ode::runge_kutta(problem, nane::ode::methods::heun, {0.0, 10.0, 1001});
}
