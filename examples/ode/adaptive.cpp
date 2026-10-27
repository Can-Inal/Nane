#include <cmath>
#include <nane/math/vector.hpp>
#include <nane/numerics/ode/runge_kutta/butcher_table.hpp>
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

    const nane::ode::butcher_table method{
        {0.0, 1.0},
        {
            {0.0, 0.0},
            {1.0, 0.0},
        },
        {0.5, 0.5},
        {1.0, 0.0},
    };

    [[maybe_unused]] const auto solution = nane::ode::runge_kutta(problem, method, {0.0, 10.0});

    return 0;
}
