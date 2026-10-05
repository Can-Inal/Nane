#include <Eigen/Core>
#include <nane/geometry/uniform_grid.hpp>
#include <nane/numerics/ode/one_step.hpp>
#include <nane/symbolic.hpp>

int main()
{
    nane::uniform_grid<1> time_grid({
        {0.0, 10.0, 1001},
    });

    const auto [t, x] = nane::symbols<2>();
    const auto derivative = nane::system(x[1], -nane::sin(x[0]) - 0.1 * x[1]);

    Eigen::VectorXd initial_value(2);
    initial_value << 1.0, 0.0;

    const auto problem = nane::ivp(derivative, initial_value);
    [[maybe_unused]] auto solution = nane::heun(problem, time_grid);

    return 0;
}
