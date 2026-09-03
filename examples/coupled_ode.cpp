#include <Eigen/Core>
#include <nane/core/symbolic.hpp>
#include <nane/geometry/uniform_grid.hpp>
#include <nane/numerics/ode/one_step.hpp>

int main()
{
    nane::uniform_grid<1> time_grid({
        {0.0, 10.0, 1001},
    });

    const auto [t, x] = nane::symbols<2>();

    const auto derivative = nane::system(x[1], -nane::sin(x[0]) - 0.1 * x[1]);

    Eigen::VectorXd initial_value(2);
    initial_value << 1.0, 0.0;

    [[maybe_unused]] auto solution = nane::heun(derivative, initial_value, time_grid);

    return 0;
}
