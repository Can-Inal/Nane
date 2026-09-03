#include <nane/core/symbolic.hpp>
#include <nane/geometry/uniform_grid.hpp>
#include <nane/numerics/ode/one_step.hpp>

int main()
{
    nane::uniform_grid<1> time_grid({
        {0.0, 2.0, 101},
    });

    const auto [t, x] = nane::symbols<2>();

    const auto derivative = nane::function(-t * t * x);

    [[maybe_unused]] auto solution = nane::heun(derivative, 5.3, time_grid);

    return 0;
}
