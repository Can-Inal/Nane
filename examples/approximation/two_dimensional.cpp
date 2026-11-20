#include <nane/geometry/sample.hpp>
#include <nane/numerics/approximation/lagrange.hpp>
#include <nane/numerics/approximation/linear.hpp>
#include <nane/numerics/approximation/spline.hpp>

int main(int argc, char* argv[])
{
    const nane::grid<2> grid{{0.0, 3.0, 4}, {0.0, 4.0, 5}};
    const nane::vector<2> position{1.5, 2.5};

    // clang-format off
    const auto values = nane::sample(grid, [](double x, double y) 
    { 
	return x * x + 2.0 * y; 
    });
    // clang-format on

    const double linear_value = nane::approximation::linear(values, position, grid);
    const double lagrange_value = nane::approximation::lagrange(values, position, grid);
    const double spline_value = nane::approximation::spline(values, position, grid);

    return 0;
}
