#include <nane/numerics/approximation/lagrange.hpp>
#include <nane/numerics/approximation/linear.hpp>
#include <nane/numerics/approximation/spline.hpp>

int main(int argc, char* argv[])
{
    const nane::grid<1> grid{0.0, 3.0, 4};
    const double position = 1.5;
    const nane::vector<> values{0.0, 1.0, 1.0, 0.0};

    const double linear_value = nane::approximation::linear(values, position, grid);
    const double lagrange_value = nane::approximation::lagrange(values, position, grid);
    const double spline_value = nane::approximation::spline(values, position, grid);

    return 0;
}
