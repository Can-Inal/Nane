#include <cmath>
#include <nane/numerics/calculus/differentiation/differentiate.hpp>
#include <nane/numerics/calculus/differentiation/stencils.hpp>
#include <numbers>

int main(int argc, char* argv[])
{
    const auto function = [](double x) { return std::sin(x); };

    const double point = std::numbers::pi / 4.0;

    [[maybe_unused]] const auto first_derivative = nane::calculus::differentiate(function, point, nane::calculus::stencils::central_first(1e-6));
    [[maybe_unused]] const auto second_derivative = nane::calculus::differentiate(function, point, nane::calculus::stencils::central_second(1e-4));

    return 0;
}
