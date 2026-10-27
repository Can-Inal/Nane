#include <cmath>
#include <nane/numerics/calculus/integration/integrate.hpp>
#include <nane/numerics/calculus/integration/quadrature_rule.hpp>
#include <numbers>

int main(int argc, char* argv[])
{
    const double node_offset = 1.0 / std::sqrt(3.0);

    const nane::calculus::quadrature_rule gauss_legendre_2{
        {0.5 * (1.0 - node_offset), 0.5},
        {0.5 * (1.0 + node_offset), 0.5},
    };

    const auto function = [](double x) { return std::sin(x); };

    const double result = nane::calculus::integrate(function, gauss_legendre_2, {0.0, std::numbers::pi, 11});

    return 0;
}
