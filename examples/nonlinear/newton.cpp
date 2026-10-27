#include <nane/math/vector.hpp>
#include <nane/numerics/nonlinear/newton.hpp>

int main(int argc, char* argv[])
{
    const auto scalar_function = [](double x) { return x * x - 2.0; };

    [[maybe_unused]] const auto scalar_root = nane::nonlinear::newton(scalar_function, 1.0);

    const auto vector_function = [](const nane::vector<2>& x)
    {
        return nane::vector<2>{
            x[0] * x[0] + x[1] * x[1] - 1.0,
            x[0] - x[1],
        };
    };

    const nane::vector<2> initial_guess{
        0.8,
        0.6,
    };

    [[maybe_unused]] const auto vector_root = nane::nonlinear::newton(vector_function, initial_guess);

    return 0;
}
