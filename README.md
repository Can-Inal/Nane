# Nane

[![Nane CI](https://github.com/Can-Inal/Nane/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/Can-Inal/Nane/actions/workflows/build.yml)

Nane is a C++23 scientific computing library focused on numerical methods,
linear algebra, ordinary differential equations, computational grids, and
scientific visualization.

Nane uses Eigen for linear algebra and VTK for visualization.

> Nane is currently under active development.

## Features

- Uniform computational grids
- Eigen-based linear algebra
- Scalar and vector-valued ODE solvers
- Generic explicit and implicit Runge-Kutta methods
- Explicit Euler
- Implicit Euler
- Heun's method
- Fixed-point iteration
- VTK-based scientific visualization

## Example

```cpp
#include <nane/geometry/uniform_grid.hpp>
#include <nane/numerics/ode/one_step.hpp>

int main()
{
    nane::uniform_grid<1> time_grid({
        {0.0, 2.0, 101},
    });

    auto derivative = [](double t, double x)
    {
        return -t * t * x;
    };

    auto solution = nane::heun(derivative, 5.3, time_grid);

    return 0;
}
