# Nane

<p align="left">
    <a href="https://github.com/Can-Inal/Nane/actions/workflows/build.yml">
        <img src="https://img.shields.io/github/actions/workflow/status/Can-Inal/Nane/build.yml?branch=main&style=for-the-badge&logo=github&logoColor=white&label=Build" alt="Build">
    </a>
    <a href="https://can-inal.github.io/Nane/">
        <img src="https://img.shields.io/badge/Documentation-Doxygen-2C4AA8?style=for-the-badge&logo=doxygen&logoColor=white" alt="Documentation">
    </a>
    <img src="https://img.shields.io/badge/C%2B%2B-23-00599C?style=for-the-badge&logo=cplusplus&logoColor=white" alt="C++23">
    <img src="https://img.shields.io/badge/CMake-3.20%2B-064F8C?style=for-the-badge&logo=cmake&logoColor=white" alt="CMake">
    <a href="LICENSE">
        <img src="https://img.shields.io/github/license/Can-Inal/Nane?style=for-the-badge" alt="License">
    </a>
</p>

Nane is a modern C++23 scientific computing library focused on numerical
methods, linear algebra, ordinary differential equations, nonlinear solvers,
and computational grids.

Nane uses Eigen directly for vector and matrix representations and provides
small, composable interfaces for numerical algorithms.

> **Nane is currently under active development.**

## Features

- Uniform computational grids
- Eigen-based linear algebra
- Scalar and vector-valued ODE solvers
- Generic explicit and implicit Runge-Kutta methods
- Explicit Euler
- Implicit Euler
- Heun's method
- Fixed-point iteration
- Optional OpenMP support
- Cross-platform builds on Linux, macOS, and Windows
- Catch2-based test suite

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
```

## Documentation

The API documentation is generated with Doxygen and published through
GitHub Pages.

[Open the Nane documentation](https://can-inal.github.io/Nane/)

## License

Nane is available under the [MIT License](LICENSE).
