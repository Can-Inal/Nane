@page getting_started Getting Started

Nane is a C++23 scientific computing library built around small numerical
objects and generic executors. This page shows the shortest path from an empty
project to a program that uses Nane.

## Requirements

You need:

- a C++23-capable compiler,
- CMake 3.20 or newer,
- Git,
- Nane and its Git submodules.

## Recommended project layout

A small application can keep Nane under `vendor/`:

```text
pendulum/
├── CMakeLists.txt
├── src/
│   └── main.cpp
└── vendor/
    └── Nane/
```

Add Nane as a Git submodule:

```bash
git submodule add https://github.com/Can-Inal/Nane.git vendor/Nane
git submodule update --init --recursive
```

The recursive update is important because Nane may itself contain vendored
submodules required by its public headers.

## Minimal CMake project

```cmake
cmake_minimum_required(VERSION 3.20)

project(pendulum LANGUAGES CXX)

add_subdirectory(vendor/Nane EXCLUDE_FROM_ALL)

add_executable(pendulum
    src/main.cpp
)

target_compile_features(pendulum PRIVATE cxx_std_23)

target_link_libraries(pendulum
    PRIVATE
        nane::nane
)
```

Configure and build with:

```bash
cmake -S . -B build
cmake --build build
```

## A first Nane program

The following example solves the nonlinear pendulum

@f[
\theta' = \omega,
\qquad
\omega' = -\frac{g}{L}\sin(\theta).
@f]

```cpp
#include <cmath>

#include <nane/math/vector.hpp>
#include <nane/numerics/ode/ivp.hpp>
#include <nane/numerics/ode/runge_kutta/methods.hpp>
#include <nane/numerics/ode/solve.hpp>

int main()
{
    const double g = 9.81;
    const double l = 1.0;

    const auto system = [&](double, const nane::vector<2>& x)
    {
        nane::vector<2> dx;

        dx[0] = x[1];
        dx[1] = -(g / l) * std::sin(x[0]);

        return dx;
    };

    const nane::vector<2> initial_value{
        1.0,
        0.0,
    };

    const auto problem = nane::ode::ivp(system, initial_value);

    [[maybe_unused]] const auto solution =
        nane::ode::solve(
            problem,
            nane::ode::methods::heun,
            {0.0, 10.0, 1001});
}
```

The three objects are intentionally visible in the call:

1. the initial value problem,
2. the numerical method,
3. the time grid.

A grid `{a, b, n}` prescribes `n` time points and therefore fixed time
stepping. A grid `{a, b}` leaves the time points unspecified and is used by
adaptive Runge–Kutta methods when the supplied Butcher table is embedded.

## Where to go next

Continue with @ref cmake_integration "Using Nane with CMake" for dependency
management, or @ref ode_tutorial "Solving ODEs" for a more complete numerical
example.
