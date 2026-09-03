# Nane {#mainpage}

Nane is a modern C++23 scientific computing library for numerical methods,
linear algebra, ordinary differential equations, computational grids, and
scientific visualization.

The library is designed around small, composable numerical APIs while using
Eigen directly for linear algebra.

## Quick Example

Consider the initial value problem

\f[
\dot{x}(t) = -t^2x(t),
\qquad
x(0) = 5.3.
\f]

It can be solved using Heun's method:

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

## Ordinary Differential Equations

Nane provides scalar and vector-valued one-step methods.

Currently available methods include:

- Explicit Euler
- Implicit Euler
- Heun's method
- Generic Runge-Kutta methods

A Runge-Kutta method is described using the Butcher coefficients

\f[
\alpha \in \mathbb{R}^s,
\qquad
B \in \mathbb{R}^{s\times s},
\qquad
\gamma \in \mathbb{R}^s.
\f]

The stages satisfy

\f[
k_i
=
f\left(
t_n + \alpha_i\tau,
x_n + \tau \sum_{j=1}^{s}\beta_{ij}k_j
\right),
\f]

followed by

\f[
x_{n+1}
=
x_n
+
\tau
\sum_{i=1}^{s}\gamma_i k_i.
\f]

Explicit and implicit Runge-Kutta methods are supported. Implicit stage systems
are currently solved using fixed-point iteration.

## Linear Algebra

Nane uses Eigen for its linear algebra representation rather than introducing
a separate matrix or vector type.

## Computational Grids

Uniform computational domains can be constructed using `nane::uniform_grid`.

```cpp
nane::uniform_grid<1> grid({
    {0.0, 1.0, 101},
});
```

## Visualization

The `nane::viz` library provides plotting functionality backed by VTK while
keeping VTK implementation details outside the public plotting interface.

## API Documentation

Use the navigation sidebar to browse namespaces, classes, functions, and files.
