@page custom_methods Defining Custom Numerical Methods

Nane is designed around a simple idea: the user can provide the numerical
scheme, while Nane provides the generic executor.

The same pattern appears throughout the library:

```text
butcher_table    -> runge_kutta / ode::solve
stencil          -> differentiate
quadrature_rule  -> integrate
```

The convenience methods shipped with Nane are intentionally small. They are
not intended to replace the underlying mathematical representations.

## Custom Runge–Kutta method

An explicit Runge–Kutta method is described by a Butcher table

@f[
\begin{array}{c|c}
\alpha & \beta \\
\hline
       & \gamma
\end{array}.
@f]

For Heun's second-order method,

@f[
\begin{array}{c|cc}
0 & 0 & 0 \\
1 & 1 & 0 \\
\hline
  & \tfrac12 & \tfrac12
\end{array}.
@f]

The corresponding Nane object is:

```cpp
const nane::ode::butcher_table method{
    {0.0, 1.0},
    {
        {0.0, 0.0},
        {1.0, 0.0},
    },
    {0.5, 0.5},
};
```

It can be passed directly to the generic ODE interface:

```cpp
const auto solution = nane::ode::solve(problem, method, {0.0, 10.0, 1001});
```

For adaptive integration, provide an embedded set of weights as part of the
Butcher table. Nane can then estimate the local error and use the tolerance
supplied to the solver.

## Custom finite-difference stencil

A stencil represents

@f[
f^{(m)}(x)
\approx
\frac{1}{h^m}
\sum_i w_i f(x+c_i h).
@f]

A centered first derivative is therefore:

```cpp
const nane::calculus::stencil central_first{
    1,
    1e-6,
    {
        {-1, -0.5},
        {1, 0.5},
    },
};
```

and can be applied with:

```cpp
const double derivative = nane::calculus::differentiate(function, point, central_first);
```

The derivative order, step size, offsets, and weights remain explicit.

## Custom quadrature rule

Nane quadrature rules are defined on the reference interval @f$[0,1]@f$.
A trapezoidal rule can be written as:

```cpp
const nane::calculus::quadrature_rule trapezoidal{
    {0.0, 0.5},
    {1.0, 0.5},
};
```

and used on any interval:

```cpp
const double integral = nane::calculus::integrate(function, trapezoidal, {0.0, 1.0});
```

For example, the two-point Gauss–Legendre rule on @f$[0,1]@f$ has nodes

@f[
c_{1,2}
=
\frac12\left(1\mp\frac{1}{\sqrt{3}}\right)
@f]

and weights

@f[
w_1=w_2=\frac12.
@f]

This can be expressed directly as another `quadrature_rule` without adding a
new integration algorithm.

## Why this matters

The executor does not need to know the name of the method. It only needs the
mathematical data required to apply it. This keeps the public API small while
allowing specialized or experimental numerical schemes to be expressed by
users of the library.
