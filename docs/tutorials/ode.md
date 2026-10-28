@page ode_tutorial Solving ODEs

This tutorial solves a nonlinear pendulum with Nane and uses the problem to
illustrate the main ODE API.

## Mathematical model

For a pendulum of length @f$L@f$ without damping or forcing,

@f[
\theta''(t) + \frac{g}{L}\sin\theta(t) = 0.
@f]

Introduce angular velocity

@f[
\omega = \theta'.
@f]

Then the second-order equation becomes the first-order system

@f[
\theta' = \omega,
\qquad
\omega' = -\frac{g}{L}\sin\theta.
@f]

The state vector is

@f[
\mathbf{x}(t)
=
\begin{bmatrix}
\theta(t) \\
\omega(t)
\end{bmatrix}.
@f]

## Define the derivative

```cpp
const double g = 9.81;
const double l = 1.0;

const auto system = [&](double, const nane::vector<2>& x)
{
    nane::vector<2> dx;

    dx[0] = x[1];
    dx[1] = -(g / l) * std::sin(x[0]);

    return dx;
};
```

The function directly mirrors the mathematical system. The first state
component is the angle and the second is the angular velocity.

## Initial value problem

For

@f[
\theta(0)=1,
\qquad
\omega(0)=0,
@f]

write:

```cpp
const nane::vector<2> initial_value{
    1.0,
    0.0,
};

const auto problem = nane::ode::ivp(system, initial_value);
```

## Solve on a fixed time grid

A grid with a point count specifies fixed time points:

```cpp
const auto solution = nane::ode::solve(problem, nane::ode::methods::heun, {0.0, 10.0, 1001});
```

Here the interval has 1001 points, hence 1000 subintervals and

@f[
h = \frac{10-0}{1000}=0.01.
@f]

The method is represented by a Butcher table and `solve` dispatches to the
Runge–Kutta implementation.

The lower-level executor remains available when the algorithm should be named
explicitly:

```cpp
nane::ode::runge_kutta(problem, method, time_grid);
```

## Adaptive time stepping

An interval without a prescribed count,

```cpp
{0.0, 10.0}
```

is intended for adaptive Runge–Kutta integration. The supplied Butcher table
must contain an embedded pair so that the method can estimate local error.

The public tolerance can then be supplied to `solve`:

```cpp
nane::ode::solve(problem, embedded_method, {0.0, 10.0}, 1e-6);
```

Nane does not hide this distinction: fixed and adaptive integration are
selected by the time-grid semantics and the numerical method supplied by the
user.

## Validate with mechanical energy

For the undamped pendulum, the mechanical energy per unit mass is

@f[
E(t)
=
\frac12 L^2\omega(t)^2
+
gL\left(1-\cos\theta(t)\right).
@f]

The exact system satisfies

@f[
\frac{dE}{dt}=0.
@f]

A numerical Runge–Kutta solution will generally not preserve this value
exactly. Monitoring

@f[
E(t)-E(0)
@f]

is therefore a useful numerical diagnostic and a good way to compare step
sizes or methods.

## Add damping

A damped pendulum satisfies

@f[
\theta'' + \gamma\theta' + \frac{g}{L}\sin\theta = 0,
@f]

or

@f[
\theta'=\omega,
\qquad
\omega'=-\frac{g}{L}\sin\theta-\gamma\omega.
@f]

Its energy obeys

@f[
\frac{dE}{dt}=-\gamma L^2\omega^2 \le 0.
@f]

This provides another useful validation case: the phase-space trajectory should
spiral toward the equilibrium and the mechanical energy should decay.
