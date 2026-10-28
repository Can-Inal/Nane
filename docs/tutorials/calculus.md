@page calculus_tutorial Numerical Calculus

Nane separates the mathematical rule from the generic operation that applies
it. Numerical differentiation is driven by a `stencil`, while numerical
integration is driven by a `quadrature_rule`.

## Differentiation

A finite-difference stencil represents

@f[
f^{(m)}(x)
\approx
\frac{1}{h^m}
\sum_i w_i f(x+c_i h),
@f]

where @f$c_i@f$ are integer offsets and @f$w_i@f$ are weights.

For a common centered first derivative, use the convenience collection:

```cpp
#include <nane/numerics/calculus/differentiation/differentiate.hpp>
#include <nane/numerics/calculus/differentiation/stencils.hpp>

const auto function = [](double x)
{
    return std::sin(x);
};

const double derivative = nane::calculus::differentiate(function, 1.0, nane::calculus::stencils::central_first(1e-6));
```

The same executor can approximate higher derivatives when the stencil has the
corresponding derivative order:

```cpp
const double second_derivative = nane::calculus::differentiate(function, 1.0, nane::calculus::stencils::central_second(1e-4));
```

For a vector-valued function

@f[
\mathbf{F}:\mathbb{R}^n\rightarrow\mathbb{R}^m,
@f]

first differentiation produces the Jacobian matrix. A first-derivative stencil
is applied independently in each coordinate direction.

## Integration

A quadrature rule on the reference interval @f$[0,1]@f$ is represented by
nodes @f$c_i@f$ and weights @f$w_i@f$:

@f[
\int_0^1 f(x)\,dx
\approx
\sum_i w_i f(c_i).
@f]

For an arbitrary interval @f$[a,b]@f$, Nane applies the affine mapping

@f[
x=a+(b-a)c_i
@f]

and therefore

@f[
\int_a^b f(x)\,dx
\approx
(b-a)\sum_i w_i f\left(a+(b-a)c_i\right).
@f]

Using a predefined rule:

```cpp
#include <nane/numerics/calculus/integration/integrate.hpp>
#include <nane/numerics/calculus/integration/rules.hpp>

const auto function = [](double x)
{
    return std::sin(x);
};

const double value = nane::calculus::integrate(function, nane::calculus::rules::gauss_legendre_2, {0.0, 3.141592653589793});
```

If the grid contains a point count, the rule is applied separately on every
subinterval and the contributions are summed:

```cpp
nane::calculus::integrate(function, nane::calculus::rules::simpson, {0.0, 3.141592653589793, 101});
```

This gives the same `grid<1>` type two natural meanings:

- `{a, b}` — apply the quadrature rule once on the interval,
- `{a, b, n}` — composite quadrature over `n - 1` subintervals.

## Convenience rules are optional

The predefined `stencils` and `rules` namespaces contain only common choices.
The core API accepts user-defined `stencil` and `quadrature_rule` objects, so
specialized numerical schemes do not require changes to the library.
