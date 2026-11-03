#pragma once

// Core
#include "nane/error.hpp"
#include "nane/version.hpp"

// Geometry
#include "nane/geometry/grid.hpp"
#include "nane/geometry/sample.hpp"

// Math
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"

// numerics/calculus
#include "nane/numerics/calculus/differentiation/differentiate.hpp"
#include "nane/numerics/calculus/differentiation/stencil.hpp"
#include "nane/numerics/calculus/differentiation/stencils.hpp"
#include "nane/numerics/calculus/integration/integrate.hpp"
#include "nane/numerics/calculus/integration/quadrature_rule.hpp"
#include "nane/numerics/calculus/integration/rules.hpp"

// numerics/nonlinear
#include "nane/numerics/nonlinear/fixed_point.hpp"
#include "nane/numerics/nonlinear/newton.hpp"

// numerics/ode
#include "nane/numerics/ode/ivp.hpp"
#include "nane/numerics/ode/runge_kutta/butcher_table.hpp"
#include "nane/numerics/ode/runge_kutta/methods.hpp"
#include "nane/numerics/ode/runge_kutta/order.hpp"
#include "nane/numerics/ode/runge_kutta/runge_kutta.hpp"
#include "nane/numerics/ode/runge_kutta/stability.hpp"
#include "nane/numerics/ode/solution.hpp"
#include "nane/numerics/ode/solve.hpp"

// numerics/pde
#include "nane/numerics/pde/assemble.hpp"
#include "nane/numerics/pde/boundary_conditions.hpp"
#include "nane/numerics/pde/finite_difference.hpp"
#include "nane/numerics/pde/ibvp.hpp"
#include "nane/numerics/pde/solution.hpp"
#include "nane/numerics/pde/solve.hpp"
