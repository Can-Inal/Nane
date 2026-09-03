#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <nane/geometry/uniform_grid.hpp>
#include <nane/linalg/vector.hpp>
#include <nane/numerics/ode/one_step.hpp>
#include <nane/viz/plotting/plot.hpp>

int main()
{
    constexpr std::size_t steps = 20;
    constexpr std::size_t exact_steps = 1000;

    // -------------------------------------------------------------------------
    // Time grids
    // -------------------------------------------------------------------------

    nane::uniform_grid<1> time_grid({
        {0.0, 2.0, steps + 1},
    });

    nane::uniform_grid<1> exact_time_grid({
        {0.0, 2.0, exact_steps + 1},
    });

    // -------------------------------------------------------------------------
    // Scalar initial value problem
    //
    // x'(t) = -t^2 x(t)
    // x(0)  = 5.3
    //
    // Exact:
    //
    // x(t) = 5.3 exp(-t^3 / 3)
    // -------------------------------------------------------------------------

    auto scalar_derivative = [](double t, double x) { return -t * t * x; };
    constexpr double scalar_initial_value = 5.3;

    auto explicit_solution = nane::explicit_euler(scalar_derivative, scalar_initial_value, time_grid);
    auto implicit_solution = nane::implicit_euler(scalar_derivative, scalar_initial_value, time_grid);
    auto heun_solution = nane::heun(scalar_derivative, scalar_initial_value, time_grid);

    // -------------------------------------------------------------------------
    // Scalar exact solution
    // -------------------------------------------------------------------------

    auto scalar_exact_solution = nane::vector(exact_steps + 1);

    for (auto i = 0; i <= (int)exact_steps; ++i)
    {
        const double t = exact_time_grid.axis(0)[i];
        scalar_exact_solution[i] = scalar_initial_value * std::exp(-(t * t * t) / 3.0);
    }

    // -------------------------------------------------------------------------
    // Scalar errors
    // -------------------------------------------------------------------------

    auto explicit_error = nane::vector(steps + 1);
    auto implicit_error = nane::vector(steps + 1);
    auto heun_error = nane::vector(steps + 1);

    for (auto i = 0; i <= (int)steps; ++i)
    {
        const double t = time_grid.axis(0)[i];
        const double exact = scalar_initial_value * std::exp(-(t * t * t) / 3.0);

        explicit_error[i] = std::abs(explicit_solution[i] - exact);
        implicit_error[i] = std::abs(implicit_solution[i] - exact);
        heun_error[i] = std::abs(heun_solution[i] - exact);
    }

    // -------------------------------------------------------------------------
    // Vector initial value problem
    //
    // x'(t) = -x(t)
    // y'(t) = -2y(t)
    //
    // x(0) = 2
    // y(0) = 1
    //
    // Exact:
    //
    // x(t) = 2 exp(-t)
    // y(t) = exp(-2t)
    // -------------------------------------------------------------------------

    auto vector_derivative = [](double, const Eigen::VectorXd& state)
    {
        auto derivative = nane::vector(2);

        derivative[0] = -state[0];
        derivative[1] = -2.0 * state[1];

        return derivative;
    };

    auto vector_initial_value = nane::vector(2);
    vector_initial_value[0] = 2.0;
    vector_initial_value[1] = 1.0;

    auto vector_solution = nane::implicit_euler(vector_derivative, vector_initial_value, time_grid);

    // -------------------------------------------------------------------------
    // Vector exact solution
    // -------------------------------------------------------------------------

    auto vector_exact_x = nane::vector(exact_steps + 1);
    auto vector_exact_y = nane::vector(exact_steps + 1);

    for (auto i = 0; i <= (int)exact_steps; ++i)
    {
        const double t = exact_time_grid.axis(0)[i];

        vector_exact_x[i] = 2.0 * std::exp(-t);
        vector_exact_y[i] = std::exp(-2.0 * t);
    }

    // -------------------------------------------------------------------------
    // Vector errors
    // -------------------------------------------------------------------------

    auto vector_error_x = nane::vector(steps + 1);
    auto vector_error_y = nane::vector(steps + 1);

    for (auto i = 0; i <= (int)steps; ++i)
    {
        const double t = time_grid.axis(0)[i];

        const double exact_x = 2.0 * std::exp(-t);
        const double exact_y = std::exp(-2.0 * t);

        vector_error_x[i] = std::abs(vector_solution(0, i) - exact_x);
        vector_error_y[i] = std::abs(vector_solution(1, i) - exact_y);
    }

    // -------------------------------------------------------------------------
    // Plot
    // -------------------------------------------------------------------------

    nane::viz::plot comparison_plot;

    // -------------------------------------------------------------------------
    // Scalar solution comparison
    // -------------------------------------------------------------------------

    comparison_plot.subplot(2, 2, 0);

    comparison_plot.add_points(time_grid.axis(0), explicit_solution, "Explicit Euler");
    comparison_plot.add_points(time_grid.axis(0), implicit_solution, "Implicit Euler");
    comparison_plot.add_points(time_grid.axis(0), heun_solution, "Heun");
    comparison_plot.add_line(exact_time_grid.axis(0), scalar_exact_solution, "Exact");

    comparison_plot.title("Scalar Solution");
    comparison_plot.x_label("t");
    comparison_plot.y_label("x(t)");

    // -------------------------------------------------------------------------
    // Scalar error comparison
    // -------------------------------------------------------------------------

    comparison_plot.subplot(2, 2, 1);

    comparison_plot.add_line(time_grid.axis(0), explicit_error, "Explicit Euler Error");
    comparison_plot.add_line(time_grid.axis(0), implicit_error, "Implicit Euler Error");
    comparison_plot.add_line(time_grid.axis(0), heun_error, "Heun Error");

    comparison_plot.title("Scalar Absolute Error");
    comparison_plot.x_label("t");
    comparison_plot.y_label("|x_n - x(t_n)|");

    // -------------------------------------------------------------------------
    // Vector solution comparison
    // -------------------------------------------------------------------------

    comparison_plot.subplot(2, 2, 2);

    comparison_plot.add_points(time_grid.axis(0), vector_solution.row(0), "Implicit Euler x");
    comparison_plot.add_points(time_grid.axis(0), vector_solution.row(1), "Implicit Euler y");

    comparison_plot.add_line(exact_time_grid.axis(0), vector_exact_x, "Exact x");
    comparison_plot.add_line(exact_time_grid.axis(0), vector_exact_y, "Exact y");

    comparison_plot.title("Vector Implicit Euler");
    comparison_plot.x_label("t");
    comparison_plot.y_label("State");

    // -------------------------------------------------------------------------
    // Vector error comparison
    // -------------------------------------------------------------------------

    comparison_plot.subplot(2, 2, 3);

    comparison_plot.add_line(time_grid.axis(0), vector_error_x, "x Error");
    comparison_plot.add_line(time_grid.axis(0), vector_error_y, "y Error");

    comparison_plot.title("Vector Absolute Error");
    comparison_plot.x_label("t");
    comparison_plot.y_label("Error");

    // -------------------------------------------------------------------------
    // Save and display
    // -------------------------------------------------------------------------

    comparison_plot.save("implicit_euler_test.png");
    comparison_plot.show();

    // -------------------------------------------------------------------------
    // Final-time errors
    // -------------------------------------------------------------------------

    std::cout << "Scalar explicit Euler error at t = 2: " << explicit_error[steps] << '\n';
    std::cout << "Scalar implicit Euler error at t = 2: " << implicit_error[steps] << '\n';
    std::cout << "Scalar Heun error at t = 2: " << heun_error[steps] << '\n';

    std::cout << '\n';

    std::cout << "Vector implicit Euler x error at t = 2: " << vector_error_x[steps] << '\n';
    std::cout << "Vector implicit Euler y error at t = 2: " << vector_error_y[steps] << '\n';

    return 0;
}
