#pragma once

#include "nane/error.hpp"
#include "nane/geometry/grid.hpp"
#include "nane/math/matrix.hpp"
#include "nane/math/vector.hpp"
#include "nane/numerics/calculus/differentiation/stencil.hpp"

#include <cmath>
#include <cstddef>

namespace nane::calculus
{
    /**
     * @ingroup calculus
     *
     * @brief Approximates a scalar derivative using a finite-difference stencil.
     *
     * For a stencil of derivative order @f$m@f$,
     *
     * @f[
     * f^{(m)}(x)
     * \approx
     * \frac{1}{h^m}
     * \sum_{i=0}^{s-1}
     * w_i f(x + c_i h),
     * @f]
     *
     * where @f$c_i@f$ are the stencil offsets, @f$w_i@f$ are the
     * corresponding weights, and @f$h@f$ is the step size.
     *
     * @param function Function to differentiate.
     * @param point Point at which the derivative is approximated.
     * @param stencil Finite-difference stencil.
     *
     * @return Approximation of the derivative.
     */
    [[nodiscard]] inline double differentiate(auto&& function, double point, const nane::calculus::stencil& stencil)
    {
        const auto& offsets = stencil.offsets();
        const auto& weights = stencil.weights();

        const double step_size = stencil.step_size();

        double sum = 0.0;

        for (std::size_t i = 0; i < stencil.size(); ++i)
        {
            const double x = point + static_cast<double>(offsets[i]) * step_size;
            sum += weights[i] * function(x);
        }

        const double scale = std::pow(step_size, static_cast<double>(stencil.derivative_order()));
        return sum / scale;
    }

    /**
     * @ingroup calculus
     *
     * @brief Approximates the derivative of a vector-valued function.
     *
     * For a function
     *
     * @f[
     * \mathbf{F}:\mathbb{R}^n \rightarrow \mathbb{R}^m,
     * @f]
     *
     * the derivative is the Jacobian matrix
     *
     * @f[
     * J_{ij}
     * =
     * \frac{\partial F_i}{\partial x_j}.
     * @f]
     *
     * Each column is approximated by applying the finite-difference
     * stencil in the corresponding coordinate direction.
     *
     * @tparam Function Type of the vector-valued function.
     * @tparam Dimension Input dimension.
     *
     * @param function Vector-valued function.
     * @param point Point at which the derivative is approximated.
     * @param stencil First-derivative finite-difference stencil.
     *
     * @return Approximation of the Jacobian matrix.
     */
    template <std::size_t Dimension>
    [[nodiscard]] auto differentiate(auto&& function, const nane::vector<Dimension>& point, const nane::calculus::stencil& stencil)
    {
        NANE_REQUIRE(stencil.derivative_order() == 1, "vector differentiation requires a first-derivative stencil");

        const auto value = function(point); // evaluate the function once to get the .size() dimension info.
        nane::matrix<> result(value.size(), point.size());

        const auto& offsets = stencil.offsets();
        const auto& weights = stencil.weights();

        const double step_size = stencil.step_size();

        // Differentiate with respect to each input component.
        // Column.
        for (std::size_t j = 0; j < point.size(); ++j)
        {
            // Stencil point.
            for (std::size_t k = 0; k < stencil.size(); ++k)
            {
                auto x = point;

                x[j] += static_cast<double>(offsets[k]) * step_size;
                const auto next = function(x);

                for (std::size_t i = 0; i < value.size(); ++i)
                    result(i, j) += weights[k] * next[i];
            }

            const double scale = std::pow(step_size, static_cast<double>(stencil.derivative_order()));
            for (std::size_t i = 0; i < value.size(); ++i)
                result(i, j) /= scale;
        }

        return result;
    }

    /**
     * @ingroup calculus
     *
     * @brief Approximates the derivative of sampled one-dimensional data.
     *
     * For sampled values
     *
     * @f[
     * u_i = u(x_i),
     * @f]
     *
     * and a stencil of derivative order @f$m@f$,
     *
     * @f[
     * u^{(m)}(x_i)
     * \approx
     * \frac{1}{h^m}
     * \sum_{k=0}^{s-1}
     * w_k u_{i+c_k},
     * @f]
     *
     * where @f$c_k@f$ are the stencil offsets, @f$w_k@f$ are the
     * corresponding weights, and @f$h@f$ is the grid spacing.
     *
     * The stencil is applied only where all required samples lie
     * inside the grid.
     *
     * @tparam Dimension Number of sampled values.
     *
     * @param values Sampled function values.
     * @param grid One-dimensional grid associated with the samples.
     * @param stencil Finite-difference stencil.
     *
     * @return Approximation of the derivative at the grid points.
     */
    template <std::size_t Dimension>
    [[nodiscard]] nane::vector<Dimension> differentiate(const nane::vector<Dimension>& values, const nane::grid<1>& grid,
                                                        const nane::calculus::stencil& stencil)
    {
        // values = [u0, u1, u2, ..., uN-1];
        // grid = [x0, x1, x2, ..., xN-1];

        NANE_REQUIRE(grid.has_count(0), "sampled differentiation requires a prescribed grid");
        const double h = *grid.spacing(0);

        const auto N = grid.axis(0).size();
        NANE_REQUIRE(values.size() == N, "number of values must match grid point count");

        const auto& offsets = stencil.offsets();
        const auto& weights = stencil.weights();
        const auto& order = stencil.derivative_order();

        nane::vector<Dimension> result(N, 0);

        // loop over grid.
        for (std::size_t i = 0; i < N; ++i)
        {
            bool valid = true;
            double sum = 0.0;

            // loop over stencil.
            for (std::size_t k = 0; k < weights.size(); ++k)
            {
                // check if offset is valid. It shouldn't be negative or greater than or equal to N.
                const auto offset = static_cast<int>(i) + offsets[k];
                valid = offset >= 0 && offset < N;

                if (!valid)
                    break;

                sum += weights[k] * values[offset];
            }

            // if the offset is valid, calculate the result.
            if (valid)
                result[i] = sum / std::pow(h, static_cast<double>(order));
        }

        return result;
    }

    /**
     * @ingroup calculus
     *
     * @brief Approximates partial derivatives of sampled two-dimensional data.
     *
     * For sampled values
     *
     * @f[
     * u_{ij} = u(x_i, y_j),
     * @f]
     *
     * and a stencil of derivative order @f$m@f$, differentiation along
     * axis 0 is approximated by
     *
     * @f[
     * \frac{\partial^m u}{\partial x^m}(x_i, y_j)
     * \approx
     * \frac{1}{h^m}
     * \sum_{k=0}^{s-1}
     * w_k u_{i+c_k,j},
     * @f]
     *
     * while differentiation along axis 1 is approximated by
     *
     * @f[
     * \frac{\partial^m u}{\partial y^m}(x_i, y_j)
     * \approx
     * \frac{1}{h^m}
     * \sum_{k=0}^{s-1}
     * w_k u_{i,j+c_k}.
     * @f]
     *
     * where @f$c_k@f$ are the stencil offsets, @f$w_k@f$ are the
     * corresponding weights, and @f$h@f$ is the grid spacing along
     * the selected axis.
     *
     * The stencil is applied only where all required samples lie
     * inside the grid.
     *
     * @tparam Rows Number of matrix rows.
     * @tparam Columns Number of matrix columns.
     *
     * @param values Sampled function values.
     * @param grid Two-dimensional grid associated with the samples.
     * @param stencil Finite-difference stencil.
     * @param axis Axis along which the derivative is approximated.
     *
     * @return Approximation of the partial derivative at the grid points.
     */
    template <std::size_t Rows, std::size_t Columns>
    [[nodiscard]] nane::matrix<Rows, Columns> differentiate(const nane::matrix<Rows, Columns>& values, const nane::grid<2>& grid,
                                                            const nane::calculus::stencil& stencil, std::size_t axis = 0)
    {
        NANE_REQUIRE(axis < 2, "matrix differentiation axis must be 0 or 1");
        NANE_REQUIRE(grid.has_count(axis), "sampled differentiation requires a prescribed grid");

        const double h = *grid.spacing(axis);

        const auto N = grid.axis(axis).size();

        if (axis == 0)
            NANE_REQUIRE(values.rows() == N, "matrix row count must match grid point count");
        else
            NANE_REQUIRE(values.columns() == N, "matrix column count must match grid point count");

        const auto& offsets = stencil.offsets();
        const auto& weights = stencil.weights();
        const auto& order = stencil.derivative_order();

        nane::matrix<Rows, Columns> result(values.rows(), values.columns(), 0.0);

        // iterate over rows.
        for (std::size_t i = 0; i < result.rows(); ++i)
        {
            // iterate over columns.
            for (std::size_t j = 0; j < result.columns(); ++j)
            {
                bool valid = true;
                double sum = 0.0;

                // iterate over stencil for each element (row, column) of matrix.
                for (std::size_t k = 0; k < weights.size(); ++k)
                {
                    // if the axis == 0 (x) we use row index (i), otherwise column index (j).
                    const auto index = static_cast<int>(axis == 0 ? i : j) + offsets[k];

                    // check that index is still inside the grid.
                    if (index < 0 || index >= N)
                    {
                        valid = false;
                        break;
                    }

                    // based on the axis we use either:
                    // u_(i + c_k, j) => for axis == 0 (x axis.)
                    // or
                    // u_(i, j + c_k) => for axis == 1 (y axis.)
                    sum += weights[k] * ((axis == 0) ? values(index, j) : values(i, index));
                }

                // final result after dividing by h^m.
                if (valid)
                    result(i, j) = sum / std::pow(h, static_cast<double>(order));
            }
        }

        return result;
    }

} // namespace nane::calculus
