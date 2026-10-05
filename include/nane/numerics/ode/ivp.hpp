#pragma once

#include <utility>

namespace nane
{
    /**
     * @ingroup ode
     *
     * @brief Represents an initial value problem for an ordinary
     * differential equation.
     *
     * An initial value problem consists of a differential equation
     *
     * @f[
     * \dot{x}(t) = f(t, x(t))
     * @f]
     *
     * together with an initial value.
     *
     * The derivative may be either a scalar-valued function or a
     * vector-valued system. Correspondingly, the initial value may be
     * either scalar-valued or vector-valued.
     *
     * The initial time is determined by the first point of the time grid
     * supplied to the numerical method.
     *
     * @tparam Derivative Type of the right-hand-side function.
     * @tparam InitialValue Type of the initial value.
     */
    template <typename Derivative, typename InitialValue>
    class ivp
    {
    public:
        using derivative_type = Derivative;
        using initial_value_type = InitialValue;

        /**
         * @brief Constructs an initial value problem.
         *
         * @param derivative Right-hand-side function of the differential
         * equation.
         * @param initial_value Initial value of the solution.
         */
        ivp(Derivative derivative, InitialValue initial_value) : derivative_(std::move(derivative)), initial_value_(std::move(initial_value))
        {
        }

        /**
         * @brief Returns the right-hand-side function.
         *
         * @return Constant reference to the derivative function.
         */
        [[nodiscard]] const Derivative& derivative() const noexcept
        {
            return derivative_;
        }

        /**
         * @brief Returns the initial value.
         *
         * @return Constant reference to the initial value.
         */
        [[nodiscard]] const InitialValue& initial_value() const noexcept
        {
            return initial_value_;
        }

    private:
        Derivative derivative_;
        InitialValue initial_value_;
    };

} // namespace nane
