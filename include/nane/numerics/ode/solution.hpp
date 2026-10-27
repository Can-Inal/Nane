#pragma once

#include <vector>

namespace nane::ode
{
    /**
     * @ingroup ode
     *
     * @brief Stores the time points and values of an ODE solution.
     *
     * @tparam Values Type storing the numerical solution values.
     */
    template <typename Values>
    struct solution
    {
        std::vector<double> time;
        Values values;
    };

} // namespace nane::ode
