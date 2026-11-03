#pragma once

#include "nane/math/matrix.hpp"

#include <cstddef>

namespace nane::pde
{
    /**
     * @ingroup pde
     *
     * @brief Represents a local finite-difference operator.
     *
     * The finite-difference operator is defined by a matrix kernel centered
     * at its middle entry. The kernel must have an odd number of rows and
     * columns so that its center is uniquely defined.
     *
     * @tparam Rows Compile-time number of kernel rows or nane::dynamic.
     * @tparam Columns Compile-time number of kernel columns or nane::dynamic.
     */
    template <std::size_t Rows, std::size_t Columns = Rows>
    class finite_difference
    {
        static_assert(Rows == dynamic || Rows % 2 == 1, "finite difference kernel row count must be odd");
        static_assert(Columns == dynamic || Columns % 2 == 1, "finite difference kernel column count must be odd");

    public:
        /**
         * @brief Constructs a finite-difference operator from a kernel.
         *
         * The kernel must have an odd number of rows and columns.
         *
         * @param kernel Local finite-difference kernel.
         */
        finite_difference(const nane::matrix<Rows, Columns>& kernel) : m_Kernel(kernel)
        {
            NANE_REQUIRE(kernel.rows() % 2 == 1, "finite difference kernel row count must be odd");
            NANE_REQUIRE(kernel.columns() % 2 == 1, "finite difference kernel column count must be odd");
        }

        /**
         * @brief Returns the finite-difference kernel.
         *
         * @return Constant reference to the kernel.
         */
        [[nodiscard]] const nane::matrix<Rows, Columns>& kernel() const noexcept
        {
            return m_Kernel;
        }

    private:
        nane::matrix<Rows, Columns> m_Kernel;
    };

} // namespace nane::pde
