#pragma once

#include "nane/error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <vector>

namespace nane
{
    inline constexpr std::size_t dynamic = std::numeric_limits<std::size_t>::max();

    /**
     * @brief Numerical vector.
     *
     * @tparam Dimension Compile-time dimension or nane::dynamic.
     */
    template <std::size_t Dimension = dynamic>
    class vector
    {
    public:
        vector() = default;

        /**
         * @brief Constructs a vector with the given dimension.
         *
         * For a fixed-size vector, the given dimension must match the
         * compile-time dimension.
         *
         * @param dimension Vector dimension.
         */
        explicit vector(std::size_t dimension)
        {
            if constexpr (s_Dynamic)
                m_Data.resize(dimension);
            else
                NANE_REQUIRE(dimension == Dimension, "vector dimension does not match");
        }

        /**
         * @brief Constructs a vector with the given dimension and value.
         *
         * @param dimension Vector dimension.
         * @param value Initial value of every component.
         */
        vector(std::size_t dimension, double value) : vector(dimension)
        {
            fill(value);
        }

        /**
         * @brief Constructs a vector from an initializer list.
         *
         * @param values Initial vector components.
         */
        vector(std::initializer_list<double> values)
        {
            if constexpr (s_Dynamic)
            {
                m_Data.assign(values.begin(), values.end());
            }
            else
            {
                NANE_REQUIRE(values.size() == Dimension, "vector initializer dimension does not match");

                std::copy(values.begin(), values.end(), m_Data.begin());
            }
        }

        /**
         * @brief Returns the vector dimension.
         */
        [[nodiscard]] constexpr std::size_t size() const noexcept
        {
            if constexpr (s_Dynamic)
                return m_Data.size();

            return Dimension;
        }

        /**
         * @brief Returns the compile-time dimension.
         */
        [[nodiscard]] static constexpr std::size_t dimension() noexcept
        {
            return Dimension;
        }

        /**
         * @brief Returns whether the vector has dynamic dimension.
         */
        [[nodiscard]] static constexpr bool is_dynamic() noexcept
        {
            return s_Dynamic;
        }

        [[nodiscard]] double& operator[](std::size_t index) noexcept
        {
            return m_Data[index];
        }

        [[nodiscard]] const double& operator[](std::size_t index) const noexcept
        {
            return m_Data[index];
        }

        /**
         * @brief Adds two vectors.
         */
        [[nodiscard]] vector operator+(const vector& other) const
        {
            NANE_REQUIRE(size() == other.size(), "vector dimensions must match");

            vector result(size());

            for (std::size_t i = 0; i < size(); ++i)
                result[i] = m_Data[i] + other[i];

            return result;
        }

        /**
         * @brief Subtracts two vectors.
         */
        [[nodiscard]] vector operator-(const vector& other) const
        {
            NANE_REQUIRE(size() == other.size(), "vector dimensions must match");

            vector result(size());

            for (std::size_t i = 0; i < size(); ++i)
                result[i] = m_Data[i] - other[i];

            return result;
        }

        /**
         * @brief Multiplies the vector by a scalar.
         */
        [[nodiscard]] vector operator*(double scalar) const
        {
            vector result(size());

            for (std::size_t i = 0; i < size(); ++i)
                result[i] = m_Data[i] * scalar;

            return result;
        }

        /**
         * @brief Divides the vector by a scalar.
         */
        [[nodiscard]] vector operator/(double scalar) const
        {
            vector result(size());

            for (std::size_t i = 0; i < size(); ++i)
                result[i] = m_Data[i] / scalar;

            return result;
        }

        /**
         * @brief Computes the dot product with another vector.
         *
         * @param other Other vector.
         *
         * @return Dot product of both vectors.
         */
        [[nodiscard]] double dot(const vector& other) const
        {
            NANE_REQUIRE(size() == other.size(), "vector dimensions must match");

            double result = 0.0;

            for (std::size_t i = 0; i < size(); ++i)
                result += m_Data[i] * other[i];

            return result;
        }

        /**
         * @brief Computes the Euclidean norm of the vector.
         */
        [[nodiscard]] double norm() const
        {
            return std::sqrt(dot(*this));
        }

        /**
         * @brief Sets every component to the same value.
         */
        void fill(double value)
        {
            std::fill(m_Data.begin(), m_Data.end(), value);
        }

        [[nodiscard]] double* data() noexcept
        {
            return m_Data.data();
        }

        [[nodiscard]] const double* data() const noexcept
        {
            return m_Data.data();
        }

        [[nodiscard]] auto begin() noexcept
        {
            return m_Data.begin();
        }

        [[nodiscard]] auto begin() const noexcept
        {
            return m_Data.begin();
        }

        [[nodiscard]] auto end() noexcept
        {
            return m_Data.end();
        }

        [[nodiscard]] auto end() const noexcept
        {
            return m_Data.end();
        }

        friend vector operator*(double scalar, const vector& value)
        {
            return value * scalar;
        }

    private:
        static constexpr bool s_Dynamic = Dimension == dynamic;
        static constexpr std::size_t s_StaticSize = s_Dynamic ? 0 : Dimension;

        using storage_type = std::conditional_t<s_Dynamic, std::vector<double>, std::array<double, s_StaticSize>>;

        storage_type m_Data{};
    };

    vector(std::initializer_list<double>) -> vector<dynamic>;
    vector(std::size_t) -> vector<dynamic>;
    vector(std::size_t, double) -> vector<dynamic>;

} // namespace nane
