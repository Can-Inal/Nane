#pragma once

#include <stdexcept>

/**
 * @brief Requires a condition imposed on user input.
 *
 * Throws std::invalid_argument if the condition is not satisfied.
 */
#define NANE_REQUIRE(condition, message)                                                                                                             \
    do                                                                                                                                               \
    {                                                                                                                                                \
        if (!(condition)) [[unlikely]]                                                                                                               \
            throw std::invalid_argument(message);                                                                                                    \
    } while (false)

/**
 * @brief Ensures a required object or internal state.
 *
 * Throws std::logic_error if the condition is not satisfied.
 */
#define NANE_ENSURE(condition, message)                                                                                                              \
    do                                                                                                                                               \
    {                                                                                                                                                \
        if (!(condition)) [[unlikely]]                                                                                                               \
            throw std::logic_error(message);                                                                                                         \
    } while (false)

/**
 * @brief Reports a runtime or numerical failure.
 *
 * Throws std::runtime_error if the condition is satisfied.
 */
#define NANE_FAIL(condition, message)                                                                                                                \
    do                                                                                                                                               \
    {                                                                                                                                                \
        if (condition) [[unlikely]]                                                                                                                  \
            throw std::runtime_error(message);                                                                                                       \
    } while (false)
