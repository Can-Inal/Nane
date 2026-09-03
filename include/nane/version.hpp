#pragma once

#include <string_view>

namespace nane
{
    /**
     * @brief Returns the current Nane version.
     *
     * @return Nane version string in semantic version form.
     */
    [[nodiscard]] std::string_view version();

} // namespace nane
