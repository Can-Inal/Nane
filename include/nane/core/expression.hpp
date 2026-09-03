#pragma once

#include <nane/core/details/expression.hpp>
#include <utility>

namespace nane
{
    template <typename Node>
    class expression : public details::expression_marker
    {
    public:
        explicit expression(Node node) : node_(std::move(node))
        {
        }

        template <typename... Arguments>
        [[nodiscard]] decltype(auto) operator()(Arguments&&... arguments) const
        {
            return node_(std::forward<Arguments>(arguments)...);
        }

    private:
        Node node_;
    };

} // namespace nane
