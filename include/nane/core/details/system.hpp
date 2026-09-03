#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <tuple>
#include <utility>

namespace nane::details
{
    template <typename... Expressions>
    class system_function
    {
    public:
        explicit system_function(Expressions... expressions) : expressions_(std::move(expressions)...)
        {
        }

        template <typename... Arguments>
        [[nodiscard]] Eigen::VectorXd operator()(Arguments&&... arguments) const
        {
            Eigen::VectorXd result(static_cast<Eigen::Index>(sizeof...(Expressions)));

            evaluate(result, std::index_sequence_for<Expressions...>{}, arguments...);

            return result;
        }

    private:
        template <std::size_t... Indices, typename... Arguments>
        void evaluate(Eigen::VectorXd& result, std::index_sequence<Indices...>, Arguments&&... arguments) const
        {
            ((result[static_cast<Eigen::Index>(Indices)] = std::get<Indices>(expressions_)(arguments...)), ...);
        }

        std::tuple<Expressions...> expressions_;
    };

} // namespace nane::details
