#include "nane/numerics/calculus/integration/exactness.hpp"
#include "nane/numerics/calculus/integration/rules.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Integration methods with their exactness", "[calculus][integration][exactness]")
{
    REQUIRE(nane::calculus::degree_of_exactness(nane::calculus::rules::trapezoidal) == 1);
    REQUIRE(nane::calculus::degree_of_exactness(nane::calculus::rules::simpson) == 3);
    REQUIRE(nane::calculus::degree_of_exactness(nane::calculus::rules::gauss_legendre_2) == 3);
    REQUIRE(nane::calculus::degree_of_exactness(nane::calculus::rules::gauss_legendre_3) == 5);
}
