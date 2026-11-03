#include "nane/numerics/pde/assemble.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Finite difference assembly with Dirichlet boundaries", "[pde][assemble]")
{
    const nane::geometry::grid<2> grid{
        {0.0, 2.0, 3},
        {0.0, 2.0, 3},
    };

    const nane::pde::finite_difference equation{
        nane::matrix<3>{
            {0.0, 1.0, 0.0},
            {1.0, -4.0, 1.0},
            {0.0, 1.0, 0.0},
        },
    };

    const nane::pde::boundaries<2> boundary_conditions{
        {nane::pde::dirichlet{0.0}, nane::pde::dirichlet{0.0}},
        {nane::pde::dirichlet{0.0}, nane::pde::dirichlet{0.0}},
    };

    const auto system = nane::pde::assemble(equation, grid, boundary_conditions);

    REQUIRE(system.matrix.rows() == 9);
    REQUIRE(system.matrix.columns() == 9);
    REQUIRE(system.rhs.size() == 9);

    nane::matrix<> expected(9, 9, 0.0);

    expected(4, 1) = 1.0;
    expected(4, 3) = 1.0;
    expected(4, 4) = -4.0;
    expected(4, 5) = 1.0;
    expected(4, 7) = 1.0;

    for (std::size_t i = 0; i < 9; ++i)
    {
        REQUIRE(system.rhs[i] == 0.0);

        for (std::size_t j = 0; j < 9; ++j)
            REQUIRE(system.matrix(i, j) == expected(i, j));
    }
}

TEST_CASE("Finite difference assembly with homogeneous Neumann boundaries", "[pde][assemble]")
{
    const nane::geometry::grid<2> grid{
        {0.0, 2.0, 3},
        {0.0, 2.0, 3},
    };

    const nane::pde::finite_difference equation{
        nane::matrix<3>{
            {0.0, 1.0, 0.0},
            {1.0, -4.0, 1.0},
            {0.0, 1.0, 0.0},
        },
    };

    const nane::pde::boundaries<2> boundary_conditions{
        {nane::pde::neumann{0.0}, nane::pde::neumann{0.0}},
        {nane::pde::neumann{0.0}, nane::pde::neumann{0.0}},
    };

    const auto system = nane::pde::assemble(equation, grid, boundary_conditions);

    // Lower-left corner.
    REQUIRE(system.matrix(0, 0) == -4.0);
    REQUIRE(system.matrix(0, 1) == 2.0);
    REQUIRE(system.matrix(0, 3) == 2.0);

    // Lower edge.
    REQUIRE(system.matrix(1, 0) == 1.0);
    REQUIRE(system.matrix(1, 1) == -4.0);
    REQUIRE(system.matrix(1, 2) == 1.0);
    REQUIRE(system.matrix(1, 4) == 2.0);

    // Left edge.
    REQUIRE(system.matrix(3, 0) == 1.0);
    REQUIRE(system.matrix(3, 3) == -4.0);
    REQUIRE(system.matrix(3, 4) == 2.0);
    REQUIRE(system.matrix(3, 6) == 1.0);

    // Interior.
    REQUIRE(system.matrix(4, 1) == 1.0);
    REQUIRE(system.matrix(4, 3) == 1.0);
    REQUIRE(system.matrix(4, 4) == -4.0);
    REQUIRE(system.matrix(4, 5) == 1.0);
    REQUIRE(system.matrix(4, 7) == 1.0);

    // Upper-right corner.
    REQUIRE(system.matrix(8, 5) == 2.0);
    REQUIRE(system.matrix(8, 7) == 2.0);
    REQUIRE(system.matrix(8, 8) == -4.0);

    for (std::size_t i = 0; i < 9; ++i)
        REQUIRE(system.rhs[i] == 0.0);
}

TEST_CASE("Finite difference assembly adds Neumann values to right-hand side", "[pde][assemble]")
{
    const nane::geometry::grid<2> grid{
        {0.0, 2.0, 3},
        {0.0, 2.0, 3},
    };

    const nane::pde::finite_difference equation{
        nane::matrix<3>{
            {0.0, 1.0, 0.0},
            {1.0, -4.0, 1.0},
            {0.0, 1.0, 0.0},
        },
    };

    const nane::pde::boundaries<2> boundary_conditions{
        {nane::pde::neumann{1.5}, nane::pde::dirichlet{0.0}},
        {nane::pde::dirichlet{0.0}, nane::pde::dirichlet{0.0}},
    };

    const auto system = nane::pde::assemble(equation, grid, boundary_conditions);

    // Point (0, 1), global index 3.
    REQUIRE(system.matrix(3, 0) == 1.0);
    REQUIRE(system.matrix(3, 3) == -4.0);
    REQUIRE(system.matrix(3, 4) == 2.0);
    REQUIRE(system.matrix(3, 6) == 1.0);

    // h = 1, q = 1.5:
    //
    // 2 h q = 3.
    REQUIRE(system.rhs[3] == 3.0);

    // Remaining boundary nodes are Dirichlet.
    REQUIRE(system.rhs[0] == 0.0);
    REQUIRE(system.rhs[2] == 0.0);
    REQUIRE(system.rhs[6] == 0.0);
    REQUIRE(system.rhs[8] == 0.0);
}
