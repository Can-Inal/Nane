# ------------------------------------------------------------------------------
# Scalar ODE
# ------------------------------------------------------------------------------

add_executable(
    nane_example_scalar_ode
    scalar_ode.cpp
)

target_link_libraries(
    nane_example_scalar_ode
    PRIVATE
        nane::nane
)

# ------------------------------------------------------------------------------
# Coupled ODE
# ------------------------------------------------------------------------------

add_executable(
    nane_example_coupled_ode
    coupled_ode.cpp
)

target_link_libraries(
    nane_example_coupled_ode
    PRIVATE
        nane::nane
)
