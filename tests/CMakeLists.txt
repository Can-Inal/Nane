# ------------------------------------------------------------------------------
# Catch2
# ------------------------------------------------------------------------------

add_subdirectory(
    "${PROJECT_SOURCE_DIR}/vendor/catch2"
    "${CMAKE_BINARY_DIR}/vendor/catch2"
    EXCLUDE_FROM_ALL
)

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

add_executable(
    nane_tests
    ode.cpp
)

target_link_libraries(
    nane_tests
    PRIVATE
        nane::nane
        Catch2::Catch2WithMain
)

# ------------------------------------------------------------------------------
# Compiler warnings
# ------------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        nane_tests
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
elseif(MSVC)
    target_compile_options(
        nane_tests
        PRIVATE
            /W4
            /permissive-
    )
endif()

# ------------------------------------------------------------------------------
# CTest
# ------------------------------------------------------------------------------

list(
    APPEND
    CMAKE_MODULE_PATH
    "${PROJECT_SOURCE_DIR}/vendor/catch2/extras"
)

include(Catch)

catch_discover_tests(nane_tests)
