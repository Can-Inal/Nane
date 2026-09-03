# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------

file(
    GLOB_RECURSE STUDIO_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# ------------------------------------------------------------------------------
# Executable
# ------------------------------------------------------------------------------

add_executable(
    Studio
    ${STUDIO_SOURCES}
)

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

target_link_libraries(
    Studio
    PRIVATE
        nane::nane
        nane::viz
)

# ------------------------------------------------------------------------------
# Compiler warnings
# ------------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        Studio
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
elseif(MSVC)
    target_compile_options(
        Studio
        PRIVATE
            /W4
            /permissive-
    )
endif()
