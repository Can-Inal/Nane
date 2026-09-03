# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------

file(
    GLOB_RECURSE NANE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

file(
    GLOB_RECURSE NANE_HEADERS
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
)

# ------------------------------------------------------------------------------
# Library
# ------------------------------------------------------------------------------

if(NANE_SOURCES)
    add_library(
        nane
        ${NANE_SOURCES}
        ${NANE_HEADERS}
    )

    set(NANE_USAGE_SCOPE PUBLIC)
else()
    add_library(nane INTERFACE)

    target_sources(
        nane
        INTERFACE
            ${NANE_HEADERS}
    )

    set(NANE_USAGE_SCOPE INTERFACE)
endif()

add_library(nane::nane ALIAS nane)

# ------------------------------------------------------------------------------
# Include directories
# ------------------------------------------------------------------------------

target_include_directories(
    nane
    ${NANE_USAGE_SCOPE}
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# ------------------------------------------------------------------------------
# C++ configuration
# ------------------------------------------------------------------------------

target_compile_features(
    nane
    ${NANE_USAGE_SCOPE}
        cxx_std_23
)

# ------------------------------------------------------------------------------
# Eigen
# ------------------------------------------------------------------------------

set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

add_subdirectory(
    vendor/eigen
    EXCLUDE_FROM_ALL
    SYSTEM
)

set(CMAKE_WARN_DEPRECATED ON CACHE BOOL "" FORCE)

target_link_libraries(
    nane
    ${NANE_USAGE_SCOPE}
        Eigen3::Eigen
)

# ------------------------------------------------------------------------------
# Compiler warnings
# ------------------------------------------------------------------------------

if(NANE_SOURCES)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(
            nane
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(
            nane
            PRIVATE
                /W4
                /permissive-
        )
    endif()
endif()

# ------------------------------------------------------------------------------
# OpenMP
# ------------------------------------------------------------------------------

option(
    NANE_ENABLE_OPENMP
    "Enable OpenMP support"
    ON
)

if(NANE_ENABLE_OPENMP)
    find_package(
        OpenMP
        REQUIRED
        COMPONENTS CXX
    )

    target_link_libraries(
        nane
        ${NANE_USAGE_SCOPE}
            OpenMP::OpenMP_CXX
    )
endif()
