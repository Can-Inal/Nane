cmake_minimum_required(VERSION 3.20)

project(
    Nane
    VERSION 0.1.0
    DESCRIPTION "A C++ scientific computing library"
    LANGUAGES CXX
)

# ------------------------------------------------------------------------------
# C++ configuration
# ------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------------------------
# Output directories
# ------------------------------------------------------------------------------

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# ------------------------------------------------------------------------------
# compile_commands.json
# ------------------------------------------------------------------------------

set(
    CMAKE_EXPORT_COMPILE_COMMANDS
    ON
    CACHE BOOL
    "Generate compile_commands.json"
    FORCE
)

if(CMAKE_GENERATOR MATCHES "Ninja|Makefiles")
    file(
        CREATE_LINK
        "${CMAKE_BINARY_DIR}/compile_commands.json"
        "${CMAKE_SOURCE_DIR}/compile_commands.json"
        SYMBOLIC
        RESULT COMPILE_COMMANDS_LINK_RESULT
    )

    if(NOT COMPILE_COMMANDS_LINK_RESULT STREQUAL "0")
        message(
            WARNING
            "Could not create compile_commands.json symbolic link: "
            "${COMPILE_COMMANDS_LINK_RESULT}"
        )
    endif()
endif()

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
    OFF
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

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    include(CTest)

    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
