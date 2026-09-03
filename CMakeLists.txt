cmake_minimum_required(VERSION 3.20)

project(
    Heron
    VERSION 0.1.0
    DESCRIPTION "A scientific computing software"
    LANGUAGES C CXX
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
# Subdirectories
# ------------------------------------------------------------------------------

add_subdirectory(lib)
add_subdirectory(viz)
add_subdirectory(studio)
