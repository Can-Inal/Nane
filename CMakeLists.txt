cmake_minimum_required(VERSION 3.20)

project(
    Nane
    VERSION 0.1.0
    DESCRIPTION "A C++ scientific computing library"
    LANGUAGES CXX
)

file(GLOB_RECURSE NANE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

if(NANE_SOURCES)
    add_library(nane ${NANE_SOURCES})
    set(NANE_SCOPE PUBLIC)
else()
    add_library(nane INTERFACE)
    set(NANE_SCOPE INTERFACE)
endif()

add_library(nane::nane ALIAS nane)

target_compile_features(nane ${NANE_SCOPE} cxx_std_23)

target_include_directories(nane ${NANE_SCOPE}
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Eigen

add_library(nane_eigen INTERFACE)

target_include_directories(nane_eigen SYSTEM INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen"
)

target_link_libraries(nane ${NANE_SCOPE} nane_eigen)

# Warnings

if(NANE_SOURCES)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(nane PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(nane PRIVATE /W4 /permissive-)
    endif()
endif()

# Development

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    include(CTest)

    if(BUILD_TESTING)
        add_subdirectory(vendor/catch2)
        list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/vendor/catch2/extras")
        add_subdirectory(tests)
    endif()

    option(NANE_BUILD_EXAMPLES "Build examples" OFF)

    if(NANE_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()
