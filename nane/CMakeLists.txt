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

add_library(nane ${NANE_SOURCES} ${NANE_HEADERS})
add_library(nane::nane ALIAS nane)

target_include_directories(nane PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")


# ------------------------------------------------------------------------------
# Eigen
# ------------------------------------------------------------------------------
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
add_subdirectory(vendor/eigen EXCLUDE_FROM_ALL SYSTEM)
set(CMAKE_WARN_DEPRECATED ON CACHE BOOL "" FORCE)


target_link_libraries(nane PUBLIC
        Eigen3::Eigen
)

# ------------------------------------------------------------------------------
# Compiler warnings
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# OpenMP
# ------------------------------------------------------------------------------

option(HERON_ENABLE_OPENMP "Enable OpenMP support" ON)

if(HERON_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(nane PRIVATE OpenMP::OpenMP_CXX)
endif() 
