# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------

file(
    GLOB_RECURSE NANE_VIZ_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

file(
    GLOB_RECURSE NANE_VIZ_HEADERS
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
)

# ------------------------------------------------------------------------------
# VTK
# ------------------------------------------------------------------------------

find_package(
    VTK REQUIRED
    COMPONENTS
        ChartsCore
        CommonCore
        CommonDataModel
        IOImage
        InteractionStyle
        RenderingContextOpenGL2
        RenderingCore
        RenderingFreeType
        RenderingOpenGL2
        ViewsContext2D
)

# ------------------------------------------------------------------------------
# Library
# ------------------------------------------------------------------------------

add_library(
    nane_viz
    ${NANE_VIZ_SOURCES}
    ${NANE_VIZ_HEADERS}
)

add_library(nane::viz ALIAS nane_viz)

target_include_directories(
    nane_viz
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_compile_features(
    nane_viz
    PUBLIC
        cxx_std_23
)

target_link_libraries(
    nane_viz
    PRIVATE
        nane::nane
        ${VTK_LIBRARIES}
)

# ------------------------------------------------------------------------------
# VTK initialization
# ------------------------------------------------------------------------------

vtk_module_autoinit(
    TARGETS
        nane_viz
    MODULES
        ${VTK_LIBRARIES}
)
