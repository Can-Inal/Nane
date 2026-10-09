# This file was generated with AI assistance and reviewed by the maintainer (Can Inal).


# Find all example source files recursively.
file(
    GLOB_RECURSE
    NANE_EXAMPLE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

foreach(source IN LISTS NANE_EXAMPLE_SOURCES)
    # Path relative to the examples directory.
    file(
        RELATIVE_PATH
        relative_source
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${source}"
    )

    # Example:
    #
    # relative_source = ode/scalar_ode.cpp
    # example_directory = ode
    # example_name = scalar_ode
    get_filename_component(
        example_directory
        "${relative_source}"
        DIRECTORY
    )

    get_filename_component(
        example_name
        "${source}"
        NAME_WE
    )

    # Create a unique CMake target name from the relative source path.
    #
    # ode/scalar_ode.cpp -> nane_example_ode_scalar_ode
    string(
        REPLACE
        ".cpp"
        ""
        target_suffix
        "${relative_source}"
    )

    string(
        REPLACE
        "/"
        "_"
        target_suffix
        "${target_suffix}"
    )

    set(target_name "nane_example_${target_suffix}")

    add_executable(
        ${target_name}
        "${source}"
    )

    target_link_libraries(
        ${target_name}
        PRIVATE
            nane::nane
    )

    # Preserve the source directory structure below build/bin/examples.
    set(
        output_directory
        "${CMAKE_BINARY_DIR}/bin/examples/${example_directory}"
    )

    set_target_properties(
        ${target_name}
        PROPERTIES
            OUTPUT_NAME
                "${example_name}"
            RUNTIME_OUTPUT_DIRECTORY
                "${output_directory}"
    )

    # Prevent multi-config generators such as Visual Studio from adding
    # Debug/, Release/, etc. below the desired output directory.
    foreach(config IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER "${config}" config_upper)

        set_target_properties(
            ${target_name}
            PROPERTIES
                "RUNTIME_OUTPUT_DIRECTORY_${config_upper}"
                    "${output_directory}"
        )
    endforeach()
endforeach()
