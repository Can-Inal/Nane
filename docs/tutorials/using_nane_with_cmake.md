@page cmake_integration Using Nane with CMake

Nane exports the CMake target `nane::nane`. Consumer projects should link that
target rather than manually adding Nane's include directories or transitive
dependencies.

## Git submodule

For an application repository, a convenient layout is:

```text
project/
├── CMakeLists.txt
├── src/
└── vendor/
    └── Nane/
```

Add Nane with:

```bash
git submodule add https://github.com/Can-Inal/Nane.git vendor/Nane
git submodule update --init --recursive
```

When cloning the application on another machine, use:

```bash
git clone --recurse-submodules <repository>
```

or initialize submodules afterwards:

```bash
git submodule update --init --recursive
```

## Link Nane into an executable

```cmake
cmake_minimum_required(VERSION 3.20)

project(my_project LANGUAGES CXX)

add_subdirectory(vendor/Nane EXCLUDE_FROM_ALL)

add_executable(my_project
    src/main.cpp
)

target_compile_features(my_project PRIVATE cxx_std_23)

target_link_libraries(my_project
    PRIVATE
        nane::nane
)
```

`add_subdirectory` makes Nane's targets available to the parent project.
`EXCLUDE_FROM_ALL` prevents unrelated Nane targets from being built as part of
the application's default build.

Because public dependencies are attached to `nane::nane`, consumer projects do
not need to duplicate Nane's dependency include paths in their own CMake files.

## Build Nane statically

If Nane contains compiled sources, CMake's `BUILD_SHARED_LIBS` variable controls
the default library type when `add_library` does not explicitly choose
`STATIC` or `SHARED`.

To request a static Nane library from the parent project, set the option before
adding Nane:

```cmake
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(vendor/Nane EXCLUDE_FROM_ALL)
```

Then link normally:

```cmake
target_link_libraries(my_project
    PRIVATE
        nane::nane
)
```

The final executable will contain the linked Nane object code. Header-only
parts of Nane are compiled directly as part of the consuming translation
units.

If the parent project uses `BUILD_SHARED_LIBS` for its own libraries and you do
not want to change that policy globally, keep the setting local to the Nane
subdirectory using a cache value only when your project structure requires it.

## Compilation database for clangd

CMake can generate `compile_commands.json`:

```cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
```

A project can additionally expose the generated file at the repository root:

```cmake
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${CMAKE_BINARY_DIR}/compile_commands.json"
            "${CMAKE_SOURCE_DIR}/compile_commands.json"
    )
endif()
```

The generated database contains the real include paths, compiler options, and
C++ standard used for each translation unit, which lets `clangd` understand
both the application and Nane headers.

## Updating Nane

From the application repository:

```bash
cd vendor/Nane
git pull
cd ../..
git add vendor/Nane
git commit -m "deps: update Nane"
```

A submodule records a specific Nane commit, so application builds remain
reproducible until that pointer is intentionally updated.

## Prefer the namespaced target

Both the underlying target and alias may exist, but documentation and consumer
projects should prefer:

```cmake
nane::nane
```

The namespaced form clearly identifies the dependency as a CMake target owned
by Nane and produces clearer errors if the target has not been created.
