file(
    GLOB_RECURSE NANE_TEST_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_executable(
    nane_tests
    ${NANE_TEST_SOURCES}
)

target_link_libraries(
    nane_tests
    PRIVATE
        nane::nane
        Catch2::Catch2WithMain
)

include(Catch)

catch_discover_tests(nane_tests)
