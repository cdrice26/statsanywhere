# Tests CMake for statsanywhere

# Unity vendor
add_library(unity STATIC ${CMAKE_SOURCE_DIR}/tests/unity/src/unity.c)
target_compile_definitions(unity PRIVATE UNITY_INCLUDE_DOUBLE)
target_include_directories(unity PUBLIC ${CMAKE_SOURCE_DIR}/tests/unity/include)

aux_source_directory(${CMAKE_SOURCE_DIR}/tests TEST_SRCS)
# Explicitly list tests to avoid picking up unrelated files
set(TEST_SOURCES
    ${CMAKE_SOURCE_DIR}/tests/test_distributions.c
    ${CMAKE_SOURCE_DIR}/tests/test_manual.c
    ${CMAKE_SOURCE_DIR}/tests/test_runner.c
)

add_executable(test_runner ${TEST_SOURCES})
# Ensure Unity double support and provide generator/reference paths
target_compile_definitions(test_runner PRIVATE UNITY_INCLUDE_DOUBLE)
# Generator script absolute path (points to SciPy-only generator)
target_compile_definitions(test_runner PRIVATE GENERATOR_SCRIPT_PATH="${CMAKE_SOURCE_DIR}/tests/generate_reference.py")
# Preferred venv python path (posix)
target_compile_definitions(test_runner PRIVATE GENERATOR_VENV_PY="${CMAKE_SOURCE_DIR}/tests/.venv/bin/python")
# Absolute path to the reference values file (so tests can find it when run from build/)
target_compile_definitions(test_runner PRIVATE REFERENCE_FILE_PATH="${CMAKE_SOURCE_DIR}/tests/reference_values.txt")

target_include_directories(test_runner PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/tests/unity/include)

target_link_libraries(test_runner PRIVATE statsanywhere unity m)

# Add a helper target to generate the reference file (requires SciPy available in Python)
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_custom_target(generate_reference
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/generate_reference.py
        DEPENDS statsanywhere
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generate tests/reference_values.txt using SciPy (requires SciPy in Python).")
endif()

add_test(NAME unit_tests COMMAND test_runner)
