cmake_minimum_required(VERSION 3.16)
project(statsanywhere C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Option to build shared library if desired
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Sources
file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.c")

# Create library (static by default; set BUILD_SHARED_LIBS=ON to create shared)
add_library(statsanywhere ${SRC_FILES})

target_include_directories(statsanywhere PUBLIC ${CMAKE_SOURCE_DIR}/src)

# Tests
enable_testing()
add_subdirectory(tests)

# Install (optional)
install(TARGETS statsanywhere
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
