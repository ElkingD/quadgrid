cmake_minimum_required(VERSION 3.10)
project(quadgrid VERSION 1.0 LANGUAGES CXX)

# Compiler settings
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Global compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -O3
  )
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include directory for public headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files for the library
file(GLOB QUADGRID_SOURCES "src/*.cpp")

# Create the static library
add_library(quadgrid STATIC ${QUADGRID_SOURCES})

# =====================
# Automatically build all demo apps in /app
# =====================
file(GLOB QUADGRID_APPS "app/*.cpp")
foreach(app_file ${QUADGRID_APPS})
    get_filename_component(app_name ${app_file} NAME_WE)
    add_executable(${app_name} ${app_file})
    target_link_libraries(${app_name} quadgrid)
endforeach()

# =====================
# Automatically build all tests in /test
# =====================
file(GLOB QUADGRID_TESTS "test/*.cpp")
foreach(test_file ${QUADGRID_TESTS})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} quadgrid)
endforeach()

# Install headers and library (optional for user installs)
install(DIRECTORY include/quadgrid DESTINATION include)
install(TARGETS quadgrid ARCHIVE DESTINATION lib)

