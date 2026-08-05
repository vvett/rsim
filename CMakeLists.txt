cmake_minimum_required(VERSION 3.15...4.2)
project(rsim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(RSIM_BUILD_PYTHON "Build the Python extension module" ON)
option(RSIM_BUILD_DEMO "Build the C++ demonstration executable" OFF)
option(RSIM_BUILD_TESTS "Build the C++ tests" OFF)

add_library(rsim_core src/MassProperties.cc)
target_include_directories(rsim_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
set_target_properties(rsim_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(RSIM_BUILD_PYTHON)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_rsim MODULE python/bindings.cpp)
    target_link_libraries(_rsim PRIVATE rsim_core)
    install(TARGETS _rsim DESTINATION rsim)
endif()

if(RSIM_BUILD_DEMO)
    add_executable(rsim_demo src/main.cc)
    target_link_libraries(rsim_demo PRIVATE rsim_core)
endif()

if(RSIM_BUILD_TESTS)
    enable_testing()
    add_executable(rsim_cpp_tests tests/test_mass_properties.cc)
    target_link_libraries(rsim_cpp_tests PRIVATE rsim_core)
    add_test(NAME rsim_cpp_tests COMMAND rsim_cpp_tests)
endif()
