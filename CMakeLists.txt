cmake_minimum_required(VERSION 3.15...4.2)
project(rsim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)

FetchContent_MakeAvailable(Eigen)

option(RSIM_BUILD_PYTHON "Build the Python extension module" ON)
option(RSIM_BUILD_DEMO "Build the C++ demonstration executable" OFF)
option(RSIM_BUILD_TESTS "Build the C++ tests" OFF)

add_library(rsim_core
    src/Frames.cc
    src/Vehicle.cc
)
target_include_directories(rsim_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
target_link_libraries(rsim_core PUBLIC Eigen3::Eigen)
set_target_properties(rsim_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(RSIM_BUILD_PYTHON)
    find_package(Python REQUIRED COMPONENTS Interpreter Development)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_rsim MODULE python/bindings.cpp)
    target_link_libraries(_rsim PRIVATE rsim_core)
    install(TARGETS _rsim DESTINATION rsim)
endif()

if(RSIM_BUILD_TESTS)
    enable_testing()
    add_executable(rsim_cpp_tests tests/test_mass_properties.cc)
    target_link_libraries(rsim_cpp_tests PRIVATE rsim_core)
    add_test(NAME rsim_cpp_tests COMMAND rsim_cpp_tests)

    add_executable(rsim_frame_tests tests/test_frames.cc)
    target_link_libraries(rsim_frame_tests PRIVATE rsim_core)
    add_test(NAME rsim_frame_tests COMMAND rsim_frame_tests)
endif()
