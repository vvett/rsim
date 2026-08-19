cmake_minimum_required(VERSION 3.15...4.2)
project(rsim VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

set(BUILD_TESTING OFF CACHE BOOL "Do not build dependency test suites" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "Do not build Eigen documentation" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "Do not register Eigen test suites" FORCE)

FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)

set(MATIO_BUILD_TESTING OFF CACHE BOOL "Do not build matio tests" FORCE)
set(MATIO_MAT73 OFF CACHE BOOL "MAT-v7.3/HDF5 support is not needed" FORCE)
set(MATIO_WITH_HDF5 OFF CACHE BOOL "Do not add the HDF5 dependency" FORCE)
set(MATIO_WITH_ZLIB OFF CACHE BOOL "Generated MAT-v5 files are uncompressed" FORCE)
set(MATIO_SHARED OFF CACHE BOOL "Link matio statically into rsim" FORCE)
set(MATIO_PIC ON CACHE BOOL "Allow matio to link into the Python extension" FORCE)
set(PACKAGE_TESTS OFF CACHE BOOL "Do not build mlinterp tests" FORCE)
FetchContent_Declare(
    matio
    GIT_REPOSITORY https://github.com/tbeu/matio.git
    GIT_TAG v1.5.30
    GIT_SHALLOW TRUE
    GIT_SUBMODULES ""
)

FetchContent_Declare(
    mlinterp
    GIT_REPOSITORY https://github.com/parsiad/mlinterp.git
    GIT_TAG e3a9a9cb20cc87cd12eea4790ac0dec90950db38
    GIT_SUBMODULES ""
)

FetchContent_MakeAvailable(Eigen matio mlinterp)

option(RSIM_BUILD_PYTHON "Build the Python extension module" ON)
option(RSIM_BUILD_DEMO "Build the C++ demonstration executable" OFF)
option(RSIM_BUILD_TESTS "Build the C++ tests" OFF)

add_library(rsim_core
    src/Actuators.cc
    src/Aerodynamics.cc
    src/DataLogger.cc
    src/EarthFrames.cc
    src/Environment.cc
    src/FlightState.cc
    src/FlightComputer.cc
    src/ForceTorque.cc
    src/Frames.cc
    src/Gravity.cc
    src/Kinematics.cc
    src/MassProperties.cc
    src/Model.cc
    src/Parachute.cc
    src/PressureController.cc
    src/Propulsion.cc
    src/Sensors.cc
    src/Simulator.cc
    src/TruthEvent.cc
    src/utility/MatFileReader.cc
    src/Vehicle.cc
)
target_include_directories(rsim_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${mlinterp_SOURCE_DIR}>
)
target_link_libraries(rsim_core
    PUBLIC Eigen3::Eigen matio ${CMAKE_DL_LIBS}
)
set_target_properties(rsim_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(rsim_vespula_flight_computer MODULE
    flight_computers/vespula/VespulaFlightComputer.cc)
target_include_directories(rsim_vespula_flight_computer PRIVATE src)
set_target_properties(rsim_vespula_flight_computer PROPERTIES
    PREFIX "" OUTPUT_NAME "vespula_flight_computer")
install(FILES src/FlightComputerAbi.h DESTINATION include/rsim)

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

    add_executable(rsim_earth_frame_tests tests/test_earth_frames.cc)
    target_link_libraries(rsim_earth_frame_tests PRIVATE rsim_core)
    add_test(NAME rsim_earth_frame_tests COMMAND rsim_earth_frame_tests)

    add_executable(rsim_kinematics_tests tests/test_kinematics.cc)
    target_link_libraries(rsim_kinematics_tests PRIVATE rsim_core)
    add_test(NAME rsim_kinematics_tests COMMAND rsim_kinematics_tests)

    add_executable(rsim_simulator_tests tests/test_simulator.cc)
    target_link_libraries(rsim_simulator_tests PRIVATE rsim_core)
    add_test(NAME rsim_simulator_tests COMMAND rsim_simulator_tests)

    add_executable(rsim_utility_tests tests/test_utility.cc)
    target_link_libraries(rsim_utility_tests PRIVATE rsim_core)
    add_test(NAME rsim_utility_tests COMMAND rsim_utility_tests)

    add_executable(rsim_propulsion_tests tests/test_propulsion.cc)
    target_link_libraries(rsim_propulsion_tests PRIVATE rsim_core)
    add_test(NAME rsim_propulsion_tests COMMAND rsim_propulsion_tests)

    add_executable(rsim_gravity_tests tests/test_gravity.cc)
    target_link_libraries(rsim_gravity_tests PRIVATE rsim_core)
    add_test(NAME rsim_gravity_tests COMMAND rsim_gravity_tests)

    add_executable(rsim_flight_state_tests tests/test_flight_state.cc)
    target_link_libraries(rsim_flight_state_tests PRIVATE rsim_core)
    add_test(NAME rsim_flight_state_tests COMMAND rsim_flight_state_tests)

    add_executable(rsim_native_runtime_tests tests/test_native_runtime.cc)
    target_link_libraries(rsim_native_runtime_tests PRIVATE rsim_core)
    add_test(NAME rsim_native_runtime_tests COMMAND rsim_native_runtime_tests)

    add_library(rsim_c_plugin_smoke MODULE
        tests/c_flight_computer_plugin.c)
    target_include_directories(rsim_c_plugin_smoke PRIVATE src)
endif()
