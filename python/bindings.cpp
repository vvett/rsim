#include "Aerodynamics.h"
#include "DataLogger.h"
#include "Actuators.h"
#include "Environment.h"
#include "FlightState.h"
#include "FlightComputer.h"
#include "ForceTorque.h"
#include "Frames.h"
#include "Gravity.h"
#include "Kinematics.h"
#include "MassProperties.h"
#include "Model.h"
#include "Parachute.h"
#include "PressureController.h"
#include "Propulsion.h"
#include "Simulator.h"
#include "Sensors.h"
#include "TruthEvent.h"
#include "Vehicle.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <limits>
#include <string>

namespace py = pybind11;

namespace {

std::vector<double> copyVector(const py::array_t<double,
                                py::array::forcecast>& array,
                               const char* label) {
    if (array.ndim() != 1) {
        throw py::value_error(std::string(label) + " must be one-dimensional");
    }
    const auto values = array.unchecked<1>();
    std::vector<double> result(static_cast<std::size_t>(values.shape(0)));
    for (py::ssize_t i = 0; i < values.shape(0); ++i) result[i] = values(i);
    return result;
}

std::vector<double> copyFortranField(
    const py::array_t<double, py::array::forcecast>& array,
    py::ssize_t pressure_count, py::ssize_t mixture_ratio_count,
    const char* label) {
    if (array.ndim() != 2 || array.shape(0) != pressure_count ||
        array.shape(1) != mixture_ratio_count) {
        throw py::value_error(std::string(label) +
                              " must have shape (len(Pc), len(MR))");
    }
    const auto values = array.unchecked<2>();
    std::vector<double> result(
        static_cast<std::size_t>(pressure_count * mixture_ratio_count));
    for (py::ssize_t j = 0; j < mixture_ratio_count; ++j) {
        for (py::ssize_t i = 0; i < pressure_count; ++i) {
            result[static_cast<std::size_t>(i + pressure_count * j)] = values(i, j);
        }
    }
    return result;
}

}  // namespace

PYBIND11_MODULE(_rsim, module) {
    module.doc() = "Python bindings for the rsim simulation framework";

    py::class_<rsim::UpdateRate>(module, "UpdateRate")
        .def_static("every_step", &rsim::UpdateRate::everyStep)
        .def_static("times_per_step", &rsim::UpdateRate::timesPerStep,
                    py::arg("count"))
        .def_static("every_n_steps", &rsim::UpdateRate::everyNSteps,
                    py::arg("count"))
        .def_property_readonly("updates_per_step",
                               &rsim::UpdateRate::updatesPerStep)
        .def_property_readonly("steps_per_update",
                               &rsim::UpdateRate::stepsPerUpdate);

    py::enum_<rsim::ModelPhase>(module, "ModelPhase")
        .value("ENVIRONMENT", rsim::ModelPhase::Environment)
        .value("SENSORS", rsim::ModelPhase::Sensors)
        .value("FLIGHT_CONTROL", rsim::ModelPhase::FlightControl)
        .value("ACTUATORS", rsim::ModelPhase::Actuators)
        .value("PLANT", rsim::ModelPhase::Plant)
        .value("MASS_PROPERTIES", rsim::ModelPhase::MassProperties)
        .value("FORCES", rsim::ModelPhase::Forces)
        .value("KINEMATICS", rsim::ModelPhase::Kinematics)
        .value("TRUTH", rsim::ModelPhase::Truth);

    // Derived Python classes expose these common properties through Model.
    py::class_<rsim::Model, std::shared_ptr<rsim::Model>>(module, "Model")
        .def_property_readonly("name", &rsim::Model::name)
        .def_property_readonly("state_namespace", &rsim::Model::stateNamespace)
        .def_property("update_rate", &rsim::Model::updateRate,
                      &rsim::Model::setUpdateRate)
        .def_property("enabled", &rsim::Model::enabled,
                      &rsim::Model::setEnabled)
        .def_property("phase", &rsim::Model::phase, &rsim::Model::setPhase)
        .def_property_readonly("simulation_time",
                               &rsim::Model::simulationTime)
        .def_property_readonly("time_step", &rsim::Model::timeStep)
        .def_property_readonly("global_step", &rsim::Model::globalStep)
        .def_property_readonly("substep", &rsim::Model::substep);

    py::enum_<rsim::InertialStatus>(module, "InertialStatus")
        .value("INERTIAL", rsim::InertialStatus::Inertial)
        .value("NON_INERTIAL", rsim::InertialStatus::NonInertial);

    py::class_<rsim::Frame, std::shared_ptr<rsim::Frame>>(module, "Frame")
        .def(py::init<std::string, rsim::InertialStatus>(),
             py::arg("name"), py::arg("inertial_status"))
        .def_property_readonly("name", &rsim::Frame::name)
        .def_property_readonly("inertial_status",
                               &rsim::Frame::inertialStatus)
        .def_property_readonly("is_inertial", &rsim::Frame::isInertial)
        .def_property_readonly("is_root", &rsim::Frame::isRoot)
        .def_property_readonly("parent", &rsim::Frame::parent);

    py::class_<rsim::AerodynamicCoefficients>(module, "AerodynamicCoefficients")
        .def_readonly("axial", &rsim::AerodynamicCoefficients::axial)
        .def_readonly("side_y", &rsim::AerodynamicCoefficients::side_y)
        .def_readonly("side_z", &rsim::AerodynamicCoefficients::side_z)
        .def_readonly("roll_moment", &rsim::AerodynamicCoefficients::roll_moment)
        .def_readonly("pitch_moment", &rsim::AerodynamicCoefficients::pitch_moment)
        .def_readonly("yaw_moment", &rsim::AerodynamicCoefficients::yaw_moment);

    py::class_<rsim::AerodynamicConditions>(module, "AerodynamicConditions")
        .def(py::init<>())
        .def_readwrite("mach", &rsim::AerodynamicConditions::mach)
        .def_readwrite("alpha_deg", &rsim::AerodynamicConditions::alpha_deg)
        .def_readwrite("phi_aero_deg", &rsim::AerodynamicConditions::phi_aero_deg)
        .def_readwrite("reynolds", &rsim::AerodynamicConditions::reynolds)
        .def_readwrite("dynamic_pressure",
                       &rsim::AerodynamicConditions::dynamic_pressure)
        .def_readwrite(
            "relative_air_velocity_in_body",
            &rsim::AerodynamicConditions::relative_air_velocity_in_body);

    py::class_<rsim::AerodynamicLoad>(module, "AerodynamicLoad")
        .def_readonly("force", &rsim::AerodynamicLoad::force)
        .def_readonly("torque", &rsim::AerodynamicLoad::torque);

    py::class_<rsim::AerodynamicDatabase>(module, "AerodynamicDatabase")
        .def(py::init<const std::string&>(), py::arg("path"))
        .def("evaluate", &rsim::AerodynamicDatabase::evaluate,
             py::arg("mach"), py::arg("alpha_deg"),
             py::arg("phi_aero_deg"), py::arg("reynolds"))
        .def_property_readonly("reference_area",
                               &rsim::AerodynamicDatabase::referenceArea)
        .def_property_readonly("reference_length",
                               &rsim::AerodynamicDatabase::referenceLength)
        .def_property_readonly("moment_reference_x",
                               &rsim::AerodynamicDatabase::momentReferenceX)
        .def_property_readonly(
            "stability_reference_length",
            &rsim::AerodynamicDatabase::stabilityReferenceLength);

    py::class_<rsim::CylindricalPropellantTank>(
        module, "CylindricalPropellantTank")
        .def(py::init<double, double, double, double, double, double>(),
             py::arg("radius"), py::arg("length"),
             py::arg("aft_position"), py::arg("propellant_density"),
             py::arg("initial_propellant_mass"),
             py::arg("mass_flow_rate") = 0.0)
        .def_property_readonly("capacity",
                               &rsim::CylindricalPropellantTank::capacity)
        .def_property_readonly(
            "propellant_mass",
            &rsim::CylindricalPropellantTank::propellantMass)
        .def_property("mass_flow_rate",
                      &rsim::CylindricalPropellantTank::massFlowRate,
                      &rsim::CylindricalPropellantTank::setMassFlowRate)
        .def_property_readonly(
            "center_of_mass",
            &rsim::CylindricalPropellantTank::centerOfMass)
        .def_property_readonly(
            "moment_of_inertia",
            &rsim::CylindricalPropellantTank::inertiaAboutCenterOfMass);

    py::class_<rsim::RigidBody, rsim::Model,
               std::shared_ptr<rsim::RigidBody>>(module, "RigidBody")
        .def(py::init<double, rsim::Matrix3d, rsim::Vector3d, std::string,
                      rsim::UpdateRate>(),
             py::arg("mass"), py::arg("moment_of_inertia"),
             py::arg("center_of_gravity") = rsim::Vector3d::Zero(),
             py::arg("name") = "rigid-body",
             py::arg("update_rate") = rsim::UpdateRate::everyStep())
        .def_property("mass", &rsim::RigidBody::mass,
                      &rsim::RigidBody::setMass)
        .def_property("moment_of_inertia",
                      &rsim::RigidBody::momentOfInertia,
                      &rsim::RigidBody::setMomentOfInertia)
        .def_property("center_of_gravity", &rsim::RigidBody::cg,
                      &rsim::RigidBody::setCg)
        .def("set_dry_mass_properties",
             &rsim::RigidBody::setDryMassProperties,
             py::arg("mass"), py::arg("moment_of_inertia"),
             py::arg("center_of_gravity"))
        .def("add_propellant_tank", &rsim::RigidBody::addPropellantTank,
             py::arg("tank"))
        .def("propellant_tank", [](rsim::RigidBody& body, std::size_t index)
             -> rsim::CylindricalPropellantTank& {
                 return body.propellantTank(index);
             }, py::arg("index"), py::return_value_policy::reference_internal)
        .def("set_propellant_mass", &rsim::RigidBody::setPropellantMass,
             py::arg("index"), py::arg("mass"))
        .def("set_tank_mass_flow_rate",
             &rsim::RigidBody::setTankMassFlowRate,
             py::arg("index"), py::arg("mass_flow_rate"))
        .def_property_readonly("propellant_tank_count",
                               &rsim::RigidBody::propellantTankCount)
        .def("add_variable_mass_source", &rsim::RigidBody::addVariableMassSource,
             py::arg("source"), py::keep_alive<1, 2>())
        .def("remove_variable_mass_source",
             &rsim::RigidBody::removeVariableMassSource, py::arg("source"))
        .def_property_readonly("variable_mass_source_count",
                               &rsim::RigidBody::variableMassSourceCount)
        .def("print_mass", &rsim::RigidBody::printMass)
        .def("__repr__", [](const rsim::RigidBody& body) {
            return "rsim.RigidBody(mass=" + std::to_string(body.mass()) + ")";
        });

    // No trampoline is provided: force-model update() calls stay entirely in C++.
    py::class_<rsim::ForceTorque, rsim::Model,
               std::shared_ptr<rsim::ForceTorque>>(module, "ForceTorque")
        .def_property_readonly("force", &rsim::ForceTorque::force)
        .def("torque_about", &rsim::ForceTorque::torqueAbout,
             py::arg("reference_point"))
        .def("force_torque_about", &rsim::ForceTorque::forceTorqueAbout,
             py::arg("reference_point"))
        .def_property_readonly("application_point",
                               &rsim::ForceTorque::applicationPoint);

    py::class_<rsim::ConstantForceTorque, rsim::ForceTorque,
               std::shared_ptr<rsim::ConstantForceTorque>>(
                   module, "ConstantForceTorque")
        .def(py::init<std::string>(), py::arg("name") = "constant-load")
        .def("set_load", &rsim::ConstantForceTorque::setLoad,
             py::arg("force"), py::arg("torque"),
             py::arg("application_point") = rsim::Vector3d::Zero());

    auto gravity_class =
        py::class_<rsim::Gravity, rsim::ForceTorque,
                   std::shared_ptr<rsim::Gravity>>(module, "Gravity")
            .def(py::init<std::string, double, rsim::Vector3d>(),
                 py::arg("name") = "gravity",
                 py::arg("gravitational_parameter") =
                     rsim::Gravity::EarthGravitationalParameter,
                 py::arg("center_in_integration_frame") =
                     rsim::Vector3d::Zero())
            .def_static(
                "earth",
                [](std::string name,
                   rsim::Vector3d center_in_integration_frame) {
                    return std::make_shared<rsim::Gravity>(
                        std::move(name),
                        rsim::Gravity::EarthGravitationalParameter,
                        std::move(center_in_integration_frame));
                },
                py::arg("name") = "earth-gravity",
                py::arg("center_in_integration_frame") =
                    rsim::Vector3d::Zero())
            .def_property_readonly(
                "gravitational_parameter",
                &rsim::Gravity::gravitationalParameter)
            .def_property_readonly(
                "center_in_integration_frame",
                &rsim::Gravity::centerInIntegrationFrame);
    gravity_class.attr("EARTH_GRAVITATIONAL_PARAMETER") =
        rsim::Gravity::EarthGravitationalParameter;

    py::class_<rsim::ForceTorqueRegistry, rsim::Model,
               std::shared_ptr<rsim::ForceTorqueRegistry>>(
                   module, "ForceTorqueRegistry")
        .def(py::init<std::string, rsim::UpdateRate>(),
             py::arg("name") = "forces",
             py::arg("update_rate") = rsim::UpdateRate::everyStep())
        .def_property_readonly("force", &rsim::ForceTorqueRegistry::getForce)
        .def("torque_about", &rsim::ForceTorqueRegistry::getTorque,
             py::arg("reference_point"))
        .def("force_torque_about", &rsim::ForceTorqueRegistry::getForceTorque,
             py::arg("reference_point"))
        .def("__len__", &rsim::ForceTorqueRegistry::size)
        .def("__getitem__", [](rsim::ForceTorqueRegistry& registry,
                                std::size_t index) -> rsim::ForceTorque& {
                 return registry.at(index);
             }, py::arg("index"), py::return_value_policy::reference_internal);

    py::class_<rsim::Aerodynamics, rsim::ForceTorque,
               std::shared_ptr<rsim::Aerodynamics>>(module, "Aerodynamics")
        .def(py::init<std::string, rsim::AerodynamicDatabase>(),
             py::arg("name"), py::arg("database"))
        .def("set_conditions", &rsim::Aerodynamics::setConditions,
             py::arg("conditions"))
        .def_property_readonly("conditions", &rsim::Aerodynamics::conditions)
        .def_property_readonly("coefficients", &rsim::Aerodynamics::coefficients)
        .def_property_readonly("load", &rsim::Aerodynamics::load)
        .def_property_readonly("center_pressure_x",
                               &rsim::Aerodynamics::centerPressureX)
        .def_property_readonly("static_stability_margin",
                               &rsim::Aerodynamics::staticStabilityMargin)
        .def("calculate_load", &rsim::Aerodynamics::calculateLoad,
             py::arg("conditions"));

    py::enum_<rsim::FluidPhase>(module, "FluidPhase")
        .value("LIQUID", rsim::FluidPhase::Liquid)
        .value("GAS", rsim::FluidPhase::Gas);
    py::enum_<rsim::InletRole>(module, "InletRole")
        .value("NONE", rsim::InletRole::None)
        .value("OXIDIZER", rsim::InletRole::Oxidizer)
        .value("FUEL", rsim::InletRole::Fuel);
    py::enum_<rsim::BlowdownMode>(module, "BlowdownMode")
        .value("ADIABATIC", rsim::BlowdownMode::Adiabatic)
        .value("ISENTROPIC", rsim::BlowdownMode::Isentropic)
        .value("POLYTROPIC", rsim::BlowdownMode::Polytropic)
        .value("ENERGY_CONSERVING", rsim::BlowdownMode::EnergyConserving);

    py::class_<rsim::VariableMassSource,
               std::shared_ptr<rsim::VariableMassSource>>(
                   module, "VariableMassSource");
    py::class_<rsim::PropulsionNode, std::shared_ptr<rsim::PropulsionNode>>(
        module, "PropulsionNode")
        .def_property_readonly("name", &rsim::PropulsionNode::name)
        .def_property_readonly("pressure", &rsim::PropulsionNode::pressure)
        .def_property_readonly("temperature", &rsim::PropulsionNode::temperature)
        .def("density", &rsim::PropulsionNode::density, py::arg("phase"))
        .def("phase_available", &rsim::PropulsionNode::phaseAvailable,
             py::arg("phase"));

    py::class_<rsim::Ambient, rsim::PropulsionNode,
               std::shared_ptr<rsim::Ambient>>(module, "Ambient")
        .def(py::init<std::string, double, double>(), py::arg("name") = "ambient",
             py::arg("pressure") = 101325.0, py::arg("temperature") = 288.15)
        .def_property("pressure", &rsim::Ambient::pressure,
                      &rsim::Ambient::setPressure)
        .def_property("temperature", &rsim::Ambient::temperature,
                      &rsim::Ambient::setTemperature);

    py::class_<rsim::Tank, rsim::PropulsionNode, rsim::VariableMassSource,
               std::shared_ptr<rsim::Tank>>(module, "Tank")
        .def(py::init<std::string, double, double, double, double, double,
                      double, double, double, double, rsim::BlowdownMode,
                      double, std::string, std::string>(),
             py::arg("name"), py::arg("radius"), py::arg("length"),
             py::arg("aft_position"), py::arg("liquid_density"),
             py::arg("propellant_mass"), py::arg("pressure"),
             py::arg("temperature") = 293.15,
             py::arg("pressurant_mass") = 0.0,
             py::arg("collapse_factor") = 1.0,
             py::arg("blowdown_mode") = rsim::BlowdownMode::Adiabatic,
             py::arg("polytropic_exponent") = 1.2,
             py::arg("propellant") = "propellant",
             py::arg("pressurant") = "Helium")
        .def_property("propellant_mass", &rsim::Tank::propellantMass,
                      &rsim::Tank::setPropellantMass)
        .def_property_readonly("pressurant_mass", &rsim::Tank::pressurantMass)
        .def_property_readonly("volume", &rsim::Tank::volume)
        .def_property_readonly("liquid_volume", &rsim::Tank::liquidVolume)
        .def_property_readonly("ullage_volume", &rsim::Tank::ullageVolume)
        .def_property_readonly("propellant_internal_energy",
                               &rsim::Tank::propellantInternalEnergy)
        .def_property_readonly("pressurant_internal_energy",
                               &rsim::Tank::pressurantInternalEnergy)
        .def_property_readonly("collapse_factor", &rsim::Tank::collapseFactor)
        .def_property_readonly("blowdown_mode", &rsim::Tank::blowdownMode)
        .def_property_readonly("propellant", &rsim::Tank::propellant)
        .def_property_readonly("pressurant", &rsim::Tank::pressurant)
        .def_property_readonly("center_of_mass",
                               &rsim::Tank::variableCenterOfMass)
        .def_property_readonly("moment_of_inertia",
                               &rsim::Tank::variableInertiaAboutCenterOfMass);

    py::class_<rsim::Junction, rsim::PropulsionNode,
               std::shared_ptr<rsim::Junction>>(module, "Junction")
        .def(py::init<std::string, double, double>(), py::arg("name"),
             py::arg("pressure") = 101325.0, py::arg("temperature") = 293.15);

    py::class_<rsim::CombustionTableData,
               std::shared_ptr<rsim::CombustionTableData>>(
                   module, "CombustionTableData")
        .def(py::init([](
                 const py::array_t<double, py::array::forcecast>& pc,
                 const py::array_t<double, py::array::forcecast>& mr,
                 const py::array_t<double, py::array::forcecast>& temperature,
                 const py::array_t<double, py::array::forcecast>& cstar,
                 const py::array_t<double, py::array::forcecast>& cf_vac,
                 const py::array_t<double, py::array::forcecast>& isp_vac,
                 double expansion_ratio, std::string oxidizer,
                 std::string fuel, bool frozen_at_throat) {
            if (pc.ndim() != 1 || mr.ndim() != 1) {
                throw py::value_error("Pc and MR axes must be one-dimensional");
            }
            return std::make_shared<rsim::CombustionTableData>(
                copyVector(pc, "Pc"), copyVector(mr, "MR"),
                copyFortranField(temperature, pc.shape(0), mr.shape(0),
                                 "chamber_temperature"),
                copyFortranField(cstar, pc.shape(0), mr.shape(0), "cstar"),
                copyFortranField(cf_vac, pc.shape(0), mr.shape(0), "cf_vac"),
                copyFortranField(isp_vac, pc.shape(0), mr.shape(0), "isp_vac"),
                expansion_ratio, std::move(oxidizer), std::move(fuel),
                frozen_at_throat);
        }), py::arg("chamber_pressure"), py::arg("mixture_ratio"),
             py::arg("chamber_temperature"), py::arg("cstar"),
             py::arg("cf_vac"), py::arg("isp_vac"),
             py::arg("expansion_ratio"), py::arg("oxidizer") = "",
             py::arg("fuel") = "", py::arg("frozen_at_throat") = false)
        .def("chamber_temperature", &rsim::CombustionTableData::chamberTemperature)
        .def("cstar", &rsim::CombustionTableData::cstar)
        .def("cf_vac", &rsim::CombustionTableData::cfVac)
        .def("isp_vac", &rsim::CombustionTableData::ispVac)
        .def_property_readonly("expansion_ratio",
                               &rsim::CombustionTableData::expansionRatio)
        .def_property_readonly("oxidizer", &rsim::CombustionTableData::oxidizer)
        .def_property_readonly("fuel", &rsim::CombustionTableData::fuel)
        .def_property_readonly("frozen_at_throat",
                               &rsim::CombustionTableData::frozenAtThroat);

    py::class_<rsim::Combustor, rsim::PropulsionNode,
               std::shared_ptr<rsim::Combustor>>(module, "Combustor")
        .def(py::init<std::string, std::shared_ptr<rsim::CombustionTableData>,
                      double, double>(), py::arg("name"), py::arg("table"),
             py::arg("initial_pressure") = 101325.0,
             py::arg("initial_mixture_ratio") = 1.0)
        .def_property_readonly("mixture_ratio", &rsim::Combustor::mixtureRatio)
        .def_property_readonly("oxidizer_mass_flow",
                               &rsim::Combustor::oxidizerMassFlow)
        .def_property_readonly("fuel_mass_flow", &rsim::Combustor::fuelMassFlow)
        .def_property_readonly("total_mass_flow", &rsim::Combustor::totalMassFlow);

    py::class_<rsim::Engine, rsim::Combustor, rsim::ForceTorque,
               std::shared_ptr<rsim::Engine>>(module, "Engine",
                                               py::multiple_inheritance())
        .def(py::init<std::string, std::shared_ptr<rsim::CombustionTableData>,
                      double, double, rsim::Vector3d, rsim::Vector3d,
                      std::shared_ptr<rsim::Ambient>, double, double>(),
             py::arg("name"), py::arg("table"), py::arg("throat_area"),
             py::arg("expansion_ratio"),
             py::arg("thrust_direction") = rsim::Vector3d::UnitX(),
             py::arg("application_point") = rsim::Vector3d::Zero(),
             py::arg("ambient") = nullptr,
             py::arg("initial_pressure") = 101325.0,
             py::arg("initial_mixture_ratio") = 1.0)
        .def_property_readonly("thrust", &rsim::Engine::thrust)
        .def_property_readonly("cstar", &rsim::Engine::cstar)
        .def_property_readonly("thrust_coefficient",
                               &rsim::Engine::thrustCoefficient)
        .def_property_readonly("ambient_pressure", &rsim::Engine::ambientPressure);

    py::class_<rsim::FlowComponent, std::shared_ptr<rsim::FlowComponent>>(
        module, "FlowComponent")
        .def("mass_flow", &rsim::FlowComponent::massFlow,
             py::arg("upstream"), py::arg("downstream"), py::arg("phase"));
    py::class_<rsim::Orifice, rsim::FlowComponent,
               std::shared_ptr<rsim::Orifice>>(module, "Orifice")
        .def(py::init<double>(), py::arg("cd_area"))
        .def(py::init<double, double>(), py::arg("area"),
             py::arg("discharge_coefficient"))
        .def_property_readonly("cd_area", &rsim::Orifice::cdArea);
    py::class_<rsim::Valve, rsim::Orifice,
               std::shared_ptr<rsim::Valve>>(module, "Valve")
        .def(py::init<double, double>(), py::arg("cd_area"),
             py::arg("opening") = 1.0)
        .def(py::init<double, double, double>(), py::arg("area"),
             py::arg("discharge_coefficient"), py::arg("opening"))
        .def_property("opening", &rsim::Valve::opening, &rsim::Valve::setOpening);
    py::class_<rsim::Line, rsim::FlowComponent,
               std::shared_ptr<rsim::Line>>(module, "Line")
        .def(py::init<double, double, double, double>(), py::arg("length"),
             py::arg("diameter"), py::arg("roughness") = 0.0,
             py::arg("dynamic_viscosity") = 1.0e-3);
    py::class_<rsim::CompoundEdge, rsim::FlowComponent,
               std::shared_ptr<rsim::CompoundEdge>>(module, "CompoundEdge")
        .def(py::init<std::vector<std::shared_ptr<rsim::FlowComponent>>>(),
             py::arg("parts"))
        .def_property_readonly("intermediate_pressures",
                               &rsim::CompoundEdge::intermediatePressures);
    py::class_<rsim::Pump, rsim::FlowComponent,
               std::shared_ptr<rsim::Pump>>(module, "Pump").def(py::init<>());
    py::class_<rsim::Turbine, rsim::FlowComponent,
               std::shared_ptr<rsim::Turbine>>(module, "Turbine").def(py::init<>());

    py::class_<rsim::ThermalContact>(module, "ThermalContact")
        .def(py::init<>())
        .def_readwrite("node", &rsim::ThermalContact::node)
        .def_readwrite("phase", &rsim::ThermalContact::phase)
        .def_readwrite("environment_temperature",
                       &rsim::ThermalContact::environment_temperature)
        .def_readwrite("heat_flux", &rsim::ThermalContact::heat_flux)
        .def_readwrite("area", &rsim::ThermalContact::area);

    py::class_<rsim::Propulsion, rsim::Model,
               std::shared_ptr<rsim::Propulsion>>(module, "Propulsion")
        .def(py::init<std::string, rsim::UpdateRate>(),
             py::arg("name") = "propulsion",
             py::arg("update_rate") = rsim::UpdateRate::everyStep())
        .def("connect", py::overload_cast<
                 std::shared_ptr<rsim::PropulsionNode>,
                 std::shared_ptr<rsim::PropulsionNode>,
                 std::shared_ptr<rsim::FlowComponent>, rsim::FluidPhase,
                 rsim::InletRole>(&rsim::Propulsion::connect),
             py::arg("upstream"), py::arg("downstream"), py::arg("edge"),
             py::arg("upstream_phase") = rsim::FluidPhase::Liquid,
             py::arg("downstream_role") = rsim::InletRole::None)
        .def("connect", py::overload_cast<
                 std::shared_ptr<rsim::PropulsionNode>,
                 std::shared_ptr<rsim::PropulsionNode>,
                 std::vector<std::shared_ptr<rsim::FlowComponent>>,
                 rsim::FluidPhase, rsim::InletRole>(&rsim::Propulsion::connect),
             py::arg("upstream"), py::arg("downstream"), py::arg("edges"),
             py::arg("upstream_phase") = rsim::FluidPhase::Liquid,
             py::arg("downstream_role") = rsim::InletRole::None)
        .def("vent", &rsim::Propulsion::vent, py::arg("upstream"),
             py::arg("ambient"), py::arg("edge"),
             py::arg("upstream_phase") = rsim::FluidPhase::Gas)
        .def("add_thermal_contact", &rsim::Propulsion::addThermalContact)
        .def("finalize", &rsim::Propulsion::finalize)
        .def_property_readonly("finalized", &rsim::Propulsion::finalized)
        .def_property_readonly("connection_count",
                               &rsim::Propulsion::connectionCount)
        .def("connection_mass_flow", &rsim::Propulsion::connectionMassFlow)
        .def_property_readonly("total_stored_mass",
                               &rsim::Propulsion::totalStoredMass)
        .def_property_readonly("last_solve_summary",
                               &rsim::Propulsion::lastSolveSummary);

    py::class_<rsim::BangBangPressureController, rsim::Model,
               std::shared_ptr<rsim::BangBangPressureController>>(
                   module, "BangBangPressureController")
        .def(py::init<
                 std::string, std::shared_ptr<rsim::Tank>,
                 std::shared_ptr<rsim::Valve>, double, double, bool,
                 rsim::UpdateRate>(),
             py::arg("name"), py::arg("controlled_tank"),
             py::arg("branch_valve"), py::arg("pressure_setpoint"),
             py::arg("deadband_half_width"),
             py::arg("initially_open") = false,
             py::arg("update_rate") = rsim::UpdateRate::everyStep())
        .def_property_readonly(
            "controlled_tank",
            &rsim::BangBangPressureController::controlledTank)
        .def_property_readonly(
            "branch_valve",
            &rsim::BangBangPressureController::branchValve)
        .def_property_readonly(
            "pressure_setpoint",
            &rsim::BangBangPressureController::pressureSetpoint)
        .def_property_readonly(
            "deadband_half_width",
            &rsim::BangBangPressureController::deadbandHalfWidth)
        .def_property_readonly(
            "lower_threshold",
            &rsim::BangBangPressureController::lowerThreshold)
        .def_property_readonly(
            "upper_threshold",
            &rsim::BangBangPressureController::upperThreshold)
        .def_property_readonly(
            "valve_open",
            &rsim::BangBangPressureController::valveOpen);

    py::class_<rsim::Kinematics, rsim::Model,
               std::shared_ptr<rsim::Kinematics>>(module, "Kinematics")
        .def(py::init<std::string, std::shared_ptr<rsim::Frame>,
                      rsim::ForceTorqueRegistry&, rsim::RigidBody&,
                      rsim::UpdateRate>(),
             py::arg("body_frame_name"), py::arg("parent_frame"),
             py::arg("force_torque"), py::arg("rigid_body"),
             py::arg("update_rate") = rsim::UpdateRate::everyStep(),
             // keep_alive retains each non-owning constructor dependency.
             py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def(py::init([](std::string body_frame_name,
                         rsim::FrameGraph& frame_graph,
                         const std::string& integration_frame_name,
                         rsim::ForceTorqueRegistry& force_torque,
                         rsim::RigidBody& rigid_body,
                         rsim::UpdateRate update_rate) {
                 const std::shared_ptr<rsim::Frame> integration_frame =
                     frame_graph.frame(integration_frame_name);
                 if (!integration_frame->isInertial() &&
                     PyErr_WarnEx(PyExc_RuntimeWarning,
                                  "selected integration frame is non-inertial; "
                                  "transport accelerations are not yet included",
                                  1) < 0) {
                     throw py::error_already_set();
                 }
                 return std::make_unique<rsim::Kinematics>(
                     std::move(body_frame_name), integration_frame,
                     force_torque, rigid_body, update_rate);
             }),
             py::arg("body_frame_name"), py::arg("frame_graph"),
             py::arg("integration_frame_name"), py::arg("force_torque"),
             py::arg("rigid_body"),
             py::arg("update_rate") = rsim::UpdateRate::everyStep(),
             py::keep_alive<1, 3>(), py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>())
        .def("set_state", [](rsim::Kinematics& kinematics,
                             rsim::Vector3d position_in_parent,
                             rsim::Vector3d velocity_in_parent,
                             const rsim::Matrix3d& parent_from_body,
                             rsim::Vector3d angular_velocity_in_body) {
                 // Python uses a 3-by-3 rotation matrix instead of Eigen's quaternion type.
                 kinematics.setState(
                     std::move(position_in_parent),
                     std::move(velocity_in_parent),
                     Eigen::Quaterniond(parent_from_body),
                     std::move(angular_velocity_in_body));
             },
             py::arg("position_in_parent"),
             py::arg("velocity_in_parent"),
             py::arg("parent_from_body"),
             py::arg("angular_velocity_in_body"))
        .def_property_readonly("position_in_parent",
                               &rsim::Kinematics::positionInParent)
        .def_property_readonly("velocity_in_parent",
                               &rsim::Kinematics::velocityInParent)
        .def_property_readonly("parent_from_body", [](const rsim::Kinematics& kinematics) {
            return kinematics.rotationIntoParent().toRotationMatrix();
        })
        .def_property_readonly("angular_velocity_in_body",
                               &rsim::Kinematics::angularVelocityInBody)
        .def_property_readonly("integration_frame",
                               &rsim::Kinematics::integrationFrame)
        .def_property_readonly("body_frame",
                               [](const rsim::Kinematics& kinematics) {
                                   return kinematics.bodyFrame;
                               });

    py::class_<rsim::FrameGraph, rsim::Model,
               std::shared_ptr<rsim::FrameGraph>>(module, "FrameGraph")
        .def(py::init<std::string, rsim::UpdateRate>(),
             py::arg("name") = "frames",
             py::arg("update_rate") = rsim::UpdateRate::everyStep())
        .def("add_frame", &rsim::FrameGraph::addFrame, py::arg("frame"))
        .def("frame", &rsim::FrameGraph::frame, py::arg("name"));

    py::enum_<rsim::VehicleState>(module, "VehicleState")
        .value("ON_RAIL", rsim::VehicleState::OnRail)
        .value("THRUSTING", rsim::VehicleState::Thrusting)
        .value("COASTING", rsim::VehicleState::Coasting)
        .value("DESCENDING", rsim::VehicleState::Descending)
        .value("DROGUE_DEPLOYED", rsim::VehicleState::DrogueDeployed)
        .value("MAIN_DEPLOYED", rsim::VehicleState::MainDeployed)
        .value("LANDED", rsim::VehicleState::Landed);

    py::class_<rsim::DeploymentTrigger>(module, "DeploymentTrigger")
        .def(py::init<>())
        .def_readwrite("maximum_altitude",
                       &rsim::DeploymentTrigger::maximum_altitude)
        .def_readwrite("maximum_vertical_velocity",
                       &rsim::DeploymentTrigger::maximum_vertical_velocity)
        .def_readwrite("earliest_simulation_time",
                       &rsim::DeploymentTrigger::earliest_simulation_time)
        .def_readwrite("deploy_at_apogee",
                       &rsim::DeploymentTrigger::deploy_at_apogee)
        .def_readwrite(
            "require_descending_altitude_crossing",
            &rsim::DeploymentTrigger::require_descending_altitude_crossing);

    py::class_<rsim::FlightStateConfiguration>(
        module, "FlightStateConfiguration")
        .def(py::init<>())
        .def_readwrite(
            "origin_in_integration_frame",
            &rsim::FlightStateConfiguration::origin_in_integration_frame)
        .def_readwrite(
            "up_direction_in_integration_frame",
            &rsim::FlightStateConfiguration::up_direction_in_integration_frame)
        .def_readwrite(
            "rail_direction_in_integration_frame",
            &rsim::FlightStateConfiguration::rail_direction_in_integration_frame)
        .def_readwrite("rail_length",
                       &rsim::FlightStateConfiguration::rail_length)
        .def_readwrite("thrust_on_threshold",
                       &rsim::FlightStateConfiguration::thrust_on_threshold)
        .def_readwrite("thrust_off_threshold",
                       &rsim::FlightStateConfiguration::thrust_off_threshold)
        .def_readwrite(
            "descending_velocity_threshold",
            &rsim::FlightStateConfiguration::descending_velocity_threshold)
        .def_readwrite(
            "apogee_velocity_deadband",
            &rsim::FlightStateConfiguration::apogee_velocity_deadband)
        .def_readwrite("drogue_stage_enabled",
                       &rsim::FlightStateConfiguration::drogue_stage_enabled)
        .def_readwrite("drogue_trigger",
                       &rsim::FlightStateConfiguration::drogue_trigger)
        .def_readwrite("main_stage_enabled",
                       &rsim::FlightStateConfiguration::main_stage_enabled)
        .def_readwrite("main_trigger",
                       &rsim::FlightStateConfiguration::main_trigger)
        .def_readwrite("ground_altitude",
                       &rsim::FlightStateConfiguration::ground_altitude)
        .def_readwrite(
            "maximum_landing_vertical_velocity",
            &rsim::FlightStateConfiguration::maximum_landing_vertical_velocity)
        .def_readwrite(
            "legacy_automatic_recovery",
            &rsim::FlightStateConfiguration::legacy_automatic_recovery);

    py::class_<rsim::TruthState, rsim::Model,
               std::shared_ptr<rsim::TruthState>>(
                   module, "TruthState")
        .def(py::init<
                 std::string, rsim::FlightStateConfiguration,
                 std::vector<std::shared_ptr<rsim::ForceTorque>>,
                 rsim::UpdateRate>(),
             py::arg("name") = "flight-state",
             py::arg("configuration") = rsim::FlightStateConfiguration{},
             py::arg("thrust_sources") =
                 std::vector<std::shared_ptr<rsim::ForceTorque>>{},
             py::arg("update_rate") = rsim::UpdateRate::everyStep())
        .def("add_thrust_source", &rsim::TruthState::addThrustSource,
             py::arg("thrust_source"))
        .def("set_recovery_parachutes", &rsim::TruthState::setRecoveryParachutes,
             py::arg("drogue"), py::arg("main"), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def("enable_console_output", &rsim::TruthState::enableConsoleOutput,
             py::arg("flight_conditions"), py::arg("aerodynamics"),
             py::arg("period") = 1.0, py::keep_alive<1, 2>(),
             py::keep_alive<1, 3>())
        .def_property_readonly("configuration",
                               &rsim::TruthState::configuration)
        .def_property_readonly("altitude", &rsim::TruthState::altitude)
        .def_property_readonly("vertical_velocity",
                               &rsim::TruthState::verticalVelocity)
        .def_property_readonly("rail_distance",
                               &rsim::TruthState::railDistance)
        .def_property_readonly("thrust_magnitude",
                               &rsim::TruthState::thrustMagnitude)
        .def_property_readonly("apogee_crossed",
                               &rsim::TruthState::apogeeCrossed);
    module.attr("FlightStateModel") = module.attr("TruthState");

    py::class_<rsim::Parachute, rsim::ForceTorque,
               std::shared_ptr<rsim::Parachute>>(module, "Parachute")
        .def(py::init<std::string, double, double, rsim::VehicleState,
                      double, rsim::Vector3d, rsim::Vector3d>(),
             py::arg("name"), py::arg("drag_coefficient"),
             py::arg("reference_area"), py::arg("deployment_state"),
             py::arg("air_density") = 1.225,
             py::arg("wind_velocity_in_integration_frame") =
                 rsim::Vector3d::Zero(),
             py::arg("application_point_in_body") = rsim::Vector3d::Zero())
        .def_static(
            "from_diameter", &rsim::Parachute::fromDiameter,
            py::arg("name"), py::arg("drag_coefficient"),
            py::arg("diameter"), py::arg("deployment_state"),
            py::arg("air_density") = 1.225,
            py::arg("wind_velocity_in_integration_frame") =
                rsim::Vector3d::Zero(),
            py::arg("application_point_in_body") = rsim::Vector3d::Zero())
        .def_static(
            "from_terminal_descent_speed",
            &rsim::Parachute::fromTerminalDescentSpeed,
            py::arg("name"), py::arg("drag_coefficient"),
            py::arg("terminal_descent_speed"), py::arg("design_mass"),
            py::arg("design_air_density"), py::arg("deployment_state"),
            py::arg("wind_velocity_in_integration_frame") =
                rsim::Vector3d::Zero(),
            py::arg("application_point_in_body") = rsim::Vector3d::Zero())
        .def("set_atmosphere", &rsim::Parachute::setAtmosphere,
             py::arg("air_density"),
             py::arg("wind_velocity_in_integration_frame") =
                 rsim::Vector3d::Zero())
        .def("deploy", &rsim::Parachute::deploy)
        .def_property_readonly("drag_coefficient",
                               &rsim::Parachute::dragCoefficient)
        .def_property_readonly("reference_area",
                               &rsim::Parachute::referenceArea)
        .def_property_readonly("equivalent_diameter",
                               &rsim::Parachute::equivalentDiameter)
        .def_property_readonly("deployment_state",
                               &rsim::Parachute::deploymentState)
        .def_property_readonly("air_density", &rsim::Parachute::airDensity)
        .def_property_readonly(
            "wind_velocity_in_integration_frame",
            &rsim::Parachute::windVelocityInIntegrationFrame)
        .def_property_readonly("deployed", &rsim::Parachute::deployed);

    py::class_<rsim::Vehicle, rsim::Model,
               std::shared_ptr<rsim::Vehicle>>(module, "Vehicle")
        .def(py::init<std::string, std::shared_ptr<rsim::Frame>, double,
                      rsim::Matrix3d, rsim::Vector3d, rsim::UpdateRate,
                      rsim::UpdateRate, rsim::UpdateRate, rsim::UpdateRate,
                      std::vector<std::shared_ptr<rsim::Model>>,
                      std::vector<std::shared_ptr<rsim::ForceTorque>>>(),
             py::arg("name"), py::arg("parent_frame"), py::arg("mass"),
             py::arg("moment_of_inertia"),
             py::arg("center_of_gravity") = rsim::Vector3d::Zero(),
             py::arg("vehicle_rate") = rsim::UpdateRate::everyStep(),
             py::arg("rigid_body_rate") = rsim::UpdateRate::everyStep(),
             py::arg("force_rate") = rsim::UpdateRate::everyStep(),
             py::arg("kinematics_rate") = rsim::UpdateRate::everyStep(),
             py::arg("models") =
                 std::vector<std::shared_ptr<rsim::Model>>{},
             py::arg("force_torques") =
                 std::vector<std::shared_ptr<rsim::ForceTorque>>{})
        .def(py::init([](std::string name,
                         rsim::FrameGraph& frame_graph,
                         const std::string& integration_frame_name,
                         double mass,
                         rsim::Matrix3d moment_of_inertia,
                         rsim::Vector3d center_of_gravity,
                         rsim::UpdateRate vehicle_rate,
                         rsim::UpdateRate rigid_body_rate,
                         rsim::UpdateRate force_rate,
                         rsim::UpdateRate kinematics_rate,
                         std::vector<std::shared_ptr<rsim::Model>> models,
                         std::vector<std::shared_ptr<rsim::ForceTorque>>
                             force_torques) {
                 const std::shared_ptr<rsim::Frame> integration_frame =
                     frame_graph.frame(integration_frame_name);
                 if (!integration_frame->isInertial() &&
                     PyErr_WarnEx(PyExc_RuntimeWarning,
                                  "selected integration frame is non-inertial; "
                                  "transport accelerations are not yet included",
                                  1) < 0) {
                     throw py::error_already_set();
                 }
                 return std::make_unique<rsim::Vehicle>(
                     std::move(name), integration_frame, mass,
                     std::move(moment_of_inertia),
                     std::move(center_of_gravity), vehicle_rate,
                     rigid_body_rate, force_rate, kinematics_rate,
                     std::move(models),
                     std::move(force_torques));
             }),
             py::arg("name"), py::arg("frame_graph"),
             py::arg("integration_frame_name"), py::arg("mass"),
             py::arg("moment_of_inertia"),
             py::arg("center_of_gravity") = rsim::Vector3d::Zero(),
             py::arg("vehicle_rate") = rsim::UpdateRate::everyStep(),
             py::arg("rigid_body_rate") = rsim::UpdateRate::everyStep(),
             py::arg("force_rate") = rsim::UpdateRate::everyStep(),
             py::arg("kinematics_rate") = rsim::UpdateRate::everyStep(),
             py::arg("models") =
                 std::vector<std::shared_ptr<rsim::Model>>{},
             py::arg("force_torques") =
                 std::vector<std::shared_ptr<rsim::ForceTorque>>{},
             py::keep_alive<1, 3>())
        .def_property("state", &rsim::Vehicle::state,
                      &rsim::Vehicle::setState)
        .def("add_force_torque", &rsim::Vehicle::addSharedForceTorque,
             py::arg("force_torque"))
        .def("remove_force_torque", &rsim::Vehicle::removeForceTorque,
             py::arg("force_torque"))
        .def("add_model", &rsim::Vehicle::addModel, py::arg("model"))
        .def("remove_model", &rsim::Vehicle::removeModel, py::arg("model"))
        // reference_internal ties each borrowed child to its owning Vehicle.
        .def_property_readonly("rigid_body", [](rsim::Vehicle& vehicle)
                               -> rsim::RigidBody& {
                                   return vehicle.rigidBody();
                               },
                               py::return_value_policy::reference_internal)
        .def_property_readonly("forces", [](rsim::Vehicle& vehicle)
                               -> rsim::ForceTorqueRegistry& {
                                   return vehicle.forces();
                               },
                               py::return_value_policy::reference_internal)
        .def_property_readonly("kinematics", [](rsim::Vehicle& vehicle)
                               -> rsim::Kinematics& {
                                   return vehicle.kinematics();
                               },
                               py::return_value_policy::reference_internal)
        .def_property_readonly("body_frame", &rsim::Vehicle::bodyFrame);

    py::class_<rsim::ScalarSensor, rsim::Model,
               std::shared_ptr<rsim::ScalarSensor>>(module, "ScalarSensor")
        .def(py::init<std::string, const rsim::Model&, std::string, rsim::UpdateRate>(),
             py::arg("name"), py::arg("source"), py::arg("alias"),
             py::arg("update_rate") = rsim::UpdateRate::everyStep(),
             py::keep_alive<1, 3>())
        .def_property_readonly("value", &rsim::ScalarSensor::value)
        .def_property_readonly("valid", &rsim::ScalarSensor::valid)
        .def_property_readonly("sample_time", &rsim::ScalarSensor::sampleTime)
        .def_property_readonly("sequence", &rsim::ScalarSensor::sequence);

    py::class_<rsim::SpecificForceSensor, rsim::Model,
               std::shared_ptr<rsim::SpecificForceSensor>>(
                   module, "SpecificForceSensor")
        .def(py::init<std::string, rsim::Vehicle&, rsim::UpdateRate>(),
             py::arg("name"), py::arg("vehicle"),
             py::arg("update_rate") = rsim::UpdateRate::everyStep(),
             py::keep_alive<1, 3>())
        .def_property_readonly(
            "specific_force_in_body",
            &rsim::SpecificForceSensor::specificForceInBody)
        .def_property_readonly("valid", &rsim::SpecificForceSensor::valid)
        .def_property_readonly(
            "sample_time", &rsim::SpecificForceSensor::sampleTime)
        .def_property_readonly("sequence", &rsim::SpecificForceSensor::sequence);

    py::class_<rsim::IdealActuator, rsim::Model,
               std::shared_ptr<rsim::IdealActuator>>(module, "IdealActuator")
        .def("command", &rsim::IdealActuator::command)
        .def_property_readonly("actual_value",
                               &rsim::IdealActuator::actualValue);

    py::class_<rsim::ValveActuator, rsim::IdealActuator,
               std::shared_ptr<rsim::ValveActuator>>(module, "ValveActuator")
        .def(py::init<std::string, rsim::Valve&>(), py::arg("name"), py::arg("valve"), py::keep_alive<1, 3>())
        .def("command", &rsim::ValveActuator::command)
        .def_property_readonly("commanded", &rsim::ValveActuator::commanded)
        .def_property_readonly("actual", &rsim::ValveActuator::actual);

    py::class_<rsim::ParachuteActuator, rsim::IdealActuator,
               std::shared_ptr<rsim::ParachuteActuator>>(module, "ParachuteActuator")
        .def(py::init<std::string, rsim::Parachute&>(), py::arg("name"), py::arg("parachute"), py::keep_alive<1, 3>())
        .def("command_deploy", &rsim::ParachuteActuator::commandDeploy)
        .def_property_readonly("commanded", &rsim::ParachuteActuator::commanded)
        .def_property_readonly("actual", &rsim::ParachuteActuator::actual);

    py::class_<rsim::ModelEnableActuator, rsim::IdealActuator,
               std::shared_ptr<rsim::ModelEnableActuator>>(
                   module, "ModelEnableActuator")
        .def(py::init<std::string, rsim::Model&>(),
             py::arg("name"), py::arg("target"), py::keep_alive<1, 3>())
        .def("command", &rsim::ModelEnableActuator::command)
        .def_property_readonly("commanded",
                               &rsim::ModelEnableActuator::commanded)
        .def_property_readonly("actual",
                               &rsim::ModelEnableActuator::actual);

    py::class_<rsim::AtmosphereTable>(module, "AtmosphereTable")
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<rsim::Vector3d>>(),
             py::arg("altitude"), py::arg("density"),
             py::arg("speed_of_sound"), py::arg("dynamic_viscosity"),
             py::arg("wind") = std::vector<rsim::Vector3d>{});

    py::class_<rsim::FlightConditions, rsim::Model,
               std::shared_ptr<rsim::FlightConditions>>(module, "FlightConditions")
        .def(py::init<std::string, rsim::Vehicle&, rsim::AtmosphereTable, rsim::Vector3d, double, rsim::Vector3d>(),
             py::arg("name"), py::arg("vehicle"), py::arg("atmosphere"),
             py::arg("up_direction") = rsim::Vector3d::UnitX(), py::arg("reference_length") = 1.0,
             py::arg("altitude_origin") = rsim::Vector3d::Zero(),
             py::keep_alive<1, 3>())
        .def("add_aerodynamics", &rsim::FlightConditions::addAerodynamics, py::keep_alive<1, 2>())
        .def("add_parachute", &rsim::FlightConditions::addParachute, py::keep_alive<1, 2>())
        .def_property_readonly("altitude", &rsim::FlightConditions::altitude);

    py::class_<rsim::FlightComputer, rsim::Model,
               std::shared_ptr<rsim::FlightComputer>>(module, "FlightComputer")
        .def(py::init<std::string, std::string,
                      std::unordered_map<std::string,
                          std::shared_ptr<rsim::ScalarSensor>>,
                      std::unordered_map<std::string,
                          std::shared_ptr<rsim::IdealActuator>>,
                      std::unordered_map<std::string, double>>(),
             py::arg("name"), py::arg("plugin_path"), py::arg("sensors"),
             py::arg("actuators"),
             py::arg("configuration") =
                 std::unordered_map<std::string, double>{})
        .def_property_readonly("plugin_path", &rsim::FlightComputer::pluginPath);

    py::class_<rsim::RunResult>(module, "RunResult")
        .def_readonly("reason", &rsim::RunResult::reason)
        .def_readonly("event_name", &rsim::RunResult::event_name)
        .def_readonly("final_time", &rsim::RunResult::final_time)
        .def_readonly("step_count", &rsim::RunResult::step_count);

    py::class_<rsim::DataLogger>(module, "DataLogger")
        .def(py::init<double, bool>(), py::arg("sample_period") = 0.0,
             py::arg("include_initial_sample") = true)
        .def("reserve_samples", &rsim::DataLogger::reserveSamples)
        .def("clear", &rsim::DataLogger::clear)
        .def_property_readonly("sample_count", &rsim::DataLogger::sampleCount)
        .def_property_readonly("channel_names", &rsim::DataLogger::channelNames)
        .def_property_readonly("channel_metadata", [](const rsim::DataLogger& logger) {
            py::dict metadata;
            for (const auto& channel : logger.metadata()) {
                py::dict item;
                item["units"] = channel.units;
                switch (channel.scalar_type) {
                    case rsim::StateScalarType::Float64:
                        item["type"] = "float64";
                        break;
                    case rsim::StateScalarType::Int64:
                        item["type"] = "int64";
                        break;
                    case rsim::StateScalarType::Boolean:
                        item["type"] = "bool";
                        break;
                }
                metadata[py::str(channel.name)] = std::move(item);
            }
            return metadata;
        })
        .def("to_dict", [](const rsim::DataLogger& logger) {
            py::dict result;
            const auto& times = logger.times();
            py::array_t<double> time_array(times.size());
            std::copy(times.begin(), times.end(), time_array.mutable_data());
            result["time"] = std::move(time_array);
            for (const std::string& name : logger.channelNames()) {
                const rsim::StateSamples& values = logger.channelValues(name);
                std::visit([&result, &name](const auto& samples) {
                    using Value = typename std::decay_t<decltype(samples)>::value_type;
                    py::array_t<Value> array(samples.size());
                    for (std::size_t index = 0; index < samples.size(); ++index) {
                        array.mutable_data()[index] = samples[index];
                    }
                    result[py::str(name)] = std::move(array);
                }, values);
            }
            return result;
        });

    py::class_<rsim::TruthEvent>(module, "TruthEvent")
        .def(py::init<std::string, std::string, bool, bool>(),
             py::arg("name"), py::arg("trigger"),
             py::arg("stop_simulation") = false,
             py::arg("one_shot") = true)
        .def_property_readonly("name", &rsim::TruthEvent::name)
        .def_property_readonly("trigger", &rsim::TruthEvent::trigger)
        .def_property_readonly("fired", &rsim::TruthEvent::fired)
        .def_property_readonly("fire_count", &rsim::TruthEvent::fireCount)
        .def("add_model_enable_command", &rsim::TruthEvent::addModelEnableCommand,
             py::arg("model"), py::arg("enabled"), py::keep_alive<1, 2>())
        .def("add_parachute_deploy_command", &rsim::TruthEvent::addParachuteDeployCommand,
             py::arg("parachute"), py::keep_alive<1, 2>())
        .def("add_actuator_command", &rsim::TruthEvent::addActuatorCommand,
             py::arg("actuator"), py::arg("value"),
             py::keep_alive<1, 2>());

    py::class_<rsim::Simulator>(module, "Simulator")
        .def(py::init<double>(), py::arg("global_time_step"))
        .def("add_model", &rsim::Simulator::addModel, py::arg("model"),
             // Simulator borrows the model, so Python must keep it alive.
             py::keep_alive<1, 2>())
        .def("add_vehicle", &rsim::Simulator::addVehicle,
             py::arg("vehicle"), py::keep_alive<1, 2>())
        .def("set_vehicle_state_time_steps",
             &rsim::Simulator::setVehicleStateTimeSteps,
             py::arg("vehicle"), py::arg("time_steps"),
             // The copied map has no Python runtime dependency; keep only the vehicle alive.
             py::keep_alive<1, 2>())
        .def("clear_vehicle_state_time_steps",
             &rsim::Simulator::clearVehicleStateTimeSteps,
             py::arg("vehicle"))
        // Release the GIL because the complete simulation loop is native C++.
        .def("step", &rsim::Simulator::step,
             py::call_guard<py::gil_scoped_release>())
        .def("run_steps", &rsim::Simulator::runSteps, py::arg("count"),
             py::call_guard<py::gil_scoped_release>())
        .def("add_data_logger", &rsim::Simulator::addDataLogger,
             py::arg("logger"), py::keep_alive<1, 2>())
        .def("add_truth_event", &rsim::Simulator::addTruthEvent,
             py::arg("event"), py::keep_alive<1, 2>())
        .def("run", &rsim::Simulator::run,
             py::arg("maximum_time") = std::numeric_limits<double>::infinity(),
             py::arg("maximum_steps") = 0U,
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("time", &rsim::Simulator::time)
        .def_property("global_time_step",
                      &rsim::Simulator::globalTimeStep,
                      &rsim::Simulator::setGlobalTimeStep)
        .def_property_readonly("last_time_step",
                               &rsim::Simulator::lastTimeStep)
        .def_property_readonly("step_number", &rsim::Simulator::stepNumber);
}
