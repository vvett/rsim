#include "MassProperties.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rsim {

namespace {

// Shift a centroidal inertia tensor to another point using the parallel-axis theorem.
Matrix3d shiftInertia(Matrix3d inertia,
                      double mass,
                      const Vector3d& displacement) {
    return inertia + mass *
        (displacement.squaredNorm() * Matrix3d::Identity() -
         displacement * displacement.transpose());
}

}  // namespace

CylindricalPropellantTank::CylindricalPropellantTank(
    double radius,
    double length,
    double aft_position,
    double propellant_density,
    double initial_propellant_mass,
    double mass_flow_rate)
    : radius_(radius),
      length_(length),
      aft_position_(aft_position),
      density_(propellant_density),
      propellant_mass_(initial_propellant_mass),
      reference_mass_(initial_propellant_mass),
      mass_flow_rate_(mass_flow_rate) {
    if (radius <= 0.0 || length <= 0.0 || propellant_density <= 0.0) {
        throw std::invalid_argument("tank dimensions and density must be positive");
    }
    if (initial_propellant_mass < 0.0 ||
        initial_propellant_mass > capacity()) {
        throw std::invalid_argument("initial propellant mass exceeds tank capacity");
    }
    setMassFlowRate(mass_flow_rate);
}

double CylindricalPropellantTank::capacity() const noexcept {
    constexpr double pi = 3.14159265358979323846;
    return density_ * pi * radius_ * radius_ * length_;
}

double CylindricalPropellantTank::propellantMass() const noexcept {
    return propellant_mass_;
}

void CylindricalPropellantTank::setPropellantMass(
    double mass, double reference_time) {
    if (mass < 0.0 || mass > capacity()) {
        throw std::invalid_argument("propellant mass must be within tank capacity");
    }
    propellant_mass_ = mass;
    reference_mass_ = mass;
    reference_time_ = reference_time;
}

double CylindricalPropellantTank::massFlowRate() const noexcept {
    return mass_flow_rate_;
}

void CylindricalPropellantTank::setMassFlowRate(double mass_flow_rate) {
    if (mass_flow_rate < 0.0) {
        throw std::invalid_argument("mass flow rate must be nonnegative");
    }
    mass_flow_rate_ = mass_flow_rate;
}

void CylindricalPropellantTank::updateMass(double simulation_time) noexcept {
    const double elapsed = std::max(0.0, simulation_time - reference_time_);
    propellant_mass_ = std::max(0.0, reference_mass_ - mass_flow_rate_ * elapsed);
}

Vector3d CylindricalPropellantTank::centerOfMass() const noexcept {
    const double fill_length = length_ * propellant_mass_ / capacity();
    return Vector3d{aft_position_ + 0.5 * fill_length, 0.0, 0.0};
}

Matrix3d CylindricalPropellantTank::inertiaAboutCenterOfMass() const noexcept {
    const double fill_length = length_ * propellant_mass_ / capacity();
    const double axial = 0.5 * propellant_mass_ * radius_ * radius_;
    const double transverse = propellant_mass_ *
        (3.0 * radius_ * radius_ + fill_length * fill_length) / 12.0;
    Matrix3d inertia = Matrix3d::Zero();
    inertia.diagonal() << axial, transverse, transverse;
    return inertia;
}

RigidBody::RigidBody(double mass,
                     Matrix3d moment_of_inertia,
                     Vector3d center_of_gravity,
                     std::string name,
                     UpdateRate update_rate)
    : Model(std::move(name), update_rate),
      moment_of_inertia_(std::move(moment_of_inertia)),
      mass_(mass),
      center_of_gravity_(std::move(center_of_gravity)),
      dry_mass_(mass),
      dry_moment_of_inertia_(moment_of_inertia_),
      dry_center_of_gravity_(center_of_gravity_) {
    setPhase(ModelPhase::MassProperties);
    if (mass <= 0.0) {
        throw std::invalid_argument("mass must be greater than zero");
    }
    addStateVariable(mass_, "mass", "kg");
}

double RigidBody::mass() const noexcept {
    return mass_;
}

void RigidBody::setMass(double mass) {
    if (mass <= 0.0) {
        throw std::invalid_argument("mass must be greater than zero");
    }
    dry_mass_ = mass;
    update();
}

Matrix3d RigidBody::momentOfInertia() const noexcept {
    return moment_of_inertia_;
}

void RigidBody::setMomentOfInertia(Matrix3d moment_of_inertia) {
    dry_moment_of_inertia_ = std::move(moment_of_inertia);
    update();
}

Vector3d RigidBody::cg() const noexcept {
    return center_of_gravity_;
}

void RigidBody::setCg(Vector3d center_of_gravity) noexcept {
    dry_center_of_gravity_ = std::move(center_of_gravity);
    update();
}

void RigidBody::setDryMassProperties(
    double mass,
    Matrix3d moment_of_inertia,
    Vector3d center_of_gravity) {
    if (mass <= 0.0) {
        throw std::invalid_argument("dry mass must be greater than zero");
    }
    dry_mass_ = mass;
    dry_moment_of_inertia_ = std::move(moment_of_inertia);
    dry_center_of_gravity_ = std::move(center_of_gravity);
    update();
}

std::size_t RigidBody::addPropellantTank(CylindricalPropellantTank tank) {
    propellant_tanks_.push_back(std::move(tank));
    update();
    return propellant_tanks_.size() - 1U;
}

CylindricalPropellantTank& RigidBody::propellantTank(std::size_t index) {
    return propellant_tanks_.at(index);
}

const CylindricalPropellantTank& RigidBody::propellantTank(
    std::size_t index) const {
    return propellant_tanks_.at(index);
}

void RigidBody::setPropellantMass(std::size_t index, double mass) {
    propellant_tanks_.at(index).setPropellantMass(mass, simulationTime());
    update();
}

void RigidBody::setTankMassFlowRate(
    std::size_t index, double mass_flow_rate) {
    CylindricalPropellantTank& tank = propellant_tanks_.at(index);
    tank.updateMass(simulationTime());
    // Restart the time-based burn law so a rate change is not applied retroactively.
    tank.setPropellantMass(tank.propellantMass(), simulationTime());
    tank.setMassFlowRate(mass_flow_rate);
}

std::size_t RigidBody::propellantTankCount() const noexcept {
    return propellant_tanks_.size();
}

void RigidBody::addVariableMassSource(VariableMassSource& source) {
    if (hasVariableMassSource(source)) {
        throw std::invalid_argument("variable mass source is already registered");
    }
    variable_mass_sources_.push_back(&source);
    update();
}

void RigidBody::removeVariableMassSource(VariableMassSource& source) {
    const auto found = std::find(variable_mass_sources_.begin(),
                                 variable_mass_sources_.end(), &source);
    if (found == variable_mass_sources_.end()) {
        throw std::invalid_argument("variable mass source is not registered");
    }
    variable_mass_sources_.erase(found);
    update();
}

bool RigidBody::hasVariableMassSource(
    const VariableMassSource& source) const noexcept {
    return std::find(variable_mass_sources_.begin(),
                     variable_mass_sources_.end(), &source) !=
           variable_mass_sources_.end();
}

std::size_t RigidBody::variableMassSourceCount() const noexcept {
    return variable_mass_sources_.size();
}

void RigidBody::update() {
    double total_mass = dry_mass_;
    Vector3d first_moment = dry_mass_ * dry_center_of_gravity_;
    for (CylindricalPropellantTank& tank : propellant_tanks_) {
        tank.updateMass(simulationTime());
        total_mass += tank.propellantMass();
        first_moment += tank.propellantMass() * tank.centerOfMass();
    }
    for (const VariableMassSource* source : variable_mass_sources_) {
        const double source_mass = source->variableMass();
        if (source_mass < 0.0 || !std::isfinite(source_mass)) {
            throw std::runtime_error(
                "variable mass source returned an invalid mass");
        }
        total_mass += source_mass;
        first_moment += source_mass * source->variableCenterOfMass();
    }

    mass_ = total_mass;
    center_of_gravity_ = first_moment / total_mass;
    moment_of_inertia_ = shiftInertia(
        dry_moment_of_inertia_, dry_mass_,
        dry_center_of_gravity_ - center_of_gravity_);
    for (const CylindricalPropellantTank& tank : propellant_tanks_) {
        moment_of_inertia_ += shiftInertia(
            tank.inertiaAboutCenterOfMass(), tank.propellantMass(),
            tank.centerOfMass() - center_of_gravity_);
    }
    for (const VariableMassSource* source : variable_mass_sources_) {
        moment_of_inertia_ += shiftInertia(
            source->variableInertiaAboutCenterOfMass(),
            source->variableMass(),
            source->variableCenterOfMass() - center_of_gravity_);
    }
}

void RigidBody::printMass() const {
    std::cout << "Body mass is " << mass_ << '\n';
}

}  // namespace rsim
