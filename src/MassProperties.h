#pragma once

#include "Model.h"
#include "defs.h"

#include <string>
#include <vector>

namespace rsim {

// Supplies one independently evolving contribution to vehicle mass properties.
class VariableMassSource {
public:
    virtual ~VariableMassSource() = default;
    [[nodiscard]] virtual double variableMass() const noexcept = 0;
    [[nodiscard]] virtual Vector3d variableCenterOfMass() const noexcept = 0;
    [[nodiscard]] virtual Matrix3d variableInertiaAboutCenterOfMass()
        const noexcept = 0;
};

// Describes liquid propellant occupying an axis-aligned cylindrical tank.
class CylindricalPropellantTank {
public:
    CylindricalPropellantTank(double radius,
                              double length,
                              double aft_position,
                              double propellant_density,
                              double initial_propellant_mass,
                              double mass_flow_rate = 0.0);

    // Return the geometric propellant capacity in kilograms.
    [[nodiscard]] double capacity() const noexcept;

    // Return the current propellant mass.
    [[nodiscard]] double propellantMass() const noexcept;

    // Replace the current mass and restart consumption at reference_time.
    void setPropellantMass(double mass, double reference_time);

    // Return or replace the nonnegative consumption rate in kilograms per second.
    [[nodiscard]] double massFlowRate() const noexcept;
    void setMassFlowRate(double mass_flow_rate);

    // Recalculate remaining mass from absolute time without cumulative drift.
    void updateMass(double simulation_time) noexcept;

    // Return the liquid center of mass in body coordinates.
    [[nodiscard]] Vector3d centerOfMass() const noexcept;

    // Return liquid inertia about its own center of mass in body axes.
    [[nodiscard]] Matrix3d inertiaAboutCenterOfMass() const noexcept;

private:
    double radius_;
    double length_;
    double aft_position_;
    double density_;
    double propellant_mass_;
    double reference_mass_;
    double reference_time_ = 0.0;
    double mass_flow_rate_;
};

// Stores the mass, center of gravity, and inertia of one rigid body.
class RigidBody : public Model {
public:
    RigidBody(double mass,
              Matrix3d moment_of_inertia,
              Vector3d center_of_gravity = Vector3d::Zero(),
              std::string name = "rigid-body",
              UpdateRate update_rate = UpdateRate::everyStep());

    [[nodiscard]] double mass() const noexcept;
    void setMass(double mass);

    [[nodiscard]] Matrix3d momentOfInertia() const noexcept;
    void setMomentOfInertia(Matrix3d moment_of_inertia);

    [[nodiscard]] Vector3d cg() const noexcept;
    void setCg(Vector3d center_of_gravity) noexcept;

    // Replace the fixed dry mass properties used during tank aggregation.
    void setDryMassProperties(double mass,
                              Matrix3d moment_of_inertia,
                              Vector3d center_of_gravity);

    // Add a liquid tank and return its stable numeric index.
    std::size_t addPropellantTank(CylindricalPropellantTank tank);

    // Return a tank by index so its mass or flow rate can be changed.
    [[nodiscard]] CylindricalPropellantTank& propellantTank(std::size_t index);
    [[nodiscard]] const CylindricalPropellantTank& propellantTank(
        std::size_t index) const;

    // Replace one tank's mass using the model's current simulation time.
    void setPropellantMass(std::size_t index, double mass);

    // Replace one tank's consumption rate.
    void setTankMassFlowRate(std::size_t index, double mass_flow_rate);

    // Return the number of configured liquid tanks.
    [[nodiscard]] std::size_t propellantTankCount() const noexcept;

    // Attach a native variable-mass provider such as a propulsion Tank.
    void addVariableMassSource(VariableMassSource& source);
    void removeVariableMassSource(VariableMassSource& source);
    [[nodiscard]] bool hasVariableMassSource(
        const VariableMassSource& source) const noexcept;
    [[nodiscard]] std::size_t variableMassSourceCount() const noexcept;

    // Recalculate total mass, center of mass, and inertia for the current time.
    void update() override;

    void printMass() const;

private:
    Matrix3d moment_of_inertia_;
    double mass_;
    Vector3d center_of_gravity_;

    double dry_mass_;
    Matrix3d dry_moment_of_inertia_;
    Vector3d dry_center_of_gravity_;
    std::vector<CylindricalPropellantTank> propellant_tanks_;
    std::vector<VariableMassSource*> variable_mass_sources_;
};

}  // namespace rsim
