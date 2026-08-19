#pragma once

#include "ForceTorque.h"
#include "MassProperties.h"
#include "Model.h"
#include "defs.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace rsim {

enum class FluidPhase { Liquid, Gas };
enum class InletRole { None, Oxidizer, Fuel };
enum class BlowdownMode {
    Adiabatic,
    Isentropic,
    Polytropic,
    EnergyConserving
};

// Base for native nodes. Nodes deliberately are not Models; Propulsion owns
// and advances all network state as one coupled model.
class PropulsionNode {
public:
    explicit PropulsionNode(std::string name);
    virtual ~PropulsionNode() = default;
    PropulsionNode(const PropulsionNode&) = delete;
    PropulsionNode& operator=(const PropulsionNode&) = delete;

    [[nodiscard]] const std::string& name() const noexcept;
    [[nodiscard]] virtual double pressure() const noexcept = 0;
    [[nodiscard]] virtual double temperature() const noexcept = 0;
    [[nodiscard]] virtual double density(FluidPhase phase) const noexcept;
    [[nodiscard]] virtual double heatCapacityRatio() const noexcept;
    [[nodiscard]] virtual double molarMass() const noexcept;
    [[nodiscard]] virtual bool phaseAvailable(FluidPhase phase) const noexcept;
    [[nodiscard]] virtual double specificEnthalpy(
        FluidPhase phase) const noexcept;
    [[nodiscard]] virtual const std::string& fluidSpecies(
        FluidPhase phase) const noexcept;

private:
    std::string name_;
};

class Ambient final : public PropulsionNode {
public:
    Ambient(std::string name = "ambient", double pressure = 101325.0,
            double temperature = 288.15);
    [[nodiscard]] double pressure() const noexcept override;
    [[nodiscard]] double temperature() const noexcept override;
    void setPressure(double pressure);
    void setTemperature(double temperature);

private:
    double pressure_;
    double temperature_;
};

// A settled cylindrical liquid tank with an ideal-gas ullage closure. Conserved
// liquid/pressurant masses and internal energies are exposed for diagnostics.
class Tank final : public PropulsionNode, public VariableMassSource {
public:
    Tank(std::string name, double radius, double length, double aft_position,
         double liquid_density, double propellant_mass, double pressure,
         double temperature = 293.15, double pressurant_mass = 0.0,
         double collapse_factor = 1.0,
         BlowdownMode blowdown_mode = BlowdownMode::Adiabatic,
         double polytropic_exponent = 1.2,
         std::string propellant = "propellant",
         std::string pressurant = "Helium");

    [[nodiscard]] double pressure() const noexcept override;
    [[nodiscard]] double temperature() const noexcept override;
    [[nodiscard]] double density(FluidPhase phase) const noexcept override;
    [[nodiscard]] bool phaseAvailable(FluidPhase phase) const noexcept override;
    [[nodiscard]] double heatCapacityRatio() const noexcept override;
    [[nodiscard]] double molarMass() const noexcept override;
    [[nodiscard]] double specificEnthalpy(
        FluidPhase phase) const noexcept override;
    [[nodiscard]] const std::string& fluidSpecies(
        FluidPhase phase) const noexcept override;

    [[nodiscard]] double volume() const noexcept;
    [[nodiscard]] double liquidVolume() const noexcept;
    [[nodiscard]] double ullageVolume() const noexcept;
    [[nodiscard]] double propellantMass() const noexcept;
    [[nodiscard]] double pressurantMass() const noexcept;
    [[nodiscard]] double propellantInternalEnergy() const noexcept;
    [[nodiscard]] double pressurantInternalEnergy() const noexcept;
    [[nodiscard]] double collapseFactor() const noexcept;
    [[nodiscard]] BlowdownMode blowdownMode() const noexcept;
    [[nodiscard]] const std::string& propellant() const noexcept;
    [[nodiscard]] const std::string& pressurant() const noexcept;
    void setPropellantMass(double mass);

    [[nodiscard]] double variableMass() const noexcept override;
    [[nodiscard]] Vector3d variableCenterOfMass() const noexcept override;
    [[nodiscard]] Matrix3d variableInertiaAboutCenterOfMass()
        const noexcept override;

private:
    friend class Propulsion;
    void validateTransfer(FluidPhase phase, double mass_delta,
                          double energy_delta) const;
    void applyTransfer(FluidPhase phase, double mass_delta,
                       double energy_delta);
    void refreshPressure() noexcept;

    double radius_;
    double length_;
    double aft_position_;
    double liquid_density_;
    double propellant_mass_;
    double pressurant_mass_;
    double pressure_;
    double temperature_;
    double propellant_energy_;
    double pressurant_energy_;
    double collapse_factor_;
    BlowdownMode blowdown_mode_;
    double polytropic_exponent_;
    std::string propellant_;
    std::string pressurant_;
    double initial_pressure_;
    double initial_temperature_;
    double initial_ullage_volume_;
    double initial_pressurant_mass_;
};

class Junction : public PropulsionNode {
public:
    Junction(std::string name, double pressure = 101325.0,
             double temperature = 293.15);
    [[nodiscard]] double pressure() const noexcept override;
    [[nodiscard]] double temperature() const noexcept override;

protected:
    friend class Propulsion;
    void setAlgebraicPressure(double pressure) noexcept;
    double pressure_;
    double temperature_;
};

// In-memory chamber-performance table. Values use MATLAB/Fortran ordering:
// flat index = pc_index + pc_count * mixture_ratio_index.
class CombustionTableData {
public:
    CombustionTableData(std::vector<double> chamber_pressure,
                        std::vector<double> mixture_ratio,
                        std::vector<double> chamber_temperature,
                        std::vector<double> cstar,
                        std::vector<double> cf_vac,
                        std::vector<double> isp_vac,
                        double expansion_ratio,
                        std::string oxidizer = {}, std::string fuel = {},
                        bool frozen_at_throat = false);

    [[nodiscard]] double chamberTemperature(double pc, double mr) const;
    [[nodiscard]] double cstar(double pc, double mr) const;
    [[nodiscard]] double cfVac(double pc, double mr) const;
    [[nodiscard]] double ispVac(double pc, double mr) const;
    [[nodiscard]] double expansionRatio() const noexcept;
    [[nodiscard]] const std::string& oxidizer() const noexcept;
    [[nodiscard]] const std::string& fuel() const noexcept;
    [[nodiscard]] bool frozenAtThroat() const noexcept;
    [[nodiscard]] const std::vector<double>& pressureAxis() const noexcept;
    [[nodiscard]] const std::vector<double>& mixtureRatioAxis() const noexcept;

private:
    [[nodiscard]] double interpolate(const std::vector<double>& values,
                                     double pc, double mr) const;
    std::vector<double> pc_;
    std::vector<double> mr_;
    std::vector<double> chamber_temperature_;
    std::vector<double> cstar_;
    std::vector<double> cf_vac_;
    std::vector<double> isp_vac_;
    double expansion_ratio_;
    std::string oxidizer_;
    std::string fuel_;
    bool frozen_at_throat_;
};

class Combustor : public PropulsionNode {
public:
    Combustor(std::string name, std::shared_ptr<CombustionTableData> table,
              double initial_pressure = 101325.0,
              double initial_mixture_ratio = 1.0);
    [[nodiscard]] double pressure() const noexcept override;
    [[nodiscard]] double temperature() const noexcept override;
    [[nodiscard]] double mixtureRatio() const noexcept;
    [[nodiscard]] double oxidizerMassFlow() const noexcept;
    [[nodiscard]] double fuelMassFlow() const noexcept;
    [[nodiscard]] double totalMassFlow() const noexcept;

protected:
    friend class Propulsion;
    virtual double throatAreaForClosure() const noexcept;
    void cacheChamber(double pressure, double oxidizer_flow, double fuel_flow,
                      double temperature) noexcept;
    std::shared_ptr<CombustionTableData> table_;

private:
    double pressure_;
    double temperature_ = 293.15;
    double mixture_ratio_;
    double oxidizer_flow_ = 0.0;
    double fuel_flow_ = 0.0;
};

class Engine final : public Combustor, public ForceTorque {
public:
    Engine(std::string name, std::shared_ptr<CombustionTableData> table,
           double throat_area, double expansion_ratio,
           Vector3d thrust_direction = Vector3d::UnitX(),
           Vector3d application_point = Vector3d::Zero(),
           std::shared_ptr<Ambient> ambient = nullptr,
           double initial_pressure = 101325.0,
           double initial_mixture_ratio = 1.0);
    void update() override;
    [[nodiscard]] double throatArea() const noexcept;
    [[nodiscard]] double expansionRatio() const noexcept;
    [[nodiscard]] double thrust() const noexcept;
    [[nodiscard]] double ambientPressure() const noexcept;
    [[nodiscard]] double cstar() const noexcept;
    [[nodiscard]] double thrustCoefficient() const noexcept;

protected:
    double throatAreaForClosure() const noexcept override;

private:
    friend class Propulsion;
    void cachePerformance(double thrust, double cstar, double cf) noexcept;
    double throat_area_;
    double expansion_ratio_;
    Vector3d direction_;
    Vector3d application_point_;
    std::shared_ptr<Ambient> ambient_;
    double cached_thrust_ = 0.0;
    double cached_cstar_ = 0.0;
    double cached_cf_ = 0.0;
};

class FlowComponent {
public:
    virtual ~FlowComponent() = default;
    [[nodiscard]] virtual double massFlow(const PropulsionNode& upstream,
                                          const PropulsionNode& downstream,
                                          FluidPhase phase) const = 0;
};

class Orifice : public FlowComponent {
public:
    explicit Orifice(double cd_area);
    Orifice(double area, double discharge_coefficient);
    [[nodiscard]] double cdArea() const noexcept;
    [[nodiscard]] double massFlow(const PropulsionNode& upstream,
                                  const PropulsionNode& downstream,
                                  FluidPhase phase) const override;

protected:
    double cd_area_;
};

class Valve final : public Orifice {
public:
    Valve(double cd_area, double opening = 1.0);
    Valve(double area, double discharge_coefficient, double opening);
    [[nodiscard]] double opening() const noexcept;
    void setOpening(double opening);
    [[nodiscard]] double massFlow(const PropulsionNode& upstream,
                                  const PropulsionNode& downstream,
                                  FluidPhase phase) const override;

private:
    double opening_;
};

class Line final : public FlowComponent {
public:
    Line(double length, double diameter, double roughness = 0.0,
         double dynamic_viscosity = 1.0e-3);
    [[nodiscard]] double massFlow(const PropulsionNode& upstream,
                                  const PropulsionNode& downstream,
                                  FluidPhase phase) const override;

private:
    double length_;
    double diameter_;
    double roughness_;
    double viscosity_;
};

class CompoundEdge final : public FlowComponent {
public:
    explicit CompoundEdge(std::vector<std::shared_ptr<FlowComponent>> parts);
    [[nodiscard]] const std::vector<double>& intermediatePressures() const noexcept;
    [[nodiscard]] bool containsUnimplementedComponent() const noexcept;
    [[nodiscard]] double massFlow(const PropulsionNode& upstream,
                                  const PropulsionNode& downstream,
                                  FluidPhase phase) const override;

private:
    std::vector<std::shared_ptr<FlowComponent>> parts_;
    mutable std::vector<double> intermediate_pressures_;
};

class Pump final : public FlowComponent {
public:
    Pump() = default;
    [[nodiscard]] double massFlow(const PropulsionNode&, const PropulsionNode&,
                                  FluidPhase) const override;
};

class Turbine final : public FlowComponent {
public:
    Turbine() = default;
    [[nodiscard]] double massFlow(const PropulsionNode&, const PropulsionNode&,
                                  FluidPhase) const override;
};

struct ThermalContact {
    std::shared_ptr<PropulsionNode> node;
    FluidPhase phase = FluidPhase::Liquid;
    double environment_temperature = 293.15;
    double heat_flux = 0.0;
    double area = 0.0;
};

class Propulsion final : public Model {
public:
    explicit Propulsion(std::string name = "propulsion",
                        UpdateRate update_rate = UpdateRate::everyStep());
    void connect(std::shared_ptr<PropulsionNode> upstream,
                 std::shared_ptr<PropulsionNode> downstream,
                 std::shared_ptr<FlowComponent> edge,
                 FluidPhase upstream_phase = FluidPhase::Liquid,
                 InletRole downstream_role = InletRole::None);
    void connect(std::shared_ptr<PropulsionNode> upstream,
                 std::shared_ptr<PropulsionNode> downstream,
                 std::vector<std::shared_ptr<FlowComponent>> edges,
                 FluidPhase upstream_phase = FluidPhase::Liquid,
                 InletRole downstream_role = InletRole::None);
    void vent(std::shared_ptr<PropulsionNode> upstream,
              std::shared_ptr<Ambient> ambient,
              std::shared_ptr<FlowComponent> edge,
              FluidPhase upstream_phase = FluidPhase::Gas);
    void addThermalContact(ThermalContact contact);
    void finalize();
    [[nodiscard]] bool finalized() const noexcept;
    void update() override;
    [[nodiscard]] std::size_t connectionCount() const noexcept;
    [[nodiscard]] double connectionMassFlow(std::size_t index) const;
    [[nodiscard]] double totalStoredMass() const noexcept;
    [[nodiscard]] const std::string& lastSolveSummary() const noexcept;

private:
    struct Connection {
        std::shared_ptr<PropulsionNode> upstream;
        std::shared_ptr<PropulsionNode> downstream;
        std::shared_ptr<FlowComponent> edge;
        FluidPhase upstream_phase;
        InletRole downstream_role;
        double mass_flow = 0.0;
    };
    void requireMutable() const;
    void ownNode(const std::shared_ptr<PropulsionNode>& node);
    void solveAlgebraicFlows();
    std::vector<std::shared_ptr<PropulsionNode>> nodes_;
    std::vector<std::shared_ptr<FlowComponent>> edges_;
    std::vector<Connection> connections_;
    std::vector<ThermalContact> thermal_contacts_;
    bool finalized_ = false;
    std::string last_solve_summary_ = "network not yet solved";
};

}  // namespace rsim
