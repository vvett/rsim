#include "Propulsion.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace rsim {
namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double universal_gas_constant = 8.31446261815324;

double requirePositive(double value, const char* label) {
    if (!(value > 0.0) || !std::isfinite(value)) {
        throw std::invalid_argument(std::string(label) + " must be positive");
    }
    return value;
}

double signOf(double value) noexcept { return value < 0.0 ? -1.0 : 1.0; }

double gasRestrictionFlow(const PropulsionNode& high,
                          const PropulsionNode& low,
                          double cd_area) {
    const double p0 = high.pressure();
    const double p1 = std::max(0.0, low.pressure());
    if (p0 <= 0.0 || p1 >= p0) {
        return 0.0;
    }
    const double gamma = std::max(1.01, high.heatCapacityRatio());
    const double specific_r = universal_gas_constant /
                              std::max(1e-6, high.molarMass());
    const double temperature = std::max(1.0, high.temperature());
    const double ratio = p1 / p0;
    const double critical = std::pow(2.0 / (gamma + 1.0),
                                     gamma / (gamma - 1.0));
    double flux = 0.0;
    if (ratio <= critical) {
        flux = p0 * std::sqrt(gamma / (specific_r * temperature)) *
               std::pow(2.0 / (gamma + 1.0),
                        (gamma + 1.0) / (2.0 * (gamma - 1.0)));
    } else {
        const double term = 2.0 * gamma / (specific_r * temperature *
                                            (gamma - 1.0)) *
            (std::pow(ratio, 2.0 / gamma) -
             std::pow(ratio, (gamma + 1.0) / gamma));
        flux = p0 * std::sqrt(std::max(0.0, term));
    }
    return cd_area * flux;
}

std::size_t lowerCell(const std::vector<double>& axis, double value) {
    if (value < axis.front() || value > axis.back()) {
        throw std::out_of_range("combustion-table query would extrapolate");
    }
    if (value == axis.back()) {
        return axis.size() - 2U;
    }
    return static_cast<std::size_t>(
        std::upper_bound(axis.begin(), axis.end(), value) - axis.begin() - 1);
}

void validateAxis(const std::vector<double>& axis, const char* label) {
    if (axis.size() < 2U) {
        throw std::invalid_argument(std::string(label) +
                                    " axis needs at least two values");
    }
    for (std::size_t i = 0; i < axis.size(); ++i) {
        if (!std::isfinite(axis[i]) ||
            (i > 0U && axis[i] <= axis[i - 1U])) {
            throw std::invalid_argument(std::string(label) +
                                        " axis must be finite and increasing");
        }
    }
}

}  // namespace

PropulsionNode::PropulsionNode(std::string name) : name_(std::move(name)) {
    if (name_.empty()) {
        throw std::invalid_argument("propulsion node name must not be empty");
    }
}

const std::string& PropulsionNode::name() const noexcept { return name_; }

double PropulsionNode::density(FluidPhase phase) const noexcept {
    if (phase == FluidPhase::Liquid) {
        return 1000.0;
    }
    return pressure() * molarMass() /
           (universal_gas_constant * std::max(1.0, temperature()));
}

double PropulsionNode::heatCapacityRatio() const noexcept { return 1.4; }
double PropulsionNode::molarMass() const noexcept { return 0.02897; }
bool PropulsionNode::phaseAvailable(FluidPhase) const noexcept { return true; }
double PropulsionNode::specificEnthalpy(FluidPhase phase) const noexcept {
    if (phase == FluidPhase::Liquid) return 2000.0 * temperature();
    const double specific_r = universal_gas_constant / molarMass();
    return heatCapacityRatio() / (heatCapacityRatio() - 1.0) *
           specific_r * temperature();
}
const std::string& PropulsionNode::fluidSpecies(FluidPhase) const noexcept {
    static const std::string unspecified;
    return unspecified;
}

Ambient::Ambient(std::string name, double pressure, double temperature)
    : PropulsionNode(std::move(name)),
      pressure_(pressure), temperature_(temperature) {
    if (pressure < 0.0 || !std::isfinite(pressure)) {
        throw std::invalid_argument("ambient pressure must be nonnegative");
    }
    requirePositive(temperature, "ambient temperature");
}
double Ambient::pressure() const noexcept { return pressure_; }
double Ambient::temperature() const noexcept { return temperature_; }
void Ambient::setPressure(double pressure) {
    if (pressure < 0.0 || !std::isfinite(pressure)) {
        throw std::invalid_argument("ambient pressure must be nonnegative");
    }
    pressure_ = pressure;
}
void Ambient::setTemperature(double temperature) {
    temperature_ = requirePositive(temperature, "ambient temperature");
}

Tank::Tank(std::string name, double radius, double length,
           double aft_position, double liquid_density,
           double propellant_mass, double pressure, double temperature,
           double pressurant_mass, double collapse_factor,
           BlowdownMode blowdown_mode, double polytropic_exponent,
           std::string propellant, std::string pressurant)
    : PropulsionNode(std::move(name)),
      radius_(requirePositive(radius, "tank radius")),
      length_(requirePositive(length, "tank length")),
      aft_position_(aft_position),
      liquid_density_(requirePositive(liquid_density, "liquid density")),
      propellant_mass_(propellant_mass),
      pressurant_mass_(pressurant_mass),
      pressure_(requirePositive(pressure, "tank pressure")),
      temperature_(requirePositive(temperature, "tank temperature")),
      propellant_energy_(propellant_mass * 2000.0 * temperature),
      pressurant_energy_(0.0),
      collapse_factor_(requirePositive(collapse_factor, "collapse factor")),
      blowdown_mode_(blowdown_mode),
      polytropic_exponent_(polytropic_exponent),
      propellant_(std::move(propellant)),
      pressurant_(std::move(pressurant)),
      initial_pressure_(pressure), initial_temperature_(temperature),
      initial_ullage_volume_(0.0),
      initial_pressurant_mass_(pressurant_mass) {
    if (propellant_mass < 0.0 || propellant_mass > liquid_density_ * volume()) {
        throw std::invalid_argument("propellant mass is outside tank capacity");
    }
    if (pressurant_mass < 0.0 || collapse_factor < 1.0) {
        throw std::invalid_argument(
            "pressurant mass must be nonnegative and collapse factor >= 1");
    }
    if (blowdown_mode == BlowdownMode::Polytropic &&
        (!(polytropic_exponent > 0.0) || !std::isfinite(polytropic_exponent))) {
        throw std::invalid_argument("polytropic exponent must be positive");
    }
    initial_ullage_volume_ = ullageVolume();
    if (initial_ullage_volume_ <= 0.0) {
        throw std::invalid_argument("tank requires nonzero initial ullage");
    }
    if (pressurant_mass_ == 0.0) {
        // Collapse factor is actual mass / ideal mass. Pressure calculations use
        // the effective ideal fraction, avoiding fictitious extra pressure.
        pressurant_mass_ = collapse_factor_ * pressure_ * initial_ullage_volume_ *
            molarMass() / (universal_gas_constant * temperature_);
        initial_pressurant_mass_ = pressurant_mass_;
    }
    pressurant_energy_ = pressurant_mass_ *
        universal_gas_constant / molarMass() /
        (heatCapacityRatio() - 1.0) * temperature_;
}

double Tank::pressure() const noexcept { return pressure_; }
double Tank::temperature() const noexcept { return temperature_; }
double Tank::density(FluidPhase phase) const noexcept {
    return phase == FluidPhase::Liquid
        ? liquid_density_
        : pressure_ * molarMass() /
          (universal_gas_constant * std::max(1.0, temperature_));
}
bool Tank::phaseAvailable(FluidPhase phase) const noexcept {
    return phase == FluidPhase::Liquid ? propellant_mass_ > 1e-12
                                       : pressurant_mass_ > 1e-12;
}
double Tank::heatCapacityRatio() const noexcept {
    return pressurant_ == "Helium" || pressurant_ == "helium" ? 5.0 / 3.0 : 1.4;
}
double Tank::molarMass() const noexcept {
    if (pressurant_ == "Helium" || pressurant_ == "helium") return 0.004002602;
    if (pressurant_ == "Nitrogen" || pressurant_ == "nitrogen") return 0.0280134;
    return 0.02897;
}
double Tank::specificEnthalpy(FluidPhase phase) const noexcept {
    const double mass = phase == FluidPhase::Liquid
        ? propellant_mass_ : pressurant_mass_;
    const double energy = phase == FluidPhase::Liquid
        ? propellant_energy_ : pressurant_energy_;
    const double specific_internal_energy = mass > 1e-12
        ? energy / mass : 0.0;
    return specific_internal_energy + pressure_ /
        std::max(1e-12, density(phase));
}
const std::string& Tank::fluidSpecies(FluidPhase phase) const noexcept {
    return phase == FluidPhase::Liquid ? propellant_ : pressurant_;
}
double Tank::volume() const noexcept { return pi * radius_ * radius_ * length_; }
double Tank::liquidVolume() const noexcept {
    return std::min(volume(), propellant_mass_ / liquid_density_);
}
double Tank::ullageVolume() const noexcept {
    return std::max(1e-12, volume() - liquidVolume());
}
double Tank::propellantMass() const noexcept { return propellant_mass_; }
double Tank::pressurantMass() const noexcept { return pressurant_mass_; }
double Tank::propellantInternalEnergy() const noexcept { return propellant_energy_; }
double Tank::pressurantInternalEnergy() const noexcept { return pressurant_energy_; }
double Tank::collapseFactor() const noexcept { return collapse_factor_; }
BlowdownMode Tank::blowdownMode() const noexcept { return blowdown_mode_; }
const std::string& Tank::propellant() const noexcept { return propellant_; }
const std::string& Tank::pressurant() const noexcept { return pressurant_; }
void Tank::setPropellantMass(double mass) {
    if (mass < 0.0 || mass > liquid_density_ * volume()) {
        throw std::invalid_argument("propellant mass is outside tank capacity");
    }
    propellant_mass_ = mass;
    propellant_energy_ = mass * 2000.0 * temperature_;
    refreshPressure();
}
double Tank::variableMass() const noexcept {
    return propellant_mass_ + pressurant_mass_;
}
Vector3d Tank::variableCenterOfMass() const noexcept {
    const double fill_length = liquidVolume() / (pi * radius_ * radius_);
    const Vector3d liquid_center(aft_position_ + 0.5 * fill_length, 0.0, 0.0);
    const Vector3d gas_center(aft_position_ + 0.5 * (length_ + fill_length),
                              0.0, 0.0);
    const double total = variableMass();
    if (total <= 0.0) return Vector3d(aft_position_ + 0.5 * length_, 0.0, 0.0);
    return (propellant_mass_ * liquid_center + pressurant_mass_ * gas_center) /
           total;
}
Matrix3d Tank::variableInertiaAboutCenterOfMass() const noexcept {
    const double mass = variableMass();
    const double fill_length = liquidVolume() / (pi * radius_ * radius_);
    Matrix3d inertia = Matrix3d::Zero();
    inertia(0, 0) = 0.5 * mass * radius_ * radius_;
    inertia(1, 1) = inertia(2, 2) = mass *
        (3.0 * radius_ * radius_ + fill_length * fill_length) / 12.0;
    return inertia;
}
void Tank::validateTransfer(FluidPhase phase, double mass_delta,
                            double energy_delta) const {
    if (!std::isfinite(mass_delta) || !std::isfinite(energy_delta)) {
        throw std::runtime_error("tank transfer must be finite");
    }
    const double new_mass = (phase == FluidPhase::Liquid
        ? propellant_mass_ : pressurant_mass_) + mass_delta;
    if (new_mass < -1e-12) {
        std::ostringstream message;
        message << "tank transfer exceeds available " << fluidSpecies(phase)
                << " mass in " << name() << ": current="
                << (phase == FluidPhase::Liquid ? propellant_mass_
                                                : pressurant_mass_)
                << " kg, requested delta=" << mass_delta << " kg";
        throw std::runtime_error(message.str());
    }
    if (phase == FluidPhase::Liquid &&
        new_mass > liquid_density_ * volume() + 1e-12) {
        throw std::runtime_error("liquid transfer exceeds capacity of " + name());
    }
}
void Tank::applyTransfer(FluidPhase phase, double mass_delta,
                         double energy_delta) {
    validateTransfer(phase, mass_delta, energy_delta);
    if (phase == FluidPhase::Liquid) {
        propellant_mass_ = std::max(0.0, propellant_mass_ + mass_delta);
        propellant_energy_ += energy_delta;
    } else {
        pressurant_mass_ = std::max(0.0, pressurant_mass_ + mass_delta);
        pressurant_energy_ += energy_delta;
    }
    refreshPressure();
}
void Tank::refreshPressure() noexcept {
    const double mass_ratio = initial_pressurant_mass_ > 0.0
        ? pressurant_mass_ / initial_pressurant_mass_ : 0.0;
    if (blowdown_mode_ == BlowdownMode::EnergyConserving) {
        // Apply dU = -h_out dm to a well-mixed ideal gas. Since
        // U = m*Cv*T and h_out = Cp*T, this is evaluated directly from
        // the conserved internal energy carried by the tank state.
        const double cv = universal_gas_constant /
            (molarMass() * (heatCapacityRatio() - 1.0));
        if (pressurant_mass_ > 1e-12 && cv > 0.0) {
            temperature_ = std::max(1.0,
                pressurant_energy_ / (pressurant_mass_ * cv));
        } else {
            temperature_ = 1.0;
        }
        const double volume_ratio = initial_ullage_volume_ / ullageVolume();
        pressure_ = std::max(1.0, initial_pressure_ * mass_ratio *
            (temperature_ / initial_temperature_) * volume_ratio);
        return;
    }
    double exponent = 1.0;
    if (blowdown_mode_ == BlowdownMode::Adiabatic) {
        exponent = heatCapacityRatio();
    } else if (blowdown_mode_ == BlowdownMode::Polytropic) {
        exponent = polytropic_exponent_;
    } else {
        exponent = heatCapacityRatio();
    }
    const double volume_ratio = initial_ullage_volume_ / ullageVolume();
    pressure_ = std::max(1.0, initial_pressure_ * mass_ratio *
                                  std::pow(volume_ratio, exponent));
    // For p V^n = constant, the ideal-gas relation gives
    // T/T0 = (V0/V)^(n-1). This also applies to the adiabatic/isentropic
    // cases with n equal to the heat-capacity ratio.
    temperature_ = std::max(1.0, initial_temperature_ *
        std::pow(std::max(1e-12, volume_ratio), exponent - 1.0));
}

Junction::Junction(std::string name, double pressure, double temperature)
    : PropulsionNode(std::move(name)), pressure_(pressure), temperature_(temperature) {
    requirePositive(pressure, "junction pressure");
    requirePositive(temperature, "junction temperature");
}
double Junction::pressure() const noexcept { return pressure_; }
double Junction::temperature() const noexcept { return temperature_; }
void Junction::setAlgebraicPressure(double pressure) noexcept {
    pressure_ = std::max(1.0, pressure);
}

CombustionTableData::CombustionTableData(
    std::vector<double> chamber_pressure, std::vector<double> mixture_ratio,
    std::vector<double> chamber_temperature, std::vector<double> cstar,
    std::vector<double> cf_vac, std::vector<double> isp_vac,
    double expansion_ratio, std::string oxidizer, std::string fuel,
    bool frozen_at_throat)
    : pc_(std::move(chamber_pressure)), mr_(std::move(mixture_ratio)),
      chamber_temperature_(std::move(chamber_temperature)),
      cstar_(std::move(cstar)), cf_vac_(std::move(cf_vac)),
      isp_vac_(std::move(isp_vac)),
      expansion_ratio_(requirePositive(expansion_ratio, "expansion ratio")),
      oxidizer_(std::move(oxidizer)), fuel_(std::move(fuel)),
      frozen_at_throat_(frozen_at_throat) {
    validateAxis(pc_, "chamber pressure");
    validateAxis(mr_, "mixture ratio");
    const std::size_t count = pc_.size() * mr_.size();
    if (chamber_temperature_.size() != count || cstar_.size() != count ||
        cf_vac_.size() != count || isp_vac_.size() != count) {
        throw std::invalid_argument(
            "combustion fields must match pressure x mixture-ratio axes");
    }
    for (const std::vector<double>* field :
         {&chamber_temperature_, &cstar_, &cf_vac_, &isp_vac_}) {
        if (std::any_of(field->begin(), field->end(),
                        [](double v) { return !std::isfinite(v); })) {
            throw std::invalid_argument("combustion fields must be finite");
        }
    }
}
double CombustionTableData::interpolate(const std::vector<double>& values,
                                         double pc, double mr) const {
    const std::size_t i = lowerCell(pc_, pc);
    const std::size_t j = lowerCell(mr_, mr);
    const double x = (pc - pc_[i]) / (pc_[i + 1U] - pc_[i]);
    const double y = (mr - mr_[j]) / (mr_[j + 1U] - mr_[j]);
    const auto at = [&](std::size_t ii, std::size_t jj) {
        return values[ii + pc_.size() * jj];
    };
    return (1.0 - x) * (1.0 - y) * at(i, j) +
           x * (1.0 - y) * at(i + 1U, j) +
           (1.0 - x) * y * at(i, j + 1U) +
           x * y * at(i + 1U, j + 1U);
}
double CombustionTableData::chamberTemperature(double pc, double mr) const {
    return interpolate(chamber_temperature_, pc, mr);
}
double CombustionTableData::cstar(double pc, double mr) const {
    return interpolate(cstar_, pc, mr);
}
double CombustionTableData::cfVac(double pc, double mr) const {
    return interpolate(cf_vac_, pc, mr);
}
double CombustionTableData::ispVac(double pc, double mr) const {
    return interpolate(isp_vac_, pc, mr);
}
double CombustionTableData::expansionRatio() const noexcept { return expansion_ratio_; }
const std::string& CombustionTableData::oxidizer() const noexcept { return oxidizer_; }
const std::string& CombustionTableData::fuel() const noexcept { return fuel_; }
bool CombustionTableData::frozenAtThroat() const noexcept { return frozen_at_throat_; }
const std::vector<double>& CombustionTableData::pressureAxis() const noexcept { return pc_; }
const std::vector<double>& CombustionTableData::mixtureRatioAxis() const noexcept { return mr_; }

Combustor::Combustor(std::string name,
                     std::shared_ptr<CombustionTableData> table,
                     double initial_pressure, double initial_mixture_ratio)
    : PropulsionNode(std::move(name)), table_(std::move(table)),
      pressure_(requirePositive(initial_pressure, "chamber pressure")),
      mixture_ratio_(requirePositive(initial_mixture_ratio, "mixture ratio")) {
    if (!table_) throw std::invalid_argument("combustor table must not be null");
}
double Combustor::pressure() const noexcept { return pressure_; }
double Combustor::temperature() const noexcept { return temperature_; }
double Combustor::mixtureRatio() const noexcept { return mixture_ratio_; }
double Combustor::oxidizerMassFlow() const noexcept { return oxidizer_flow_; }
double Combustor::fuelMassFlow() const noexcept { return fuel_flow_; }
double Combustor::totalMassFlow() const noexcept {
    return oxidizer_flow_ + fuel_flow_;
}
double Combustor::throatAreaForClosure() const noexcept { return 0.0; }
void Combustor::cacheChamber(double pressure, double oxidizer_flow,
                            double fuel_flow, double temperature) noexcept {
    pressure_ = std::max(1.0, pressure);
    oxidizer_flow_ = std::max(0.0, oxidizer_flow);
    fuel_flow_ = std::max(0.0, fuel_flow);
    if (fuel_flow_ > 0.0) mixture_ratio_ = oxidizer_flow_ / fuel_flow_;
    temperature_ = std::max(1.0, temperature);
}

Engine::Engine(std::string name, std::shared_ptr<CombustionTableData> table,
               double throat_area, double expansion_ratio,
               Vector3d thrust_direction, Vector3d application_point,
               std::shared_ptr<Ambient> ambient, double initial_pressure,
               double initial_mixture_ratio)
    : Combustor(name, std::move(table), initial_pressure, initial_mixture_ratio),
      ForceTorque(name),
      throat_area_(requirePositive(throat_area, "throat area")),
      expansion_ratio_(requirePositive(expansion_ratio, "expansion ratio")),
      direction_(std::move(thrust_direction)),
      application_point_(std::move(application_point)),
      ambient_(std::move(ambient)) {
    if (direction_.norm() <= 0.0) {
        throw std::invalid_argument("thrust direction must be nonzero");
    }
    direction_.normalize();
    if (std::abs(expansion_ratio_ - table_->expansionRatio()) >
        1e-9 * std::max(1.0, expansion_ratio_)) {
        throw std::invalid_argument(
            "engine expansion ratio does not match combustion table metadata");
    }
    addComputedStateVariable([this]() -> StateValue { return pressure(); }, "chamberPressure", "Pa");
    addComputedStateVariable([this]() -> StateValue { return mixtureRatio(); }, "mixtureRatio");
    addComputedStateVariable([this]() -> StateValue { return oxidizerMassFlow(); }, "oxidizerMassFlow", "kg/s");
    addComputedStateVariable([this]() -> StateValue { return fuelMassFlow(); }, "fuelMassFlow", "kg/s");
    addComputedStateVariable(
        [this]() -> StateValue { return thrust(); }, "thrust", "N");
}
void Engine::update() {
    setLoad(cached_thrust_ * direction_, Vector3d::Zero(), application_point_);
}
double Engine::throatArea() const noexcept { return throat_area_; }
double Engine::expansionRatio() const noexcept { return expansion_ratio_; }
double Engine::thrust() const noexcept {
    return enabled() ? cached_thrust_ : 0.0;
}
double Engine::ambientPressure() const noexcept {
    return ambient_ ? ambient_->pressure() : 101325.0;
}
double Engine::cstar() const noexcept { return cached_cstar_; }
double Engine::thrustCoefficient() const noexcept { return cached_cf_; }
double Engine::throatAreaForClosure() const noexcept { return throat_area_; }
void Engine::cachePerformance(double thrust, double cstar, double cf) noexcept {
    cached_thrust_ = std::max(0.0, thrust);
    cached_cstar_ = cstar;
    cached_cf_ = cf;
}

Orifice::Orifice(double cd_area)
    : cd_area_(requirePositive(cd_area, "CdA")) {}
Orifice::Orifice(double area, double discharge_coefficient)
    : Orifice(requirePositive(area, "orifice area") *
              requirePositive(discharge_coefficient, "discharge coefficient")) {}
double Orifice::cdArea() const noexcept { return cd_area_; }
double Orifice::massFlow(const PropulsionNode& upstream,
                         const PropulsionNode& downstream,
                         FluidPhase phase) const {
    const double dp = upstream.pressure() - downstream.pressure();
    if (std::abs(dp) < 1e-12) return 0.0;
    const PropulsionNode& high = dp > 0.0 ? upstream : downstream;
    const PropulsionNode& low = dp > 0.0 ? downstream : upstream;
    if (!high.phaseAvailable(phase)) return 0.0;
    if (phase == FluidPhase::Gas) {
        return signOf(dp) * gasRestrictionFlow(high, low, cd_area_);
    }
    const double rho = high.density(phase);
    if (!(rho > 0.0) || !std::isfinite(rho)) {
        throw std::runtime_error("invalid liquid density in orifice flow");
    }
    return signOf(dp) * cd_area_ * std::sqrt(2.0 * rho * std::abs(dp));
}

Valve::Valve(double cd_area, double opening)
    : Orifice(cd_area), opening_(0.0) { setOpening(opening); }
Valve::Valve(double area, double discharge_coefficient, double opening)
    : Orifice(area, discharge_coefficient), opening_(0.0) { setOpening(opening); }
double Valve::opening() const noexcept { return opening_; }
void Valve::setOpening(double opening) {
    if (opening < 0.0 || opening > 1.0 || !std::isfinite(opening)) {
        throw std::invalid_argument("valve opening must be in [0, 1]");
    }
    opening_ = opening;
}
double Valve::massFlow(const PropulsionNode& upstream,
                       const PropulsionNode& downstream,
                       FluidPhase phase) const {
    return opening_ * Orifice::massFlow(upstream, downstream, phase);
}

Line::Line(double length, double diameter, double roughness,
           double dynamic_viscosity)
    : length_(requirePositive(length, "line length")),
      diameter_(requirePositive(diameter, "line diameter")),
      roughness_(roughness),
      viscosity_(requirePositive(dynamic_viscosity, "dynamic viscosity")) {
    if (roughness < 0.0 || !std::isfinite(roughness)) {
        throw std::invalid_argument("line roughness must be nonnegative");
    }
}
double Line::massFlow(const PropulsionNode& upstream,
                      const PropulsionNode& downstream,
                      FluidPhase phase) const {
    const double dp = upstream.pressure() - downstream.pressure();
    if (std::abs(dp) < 1e-12) return 0.0;
    const PropulsionNode& source = dp > 0.0 ? upstream : downstream;
    if (!source.phaseAvailable(phase)) return 0.0;
    double rho = source.density(phase);
    if (phase == FluidPhase::Gas) {
        rho = 0.5 * (upstream.density(phase) + downstream.density(phase));
    }
    const double area = 0.25 * pi * diameter_ * diameter_;
    double friction = 0.02;
    double mdot = 0.0;
    for (int iteration = 0; iteration < 30; ++iteration) {
        const double velocity = std::sqrt(
            2.0 * std::abs(dp) /
            std::max(1e-30, rho * friction * length_ / diameter_));
        mdot = rho * area * velocity;
        const double reynolds = rho * velocity * diameter_ / viscosity_;
        const double laminar = 64.0 / std::max(1.0, reynolds);
        const double turbulent = 0.25 / std::pow(std::log10(
            roughness_ / (3.7 * diameter_) +
            5.74 / std::pow(std::max(1.0, reynolds), 0.9)), 2.0);
        // Smoothstep avoids a discontinuous residual around transition.
        const double blend = std::clamp((reynolds - 2000.0) / 2000.0, 0.0, 1.0);
        const double smooth = blend * blend * (3.0 - 2.0 * blend);
        const double next = (1.0 - smooth) * laminar + smooth * turbulent;
        if (std::abs(next - friction) < 1e-10) break;
        friction = 0.5 * (friction + next);
    }
    return signOf(dp) * mdot;
}

CompoundEdge::CompoundEdge(std::vector<std::shared_ptr<FlowComponent>> parts)
    : parts_(std::move(parts)) {
    if (parts_.empty() || std::any_of(parts_.begin(), parts_.end(),
                                      [](const auto& p) { return !p; })) {
        throw std::invalid_argument("compound edge needs non-null components");
    }
}
const std::vector<double>& CompoundEdge::intermediatePressures() const noexcept {
    return intermediate_pressures_;
}
bool CompoundEdge::containsUnimplementedComponent() const noexcept {
    for (const auto& part : parts_) {
        if (dynamic_cast<const Pump*>(part.get()) ||
            dynamic_cast<const Turbine*>(part.get())) {
            return true;
        }
        const auto* compound = dynamic_cast<const CompoundEdge*>(part.get());
        if (compound && compound->containsUnimplementedComponent()) return true;
    }
    return false;
}
double CompoundEdge::massFlow(const PropulsionNode& upstream,
                              const PropulsionNode& downstream,
                              FluidPhase phase) const {
    class Proxy final : public PropulsionNode {
    public:
        Proxy(double p, const PropulsionNode& source)
            : PropulsionNode("compound-intermediate"), p_(p), source_(source) {}
        double pressure() const noexcept override { return p_; }
        double temperature() const noexcept override { return source_.temperature(); }
        double density(FluidPhase p) const noexcept override { return source_.density(p); }
        double heatCapacityRatio() const noexcept override {
            return source_.heatCapacityRatio();
        }
        double molarMass() const noexcept override { return source_.molarMass(); }
    private:
        double p_;
        const PropulsionNode& source_;
    };
    intermediate_pressures_.clear();
    double common = std::numeric_limits<double>::infinity();
    const double step = (downstream.pressure() - upstream.pressure()) /
                        static_cast<double>(parts_.size());
    for (std::size_t i = 0; i < parts_.size(); ++i) {
        const double pa = upstream.pressure() + step * i;
        const double pb = upstream.pressure() + step * (i + 1U);
        Proxy a(pa, upstream), b(pb, upstream);
        const double candidate = parts_[i]->massFlow(a, b, phase);
        common = std::min(common, std::abs(candidate));
        if (i + 1U < parts_.size()) intermediate_pressures_.push_back(pb);
    }
    if (!std::isfinite(common)) return 0.0;
    return signOf(upstream.pressure() - downstream.pressure()) * common;
}

double Pump::massFlow(const PropulsionNode&, const PropulsionNode&,
                      FluidPhase) const {
    throw std::logic_error("Pump is a configuration skeleton and is unimplemented");
}
double Turbine::massFlow(const PropulsionNode&, const PropulsionNode&,
                         FluidPhase) const {
    throw std::logic_error("Turbine is a configuration skeleton and is unimplemented");
}

Propulsion::Propulsion(std::string name, UpdateRate update_rate)
    : Model(std::move(name), update_rate) {}
void Propulsion::requireMutable() const {
    if (finalized_) {
        throw std::logic_error("propulsion topology is immutable after finalize()");
    }
}
void Propulsion::ownNode(const std::shared_ptr<PropulsionNode>& node) {
    if (!node) throw std::invalid_argument("propulsion node must not be null");
    if (std::none_of(nodes_.begin(), nodes_.end(),
                     [&node](const auto& n) { return n.get() == node.get(); })) {
        nodes_.push_back(node);
        const std::string prefix =
            Model::makeStateIdentifier(node->name()) + ".";
        addComputedStateVariable([node]() -> StateValue { return node->pressure(); }, prefix + "pressure", "Pa");
        addComputedStateVariable([node]() -> StateValue { return node->temperature(); }, prefix + "temperature", "K");
        if (const auto tank = std::dynamic_pointer_cast<Tank>(node)) {
            addComputedStateVariable([tank]() -> StateValue { return tank->propellantMass(); }, prefix + "propellantMass", "kg");
            addComputedStateVariable([tank]() -> StateValue { return tank->pressurantMass(); }, prefix + "pressurantMass", "kg");
        }
    }
}
void Propulsion::connect(std::shared_ptr<PropulsionNode> upstream,
                         std::shared_ptr<PropulsionNode> downstream,
                         std::shared_ptr<FlowComponent> edge,
                         FluidPhase upstream_phase,
                         InletRole downstream_role) {
    requireMutable();
    if (!upstream || !downstream || !edge) {
        throw std::invalid_argument("connect arguments must not be null");
    }
    if (upstream.get() == downstream.get()) {
        throw std::invalid_argument("self-connected propulsion edge is invalid");
    }
    if (dynamic_cast<Combustor*>(downstream.get()) == nullptr &&
        downstream_role != InletRole::None) {
        throw std::invalid_argument("inlet roles apply only to combustors");
    }
    const auto* upstream_tank = dynamic_cast<const Tank*>(upstream.get());
    const auto* downstream_tank = dynamic_cast<const Tank*>(downstream.get());
    if (upstream_tank && downstream_tank &&
        upstream_tank->fluidSpecies(upstream_phase) !=
            downstream_tank->fluidSpecies(upstream_phase)) {
        throw std::invalid_argument(
            "incompatible " + upstream_tank->fluidSpecies(upstream_phase) +
            " -> " + downstream_tank->fluidSpecies(upstream_phase) +
            " tank connection");
    }
    ownNode(upstream);
    ownNode(downstream);
    edges_.push_back(edge);
    connections_.push_back({std::move(upstream), std::move(downstream),
                            std::move(edge), upstream_phase,
                            downstream_role, 0.0});
}
void Propulsion::connect(std::shared_ptr<PropulsionNode> upstream,
                         std::shared_ptr<PropulsionNode> downstream,
                         std::vector<std::shared_ptr<FlowComponent>> edges,
                         FluidPhase upstream_phase,
                         InletRole downstream_role) {
    connect(std::move(upstream), std::move(downstream),
            std::make_shared<CompoundEdge>(std::move(edges)),
            upstream_phase, downstream_role);
}
void Propulsion::vent(std::shared_ptr<PropulsionNode> upstream,
                      std::shared_ptr<Ambient> ambient,
                      std::shared_ptr<FlowComponent> edge,
                      FluidPhase upstream_phase) {
    connect(std::move(upstream), std::move(ambient), std::move(edge),
            upstream_phase, InletRole::None);
}
void Propulsion::addThermalContact(ThermalContact contact) {
    requireMutable();
    if (!contact.node || contact.area < 0.0 || !std::isfinite(contact.area) ||
        !std::isfinite(contact.heat_flux)) {
        throw std::invalid_argument("invalid thermal-contact configuration");
    }
    ownNode(contact.node);
    thermal_contacts_.push_back(std::move(contact));
}
void Propulsion::finalize() {
    requireMutable();
    if (connections_.empty()) {
        throw std::logic_error("cannot finalize an empty propulsion network");
    }
    for (const auto& edge : edges_) {
        const auto* compound = dynamic_cast<CompoundEdge*>(edge.get());
        if (dynamic_cast<Pump*>(edge.get()) || dynamic_cast<Turbine*>(edge.get()) ||
            (compound && compound->containsUnimplementedComponent())) {
            throw std::logic_error(
                "Pump/Turbine are unimplemented and cannot be finalized");
        }
    }
    std::unordered_set<const PropulsionNode*> connected;
    for (const Connection& connection : connections_) {
        connected.insert(connection.upstream.get());
        connected.insert(connection.downstream.get());
        if (dynamic_cast<Combustor*>(connection.downstream.get()) &&
            connection.downstream_role == InletRole::None) {
            throw std::logic_error(
                "combustor connections require OXIDIZER or FUEL inlet role");
        }
    }
    for (const auto& node : nodes_) {
        if (!connected.count(node.get())) {
            throw std::logic_error("unconnected propulsion node: " + node->name());
        }
    }
    bool inactive_thermal = false;
    for (const ThermalContact& contact : thermal_contacts_) {
        inactive_thermal = inactive_thermal ||
            (contact.area > 0.0 && contact.heat_flux != 0.0);
    }
    if (inactive_thermal) {
        std::cerr << "rsim warning: configured propulsion ThermalContact is "
                     "stored but thermal coupling is not active\n";
    }
    finalized_ = true;
}
bool Propulsion::finalized() const noexcept { return finalized_; }

void Propulsion::solveAlgebraicFlows() {
    std::vector<PropulsionNode*> algebraic_nodes;
    for (const auto& node : nodes_) {
        if (dynamic_cast<Junction*>(node.get()) ||
            dynamic_cast<Combustor*>(node.get())) {
            algebraic_nodes.push_back(node.get());
        }
    }

    Eigen::VectorXd state(static_cast<Eigen::Index>(algebraic_nodes.size()));
    for (std::size_t i = 0; i < algebraic_nodes.size(); ++i) {
        state[static_cast<Eigen::Index>(i)] = algebraic_nodes[i]->pressure();
    }

    auto residual = [&](const Eigen::VectorXd& candidate,
                        bool publish) -> Eigen::VectorXd {
        for (std::size_t i = 0; i < algebraic_nodes.size(); ++i) {
            const double pressure = std::max(1.0,
                candidate[static_cast<Eigen::Index>(i)]);
            if (auto* junction = dynamic_cast<Junction*>(algebraic_nodes[i])) {
                junction->setAlgebraicPressure(pressure);
            } else if (auto* combustor =
                       dynamic_cast<Combustor*>(algebraic_nodes[i])) {
                combustor->cacheChamber(
                    pressure, combustor->oxidizerMassFlow(),
                    combustor->fuelMassFlow(), combustor->temperature());
            }
        }
        for (Connection& connection : connections_) {
            connection.mass_flow = connection.edge->massFlow(
                *connection.upstream, *connection.downstream,
                connection.upstream_phase);
            if (!std::isfinite(connection.mass_flow)) {
                throw std::runtime_error("non-finite flow on edge into " +
                                         connection.downstream->name());
            }
        }

        Eigen::VectorXd result = Eigen::VectorXd::Zero(candidate.size());
        for (std::size_t i = 0; i < algebraic_nodes.size(); ++i) {
            PropulsionNode* algebraic = algebraic_nodes[i];
            double net_inflow = 0.0;
            double oxidizer = 0.0;
            double fuel = 0.0;
            for (const Connection& connection : connections_) {
                if (connection.downstream.get() == algebraic) {
                    net_inflow += connection.mass_flow;
                    const double flow = std::max(0.0, connection.mass_flow);
                    if (connection.downstream_role == InletRole::Oxidizer) {
                        oxidizer += flow;
                    }
                    if (connection.downstream_role == InletRole::Fuel) {
                        fuel += flow;
                    }
                }
                if (connection.upstream.get() == algebraic) {
                    net_inflow -= connection.mass_flow;
                }
            }
            auto* combustor = dynamic_cast<Combustor*>(algebraic);
            if (combustor && combustor->throatAreaForClosure() > 0.0) {
                const double mr = fuel > 1e-14 ? oxidizer / fuel
                                               : combustor->mixtureRatio();
                const double pc_query = std::clamp(
                    combustor->pressure(), combustor->table_->pressureAxis().front(),
                    combustor->table_->pressureAxis().back());
                const double mr_query = std::clamp(
                    mr, combustor->table_->mixtureRatioAxis().front(),
                    combustor->table_->mixtureRatioAxis().back());
                const double nozzle_flow = combustor->pressure() *
                    combustor->throatAreaForClosure() /
                    combustor->table_->cstar(pc_query, mr_query);
                result[static_cast<Eigen::Index>(i)] = net_inflow - nozzle_flow;
                if (publish) {
                    combustor->cacheChamber(
                        combustor->pressure(), oxidizer, fuel,
                        combustor->table_->chamberTemperature(pc_query, mr_query));
                }
            } else {
                result[static_cast<Eigen::Index>(i)] = net_inflow;
            }
        }
        return result;
    };

    if (state.size() == 0) {
        static_cast<void>(residual(state, true));
        return;
    }

    Eigen::VectorXd values = residual(state, false);
    for (int iteration = 0; iteration < 60; ++iteration) {
        if (values.lpNorm<Eigen::Infinity>() < 1e-8) {
            static_cast<void>(residual(state, true));
            return;
        }
        Eigen::MatrixXd jacobian(values.size(), state.size());
        for (Eigen::Index column = 0; column < state.size(); ++column) {
            Eigen::VectorXd perturbed = state;
            const double delta = std::max(1.0, std::abs(state[column]) * 1e-6);
            perturbed[column] += delta;
            jacobian.col(column) = (residual(perturbed, false) - values) / delta;
        }
        const Eigen::VectorXd correction =
            jacobian.colPivHouseholderQr().solve(-values);
        if (!correction.allFinite()) {
            throw std::runtime_error(
                "propulsion numeric-difference dense QR produced a non-finite step");
        }
        double scale = 1.0;
        bool accepted = false;
        for (int line_search = 0; line_search < 16; ++line_search) {
            Eigen::VectorXd trial = state + scale * correction;
            for (Eigen::Index i = 0; i < trial.size(); ++i) {
                trial[i] = std::max(1.0, trial[i]);
            }
            Eigen::VectorXd trial_values = residual(trial, false);
            if (trial_values.squaredNorm() < values.squaredNorm()) {
                state = std::move(trial);
                values = std::move(trial_values);
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if (!accepted) break;
    }
    std::ostringstream message;
    message << "propulsion monolithic algebraic solve failed: residual infinity "
            << "norm=" << values.lpNorm<Eigen::Infinity>()
            << ", unknown pressures=" << state.size();
    throw std::runtime_error(message.str());
}

void Propulsion::update() {
    if (!finalized_) {
        throw std::logic_error("Propulsion::finalize() must precede update()");
    }
    solveAlgebraicFlows();

    const double dt = timeStep();
    if (!(dt >= 0.0) || !std::isfinite(dt)) {
        throw std::runtime_error("propulsion received invalid scheduler time step");
    }

    // Limit every group of tank outflows to the inventory available this step.
    // This lets the final partial flow empty a tank exactly instead of asking the
    // tank to provide a negative mass and aborting the simulation.
    std::vector<double> flow_scales(connections_.size(), 1.0);
    for (std::size_t index = 0; index < connections_.size(); ++index) {
        const Connection& connection = connections_[index];
        PropulsionNode* source = connection.mass_flow >= 0.0
            ? connection.upstream.get() : connection.downstream.get();
        auto* source_tank = dynamic_cast<Tank*>(source);
        if (!source_tank || connection.mass_flow == 0.0 || dt == 0.0) continue;

        double requested_mass = 0.0;
        for (const Connection& other : connections_) {
            PropulsionNode* other_source = other.mass_flow >= 0.0
                ? other.upstream.get() : other.downstream.get();
            if (other_source == source_tank &&
                other.upstream_phase == connection.upstream_phase) {
                requested_mass += std::abs(other.mass_flow) * dt;
            }
        }
        const double available_mass = connection.upstream_phase == FluidPhase::Liquid
            ? source_tank->propellantMass() : source_tank->pressurantMass();
        if (requested_mass > available_mass && requested_mass > 0.0) {
            flow_scales[index] = available_mass / requested_mass;
        }
    }
    for (std::size_t index = 0; index < connections_.size(); ++index) {
        connections_[index].mass_flow *= flow_scales[index];
    }

    // Keep all inlets to one combustor at the same limiting fraction. This
    // preserves mixture ratio during the final partial burn instead of dumping
    // the still-plentiful propellant through an extinguished chamber.
    for (const auto& node : nodes_) {
        auto* combustor = dynamic_cast<Combustor*>(node.get());
        if (!combustor) continue;
        double limiting_scale = 1.0;
        for (std::size_t index = 0; index < connections_.size(); ++index) {
            if (connections_[index].downstream.get() == combustor &&
                connections_[index].mass_flow > 0.0) {
                limiting_scale = std::min(limiting_scale, flow_scales[index]);
            }
        }
        for (std::size_t index = 0; index < connections_.size(); ++index) {
            if (connections_[index].downstream.get() == combustor &&
                connections_[index].mass_flow > 0.0 &&
                flow_scales[index] > limiting_scale) {
                connections_[index].mass_flow *=
                    limiting_scale / flow_scales[index];
                flow_scales[index] = limiting_scale;
            }
        }
    }

    // Publish the inventory-limited inlet flows. A bipropellant combustor
    // flames out when either propellant stream is absent; its chamber pressure
    // and flows then fall to zero-equivalent values and an Engine produces no thrust.
    for (const auto& node : nodes_) {
        auto* combustor = dynamic_cast<Combustor*>(node.get());
        if (!combustor) continue;
        double oxidizer = 0.0;
        double fuel = 0.0;
        bool inventory_limited = false;
        for (std::size_t index = 0; index < connections_.size(); ++index) {
            const Connection& connection = connections_[index];
            if (connection.downstream.get() != combustor) continue;
            inventory_limited = inventory_limited || flow_scales[index] < 1.0;
            const double inlet_flow = std::max(0.0, connection.mass_flow);
            if (connection.downstream_role == InletRole::Oxidizer) {
                oxidizer += inlet_flow;
            } else if (connection.downstream_role == InletRole::Fuel) {
                fuel += inlet_flow;
            }
        }
        if (inventory_limited || oxidizer <= 1e-14 || fuel <= 1e-14) {
            combustor->cacheChamber(1.0, 0.0, 0.0, 1.0);
        } else {
            combustor->cacheChamber(combustor->pressure(), oxidizer, fuel,
                                    combustor->temperature());
        }
    }

    struct TankDelta {
        Tank* tank;
        FluidPhase phase;
        double mass = 0.0;
        double energy = 0.0;
    };
    std::vector<TankDelta> deltas;
    const auto add_delta = [&](Tank* tank, FluidPhase phase,
                               double mass, double energy) {
        if (!tank) return;
        auto found = std::find_if(deltas.begin(), deltas.end(),
            [tank, phase](const TankDelta& delta) {
                return delta.tank == tank && delta.phase == phase;
            });
        if (found == deltas.end()) {
            deltas.push_back({tank, phase, mass, energy});
        } else {
            found->mass += mass;
            found->energy += energy;
        }
    };
    for (Connection& connection : connections_) {
        if (connection.mass_flow == 0.0) continue;
        PropulsionNode* source = connection.mass_flow > 0.0
            ? connection.upstream.get() : connection.downstream.get();
        PropulsionNode* destination = connection.mass_flow > 0.0
            ? connection.downstream.get() : connection.upstream.get();
        auto* source_tank = dynamic_cast<Tank*>(source);
        auto* destination_tank = dynamic_cast<Tank*>(destination);
        if (destination_tank && !source_tank) {
            throw std::runtime_error(
                "cannot determine fluid species entering tank " +
                destination_tank->name() + " from " + source->name());
        }
        if (source_tank && destination_tank &&
            source_tank->fluidSpecies(connection.upstream_phase) !=
                destination_tank->fluidSpecies(connection.upstream_phase)) {
            throw std::runtime_error(
                "incompatible tank species during flow from " + source->name() +
                " to " + destination->name());
        }
        const double transported_mass = std::abs(connection.mass_flow) * dt;
        const double source_enthalpy = source->specificEnthalpy(
            connection.upstream_phase);
        const double transported_energy = transported_mass * source_enthalpy;
        add_delta(source_tank, connection.upstream_phase,
                  -transported_mass, -transported_energy);
        add_delta(destination_tank, connection.upstream_phase,
                  transported_mass, transported_energy);
    }
    for (const TankDelta& delta : deltas) {
        delta.tank->validateTransfer(delta.phase, delta.mass, delta.energy);
    }
    for (const TankDelta& delta : deltas) {
        delta.tank->applyTransfer(delta.phase, delta.mass, delta.energy);
    }

    for (const auto& node : nodes_) {
        auto* engine = dynamic_cast<Engine*>(node.get());
        if (!engine) continue;
        if (engine->totalMassFlow() <= 0.0) {
            engine->cachePerformance(0.0, 0.0, 0.0);
            continue;
        }
        try {
            const double cstar = engine->table_->cstar(
                engine->pressure(), engine->mixtureRatio());
            const double cf_vac = engine->table_->cfVac(
                engine->pressure(), engine->mixtureRatio());
            const double cf_ambient = cf_vac - engine->expansionRatio() *
                engine->ambientPressure() / engine->pressure();
            engine->cachePerformance(
                std::max(0.0, cf_ambient * engine->pressure() *
                                  engine->throatArea()),
                cstar, cf_ambient);
        } catch (const std::out_of_range&) {
            throw std::runtime_error(
                "engine state is outside its combustion table: " +
                static_cast<const PropulsionNode&>(*engine).name());
        }
    }
    std::ostringstream summary;
    summary << "converged native network solve: " << connections_.size()
            << " flows, stored mass=" << totalStoredMass() << " kg";
    last_solve_summary_ = summary.str();
}

std::size_t Propulsion::connectionCount() const noexcept {
    return connections_.size();
}
double Propulsion::connectionMassFlow(std::size_t index) const {
    return connections_.at(index).mass_flow;
}
double Propulsion::totalStoredMass() const noexcept {
    double mass = 0.0;
    for (const auto& node : nodes_) {
        if (const auto* tank = dynamic_cast<const Tank*>(node.get())) {
            mass += tank->variableMass();
        }
    }
    return mass;
}
const std::string& Propulsion::lastSolveSummary() const noexcept {
    return last_solve_summary_;
}

}  // namespace rsim
