"""Run a short native liquid sounding-rocket simulation.

The bundled combustion table is deliberately small and illustrative.  It is
useful for a deterministic software example, but not for engine design.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import rsim


HERE = Path(__file__).resolve().parent
TABLE_PATH = HERE / "data" / "n2o_ethanol_illustrative.csv"
EXPANSION_RATIO = 7.0


def load_performance_table(path: Path = TABLE_PATH) -> rsim.CombustionTableData:
    """Load the CSV in Python, then copy its arrays into native C++ storage."""

    rows = np.genfromtxt(path, delimiter=",", names=True, comments="#")
    pc = np.unique(rows["chamber_pressure_pa"])
    mr = np.unique(rows["mixture_ratio"])
    shape = (pc.size, mr.size)
    fields = {
        name: np.empty(shape, dtype=np.float64, order="F")
        for name in (
            "chamber_temperature_k",
            "cstar_m_per_s",
            "cf_vac",
            "isp_vac_s",
        )
    }
    for row in rows:
        i = int(np.searchsorted(pc, row["chamber_pressure_pa"]))
        j = int(np.searchsorted(mr, row["mixture_ratio"]))
        for name, values in fields.items():
            values[i, j] = row[name]

    # .native constructs CombustionTableData and immediately copies every
    # NumPy array. The simulator therefore makes no Python callbacks.
    return rsim.CombustionTable(
        chamber_pressure=pc,
        mixture_ratio=mr,
        chamber_temperature=fields["chamber_temperature_k"],
        cstar=fields["cstar_m_per_s"],
        cf_vac=fields["cf_vac"],
        isp_vac=fields["isp_vac_s"],
        expansion_ratio=EXPANSION_RATIO,
        oxidizer="N2O",
        fuel="Ethanol",
    ).native


def make_native_tank(
    definition: rsim.CylindricalTank,
    *,
    name: str,
    pressure: float,
    propellant: str,
) -> rsim.Tank:
    """Reuse one geometry tank's interior for native mass and flow state."""

    interior = definition.cpp_definition()
    return rsim.Tank(
        name=name,
        radius=interior["radius"],
        length=interior["length"],
        aft_position=interior["aft_position"],
        liquid_density=interior["propellant_density"],
        propellant_mass=interior["initial_propellant_mass"],
        pressure=pressure,
        temperature=293.15,
        collapse_factor=1.15,
        blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
        polytropic_exponent=1.25,
        propellant=propellant,
        pressurant="Helium",
    )


def build_vehicle():
    aluminum = rsim.Material("6061-T6 aluminum", 2700.0)
    equipment = rsim.Material("equivalent installed equipment", 700.0)
    geometry = rsim.LiquidSoundingRocketGeometry()
    geometry.add_component(rsim.CylindricalShell(0.16, 0.002, 2.40, 1.20, aluminum))
    geometry.add_component(rsim.ConicalShell(0.16, 0.55, 0.002, 2.40, aluminum))
    geometry.add_component(rsim.SolidCylinder(0.12, 0.25, 0.13, equipment))

    oxidizer_geometry = rsim.CylindricalTank(
        outer_radius=0.15,
        length=0.70,
        wall_thickness=0.004,
        center_x=0.85,
        material=aluminum,
        propellant_density=750.0,
        initial_fill_fraction=0.80,
    )
    fuel_geometry = rsim.CylindricalTank(
        outer_radius=0.15,
        length=0.28,
        wall_thickness=0.004,
        center_x=1.48,
        material=aluminum,
        propellant_density=789.0,
        initial_fill_fraction=0.90,
    )
    geometry.add_tank(oxidizer_geometry)
    geometry.add_tank(fuel_geometry)
    dry = geometry.dry_mass_properties()

    oxidizer_tank = make_native_tank(
        oxidizer_geometry, name="oxidizer-tank", pressure=4.5e6, propellant="N2O"
    )
    fuel_tank = make_native_tank(
        fuel_geometry, name="fuel-tank", pressure=4.5e6, propellant="Ethanol"
    )
    ambient = rsim.Ambient("12-km ambient", pressure=19_330.0, temperature=216.7)
    engine = rsim.Engine(
        "pressure-fed engine",
        load_performance_table(),
        throat_area=2.0e-3,
        expansion_ratio=EXPANSION_RATIO,
        thrust_direction=np.array([1.0, 0.0, 0.0]),
        application_point=np.array([0.05, 0.0, 0.0]),
        ambient=ambient,
        initial_pressure=2.0e6,
        initial_mixture_ratio=5.0,
    )

    propulsion = rsim.Propulsion("pressure-fed N2O/ethanol")
    propulsion.connect(
        oxidizer_tank,
        engine,
        [
            rsim.Line(0.55, 0.028, roughness=1.5e-5, dynamic_viscosity=3.0e-4),
            rsim.Valve(1.6e-4, opening=0.95),
            rsim.Orifice(7.7e-5),
        ],
        upstream_phase=rsim.FluidPhase.LIQUID,
        downstream_role=rsim.InletRole.OXIDIZER,
    )
    propulsion.connect(
        fuel_tank,
        engine,
        [
            rsim.Line(0.45, 0.016, roughness=1.5e-5, dynamic_viscosity=1.2e-3),
            rsim.Valve(3.2e-5, opening=0.95),
            rsim.Orifice(1.5e-5),
        ],
        upstream_phase=rsim.FluidPhase.LIQUID,
        downstream_role=rsim.InletRole.FUEL,
    )
    propulsion.finalize()

    world = rsim.Frame("inertial-world", rsim.InertialStatus.INERTIAL)
    vehicle = rsim.Vehicle(
        "liquid-sounding-rocket",
        world,
        dry.mass,
        dry.moment_of_inertia,
        center_of_gravity=dry.center_of_mass,
        models=[propulsion],
        force_torques=[engine],
    )
    vehicle.rigid_body.add_variable_mass_source(oxidizer_tank)
    vehicle.rigid_body.add_variable_mass_source(fuel_tank)
    vehicle.kinematics.set_state(
        position_in_parent=np.zeros(3),
        velocity_in_parent=np.zeros(3),
        parent_from_body=np.eye(3),
        angular_velocity_in_body=np.zeros(3),
    )
    return vehicle, propulsion, engine, oxidizer_tank, fuel_tank, dry


def main() -> None:
    vehicle, propulsion, engine, oxidizer, fuel, dry = build_vehicle()
    initial_oxidizer = oxidizer.propellant_mass
    initial_fuel = fuel.propellant_mass
    initial_total_mass = dry.mass + propulsion.total_stored_mass

    simulator = rsim.Simulator(0.01)
    simulator.add_vehicle(vehicle)
    simulator.run_steps(50)

    final_mass = vehicle.rigid_body.mass
    final_cg = vehicle.rigid_body.center_of_gravity
    final_position = vehicle.kinematics.position_in_parent
    final_velocity = vehicle.kinematics.velocity_in_parent

    assert oxidizer.propellant_mass < initial_oxidizer
    assert fuel.propellant_mass < initial_fuel
    assert 4.0 <= engine.mixture_ratio <= 7.0
    assert 4.0e5 <= engine.pressure <= 4.0e6
    assert engine.temperature > 2500.0
    assert engine.cstar > 1400.0
    assert engine.thrust > 1000.0
    assert final_mass < initial_total_mass
    assert final_velocity[0] > 0.0 and final_position[0] > 0.0
    assert "converged" in propulsion.last_solve_summary

    print("Illustrative table: N2O/ethanol values are not design-quality")
    print(f"Simulated time:       {simulator.time:.3f} s")
    print(f"Oxidizer depletion:   {initial_oxidizer - oxidizer.propellant_mass:.6f} kg")
    print(f"Fuel depletion:       {initial_fuel - fuel.propellant_mass:.6f} kg")
    print(f"Chamber pressure:     {engine.pressure / 1e6:.6f} MPa")
    print(f"Mixture ratio (O/F):  {engine.mixture_ratio:.6f}")
    print(f"Chamber temperature:  {engine.temperature:.3f} K")
    print(f"Characteristic vel.:  {engine.cstar:.3f} m/s")
    print(f"Thrust:               {engine.thrust:.3f} N")
    print(f"Total vehicle mass:   {final_mass:.6f} kg")
    print(f"Vehicle CG [m]:       {np.array2string(final_cg, precision=6)}")
    print(f"Final position [m]:   {np.array2string(final_position, precision=6)}")
    print(f"Final velocity [m/s]: {np.array2string(final_velocity, precision=6)}")


if __name__ == "__main__":
    main()
