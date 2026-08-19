# %%
"""Run a Vespula-inspired ascent plus targeted recovery-transition check.

Authoritative mockup inputs are separated from the explicit engineering
assumptions in examples/README.md. The combustion table and external geometry
are illustrative software inputs, not flight- or engine-design data.
"""

from __future__ import annotations

from math import pi, sqrt
from pathlib import Path
import tempfile

import numpy as np

import rsim


HERE = Path(__file__).resolve().parent
TABLE_PATH = HERE / "data" / "lox_rp1_illustrative.csv"
PSI_TO_PA = 6894.757293168
LITERS_TO_M3 = 1.0e-3

# Authoritative Vespula mockup inputs, converted to SI once here.
DESIGN_PC = 300.0 * PSI_TO_PA
INJECTOR_PRESSURE_DROP = 54.0 * PSI_TO_PA
TARGET_MIXTURE_RATIO = 2.0
TARGET_LOX_FLOW = 3.2
TARGET_RP1_FLOW = 1.6
LOX_INTERNAL_VOLUME = 75.0 * LITERS_TO_M3
INITIAL_FILL_FRACTION = 0.95
COPV_INTERNAL_VOLUME = 30.0 * LITERS_TO_M3
COPV_INITIAL_PRESSURE = 5400.0 * PSI_TO_PA

# Explicit example assumptions.
LOX_DENSITY = 1141.0
RP1_DENSITY = 810.0
DESIGN_BURN_DURATION = 18.0
POWERED_DURATION = 16.0
MAX_APOGEE_TIME = 150.0
EXPANSION_RATIO = 6.0
TIME_STEP = 0.01
SAMPLE_PERIOD = 0.10
PRESSURE_DEADBAND_HALF_WIDTH = 5.0 * PSI_TO_PA
PRESSURANT_VALVE_CAPACITY_FACTOR = 4.0
APOGEE_VELOCITY_DEADBAND = 0.1
MAIN_DEPLOY_ALTITUDE = 609.6
RECOVERY_DESIGN_MASS = 125.0
RECOVERY_DESIGN_DENSITY = 1.0
DROGUE_DRAG_COEFFICIENT = 0.8
MAIN_DRAG_COEFFICIENT = 1.5
DROGUE_TARGET_SPEED = 213.36
MAIN_TARGET_SPEED = 1.524
BODY_LENGTH = 7.00
BODY_DIAMETER = 0.254
EARTH_RADIUS = 6_378_137.0
SEA_LEVEL_DENSITY = 1.225
SEA_LEVEL_SPEED_OF_SOUND = 340.3
SEA_LEVEL_DYNAMIC_VISCOSITY = 1.789e-5
TABLE_PC_BOUNDS = (0.15e6, 5.0e6)
TABLE_MR_BOUNDS = (1.5, 2.5)

# Fuel volume follows directly from 18 s * 1.6 kg/s and 95% fill.
RP1_PROPELLANT_MASS = DESIGN_BURN_DURATION * TARGET_RP1_FLOW
RP1_INTERNAL_VOLUME = (
    RP1_PROPELLANT_MASS / (RP1_DENSITY * INITIAL_FILL_FRACTION)
)


def load_performance_table(path: Path = TABLE_PATH) -> rsim.CombustionTableData:
    """Load the CSV in Python and copy all arrays into native C++ storage."""

    rows = np.genfromtxt(path, delimiter=",", names=True)
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

    # .native immediately copies every setup-time NumPy array into C++.
    return rsim.CombustionTable(
        chamber_pressure=pc,
        mixture_ratio=mr,
        chamber_temperature=fields["chamber_temperature_k"],
        cstar=fields["cstar_m_per_s"],
        cf_vac=fields["cf_vac"],
        isp_vac=fields["isp_vac_s"],
        expansion_ratio=EXPANSION_RATIO,
        oxidizer="LOX",
        fuel="RP-1",
    ).native


def tank_outer_length(internal_volume: float) -> float:
    inner_radius = 0.120 - 0.006
    return internal_volume / (pi * inner_radius**2) + 2.0 * 0.006


def make_native_tank(
    definition: rsim.CylindricalTank,
    *,
    name: str,
    pressure: float,
    temperature: float,
    propellant: str,
) -> rsim.Tank:
    interior = definition.cpp_definition()
    return rsim.Tank(
        name=name,
        radius=interior["radius"],
        length=interior["length"],
        aft_position=interior["aft_position"],
        liquid_density=interior["propellant_density"],
        propellant_mass=interior["initial_propellant_mass"],
        pressure=pressure,
        temperature=temperature,
        collapse_factor=1.0,
        blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
        polytropic_exponent=1.0,
        propellant=propellant,
        pressurant="Nitrogen",
    )


def liquid_cd_area(target_flow: float, density: float) -> float:
    """Size effective CdA for the authoritative 54 psi injector drop."""

    return target_flow / sqrt(2.0 * density * INJECTOR_PRESSURE_DROP)


def size_line_diameter(
    tank: rsim.Tank,
    target_flow: float,
    *,
    length: float,
    roughness: float,
    viscosity: float,
) -> float:
    """Use the public Line law to size one separate 54 psi design loss."""

    downstream = rsim.Ambient(
        "line-sizing endpoint",
        tank.pressure - INJECTOR_PRESSURE_DROP,
        tank.temperature,
    )
    low, high = 1.0e-3, 0.10
    for _ in range(70):
        diameter = 0.5 * (low + high)
        trial = rsim.Line(length, diameter, roughness, viscosity)
        if trial.mass_flow(tank, downstream, rsim.FluidPhase.LIQUID) < target_flow:
            low = diameter
        else:
            high = diameter
    return 0.5 * (low + high)


def size_gas_orifice(
    source: rsim.Tank,
    destination: rsim.Tank,
    target_flow: float,
) -> float:
    unit_flow = rsim.Orifice(1.0).mass_flow(
        source, destination, rsim.FluidPhase.GAS
    )
    return target_flow / unit_flow


def make_aerodynamics() -> tuple[rsim.Aerodynamics, rsim.aero.RocketGeometry]:
    """Generate a Vespula-like table in Python, then load it into native C++."""

    geometry = rsim.aero.RocketGeometry(
        nose=rsim.aero.TangentOgiveNose(0.80, BODY_DIAMETER),
        body=rsim.aero.CylindricalBody(6.20, BODY_DIAMETER),
        fins=rsim.aero.TrapezoidalFinSet(
            count=4,
            root_chord=0.75,
            tip_chord=0.32,
            semi_span=0.28,
            tip_leading_edge_offset=0.25,
            thickness=0.012,
            root_leading_edge_x=6.15,
        ),
        moment_reference_x=4.40,
    )
    grid = rsim.aero.AeroGrid(
        mach=rsim.aero.GridAxis([0.05, 0.3, 0.8, 1.0, 1.3, 2.0, 3.0]),
        # Use the complete validated envelope of the implementation-guide methods.
        alpha_deg=rsim.aero.GridAxis([-10.0, -5.0, 0.0, 5.0, 10.0]),
        phi_aero_deg=rsim.aero.GridAxis([0.0, 45.0, 90.0]),
        reynolds=rsim.aero.GridAxis([1.0e5, 1.0e6, 1.0e7, 3.0e7]),
    )
    table = rsim.aero.generate_aero_table(
        geometry,
        grid,
        rsim.aero.AerodynamicOptions(
            body_station_count=51,
            body_azimuth_count=20,
        ),
    )
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "vespula_inspired_aero.mat"
        table.save_mat(path)
        database = rsim.AerodynamicDatabase(str(path))
    return rsim.Aerodynamics("Vespula-like aerodynamics", database), geometry


def build_vehicle():
    """Build the vehicle exclusively from public ``rsim`` Python APIs."""

    aluminum = rsim.Material("6061-T6 aluminum", 2700.0)
    installed_hardware = rsim.Material("installed hardware equivalent", 900.0)
    geometry = rsim.LiquidSoundingRocketGeometry()
    geometry.add_component(rsim.CylindricalShell(0.127, 0.002, 6.30, 3.15,
                                                  aluminum))
    geometry.add_component(rsim.ConicalShell(0.127, 0.70, 0.002, 6.30,
                                              aluminum))
    geometry.add_component(rsim.SolidCylinder(0.115, 0.65, 0.33,
                                               installed_hardware))
    geometry.add_component(rsim.SolidCylinder(0.105, 0.50, 5.70,
                                               installed_hardware))

    lox_geometry = rsim.CylindricalTank(
        0.120, tank_outer_length(LOX_INTERNAL_VOLUME), 0.006, 1.45,
        aluminum, LOX_DENSITY, INITIAL_FILL_FRACTION,
    )
    rp1_geometry = rsim.CylindricalTank(
        0.120, tank_outer_length(RP1_INTERNAL_VOLUME), 0.006, 2.75,
        aluminum, RP1_DENSITY, INITIAL_FILL_FRACTION,
    )
    geometry.add_tank(lox_geometry)
    geometry.add_tank(rp1_geometry)
    dry = geometry.dry_mass_properties()

    # CompoundEdge currently gives Line, Valve, and Orifice equal shares of
    # total pressure drop. Each is therefore sized independently at 54 psi,
    # making tank pressure Pc + 3*54 psi at the design point.
    tank_pressure = DESIGN_PC + 3.0 * INJECTOR_PRESSURE_DROP
    lox_tank = make_native_tank(
        lox_geometry, name="LOX tank", pressure=tank_pressure,
        temperature=90.0, propellant="LOX",
    )
    rp1_tank = make_native_tank(
        rp1_geometry, name="RP-1 tank", pressure=tank_pressure,
        temperature=293.15, propellant="RP-1",
    )
    copv_radius = 0.10
    copv = rsim.Tank(
        "shared 20 L N2 COPV",
        copv_radius,
        COPV_INTERNAL_VOLUME / (pi * copv_radius**2),
        3.65,
        1000.0,
        0.0,
        COPV_INITIAL_PRESSURE,
        temperature=293.15,
        collapse_factor=1.0,
        blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
        polytropic_exponent=1.0,
        propellant="none",
        pressurant="Nitrogen",
    )

    performance = load_performance_table()
    design_cstar = performance.cstar(DESIGN_PC, TARGET_MIXTURE_RATIO)
    throat_area = (
        (TARGET_LOX_FLOW + TARGET_RP1_FLOW) * design_cstar / DESIGN_PC
    )
    ambient = rsim.Ambient("sea-level reference", 101_325.0, 288.15)
    engine = rsim.Engine(
        "Vespula-class illustrative engine",
        performance,
        throat_area=throat_area,
        expansion_ratio=EXPANSION_RATIO,
        thrust_direction=np.array([1.0, 0.0, 0.0]),
        application_point=np.array([0.05, 0.0, 0.0]),
        ambient=ambient,
        initial_pressure=DESIGN_PC,
        initial_mixture_ratio=TARGET_MIXTURE_RATIO,
    )

    lox_line = size_line_diameter(
        lox_tank, TARGET_LOX_FLOW, length=1.15,
        roughness=1.5e-5, viscosity=1.9e-4,
    )
    fuel_line = size_line_diameter(
        rp1_tank, TARGET_RP1_FLOW, length=1.35,
        roughness=1.5e-5, viscosity=1.7e-3,
    )
    lox_cda = liquid_cd_area(TARGET_LOX_FLOW, LOX_DENSITY)
    fuel_cda = liquid_cd_area(TARGET_RP1_FLOW, RP1_DENSITY)

    propulsion = rsim.Propulsion("shared-N2-pressure-fed LOX/RP-1")
    propulsion.connect(
        lox_tank,
        engine,
        [
            rsim.Line(1.15, lox_line, 1.5e-5, 1.9e-4),
            rsim.Valve(lox_cda / 0.95, opening=0.95),
            rsim.Orifice(lox_cda),
        ],
        upstream_phase=rsim.FluidPhase.LIQUID,
        downstream_role=rsim.InletRole.OXIDIZER,
    )
    propulsion.connect(
        rp1_tank,
        engine,
        [
            rsim.Line(1.35, fuel_line, 1.5e-5, 1.7e-3),
            rsim.Valve(fuel_cda / 0.95, opening=0.95),
            rsim.Orifice(fuel_cda),
        ],
        upstream_phase=rsim.FluidPhase.LIQUID,
        downstream_role=rsim.InletRole.FUEL,
    )

    # Size initial gas replenishment for the liquid volume displaced at the
    # authoritative design flows. One native COPV branches into both ullages.
    lox_makeup_flow = lox_tank.density(rsim.FluidPhase.GAS) * (
        TARGET_LOX_FLOW / LOX_DENSITY
    )
    fuel_makeup_flow = rp1_tank.density(rsim.FluidPhase.GAS) * (
        TARGET_RP1_FLOW / RP1_DENSITY
    )
    lox_pressurant_valve = rsim.Valve(
        PRESSURANT_VALVE_CAPACITY_FACTOR *
        size_gas_orifice(copv, lox_tank, lox_makeup_flow),
        opening=0.0,
    )
    rp1_pressurant_valve = rsim.Valve(
        PRESSURANT_VALVE_CAPACITY_FACTOR *
        size_gas_orifice(copv, rp1_tank, fuel_makeup_flow),
        opening=0.0,
    )
    propulsion.connect(
        copv,
        lox_tank,
        lox_pressurant_valve,
        upstream_phase=rsim.FluidPhase.GAS,
    )
    propulsion.connect(
        copv,
        rp1_tank,
        rp1_pressurant_valve,
        upstream_phase=rsim.FluidPhase.GAS,
    )
    propulsion.finalize()

    lox_pressure_controller = rsim.BangBangPressureController(
        "LOX tank pressure control",
        lox_tank,
        lox_pressurant_valve,
        pressure_setpoint=tank_pressure,
        deadband_half_width=PRESSURE_DEADBAND_HALF_WIDTH,
        initially_open=False,
    )
    rp1_pressure_controller = rsim.BangBangPressureController(
        "RP-1 tank pressure control",
        rp1_tank,
        rp1_pressurant_valve,
        pressure_setpoint=tank_pressure,
        deadband_half_width=PRESSURE_DEADBAND_HALF_WIDTH,
        initially_open=False,
    )

    aerodynamics, aero_geometry = make_aerodynamics()
    gravity = rsim.Gravity.earth()
    flight_configuration = rsim.FlightStateConfiguration()
    flight_configuration.origin_in_integration_frame = np.array(
        [EARTH_RADIUS, 0.0, 0.0]
    )
    flight_configuration.up_direction_in_integration_frame = np.array(
        [1.0, 0.0, 0.0]
    )
    flight_configuration.rail_direction_in_integration_frame = np.array(
        [1.0, 0.0, 0.0]
    )
    flight_configuration.rail_length = 6.0
    flight_configuration.thrust_on_threshold = 1000.0
    flight_configuration.thrust_off_threshold = 500.0
    flight_configuration.apogee_velocity_deadband = APOGEE_VELOCITY_DEADBAND
    flight_configuration.drogue_stage_enabled = True
    flight_configuration.drogue_trigger.deploy_at_apogee = True
    flight_configuration.main_stage_enabled = True
    flight_configuration.main_trigger.maximum_altitude = MAIN_DEPLOY_ALTITUDE
    flight_configuration.main_trigger.require_descending_altitude_crossing = True
    flight_state = rsim.TruthState(
        "ascent and recovery state",
        flight_configuration,
        [engine],
    )
    drogue = rsim.Parachute.from_terminal_descent_speed(
        "700 ft/s illustrative drogue",
        drag_coefficient=DROGUE_DRAG_COEFFICIENT,
        terminal_descent_speed=DROGUE_TARGET_SPEED,
        design_mass=RECOVERY_DESIGN_MASS,
        design_air_density=RECOVERY_DESIGN_DENSITY,
        deployment_state=rsim.VehicleState.DROGUE_DEPLOYED,
    )
    main = rsim.Parachute.from_terminal_descent_speed(
        "5 ft/s illustrative main",
        drag_coefficient=MAIN_DRAG_COEFFICIENT,
        terminal_descent_speed=MAIN_TARGET_SPEED,
        design_mass=RECOVERY_DESIGN_MASS,
        design_air_density=RECOVERY_DESIGN_DENSITY,
        deployment_state=rsim.VehicleState.MAIN_DEPLOYED,
    )

    world = rsim.Frame("J2000-like inertial world", rsim.InertialStatus.INERTIAL)
    vehicle = rsim.Vehicle(
        "Vespula-inspired LOX/RP-1 sounding rocket",
        world,
        dry.mass,
        dry.moment_of_inertia,
        center_of_gravity=dry.center_of_mass,
        models=[
            lox_pressure_controller,
            rp1_pressure_controller,
            propulsion,
            flight_state,
        ],
        force_torques=[engine, aerodynamics, gravity, drogue, main],
    )
    vehicle.rigid_body.add_variable_mass_source(lox_tank)
    vehicle.rigid_body.add_variable_mass_source(rp1_tank)
    vehicle.rigid_body.add_variable_mass_source(copv)
    vehicle.kinematics.set_state(
        position_in_parent=np.array([EARTH_RADIUS, 0.0, 0.0]),
        velocity_in_parent=np.zeros(3),
        parent_from_body=np.eye(3),
        angular_velocity_in_body=np.zeros(3),
    )
    sizing = {
        "tank_pressure": tank_pressure,
        "throat_area": throat_area,
        "lox_cda": lox_cda,
        "fuel_cda": fuel_cda,
        "lox_line_diameter": lox_line,
        "fuel_line_diameter": fuel_line,
        "lox_makeup_flow": lox_makeup_flow,
        "fuel_makeup_flow": fuel_makeup_flow,
    }
    return (
        vehicle, propulsion, engine, aerodynamics, gravity, flight_state,
        drogue, main,
        lox_pressure_controller, rp1_pressure_controller,
        copv, lox_tank, rp1_tank, dry, aero_geometry, sizing,
    )


def main() -> None:
    """Run the complete mission in C++ and read native telemetry afterward."""
    (
        vehicle, propulsion, engine, aerodynamics, gravity, truth_state,
        drogue, main_parachute,
        lox_controller, rp1_controller,
        _copv, _lox, _rp1, _dry, aero_geometry, _sizing,
    ) = build_vehicle()

    # Recovery state now follows actual flight-computer actuator results.
    truth_state.set_recovery_parachutes(drogue, main_parachute)
    truth_state.phase = rsim.ModelPhase.TRUTH

    # Python creates the atmosphere once; all interpolation occurs in C++.
    altitude_axis = np.linspace(0.0, 100_000.0, 501)
    temperature = np.maximum(216.65, 288.15 - 0.0065 * altitude_axis)
    density = 1.225 * np.exp(-altitude_axis / 8_500.0)
    speed_of_sound = np.sqrt(1.4 * 287.05287 * temperature)
    viscosity = 1.7894e-5 * (temperature / 288.15) ** 0.76
    atmosphere = rsim.AtmosphereTable(
        altitude_axis, density, speed_of_sound, viscosity
    )
    flight_conditions = rsim.FlightConditions(
        "native flight conditions", vehicle, atmosphere,
        up_direction=np.array([1.0, 0.0, 0.0]),
        reference_length=aero_geometry.length,
        altitude_origin=np.array([EARTH_RADIUS, 0.0, 0.0]),
    )
    flight_conditions.add_aerodynamics(aerodynamics)
    flight_conditions.add_parachute(drogue)
    flight_conditions.add_parachute(main_parachute)
    truth_state.enable_console_output(
        flight_conditions, aerodynamics, period=1.0
    )

    altimeter = rsim.ScalarSensor(
        "altimeter", truth_state, "truthAltitude"
    )
    vertical_velocity_sensor = rsim.ScalarSensor(
        "vertical velocity sensor", truth_state, "verticalVelocity"
    )
    drogue_actuator = rsim.ParachuteActuator("drogue actuator", drogue)
    main_actuator = rsim.ParachuteActuator("main actuator", main_parachute)
    plugin_path = Path(__file__).resolve().parents[1] / "build" / "vespula_flight_computer.so"
    flight_computer = rsim.FlightComputer(
        "Vespula flight computer", str(plugin_path),
        {
            "altitude": altimeter,
            "verticalVelocity": vertical_velocity_sensor,
        },
        {
            "drogueDeploy": drogue_actuator,
            "mainDeploy": main_actuator,
        },
        {
            "mainDeployAltitudeM": MAIN_DEPLOY_ALTITUDE,
            "apogeeVelocityDeadbandMps": APOGEE_VELOCITY_DEADBAND,
        },
    )

    for model in (
        flight_conditions, altimeter, vertical_velocity_sensor,
        flight_computer, drogue_actuator, main_actuator,
    ):
        vehicle.add_model(model)

    # Burnout and termination are native truth events, not Python loop branches.
    burnout = rsim.TruthEvent(
        "powered flight complete",
        f"simulationTime greater_than_or_equal {POWERED_DURATION} and "
        f"{engine.state_namespace}.chamberPressure greater_than 0.0",
    )
    for model in (lox_controller, rp1_controller, propulsion, engine):
        burnout.add_model_enable_command(model, False)
    stop_at_main = rsim.TruthEvent(
        "main parachute deployed",
        "vehicleState equal_to MAIN_DEPLOYED",
        stop_simulation=True,
    )

    logger = rsim.DataLogger(sample_period=SAMPLE_PERIOD)
    logger.reserve_samples(10_000)
    simulator = rsim.Simulator(TIME_STEP)
    simulator.set_vehicle_state_time_steps(vehicle, {
        rsim.VehicleState.ON_RAIL: 0.005,
        rsim.VehicleState.THRUSTING: 0.01,
        rsim.VehicleState.COASTING: 0.05,
        rsim.VehicleState.DESCENDING: 0.02,
        rsim.VehicleState.DROGUE_DEPLOYED: 0.01,
        rsim.VehicleState.MAIN_DEPLOYED: 0.01,
    })
    simulator.add_vehicle(vehicle)
    simulator.add_truth_event(burnout)
    simulator.add_truth_event(stop_at_main)
    simulator.add_data_logger(logger)

    # The GIL is released for this single call; no Python runs between steps.
    result = simulator.run(maximum_time=700.0, maximum_steps=2_000_000)
    data = logger.to_dict()
    assert result.reason == "truth_event"
    assert drogue.deployed and main_parachute.deployed
    assert len(data["time"]) == logger.sample_count
    # Force-registry children are recursively discovered by the logger.
    assert f"{engine.state_namespace}.chamberPressure" in data
    assert f"{gravity.state_namespace}.forceX" in data
    print("Native Vespula run completed")
    print(f"stop={result.reason}:{result.event_name}, t={result.final_time:.2f} s")
    print(f"samples={logger.sample_count}, altitude={truth_state.altitude:.1f} m")
    print(f"drogue={drogue.deployed}, main={main_parachute.deployed}")

if __name__ == "__main__":
    main()

# %%
