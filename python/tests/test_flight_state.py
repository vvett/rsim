"""Exercise native flight-state transitions and parachute setup from Python."""

import unittest

import numpy as np

import rsim


class FlightStateTests(unittest.TestCase):
    def test_truth_state_is_primary_and_legacy_recovery_is_opt_in(self):
        self.assertIs(rsim.TruthState, rsim.FlightStateModel)
        self.assertFalse(
            rsim.FlightStateConfiguration().legacy_automatic_recovery
        )

    def test_deployment_transition_updates_parachute_in_same_native_step(self):
        world = rsim.Frame("world-flight", rsim.InertialStatus.INERTIAL)
        config = rsim.FlightStateConfiguration()
        config.legacy_automatic_recovery = True
        config.rail_length = 1.0
        config.drogue_trigger.maximum_altitude = 100.0
        config.drogue_trigger.maximum_vertical_velocity = 0.0
        config.main_stage_enabled = False

        flight_state = rsim.FlightStateModel(
            configuration=config,
        )
        drogue = rsim.Parachute(
            "drogue",
            drag_coefficient=1.0,
            reference_area=2.0,
            deployment_state=rsim.VehicleState.DROGUE_DEPLOYED,
            air_density=1.0,
        )
        vehicle = rsim.Vehicle(
            "rocket-flight",
            world,
            20.0,
            np.eye(3),
            models=[flight_state],
            force_torques=[drogue],
        )
        vehicle.state = rsim.VehicleState.DESCENDING
        vehicle.kinematics.set_state(
            np.array([50.0, 0.0, 0.0]),
            np.array([-2.0, 0.0, 0.0]),
            np.eye(3),
            np.zeros(3),
        )

        # Vehicle owns native shared pointers after Python setup references go away.
        del flight_state
        del drogue
        simulator = rsim.Simulator(0.01)
        simulator.add_vehicle(vehicle)
        simulator.step()

        self.assertEqual(vehicle.state, rsim.VehicleState.DROGUE_DEPLOYED)
        self.assertTrue(vehicle.forces[0].deployed)
        # -0.5*rho*Cd*A*|v|*v = +4 N for v = -2 m/s.
        np.testing.assert_allclose(vehicle.forces.force, [4.0, 0.0, 0.0])

    def test_diameter_constructor_and_invalid_deployment_state(self):
        chute = rsim.Parachute.from_diameter(
            "main", 1.5, 2.0, rsim.VehicleState.MAIN_DEPLOYED
        )
        self.assertAlmostEqual(chute.reference_area, np.pi)

        with self.assertRaises(ValueError):
            rsim.Parachute(
                "invalid", 1.0, 1.0, rsim.VehicleState.COASTING
            )

    def test_terminal_speed_sizing_and_recovery_trigger_configuration(self):
        feet_to_meter = 0.3048
        self.assertEqual(700.0 * feet_to_meter, 213.36)
        self.assertEqual(5.0 * feet_to_meter, 1.524)
        self.assertEqual(2000.0 * feet_to_meter, 609.6)

        mass = 80.0
        density = 1.1
        cd = 1.4
        terminal_speed = 1.524  # 5 ft/s
        chute = rsim.Parachute.from_terminal_descent_speed(
            "main",
            drag_coefficient=cd,
            terminal_descent_speed=terminal_speed,
            design_mass=mass,
            design_air_density=density,
            deployment_state=rsim.VehicleState.MAIN_DEPLOYED,
        )
        expected_area = (
            2.0 * mass * 9.80665 /
            (density * cd * terminal_speed**2)
        )
        self.assertAlmostEqual(chute.reference_area, expected_area)
        self.assertAlmostEqual(
            chute.equivalent_diameter,
            np.sqrt(4.0 * expected_area / np.pi),
        )
        drogue = rsim.Parachute.from_terminal_descent_speed(
            "drogue",
            drag_coefficient=0.8,
            terminal_descent_speed=213.36,
            design_mass=mass,
            design_air_density=density,
            deployment_state=rsim.VehicleState.DROGUE_DEPLOYED,
        )
        self.assertNotEqual(drogue.drag_coefficient,
                            chute.drag_coefficient)

        config = rsim.FlightStateConfiguration()
        config.legacy_automatic_recovery = True
        config.drogue_trigger.deploy_at_apogee = True
        config.apogee_velocity_deadband = 0.1
        config.main_trigger.maximum_altitude = 609.6
        config.main_trigger.require_descending_altitude_crossing = True
        model = rsim.FlightStateModel(configuration=config)

        self.assertTrue(model.configuration.drogue_trigger.deploy_at_apogee)
        self.assertEqual(model.configuration.main_trigger.maximum_altitude,
                         609.6)
        self.assertTrue(
            model.configuration.main_trigger
            .require_descending_altitude_crossing
        )


if __name__ == "__main__":
    unittest.main()
