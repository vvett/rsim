import pathlib
import unittest

import numpy as np

import rsim


class NativeRuntimeTests(unittest.TestCase):
    def make_vehicle(self):
        world = rsim.Frame("native-test-world", rsim.InertialStatus.INERTIAL)
        truth = rsim.TruthState("truth")
        truth.phase = rsim.ModelPhase.TRUTH
        vehicle = rsim.Vehicle(
            "native-test-vehicle", world, 1.0, np.eye(3), models=[truth]
        )
        return vehicle, truth

    def test_logger_expression_and_native_run(self):
        vehicle, _truth = self.make_vehicle()
        logger = rsim.DataLogger(sample_period=0.1)
        stop = rsim.TruthEvent(
            "time limit",
            "simulationTime greater_than_or_equal 0.3 and globalStep greater_than 1",
            stop_simulation=True,
        )
        simulator = rsim.Simulator(0.1)
        simulator.add_vehicle(vehicle)
        simulator.add_truth_event(stop)
        simulator.add_data_logger(logger)
        result = simulator.run(maximum_time=1.0, maximum_steps=20)
        data = logger.to_dict()
        self.assertEqual(result.reason, "truth_event")
        self.assertTrue(stop.fired)
        self.assertIn("truth.truthAltitude", data)
        self.assertEqual(
            data["time"].shape, data["truth.truthAltitude"].shape
        )
        self.assertEqual(logger.sample_count, 4)
        self.assertEqual(data["truth.vehicleState"].dtype, np.int64)
        self.assertEqual(data["truth.apogeeCrossed"].dtype, np.bool_)
        self.assertEqual(data["truth.truthAltitude"].dtype, np.float64)
        self.assertEqual(
            logger.channel_metadata["truth.truthAltitude"],
            {"units": "m", "type": "float64"},
        )
        self.assertEqual(
            logger.channel_metadata["truth.apogeeCrossed"]["type"], "bool"
        )

    def test_plugin_api_and_cross_computer_actuator_ownership(self):
        vehicle, truth = self.make_vehicle()
        drogue = rsim.Parachute(
            "test drogue", 1.0, 1.0,
            rsim.VehicleState.DROGUE_DEPLOYED,
        )
        actuator = rsim.ParachuteActuator("test actuator", drogue)
        altitude = rsim.ScalarSensor("altimeter", truth, "truthAltitude")
        velocity = rsim.ScalarSensor(
            "vertical velocity", truth, "verticalVelocity"
        )
        plugin = pathlib.Path(__file__).resolve().parents[2] / "build" / "vespula_flight_computer.so"
        sensors = {
            "altitude": altitude,
            "verticalVelocity": velocity,
        }
        actuators = {
            "drogueDeploy": actuator,
            "mainDeploy": rsim.ParachuteActuator("main actuator", drogue),
        }
        configuration = {
            "mainDeployAltitudeM": 609.6,
            "apogeeVelocityDeadbandMps": 0.1,
        }
        first = rsim.FlightComputer(
            "first computer", str(plugin), sensors, actuators, configuration
        )
        second = rsim.FlightComputer(
            "second computer", str(plugin), sensors, actuators, configuration
        )
        for model in (
            altitude, velocity, actuator, actuators["mainDeploy"], first, second
        ):
            vehicle.add_model(model)
        simulator = rsim.Simulator(0.1)
        simulator.add_vehicle(vehicle)
        with self.assertRaisesRegex(ValueError, "multiple flight computers"):
            simulator.run(maximum_time=0.2, maximum_steps=2)

        with self.assertRaisesRegex(ValueError, "sensor count"):
            rsim.FlightComputer(
                "bad computer", str(plugin),
                {"altitude": altitude}, actuators, configuration,
            )
        with self.assertRaisesRegex(ValueError, "missing flight-computer actuator"):
            rsim.FlightComputer(
                "bad actuator computer", str(plugin), sensors,
                {"wrong": actuator,
                 "mainDeploy": actuators["mainDeploy"]},
                configuration,
            )
        with self.assertRaisesRegex(ValueError, "missing flight-computer configuration"):
            rsim.FlightComputer(
                "bad configuration computer", str(plugin), sensors,
                actuators,
                {"wrong": 1.0,
                 "apogeeVelocityDeadbandMps": 0.1},
            )

    def test_native_wind_angles_rotation_and_initial_sample(self):
        vehicle, _truth = self.make_vehicle()
        vehicle.kinematics.set_state(
            np.array([20.0, 0.0, 0.0]),
            np.array([10.0, 0.0, 0.0]),
            np.array([[0.0, -1.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0]]),
            np.zeros(3),
        )
        atmosphere = rsim.AtmosphereTable(
            [0.0, 10.0], [1.2, 1.0], [340.0, 335.0],
            [1.8e-5, 1.8e-5],
            [np.array([1.0, 2.0, 0.0]), np.array([1.0, 2.0, 0.0])],
        )
        conditions = rsim.FlightConditions(
            "flight conditions", vehicle, atmosphere,
            altitude_origin=np.zeros(3),
        )
        vehicle.add_model(conditions)
        logger = rsim.DataLogger()
        simulator = rsim.Simulator(0.1)
        simulator.add_vehicle(vehicle)
        simulator.add_data_logger(logger)
        simulator.run(maximum_time=0.1, maximum_steps=1)
        data = logger.to_dict()
        # The t=0 row is initialized from truth/environment without integration.
        self.assertEqual(data["flight_conditions.altitude"][0], 20.0)
        self.assertTrue(data["flight_conditions.atmosphereClamped"][0])
        self.assertGreater(data["flight_conditions.alphaDeg"][0], 0.0)
        self.assertNotEqual(
            data["flight_conditions.relativeAirVelocityBodyY"][0], 0.0
        )
        # phi is atan2(w, v): +Y is zero degrees and +Z is 90 degrees.
        self.assertAlmostEqual(
            abs(data["flight_conditions.phiAeroDeg"][0]), 180.0
        )

        zero_vehicle, _ = self.make_vehicle()
        zero_atmosphere = rsim.AtmosphereTable(
            [0.0, 10.0], [1.2, 1.0], [340.0, 335.0],
            [1.8e-5, 1.8e-5],
        )
        zero_conditions = rsim.FlightConditions(
            "zero wind conditions", zero_vehicle, zero_atmosphere
        )
        zero_vehicle.add_model(zero_conditions)
        zero_logger = rsim.DataLogger()
        zero_simulator = rsim.Simulator(0.1)
        zero_simulator.add_vehicle(zero_vehicle)
        zero_simulator.add_data_logger(zero_logger)
        zero_simulator.run(maximum_time=0.1, maximum_steps=1)
        zero_data = zero_logger.to_dict()
        self.assertEqual(
            zero_data["zero_wind_conditions.relativeAirVelocityBodyX"][0],
            0.0,
        )
        self.assertEqual(zero_data["zero_wind_conditions.alphaDeg"][0], 0.0)

    def test_common_valve_and_model_enable_actuators(self):
        vehicle, truth = self.make_vehicle()
        valve = rsim.Valve(1.0e-5, opening=1.0)
        valve_actuator = rsim.ValveActuator("valve actuator", valve)
        enable_actuator = rsim.ModelEnableActuator(
            "truth enable actuator", truth
        )
        valve_actuator.command(0.25)
        enable_actuator.command(0.0)
        vehicle.add_model(valve_actuator)
        vehicle.add_model(enable_actuator)
        simulator = rsim.Simulator(0.1)
        simulator.add_vehicle(vehicle)
        simulator.run(maximum_time=0.1, maximum_steps=1)
        self.assertAlmostEqual(valve_actuator.actual, 0.25)
        self.assertFalse(enable_actuator.actual)
        self.assertFalse(truth.enabled)

    def test_truth_event_commands_through_ideal_actuator(self):
        valve = rsim.Valve(1.0e-5, opening=0.0)
        actuator = rsim.ValveActuator("event valve actuator", valve)
        event = rsim.TruthEvent(
            "open valve",
            "simulationTime greater_than_or_equal 0.1",
        )
        event.add_actuator_command(actuator, 0.6)
        simulator = rsim.Simulator(0.1)
        simulator.add_model(actuator)
        simulator.add_truth_event(event)
        simulator.run(maximum_time=0.2, maximum_steps=2)
        self.assertAlmostEqual(actuator.actual, 0.6)

    def test_specific_force_excludes_gravity_and_includes_applied_force(self):
        world = rsim.Frame("specific-force-world", rsim.InertialStatus.INERTIAL)
        gravity = rsim.Gravity.earth()
        applied = rsim.ConstantForceTorque("test thrust")
        applied.set_load(np.array([20.0, 0.0, 0.0]), np.zeros(3))
        vehicle = rsim.Vehicle(
            "specific-force vehicle", world, 10.0, np.eye(3),
            force_torques=[gravity, applied],
        )
        vehicle.kinematics.set_state(
            np.array([6_378_137.0, 0.0, 0.0]), np.zeros(3),
            np.eye(3), np.zeros(3),
        )
        sensor = rsim.SpecificForceSensor("accelerometer", vehicle)
        vehicle.add_model(sensor)
        simulator = rsim.Simulator(0.01)
        simulator.add_vehicle(vehicle)
        simulator.run(maximum_time=0.01, maximum_steps=1)
        np.testing.assert_allclose(
            sensor.specific_force_in_body, np.array([2.0, 0.0, 0.0])
        )
        gravity_only = rsim.Gravity.earth()
        gravity_only_vehicle = rsim.Vehicle(
            "gravity-only vehicle", world, 10.0, np.eye(3),
            force_torques=[gravity_only],
        )
        gravity_only_vehicle.kinematics.set_state(
            np.array([6_378_137.0, 0.0, 0.0]), np.zeros(3),
            np.eye(3), np.zeros(3),
        )
        gravity_only_sensor = rsim.SpecificForceSensor(
            "gravity-only accelerometer", gravity_only_vehicle
        )
        gravity_only_vehicle.add_model(gravity_only_sensor)
        gravity_simulator = rsim.Simulator(0.01)
        gravity_simulator.add_vehicle(gravity_only_vehicle)
        gravity_simulator.run(maximum_time=0.01, maximum_steps=1)
        np.testing.assert_allclose(
            gravity_only_sensor.specific_force_in_body, np.zeros(3)
        )


if __name__ == "__main__":
    unittest.main()
