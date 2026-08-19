"""Public Python API checks for the native central-gravity model."""

import gc
import math
import unittest

import numpy as np

import rsim


class GravityTests(unittest.TestCase):
    def test_earth_factory_registers_and_survives_python_reference_release(self):
        radius = 6_378_137.0
        mass = 15.0
        frame = rsim.Frame("J2000", rsim.InertialStatus.INERTIAL)
        gravity = rsim.Gravity.earth()

        vehicle = rsim.Vehicle(
            "rocket",
            frame,
            mass,
            np.eye(3),
            force_torques=[gravity],
        )
        vehicle.kinematics.set_state(
            np.array([radius, 0.0, 0.0]),
            np.zeros(3),
            np.eye(3),
            np.zeros(3),
        )

        del gravity
        gc.collect()

        simulator = rsim.Simulator(0.01)
        simulator.add_vehicle(vehicle)
        simulator.step()

        expected = mass * rsim.Gravity.EARTH_GRAVITATIONAL_PARAMETER / radius**2
        self.assertTrue(
            np.allclose(vehicle.forces.force, [-expected, 0.0, 0.0])
        )
        self.assertTrue(
            np.allclose(
                vehicle.forces[0].torque_about(vehicle.rigid_body.center_of_gravity),
                np.zeros(3),
            )
        )

    def test_custom_center_and_invalid_inputs(self):
        gravity = rsim.Gravity(
            "moon",
            4.9048695e12,
            np.array([100.0, -20.0, 4.0]),
        )
        self.assertEqual(gravity.name, "moon")
        self.assertTrue(
            np.allclose(gravity.center_in_integration_frame, [100.0, -20.0, 4.0])
        )
        self.assertTrue(math.isclose(gravity.gravitational_parameter, 4.9048695e12))

        with self.assertRaises(ValueError):
            rsim.Gravity("invalid", -1.0)


if __name__ == "__main__":
    unittest.main()
