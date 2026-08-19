"""Exercise the Python-facing model ownership and scheduling API."""

import unittest
import warnings
import math

import numpy as np

import rsim


class SimulatorTests(unittest.TestCase):
    def test_vehicle_children_are_exposed_and_scheduled(self):
        world = rsim.Frame("world", rsim.InertialStatus.INERTIAL)
        vehicle = rsim.Vehicle("rocket", world, 10.0, np.eye(3))
        simulator = rsim.Simulator(0.01)

        simulator.add_vehicle(vehicle)
        simulator.step()

        self.assertEqual(vehicle.rigid_body.name, "rocket-rigid-body")
        self.assertEqual(vehicle.forces.name, "rocket-forces")
        self.assertEqual(vehicle.kinematics.global_step, 0)
        self.assertEqual(vehicle.body_frame.parent.name, "world")
        self.assertAlmostEqual(simulator.time, 0.01)

    def test_kinematics_integrates_translation_and_orientation(self):
        world = rsim.Frame("world", rsim.InertialStatus.INERTIAL)
        body = rsim.RigidBody(2.0, np.eye(3))
        load = rsim.ConstantForceTorque()
        load.set_load(np.array([2.0, 0.0, 0.0]), np.zeros(3))
        vehicle = rsim.Vehicle("rocket-list", world, 2.0, np.eye(3),
                               force_torques=[load])
        # Vehicle retains the native C++ model after the setup reference is gone.
        del load
        motion = vehicle.kinematics
        motion.set_state(np.zeros(3), np.zeros(3), np.eye(3),
                         np.array([0.0, 0.0, 1.0]))
        simulator = rsim.Simulator(0.1)
        simulator.add_vehicle(vehicle)

        simulator.step()

        np.testing.assert_allclose(motion.velocity_in_parent,
                                   [0.0998334, 0.00499583, 0.0], rtol=1e-5)
        np.testing.assert_allclose(motion.position_in_parent,
                                   [0.00499583, 0.000166583, 0.0],
                                   rtol=1e-5, atol=1e-8)
        self.assertAlmostEqual(motion.angular_velocity_in_body[2], 1.0)

    def test_exact_update_rate_factories(self):
        fast = rsim.UpdateRate.times_per_step(4)
        slow = rsim.UpdateRate.every_n_steps(3)

        self.assertEqual(fast.updates_per_step, 4)
        self.assertEqual(fast.steps_per_update, 1)
        self.assertEqual(slow.updates_per_step, 1)
        self.assertEqual(slow.steps_per_update, 3)

    def test_force_models_can_be_added_and_removed_after_construction(self):
        world = rsim.Frame("world-add", rsim.InertialStatus.INERTIAL)
        vehicle = rsim.Vehicle("rocket-add", world, 1.0, np.eye(3))
        load = rsim.ConstantForceTorque("external")
        load.set_load(np.array([3.0, 0.0, 0.0]), np.zeros(3))

        vehicle.add_force_torque(load)
        self.assertEqual(len(vehicle.forces), 1)
        np.testing.assert_allclose(vehicle.forces.force, [3.0, 0.0, 0.0])

        vehicle.remove_force_torque(load)
        self.assertEqual(len(vehicle.forces), 0)

    def test_integration_frame_can_be_selected_by_name(self):
        frames = rsim.FrameGraph()
        inertial = rsim.Frame("j2000", rsim.InertialStatus.INERTIAL)
        rotating = rsim.Frame("ecef", rsim.InertialStatus.NON_INERTIAL)
        frames.add_frame(inertial)
        frames.add_frame(rotating)
        body = rsim.RigidBody(1.0, np.eye(3))
        loads = rsim.ForceTorqueRegistry()

        motion = rsim.Kinematics("body", frames, "j2000", loads, body)
        self.assertEqual(motion.integration_frame.name, "j2000")

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            non_inertial_motion = rsim.Kinematics(
                "body-ecef", frames, "ecef", loads, body
            )
        self.assertEqual(non_inertial_motion.integration_frame.name, "ecef")
        self.assertTrue(any(item.category is RuntimeWarning for item in caught))

    def test_vehicle_state_time_steps_are_copied_and_can_be_cleared(self):
        world = rsim.Frame("state-world", rsim.InertialStatus.INERTIAL)
        vehicle = rsim.Vehicle("state-rocket", world, 1.0, np.eye(3))
        simulator = rsim.Simulator(0.1)
        policy = {
            rsim.VehicleState.ON_RAIL: 0.02,
            rsim.VehicleState.THRUSTING: 0.005,
            rsim.VehicleState.MAIN_DEPLOYED: 0.2,
        }
        simulator.set_vehicle_state_time_steps(vehicle, policy)
        # Mutating the Python dict cannot change the copied native policy.
        policy[rsim.VehicleState.ON_RAIL] = 9.0
        simulator.add_vehicle(vehicle)

        simulator.step()
        self.assertAlmostEqual(simulator.time, 0.02)
        self.assertAlmostEqual(simulator.last_time_step, 0.02)

        vehicle.state = rsim.VehicleState.THRUSTING
        simulator.step()
        self.assertAlmostEqual(simulator.time, 0.025)

        vehicle.state = rsim.VehicleState.COASTING
        simulator.step()
        self.assertAlmostEqual(simulator.time, 0.125)

        vehicle.state = rsim.VehicleState.MAIN_DEPLOYED
        simulator.clear_vehicle_state_time_steps(vehicle)
        simulator.step()
        self.assertAlmostEqual(simulator.time, 0.225)

        simulator.global_time_step = 0.15
        self.assertAlmostEqual(simulator.global_time_step, 0.15)
        simulator.step()
        self.assertAlmostEqual(simulator.time, 0.375)

    def test_multiple_vehicle_policies_choose_the_smallest_request(self):
        world = rsim.Frame("multi-world", rsim.InertialStatus.INERTIAL)
        first = rsim.Vehicle("first", world, 1.0, np.eye(3))
        second = rsim.Vehicle("second", world, 1.0, np.eye(3))
        simulator = rsim.Simulator(0.1)
        simulator.add_vehicle(first)
        simulator.add_vehicle(second)
        simulator.set_vehicle_state_time_steps(
            first, {rsim.VehicleState.ON_RAIL: 0.04}
        )
        simulator.set_vehicle_state_time_steps(
            second, {rsim.VehicleState.ON_RAIL: 0.015}
        )

        simulator.step()
        self.assertAlmostEqual(simulator.last_time_step, 0.015)
        self.assertAlmostEqual(simulator.time, 0.015)

        with self.assertRaises(ValueError):
            simulator.set_vehicle_state_time_steps(
                first, {rsim.VehicleState.ON_RAIL: math.inf}
            )
        with self.assertRaises(ValueError):
            simulator.global_time_step = 0.0


if __name__ == "__main__":
    unittest.main()
