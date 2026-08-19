"""Validate dry geometry and its C++ variable-propellant handoff."""

import unittest
from math import pi

import numpy as np

import rsim


class GeometryTests(unittest.TestCase):
    def test_dry_components_are_combined_about_the_vehicle_cg(self):
        aluminum = rsim.Material("aluminum", 2700.0)
        rocket = rsim.LiquidSoundingRocketGeometry()
        rocket.add_component(rsim.SolidCylinder(0.1, 0.2, 0.0, aluminum))
        rocket.add_component(rsim.SolidCylinder(0.1, 0.2, 1.0, aluminum))

        properties = rocket.dry_mass_properties()

        expected_mass = 2.0 * 2700.0 * pi * 0.1**2 * 0.2
        self.assertAlmostEqual(properties.mass, expected_mass)
        np.testing.assert_allclose(properties.center_of_mass, [0.5, 0.0, 0.0])
        self.assertGreater(properties.moment_of_inertia[1, 1],
                           properties.moment_of_inertia[0, 0])

    def test_tank_walls_become_dry_mass_and_liquid_updates_in_cpp(self):
        aluminum = rsim.Material("aluminum", 2700.0)
        rocket = rsim.LiquidSoundingRocketGeometry()
        rocket.add_component(rsim.SolidCylinder(0.08, 0.2, -0.2, aluminum))
        tank = rsim.CylindricalTank(
            outer_radius=0.15,
            length=1.0,
            wall_thickness=0.005,
            center_x=0.5,
            material=aluminum,
            propellant_density=1000.0,
            initial_fill_fraction=0.5,
            mass_flow_rate=1.0,
        )
        rocket.add_tank(tank)

        body = rocket.create_rigid_body("rocket-mass")
        dry = rocket.dry_mass_properties()

        self.assertEqual(body.propellant_tank_count, 1)
        self.assertAlmostEqual(
            body.mass, dry.mass + 0.5 * tank.propellant_capacity
        )
        simulator = rsim.Simulator(1.0)
        simulator.add_model(body)
        simulator.step()
        simulator.step()
        self.assertAlmostEqual(
            body.mass, dry.mass + 0.5 * tank.propellant_capacity - 1.0
        )


if __name__ == "__main__":
    unittest.main()
