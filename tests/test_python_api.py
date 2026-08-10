import unittest

import numpy as np

from rsim import RigidBody


class RigidBodyTests(unittest.TestCase):
    def setUp(self):
        self.inertia = ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))

    def test_properties_round_trip(self):
        body = RigidBody(10.0, self.inertia)
        self.assertEqual(body.mass, 10.0)
        np.testing.assert_allclose(
            body.moment_of_inertia,
            [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
        body.mass = 12.0
        self.assertEqual(body.mass, 12.0)

    def test_invalid_mass_is_value_error(self):
        with self.assertRaises(ValueError):
            RigidBody(0.0, self.inertia)

    def test_inertia_requires_three_by_three_values(self):
        with self.assertRaises((TypeError, ValueError)):
            RigidBody(10.0, ((1.0, 0.0), (0.0, 1.0)))


if __name__ == "__main__":
    unittest.main()
