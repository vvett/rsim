import tempfile
import unittest
from pathlib import Path

import numpy as np
import rsim

from rsim.aero import (
    AeroCoefficientTable,
    AeroGrid,
    AerodynamicOptions,
    CylindricalBody,
    GridAxis,
    RocketGeometry,
    TangentOgiveNose,
    TrapezoidalFinSet,
    generate_aero_table,
)


def example_geometry() -> RocketGeometry:
    return RocketGeometry(
        nose=TangentOgiveNose(length=0.3, base_diameter=0.1),
        body=CylindricalBody(length=3.7, diameter=0.1),
        fins=TrapezoidalFinSet(
            count=4,
            root_chord=0.4,
            tip_chord=0.2,
            semi_span=0.15,
            tip_leading_edge_offset=0.12,
            thickness=0.008,
            root_leading_edge_x=3.4,
        ),
        moment_reference_x=2.4,
    )


class AeroGridTests(unittest.TestCase):
    def test_even_spacing_is_inclusive(self):
        axis = GridAxis.evenly_spaced(-2.0, 2.0, 1.0)
        np.testing.assert_allclose(axis.array(), [-2, -1, 0, 1, 2])

    def test_nonintegral_spacing_is_rejected(self):
        with self.assertRaises(ValueError):
            GridAxis.evenly_spaced(0.0, 1.0, 0.3)


class AeroDatabaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = AeroGrid(
            mach=GridAxis([0.5, 1.0, 1.3, 2.0, 5.0]),
            alpha_deg=GridAxis([-2.0, 0.0, 2.0]),
            phi_aero_deg=GridAxis([0.0, 45.0, 90.0]),
            reynolds=GridAxis([1.0e6, 2.0e6]),
        )
        cls.table = generate_aero_table(
            example_geometry(),
            cls.grid,
            AerodynamicOptions(body_station_count=41, body_azimuth_count=16),
        )

    def test_table_uses_ndgrid_shape_and_order(self):
        self.assertEqual(self.table.coefficients["CA"].shape, (5, 3, 3, 2))
        self.assertEqual(
            self.table.metadata["dimension_order"],
            ["mach", "alpha_deg", "phi_aero_deg", "reynolds"],
        )
        np.testing.assert_allclose(
            self.table.meshes["mach"][:, 0, 0, 0], self.grid.mach.array()
        )
        np.testing.assert_allclose(
            self.table.meshes["alpha_deg"][0, :, 0, 0], self.grid.alpha_deg.array()
        )
        np.testing.assert_allclose(
            self.table.meshes["phi_aero_deg"][0, 0, :, 0],
            self.grid.phi_aero_deg.array(),
        )
        np.testing.assert_allclose(
            self.table.meshes["reynolds"][0, 0, 0, :], self.grid.reynolds.array()
        )

    def test_symmetric_rocket_has_expected_alpha_symmetry(self):
        np.testing.assert_allclose(
            self.table.coefficients["CN"][:, 0, :, :],
            -self.table.coefficients["CN"][:, 2, :, :],
            rtol=2.0e-10,
            atol=2.0e-10,
        )
        np.testing.assert_allclose(
            self.table.coefficients["Cm"][:, 0, :, :],
            -self.table.coefficients["Cm"][:, 2, :, :],
            rtol=2.0e-10,
            atol=2.0e-10,
        )
        np.testing.assert_allclose(
            self.table.coefficients["CA"][:, 0, :, :],
            self.table.coefficients["CA"][:, 2, :, :],
            rtol=2.0e-10,
            atol=2.0e-10,
        )

    def test_zero_angle_and_reynolds_behavior(self):
        np.testing.assert_allclose(
            self.table.coefficients["CN"][:, 1, :, :], 0.0, atol=1e-12
        )
        np.testing.assert_allclose(
            self.table.coefficients["Cm"][:, 1, :, :], 0.0, atol=1e-12
        )
        self.assertTrue(
            np.all(np.isfinite(self.table.coefficients["xCP"][:, 1, :, :]))
        )
        self.assertTrue(
            np.all(
                self.table.components["ca_body_friction"][:, :, :, 1]
                < self.table.components["ca_body_friction"][:, :, :, 0]
            )
        )

    def test_router_crosses_all_primary_body_ranges(self):
        flags = self.table.method_flags["body_pressure"][:, 1, 0, 0].tolist()
        self.assertIn("CHART_NWL_TR2796_FIG3", flags)
        self.assertIn("FALLBACK_SOSE_LOCAL_PANELS", flags)
        self.assertIn("FALLBACK_MODIFIED_NEWTONIAN_PANELS", flags)

    def test_drag_is_positive_and_moment_reconstructs_cp(self):
        self.assertTrue(np.all(self.table.coefficients["CD"] > 0.0))
        nonzero_alpha = 0
        reconstructed = (
            example_geometry().moment_reference_x
            - example_geometry().length
            * self.table.coefficients["Cm"][:, nonzero_alpha, :, :]
            / self.table.coefficients["CN"][:, nonzero_alpha, :, :]
        )
        np.testing.assert_allclose(
            self.table.coefficients["xCP"][:, nonzero_alpha, :, :], reconstructed
        )

    def test_phi_rotates_body_axes_and_preserves_four_fin_periodicity(self):
        # A 90-degree roll repeats a four-fin configuration and rotates Y into Z.
        np.testing.assert_allclose(
            self.table.coefficients["CN"][:, :, 0, :],
            self.table.coefficients["CN"][:, :, 2, :],
            rtol=2e-10, atol=2e-10,
        )
        np.testing.assert_allclose(
            self.table.coefficients["CY"][:, :, 0, :],
            self.table.coefficients["CZ"][:, :, 2, :],
            rtol=2e-10, atol=2e-10,
        )
        self.assertGreater(
            np.max(np.abs(
                self.table.coefficients["CN"][:, 2, 0, :]
                - self.table.coefficients["CN"][:, 2, 1, :]
            )),
            1e-8,
        )
        self.assertTrue(np.all(np.abs(self.table.coefficients["CS"]) < 1e-10))

    def test_npz_round_trip_uses_plain_arrays(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "aero.npz"
            self.table.save_npz(path)
            loaded = AeroCoefficientTable.load_npz(path)
        self.assertEqual(loaded.metadata, self.table.metadata)
        for name, expected in self.table.coefficients.items():
            np.testing.assert_equal(loaded.coefficients[name], expected)
        for name, expected in self.table.method_flags.items():
            np.testing.assert_equal(loaded.method_flags[name], expected)

    def test_mat_export_contains_axes_coefficients_and_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "aero.mat"
            self.table.save_mat(path)
            contents = path.read_bytes()
        self.assertEqual(len(contents[:128]), 128)
        self.assertIn(b"MATLAB 5.0 MAT-file", contents[:116])
        self.assertIn(b"axis__phi_aero_deg", contents)
        self.assertIn(b"coefficient__CA", contents)
        self.assertIn(b"metadata_json_utf8", contents)

    def test_cpp_reader_interpolates_and_dimensionalizes_generated_database(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "aero.mat"
            self.table.save_mat(path)
            database = rsim.AerodynamicDatabase(str(path))

            coefficients = database.evaluate(1.3, 2.0, 45.0, 2.0e6)
            expected_index = (2, 2, 1, 1)
            self.assertAlmostEqual(
                coefficients.axial,
                self.table.coefficients["CA"][expected_index],
            )
            self.assertAlmostEqual(
                coefficients.side_y,
                self.table.coefficients["CY"][expected_index],
            )

            midpoint = database.evaluate(1.15, 1.0, 22.5, 1.5e6)
            expected_midpoint = np.mean(
                self.table.coefficients["CA"][1:3, 1:3, 0:2, 0:2]
            )
            self.assertAlmostEqual(midpoint.axial, expected_midpoint)

            # Scalar magnitudes repeat every 90 degrees, while vector
            # components rotate into the requested body-axis quadrant.
            wrapped = database.evaluate(1.3, 2.0, 135.0, 2.0e6)
            self.assertAlmostEqual(wrapped.axial, coefficients.axial)
            self.assertAlmostEqual(wrapped.side_y, -coefficients.side_z)
            self.assertAlmostEqual(wrapped.side_z, coefficients.side_y)

            opposite = database.evaluate(1.3, 2.0, 225.0, 2.0e6)
            self.assertAlmostEqual(opposite.side_y, -coefficients.side_y)
            self.assertAlmostEqual(opposite.side_z, -coefficients.side_z)

            aerodynamics = rsim.Aerodynamics("aero", database)
            conditions = rsim.AerodynamicConditions()
            conditions.mach = 1.3
            conditions.alpha_deg = 2.0
            conditions.phi_aero_deg = 45.0
            conditions.reynolds = 2.0e6
            conditions.dynamic_pressure = 5000.0
            aerodynamics.set_conditions(conditions)
            simulator = rsim.Simulator(0.01)
            simulator.add_model(aerodynamics)
            simulator.step()

            force_scale = 5000.0 * example_geometry().reference_area
            np.testing.assert_allclose(
                aerodynamics.load.force,
                force_scale * np.array([
                    -coefficients.axial,
                    -coefficients.side_y,
                    -coefficients.side_z,
                ]),
            )
            np.testing.assert_allclose(aerodynamics.force,
                                       aerodynamics.load.force)

    def test_vehicle_schedules_configured_aerodynamics_before_kinematics(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "aero.mat"
            self.table.save_mat(path)
            world = rsim.Frame("world-aero", rsim.InertialStatus.INERTIAL)
            aerodynamics = rsim.Aerodynamics(
                "aero", rsim.AerodynamicDatabase(str(path))
            )
            vehicle = rsim.Vehicle(
                "aero-rocket", world, 10.0, np.eye(3),
                center_of_gravity=np.array([2.5, 0.0, 0.0]),
                force_torques=[aerodynamics]
            )
            conditions = rsim.AerodynamicConditions()
            conditions.mach = 1.0
            conditions.alpha_deg = 2.0
            conditions.phi_aero_deg = 0.0
            conditions.reynolds = 1.0e6
            conditions.dynamic_pressure = 1000.0
            aerodynamics.set_conditions(conditions)
            simulator = rsim.Simulator(0.01)
            simulator.add_vehicle(vehicle)
            simulator.step()

            np.testing.assert_allclose(
                vehicle.forces.force, aerodynamics.load.force
            )
            database_cp = self.table.coefficients["xCP"][1, 2, 0, 0]
            expected_body_cp = example_geometry().length - database_cp
            self.assertAlmostEqual(
                aerodynamics.center_pressure_x, expected_body_cp
            )
            self.assertAlmostEqual(
                aerodynamics.static_stability_margin,
                (2.5 - expected_body_cp) / example_geometry().reference_diameter,
            )
            self.assertLess(aerodynamics.force[1], 0.0)
            self.assertGreater(
                vehicle.forces.torque_about(
                    vehicle.rigid_body.center_of_gravity
                )[2],
                0.0,
            )
            self.assertEqual(vehicle.kinematics.global_step, 0)


if __name__ == "__main__":
    unittest.main()
