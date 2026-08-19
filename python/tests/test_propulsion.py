"""Native propulsion integration, ownership, and table-copy tests."""

import gc
import tempfile
import unittest
import weakref

import numpy as np

import rsim


class PropulsionTests(unittest.TestCase):
    def make_table(self):
        pc = np.array([1.0e4, 5.0e6])
        mr = np.array([0.5, 4.0])
        shape = (pc.size, mr.size)
        temperature = np.full(shape, 3400.0, order="F")
        cstar = np.full(shape, 1500.0, order="F")
        cf = np.full(shape, 1.7, order="F")
        isp = np.full(shape, 260.0, order="F")
        reference = weakref.ref(cstar)
        table = rsim.CombustionTableData(
            pc, mr, temperature, cstar, cf, isp, 5.0, "N2O", "ethanol", False
        )
        del pc, mr, temperature, cstar, cf, isp
        gc.collect()
        self.assertIsNone(reference())
        self.assertAlmostEqual(table.cstar(1.0e6, 2.0), 1500.0)
        return table

    def test_network_vehicle_schedule_and_mass_coupling(self):
        table = self.make_table()
        ambient = rsim.Ambient(pressure=0.0)
        ox = rsim.Tank("ox", 0.2, 1.0, 0.0, 900.0, 80.0, 3.0e6)
        fuel = rsim.Tank("fuel", 0.2, 1.0, 1.2, 800.0, 70.0, 3.0e6)
        engine = rsim.Engine("engine", table, 1.0e-3, 5.0,
                             ambient=ambient,
                             application_point=np.array([0.0, 0.1, 0.0]),
                             initial_pressure=5.0e5)
        propulsion = rsim.Propulsion()
        propulsion.connect(ox, engine, rsim.Valve(1.0e-7),
                           upstream_phase=rsim.FluidPhase.LIQUID,
                           downstream_role=rsim.InletRole.OXIDIZER)
        propulsion.connect(fuel, engine, rsim.Orifice(1.0e-7),
                           upstream_phase=rsim.FluidPhase.LIQUID,
                           downstream_role=rsim.InletRole.FUEL)
        propulsion.finalize()

        world = rsim.Frame("world-prop", rsim.InertialStatus.INERTIAL)
        vehicle = rsim.Vehicle("rocket-prop", world, 10.0, np.eye(3),
                               models=[propulsion], force_torques=[engine])
        vehicle.rigid_body.add_variable_mass_source(ox)
        vehicle.rigid_body.add_variable_mass_source(fuel)
        stored = propulsion.total_stored_mass
        simulator = rsim.Simulator(0.01)
        simulator.add_vehicle(vehicle)
        simulator.step()

        self.assertLess(propulsion.total_stored_mass, stored)
        self.assertGreater(engine.thrust, 0.0)
        self.assertGreater(engine.force[0], 0.0)
        self.assertAlmostEqual(vehicle.rigid_body.mass,
                               10.0 + propulsion.total_stored_mass)
        self.assertIn("converged", propulsion.last_solve_summary)

    def test_topology_is_immutable_and_skeleton_fails(self):
        tank = rsim.Tank("tank", 0.2, 1.0, 0.0, 900.0, 80.0, 3.0e6)
        ambient = rsim.Ambient()
        network = rsim.Propulsion()
        network.vent(tank, ambient, rsim.Orifice(1.0e-8))
        network.finalize()
        with self.assertRaises(RuntimeError):
            network.vent(tank, ambient, rsim.Orifice(1.0e-8))

        invalid = rsim.Propulsion("invalid")
        invalid.connect(tank, rsim.Junction("j"), rsim.Pump())
        with self.assertRaises(RuntimeError):
            invalid.finalize()

    def test_mat_round_trip_and_flow_laws(self):
        pc = np.array([1.0e5, 2.0e6])
        mr = np.array([1.0, 3.0])
        values = np.arange(4.0).reshape((2, 2), order="F") + 1000.0
        setup = rsim.CombustionTable(
            pc, mr, values + 2000.0, values, np.full((2, 2), 1.6),
            np.full((2, 2), 250.0), 4.0, "oxidizer", "fuel", True,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = f"{directory}/combustion.mat"
            setup.save_mat(path)
            loaded = rsim.CombustionTable.load_mat(path)
            self.assertEqual(loaded.oxidizer, "oxidizer")
            self.assertTrue(loaded.frozen_at_throat)
            self.assertAlmostEqual(loaded.native.cstar(pc[1], mr[0]),
                                   values[1, 0])

        high = rsim.Ambient("high", 5.0e5, 300.0)
        vacuum = rsim.Ambient("vacuum", 0.0, 300.0)
        medium = rsim.Ambient("medium", 1.0e5, 300.0)
        restriction = rsim.Orifice(2.0e-6)
        choked_vacuum = restriction.mass_flow(high, vacuum, rsim.FluidPhase.GAS)
        choked_medium = restriction.mass_flow(high, medium, rsim.FluidPhase.GAS)
        self.assertGreater(choked_vacuum, 0.0)
        self.assertAlmostEqual(choked_vacuum, choked_medium, places=12)
        self.assertLess(restriction.mass_flow(vacuum, high,
                                              rsim.FluidPhase.GAS), 0.0)

        liquid = rsim.Tank("liquid", 0.2, 1.0, 0.0, 1000.0, 100.0,
                           5.0e5)
        line = rsim.Line(1.0, 0.01, 1.0e-5, 1.0e-3)
        self.assertGreater(line.mass_flow(liquid, medium,
                                          rsim.FluidPhase.LIQUID), 0.0)

    def test_shared_pressurant_branches_conserve_mass_and_energy(self):
        def tank(name, pressure, pressurant="Nitrogen"):
            return rsim.Tank(
                name, 0.2, 1.0, 0.0, 1000.0, 60.0, pressure,
                temperature=293.15,
                blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
                polytropic_exponent=1.0,
                propellant="water",
                pressurant=pressurant,
            )

        source = rsim.Tank(
            "shared-N2", 0.2, 1.0, 0.0, 1000.0, 0.0, 5.0e6,
            temperature=293.15,
            blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
            polytropic_exponent=1.0,
            propellant="none",
            pressurant="Nitrogen",
        )
        targets = [tank("a", 1.0e6), tank("b", 1.5e6), tank("closed", 2.0e6)]
        network = rsim.Propulsion("branched-N2")
        network.connect(source, targets[0], rsim.Orifice(2.0e-7),
                        upstream_phase=rsim.FluidPhase.GAS)
        network.connect(source, targets[1], rsim.Orifice(7.0e-8),
                        upstream_phase=rsim.FluidPhase.GAS)
        network.connect(source, targets[2], rsim.Valve(1.0e-6, opening=0.0),
                        upstream_phase=rsim.FluidPhase.GAS)
        network.finalize()

        tanks = [source, *targets]
        masses_before = np.array([item.pressurant_mass for item in tanks])
        energy_before = sum(item.pressurant_internal_energy for item in tanks)
        total_before = network.total_stored_mass
        simulator = rsim.Simulator(0.01)
        simulator.add_model(network)
        simulator.step()

        flow_a = network.connection_mass_flow(0)
        flow_b = network.connection_mass_flow(1)
        self.assertGreater(flow_a, flow_b)
        self.assertGreater(flow_b, 0.0)
        self.assertEqual(network.connection_mass_flow(2), 0.0)
        self.assertAlmostEqual(targets[0].pressurant_mass - masses_before[1],
                               flow_a * 0.01, places=12)
        self.assertAlmostEqual(targets[1].pressurant_mass - masses_before[2],
                               flow_b * 0.01, places=12)
        self.assertAlmostEqual(source.pressurant_mass - masses_before[0],
                               -(flow_a + flow_b) * 0.01, places=12)
        self.assertEqual(targets[2].pressurant_mass, masses_before[3])
        self.assertAlmostEqual(network.total_stored_mass, total_before, places=12)
        self.assertAlmostEqual(
            sum(item.pressurant_internal_energy for item in tanks),
            energy_before,
            places=7,
        )

        reverse_low = tank("reverse-low", 1.0e6)
        reverse_high = tank("reverse-high", 4.0e6)
        reverse = rsim.Propulsion("reverse")
        reverse.connect(reverse_low, reverse_high, rsim.Orifice(1.0e-7),
                        upstream_phase=rsim.FluidPhase.GAS)
        reverse.finalize()
        reverse_total = reverse.total_stored_mass
        low_before = reverse_low.pressurant_mass
        reverse_simulator = rsim.Simulator(0.01)
        reverse_simulator.add_model(reverse)
        reverse_simulator.step()
        self.assertLess(reverse.connection_mass_flow(0), 0.0)
        self.assertGreater(reverse_low.pressurant_mass, low_before)
        self.assertAlmostEqual(reverse.total_stored_mass, reverse_total, places=12)

        with self.assertRaises(ValueError):
            incompatible = rsim.Propulsion("incompatible")
            incompatible.connect(source, tank("helium", 1.0e6, "Helium"),
                                 rsim.Orifice(1.0e-7),
                                 upstream_phase=rsim.FluidPhase.GAS)

        liquid_source = tank("liquid-source", 4.0e6)
        liquid_receiver = rsim.Tank(
            "liquid-receiver", 0.2, 1.0, 0.0, 1000.0, 10.0, 1.0e6,
            temperature=293.15,
            blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
            polytropic_exponent=1.0,
            propellant="water",
            pressurant="Nitrogen",
        )
        liquid_network = rsim.Propulsion("liquid-transfer")
        liquid_network.connect(
            liquid_source, liquid_receiver, rsim.Orifice(1.0e-7),
            upstream_phase=rsim.FluidPhase.LIQUID,
        )
        liquid_network.finalize()
        liquid_total = liquid_network.total_stored_mass
        receiver_before = liquid_receiver.propellant_mass
        liquid_simulator = rsim.Simulator(0.01)
        liquid_simulator.add_model(liquid_network)
        liquid_simulator.step()
        self.assertGreater(liquid_receiver.propellant_mass, receiver_before)
        self.assertAlmostEqual(liquid_network.total_stored_mass,
                               liquid_total, places=12)

    def test_bang_bang_pressure_control_and_independent_branches(self):
        def tank(name, pressure):
            return rsim.Tank(
                name, 0.2, 1.0, 0.0, 1000.0, 60.0, pressure,
                temperature=293.15,
                blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
                polytropic_exponent=1.0,
                propellant="water",
                pressurant="Nitrogen",
            )

        low, high, middle = (
            tank("low", 0.8e6), tank("high", 1.2e6), tank("middle", 1.0e6)
        )
        low_valve = rsim.Valve(1.0e-7, opening=0.0)
        high_valve = rsim.Valve(1.0e-7, opening=1.0)
        memory_open_valve = rsim.Valve(1.0e-7, opening=0.0)
        memory_closed_valve = rsim.Valve(1.0e-7, opening=1.0)
        controllers = [
            rsim.BangBangPressureController(
                "open-low", low, low_valve, 1.0e6, 5.0e4,
                initially_open=False,
            ),
            rsim.BangBangPressureController(
                "close-high", high, high_valve, 1.0e6, 5.0e4,
                initially_open=True,
            ),
            rsim.BangBangPressureController(
                "remember-open", middle, memory_open_valve, 1.0e6, 5.0e4,
                initially_open=True,
            ),
            rsim.BangBangPressureController(
                "remember-closed", middle, memory_closed_valve, 1.0e6, 5.0e4,
                initially_open=False,
            ),
        ]
        simulator = rsim.Simulator(0.01)
        for controller in controllers:
            simulator.add_model(controller)
        simulator.step()
        self.assertTrue(controllers[0].valve_open)
        self.assertFalse(controllers[1].valve_open)
        self.assertTrue(controllers[2].valve_open)
        self.assertFalse(controllers[3].valve_open)

        source = rsim.Tank(
            "shared-source", 0.2, 1.0, 0.0, 1000.0, 0.0, 5.0e6,
            temperature=293.15,
            blowdown_mode=rsim.BlowdownMode.POLYTROPIC,
            polytropic_exponent=1.0,
            propellant="none",
            pressurant="Nitrogen",
        )
        branch_network = rsim.Propulsion("controlled-branches")
        branch_network.connect(source, low, low_valve,
                               upstream_phase=rsim.FluidPhase.GAS)
        branch_network.connect(source, high, high_valve,
                               upstream_phase=rsim.FluidPhase.GAS)
        branch_network.finalize()
        source_before = source.pressurant_mass
        branch_simulator = rsim.Simulator(0.01)
        branch_simulator.add_model(controllers[0])
        branch_simulator.add_model(controllers[1])
        branch_simulator.add_model(branch_network)
        branch_simulator.step()
        self.assertGreater(branch_network.connection_mass_flow(0), 0.0)
        self.assertEqual(branch_network.connection_mass_flow(1), 0.0)
        self.assertAlmostEqual(
            source_before - source.pressurant_mass,
            branch_network.connection_mass_flow(0) * 0.01,
            places=12,
        )


if __name__ == "__main__":
    unittest.main()
