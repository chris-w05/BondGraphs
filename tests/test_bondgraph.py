import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
import bondgraph as bg


class BondGraphTests(unittest.TestCase):
    def test_rc_has_capacitor_state(self):
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "rc.bond"))
        c = bg.causality(m)
        states, dots, _, solution = bg.derive(m, c)
        self.assertEqual([str(s) for s in states], ["q_C1"])
        self.assertIsNotNone(solution)
        self.assertEqual(str(dots[0]), "u(t) - q_C1/(C*R)")

    def test_bond_signs_do_not_change_physics(self):
        path = Path(__file__).parent / "flipped.bond"
        path.write_text("junction node 0: Is+, C1-, R1-\n"
                        "element Is Sf value=u(t)\nelement C1 C C=C\nelement R1 R R=R\n")
        self.addCleanup(path.unlink)
        m = bg.parse(str(path))
        states, dots, _, _ = bg.derive(m, bg.causality(m))
        self.assertEqual(str(dots[0]), "u(t) - q_C1/(C*R)")

    def test_gyrator_motor_conserves_power(self):
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "dc_motor.bond"))
        states, dots, _, _ = bg.derive(m, bg.causality(m))
        eqs = dict(zip(map(str, states), map(str, dots)))
        self.assertEqual(eqs["p_La"], "v(t) - R_a*p_La/L_a - K_t*p_J/J")
        self.assertEqual(eqs["p_J"], "K_t*p_La/L_a - b*p_J/J")  # torque = +K_t*i

    def test_invalid_graph_is_rejected(self):
        path = Path(__file__).parent / "bad.bond"
        path.write_text("junction j 0: a+, b+\nelement a Se value=1\nelement b Se value=2\n")
        self.addCleanup(path.unlink)
        with self.assertRaises(ValueError):
            bg.causality(bg.parse(str(path)))


if __name__ == "__main__":
    unittest.main()
