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

    def derive_text(self, text):
        path = Path(__file__).parent / "scratch.bond"
        path.write_text(text)
        self.addCleanup(path.unlink, missing_ok=True)
        m = bg.parse(str(path))
        states, dots, _, solution = bg.derive(m, bg.causality(m))
        return dict(zip(map(str, states), map(str, dots))), solution

    def test_resistor_law_in_either_direction(self):
        for law in ("R=R", "e=R*f", "f=e/R"):
            for signs in ("Is-, C1+, R1+", "Is+, C1-, R1-"):
                eqs, _ = self.derive_text(f"junction node 0: {signs}\nelement Is Sf value=u(t)\n"
                                          f"element C1 C C=C\nelement R1 R {law}\n")
                self.assertEqual(eqs["q_C1"], "u(t) - q_C1/(C*R)", f"{law} with {signs}")

    def test_diode_written_as_flow_of_effort_is_solved(self):
        eqs, solution = self.derive_text("junction l 1: V-, D+, C1+\nelement V Se value=v(t)\n"
                                         "element D R f=Piecewise((e/R_on,e>0),(0,True))\nelement C1 C C=C\n")
        self.assertIsNotNone(solution)
        self.assertEqual(eqs["q_C1"], "Piecewise(((C*v(t) - q_C1)/(C*R_on), v(t) - q_C1/C > 0), (0, True))")

    def test_storage_laws_as_functions(self):
        eqs, _ = self.derive_text("junction v 1: F-, m+, k+\nelement F Se value=F(t)\n"
                                  "element m I f=p/m\nelement k C e=k*q**3\n")
        self.assertEqual(eqs, {"p_m": "-k*q_k**3 + F(t)", "q_k": "p_m/m"})
        eqs, _ = self.derive_text("junction v 1: F-, m+, k+\nelement F Se value=F(t)\n"
                                  "element m I p=m*f\nelement k C q=e/k\n")
        self.assertEqual(eqs, {"p_m": "-k*q_k + F(t)", "q_k": "p_m/m"})

    def test_inverted_law_keeps_the_real_branch(self):
        eqs, _ = self.derive_text("junction a 1: V-, L+, R1+\nelement V Se value=v(t)\n"
                                  "element L I p=Psi*tanh(f/i0)\nelement R1 R R=R\n")
        self.assertEqual(eqs["p_L"], "-R*i0*log((Psi + p_L)/(Psi - p_L))/2 + v(t)")

    def test_law_that_cannot_be_inverted_is_reported(self):
        eqs, solution = self.derive_text("junction n 0: Q-, T+, O+\nelement Q Sf value=Q(t)\n"
                                         "element T C C=A\nelement O R e=a*f*Abs(f)\n")
        self.assertIsNone(solution)
        self.assertIn("f_O", eqs["q_T"])

    def test_law_keys_are_checked(self):
        with self.assertRaisesRegex(ValueError, r"two laws \(R= and f=\)"):
            self.derive_text("junction n 0: R1+\nelement R1 R R=1 f=e\n")
        with self.assertRaisesRegex(ValueError, "C takes C=, e=, q= not f="):
            self.derive_text("junction n 0: C1+\nelement C1 C f=e\n")

    def test_latex_separate_and_matrix(self):
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "dc_motor.bond"))
        out = bg.latex_report(m, bg.causality(m))
        self.assertTrue(out["solved"] and out["linear"])
        self.assertEqual(out["equations"][0], (r"\dot{p}_{\mathrm{La}}",
                                               r"v{\left(t \right)} - \frac{R_{a} p_{\mathrm{La}}}{L_{a}} - \frac{K_{t} p_{J}}{J}"))
        self.assertEqual(out["matrix"], (
            r"\left[\begin{matrix}\dot{p}_{\mathrm{La}} \\ \dot{p}_{J}\end{matrix}\right] = "
            r"\left[\begin{matrix}- \frac{R_{a}}{L_{a}} & - \frac{K_{t}}{J}\\\frac{K_{t}}{L_{a}} & - \frac{b}{J}\end{matrix}\right] "
            r"\left[\begin{matrix}p_{\mathrm{La}} \\ p_{J}\end{matrix}\right] + "
            r"\left[\begin{matrix}1\\0\end{matrix}\right] \left[\begin{matrix}v{\left(t \right)}\end{matrix}\right]"))

    def test_matrix_form_names_expression_inputs(self):
        m = bg.parse_text("junction v 1: F-, m+, k+, b+\nelement F Se value=F0*Heaviside(t-1)\n"
                          "element m I I=m\nelement k C C=1/k\nelement b R R=b\n")
        out = bg.latex_report(m, bg.causality(m))
        self.assertEqual(out["inputs"], [r"u_{F} = F_{0} \theta\left(t - 1\right)"])
        self.assertIn(r"\left[\begin{matrix}u_{F}\end{matrix}\right]", out["matrix"])
        self.assertIn(r"\text{where } u_{F} = ", out["matrix"])

    def test_nonlinear_matrix_form_is_a_vector(self):
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "piecewise_leak.bond"))
        out = bg.latex_report(m, bg.causality(m))
        self.assertTrue(out["solved"])
        self.assertFalse(out["linear"])
        self.assertIn(r"\begin{cases}", out["matrix"])

    def test_latex_names(self):
        self.assertEqual(bg.tex_name("p_I_Electrical"), r"p_{\mathrm{I\_Electrical}}")
        self.assertEqual(bg.tex_name("q_C1", dot=True), r"\dot{q}_{\mathrm{C1}}")
        self.assertEqual([bg.tex_name(n) for n in ("R_a", "F0", "rho", "omega_n", "eps")],
                         ["R_{a}", "F_{0}", r"\rho", r"\omega_{n}", r"\mathit{eps}"])

    def test_editor_causality_picks_the_states(self):
        m = bg.parse_text("junction n 0: S-, C1+, C2+, R1+\nelement S Sf value=u(t)\n"
                          "element C1 C C=C1\nelement C2 C C=C2\nelement R1 R R=R\n")
        states, _, _, _ = bg.derive(m, bg.storage_causality(m, {"C1": False, "C2": True}))
        self.assertEqual([str(s) for s in states], ["q_C2"])

    def test_invalid_graph_is_rejected(self):
        path = Path(__file__).parent / "bad.bond"
        path.write_text("junction j 0: a+, b+\nelement a Se value=1\nelement b Se value=2\n")
        self.addCleanup(path.unlink)
        with self.assertRaises(ValueError):
            bg.causality(bg.parse(str(path)))


if __name__ == "__main__":
    unittest.main()
