import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
import bondgraph as bg
import sympy as sp


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

    def reduce_text(self, text):
        m = bg.parse_text(text)
        reduction = {}
        states, dots, _, solution = bg.derive(m, bg.causality(m), reduction=reduction)
        self.assertIsNotNone(solution)
        return dict(zip(states, dots)), reduction

    def assertSame(self, a, b):
        self.assertEqual(sp.simplify(sp.sympify(a) - sp.sympify(b, locals={"P": sp.Function("P"), "v": sp.Function("v")})), 0, f"{a} != {b}")

    def test_inertias_across_a_transformer_are_combined(self):
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "fluid_piston.bond"))
        reduction = {}
        states, dots, _, solution = bg.derive(m, bg.causality(m), reduction=reduction)
        self.assertIsNotNone(solution)
        eqs = dict(zip(map(str, states), dots))
        self.assertEqual(list(eqs), ["q_C", "p_If"])  # the original states; m follows If
        self.assertSame(eqs["q_C"], "A*p_If/If")
        self.assertSame(eqs["p_If"], "If*(P(t) - (Rf + A**2*b)*p_If/If - A*k*q_C)/(If + A**2*m)")
        self.assertEqual([text for text, _ in reduction["steps"]],
                         ["reflected b, C and m through TF1 onto J1", "combined m and If on J1 into one inertia"])
        [(x, value)] = reduction["dependent"]
        self.assertEqual(str(x), "p_m")
        self.assertSame(value, "A*m*p_If/If")

    def test_capacitor_across_a_gyrator_becomes_an_inertia(self):
        eqs, reduction = self.reduce_text("junction a 1: V-, L+, R1+, G.1+\njunction b 0: G.2-, Cc+, R2+\n"
                                          "element V Se value=v(t)\nelement G GY n=K\nelement L I I=L\n"
                                          "element R1 R R=R1\nelement Cc C C=Cc\nelement R2 R R=R2\n")
        self.assertEqual(list(map(str, eqs)), ["p_L"])
        self.assertSame(list(eqs.values())[0], "L*(v(t) - (R1 + K**2/R2)*p_L/L)/(L + K**2*Cc)")
        self.assertSame(reduction["dependent"][0][1], "Cc*K*p_L/L")

    def test_parallel_capacitors_are_combined(self):
        eqs, reduction = self.reduce_text("junction n 0: S-, C1+, C2+, R1+\nelement S Sf value=u(t)\n"
                                          "element C1 C C=C1\nelement C2 C C=C2\nelement R1 R R=R\n")
        self.assertSame(eqs[sp.Symbol("q_C1")], "C1*(u(t) - q_C1/(C1*R))/(C1 + C2)")
        self.assertSame(reduction["dependent"][0][1], "C2*q_C1/C1")

    def test_derivative_causality_a_source_forces_is_left_alone(self):
        m = bg.parse_text("junction j 1: S-, L+, R1+\nelement S Sf value=u(t)\nelement L I I=L\nelement R1 R R=R\n")
        self.assertIsNone(bg.reduce_graph(m, bg.causality(m)))

    def test_invalid_graph_is_rejected(self):
        path = Path(__file__).parent / "bad.bond"
        path.write_text("junction j 0: a+, b+\nelement a Se value=1\nelement b Se value=2\n")
        self.addCleanup(path.unlink)
        with self.assertRaises(ValueError):
            bg.causality(bg.parse(str(path)))

    SERIES_RC = "junction a 1: V-, R1+, C1+\nelement V Se value=v(t)\nelement R1 R R=R\nelement C1 C C=C\n"

    def test_outputs_are_parsed(self):
        m = bg.parse_text(self.SERIES_RC + "output i = f_R1\noutput v_C=e_C1 / 1\n")
        self.assertEqual(m.outputs, {"i": "f_R1", "v_C": "e_C1 / 1"})
        with self.assertRaisesRegex(ValueError, "line 2: duplicate output y"):
            bg.parse_text("output y = 1\noutput y = 2\n")

    def test_outputs_in_terms_of_states_and_inputs(self):
        m = bg.parse_text(self.SERIES_RC + "output i = f_R1\noutput v_C = e_C1\noutput P = e_R1*f_R1\n")
        outputs = {}
        states, _, _, _ = bg.derive(m, bg.causality(m), outputs=outputs)
        self.assertSame(outputs["i"], "(v(t) - q_C1/C)/R")
        self.assertSame(outputs["v_C"], "q_C1/C")
        self.assertSame(outputs["P"], "(v(t) - q_C1/C)**2/R")
        self.assertEqual(bg.unsolved(m, states, outputs), [])

    def test_outputs_of_a_reduced_graph_use_its_own_bond_variables(self):
        text = (Path(__file__).parents[1] / "examples" / "fluid_piston.bond").read_text()
        m = bg.parse_text(text + "output v_m = f_m\noutput F_m = e_m\noutput p = p_m\n")
        outputs = {}
        bg.derive(m, bg.causality(m), outputs=outputs)
        self.assertSame(outputs["v_m"], "A*p_If/If")
        # m is in derivative causality: the force on it is m * dv/dt
        self.assertSame(outputs["F_m"], "A*m*(P(t) - (Rf + A**2*b)*p_If/If - A*k*q_C)/(If + A**2*m)")
        self.assertSame(outputs["p"], "A*m*p_If/If")

    def test_latex_outputs(self):
        m = bg.parse_text(self.SERIES_RC + "output i = f_R1\noutput v_C = e_C1\n")
        out = bg.latex_report(m, bg.causality(m))
        self.assertEqual(out["outputs"], [("i", r"\frac{C v{\left(t \right)} - q_{\mathrm{C1}}}{C R}"),
                                          ("v_{C}", r"\frac{q_{\mathrm{C1}}}{C}")])
        self.assertIn(r" \\[1ex]" + "\ni &= ", out["separate"])
        self.assertTrue(out["linear"] and out["outputs_linear"])
        self.assertIn(r"\left[\begin{matrix}i \\ v_{C}\end{matrix}\right] = "
                      r"\left[\begin{matrix}- \frac{1}{C R}\\\frac{1}{C}\end{matrix}\right] \left[\begin{matrix}q_{\mathrm{C1}}\end{matrix}\right] + "
                      r"\left[\begin{matrix}\frac{1}{R}\\0\end{matrix}\right] \left[\begin{matrix}v{\left(t \right)}\end{matrix}\right]", out["matrix"])

    def test_output_with_no_input_term_drops_d(self):
        m = bg.parse_text(self.SERIES_RC + "output v_C = e_C1\n")
        matrix = bg.latex_report(m, bg.causality(m))["matrix"]
        self.assertTrue(matrix.endswith(r"\left[\begin{matrix}v_{C}\end{matrix}\right] = \left[\begin{matrix}\frac{1}{C}\end{matrix}\right] "
                                        r"\left[\begin{matrix}q_{\mathrm{C1}}\end{matrix}\right]" + "\n\\end{gathered}"), matrix)

    def test_unsolved_output_is_reported(self):
        m = bg.parse_text("junction n 0: Q-, T+, O+\nelement Q Sf value=Q(t)\nelement T C C=A\nelement O R e=a*f*Abs(f)\n"
                          "output Qo = f_O\noutput h = q_T/A\n")
        out = bg.latex_report(m, bg.causality(m))
        self.assertEqual(out["unsolved_outputs"], ["Qo"])
        self.assertIsNone(out["matrix"])

    def test_integral_output_is_a_new_state(self):
        m = bg.parse_text("junction a 1: V-, L+, R1+\nelement V Se value=v(t)\nelement L I I=L\nelement R1 R R=R\n"
                          "output q_L = integral(f_L)\noutput pw = q_L*e_R1\noutput w = integral(pw + f_L)\n")
        outputs = {}
        states, dots, _, solution = bg.derive(m, bg.causality(m), outputs=outputs)
        self.assertIsNotNone(solution)
        self.assertEqual([str(s) for s in states], ["p_L", "q_L", "w"])
        self.assertSame(dots[1], "p_L/L")
        self.assertSame(dots[2], "R*p_L*q_L/L + p_L/L")  # pw is substituted; q_L is the new state
        self.assertEqual(list(outputs), ["pw"])
        self.assertSame(outputs["pw"], "R*p_L*q_L/L")

    def test_integral_output_in_latex(self):
        m = bg.parse_text(self.SERIES_RC + "output q = integral(f_C1)\n")
        out = bg.latex_report(m, bg.causality(m))
        self.assertEqual([p[0] for p in out["equations"]], [r"\dot{q}_{\mathrm{C1}}", r"\dot{q}"])
        self.assertTrue(out["linear"] and out["outputs"] == [])

    def test_output_mistakes_are_rejected(self):
        for text, why in (("output a = b\noutput b = a\n", "depend on each other"),
                          ("output a = 2*integral(f_R1)\n", "whole output on its own")):
            m = bg.parse_text("junction a 1: V-, R1+\nelement V Se value=1\nelement R1 R R=1\n" + text)
            with self.assertRaisesRegex(ValueError, why):
                bg.derive(m, bg.causality(m), outputs={})


if __name__ == "__main__":
    unittest.main()
