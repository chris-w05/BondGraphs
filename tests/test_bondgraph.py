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

    TWO_MASSES = "junction v1 1: F-, m1+, S.1+\njunction v2 1: m2+, S.2+\nelement F Se value=F(t)\nelement m1 I I=m1\nelement m2 I I=m2\n"

    def field_equations(self, text):
        m = bg.parse_text(text)
        states, dots, _, solution = bg.derive(m, bg.causality(m))
        self.assertIsNotNone(solution)
        return dict(zip(map(str, states), dots))

    def test_c_field_couples_its_ports(self):
        eqs = self.field_equations(self.TWO_MASSES + "element S Cfield e1=k1*q1+k12*(q1-q2) e2=k2*q2+k12*(q2-q1)\n")
        self.assertEqual(list(eqs), ["p_m1", "p_m2", "q_S", "q_S_2"])
        self.assertSame(eqs["p_m1"], "F(t) - k1*q_S - k12*(q_S - q_S_2)")
        self.assertSame(eqs["p_m2"], "-k2*q_S_2 - k12*(q_S_2 - q_S)")
        self.assertSame(eqs["q_S_2"], "p_m2/m2")
        # Written the other way round, q = C e, the field is inverted to the same equations
        eqs = self.field_equations(self.TWO_MASSES + "element S Cfield q1=(k2*e1+k12*(e1+e2))/(k1*k2+k12*(k1+k2)) "
                                   "q2=(k1*e2+k12*(e1+e2))/(k1*k2+k12*(k1+k2))\n")
        self.assertSame(eqs["p_m1"], "F(t) - k1*q_S - k12*(q_S - q_S_2)")

    def test_i_field_given_as_p_of_f_is_inverted(self):
        eqs = self.field_equations("junction a 1: F-, B.1+, k1+\njunction b 1: B.2+, k2+\nelement F Se value=F(t)\n"
                                   "element B Ifield p1=M11*f1+M12*f2 p2=M12*f1+M22*f2\nelement k1 C C=1/k1\nelement k2 C C=1/k2\n")
        M = sp.Matrix(sp.symbols("M11 M12 M12 M22")).reshape(2, 2)
        v = M.inv() * sp.Matrix(sp.symbols("p_B p_B_2"))
        self.assertEqual(sp.simplify(eqs["q_k1"] - v[0]), 0)
        self.assertEqual(sp.simplify(eqs["q_k2"] - v[1]), 0)
        self.assertSame(eqs["p_B"], "F(t) - k1*q_k1")

    def test_r_field_in_either_form(self):
        caps = "junction n1 0: C1+, G.1+\njunction n2 0: C2+, G.2+\nelement C1 C C=C1\nelement C2 C C=C2\n"
        eqs = self.field_equations(caps + "element G Rfield f1=G11*e1+G12*e2 f2=G12*e1+G22*e2\n")
        self.assertSame(eqs["q_C1"], "-(G11*q_C1/C1 + G12*q_C2/C2)")
        eqs = self.field_equations(caps + "element G Rfield e1=R1*f1+Rm*f2 e2=Rm*f1+R2*f2\n")
        flows = sp.Matrix([[sp.Symbol("R1"), sp.Symbol("Rm")], [sp.Symbol("Rm"), sp.Symbol("R2")]]).inv() * sp.Matrix(sp.sympify(["q_C1/C1", "q_C2/C2"]))
        self.assertEqual(sp.simplify(eqs["q_C1"] + flows[0]), 0)

    def test_ic_field_ports_take_their_type_from_their_laws(self):
        m = bg.parse_text("junction j1 1: S1-, X.1+\njunction j2 1: S2-, X.2+\nelement S1 Sf value=u1(t)\nelement S2 Se value=u2(t)\n"
                          "element X ICfield e1=k*q1+a*p2 f2=p2/m+a*q1\noutput e = e_X\n")
        outputs = {}
        states, dots, _, _ = bg.derive(m, bg.causality(m), outputs=outputs)
        self.assertEqual([str(s) for s in states], ["q_X", "p_X_2"])
        self.assertEqual([str(d) for d in dots], ["u1(t)", "u2(t)"])
        self.assertSame(outputs["e"], "k*q_X + a*p_X_2")

    def test_editor_causality_covers_field_ports(self):
        m = bg.parse_text(self.TWO_MASSES + "element S Cfield e1=k*q1 e2=k*q2\n")
        states, _, _, _ = bg.derive(m, bg.storage_causality(m, {"S.2": False}))
        self.assertEqual([str(s) for s in states], ["p_m1", "p_m2", "q_S"])

    def test_storage_a_source_forces_follows_it(self):
        # Port 1 of the field and the inductor are fixed by flow sources: each is
        # solved for, and reported as following, rather than leaving the rest unsolved.
        m = bg.parse_text("junction a 1: S-, X.1+\njunction b 1: V-, X.2+\nelement S Sf value=sin(t)\nelement V Se value=1\n"
                          "element X Ifield p1=m1*f1 f2=p2/m2\n")
        reduction = {}
        states, dots, _, solution = bg.derive(m, bg.causality(m), reduction=reduction)
        self.assertIsNotNone(solution)
        self.assertEqual((list(map(str, states)), dots), (["p_X_2"], [1]))
        self.assertEqual([(str(x), str(y)) for x, y in reduction["dependent"]], [("p_X", "m1*sin(t)")])
        m = bg.parse_text("junction j 1: S-, L+, R1+\nelement S Sf value=u(t)\nelement L I I=L\nelement R1 R R=R\noutput v = e_R1\n")
        outputs = {}
        bg.derive(m, bg.causality(m), outputs=outputs)
        self.assertSame(outputs["v"], "R*u(t)")
        self.assertEqual(bg.latex_report(m, bg.causality(m))["dependent"], [r"p_{L} = L u{\left(t \right)}"])

    HARD_SPRING = ("junction v 1: F-, m+, k+, b+\nelement F Se value=F(t)\nelement m I I=m\nelement k C e=k*q+k3*q**3\n"
                   "element b R R=b\noutput x = q_k\noutput v = f_m\n")

    def test_power_variables_replace_the_states(self):
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "dc_motor.bond"))
        v = bg.view(m, bg.causality(m), "power")
        self.assertEqual([str(x) for x in v["states"]], ["f_La", "f_J"])
        u = v["columns"][0]["slot"]
        self.assertSame(v["dots"][0].xreplace({u: sp.Symbol("V")}), "(V - R_a*f_La - K_t*f_J)/L_a")
        self.assertSame(v["dots"][1], "(K_t*f_La - b*f_J)/J")
        # After a reduction, storage that follows the states does so in power variables too
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "fluid_piston.bond"))
        v = bg.view(m, bg.causality(m), "power")
        self.assertEqual([(str(x), str(y)) for x, y in v["dependent"]], [("f_m", "A*f_If")])
        self.assertSame(v["dots"][0], "A*k*f_If")

    def test_lagrange_variables(self):
        # Displacements q stay and each I's flow replaces its p: m x'' = F - b x' - k x - k3 x^3,
        # with no storage law inverted.
        m = bg.parse_text(self.HARD_SPRING)
        c = bg.causality(m)
        v = bg.view(m, c, "lagrange")
        self.assertEqual(([str(x) for x in v["states"]], v["note"], v["hint"]), (["f_m", "q_k"], "", ""))
        inputs = {col["slot"]: col["value"] for col in v["columns"]}
        self.assertSame(v["dots"][0].xreplace(inputs), "(F(t) - b*f_m - k*q_k - k3*q_k**3)/m")
        self.assertEqual(str(v["dots"][1]), "f_m")
        self.assertEqual([(name, str(y)) for name, y in v["outputs"]], [("x", "q_k"), ("v", "f_m")])
        out = bg.latex_report(m, c, "lagrange")
        self.assertIn(r"\dot{f}_{m} &= ", out["separate"])
        # Linearized: dδf_m/dt = (δF - b δf_m - (k + 3 k3 q_k‾^2) δq_k)/m
        v = bg.view(m, c, "lagrange", linearize=True)
        self.assertEqual([str(x) for x in v["states"]], ["δf_m", "δq_k"])
        d = {sp.Symbol("δF"): sp.Symbol("dF"), sp.Symbol("δf_m"): sp.Symbol("df"), sp.Symbol("δq_k"): sp.Symbol("dq"), sp.Symbol("q_k‾"): sp.Symbol("qb")}
        self.assertSame(v["dots"][0].xreplace(d), "(dF - b*df - (k + 3*k3*qb**2)*dq)/m")

    def test_linear_approximation(self):
        m = bg.parse_text(self.HARD_SPRING)
        c = bg.causality(m)
        self.assertTrue(bg.view(m, c)["nonlinear"])
        v = bg.view(m, c, linearize=True)
        self.assertEqual([str(x) for x in v["states"]], ["δp_m", "δq_k"])
        bar = {sp.Symbol("q_k‾"): sp.Symbol("qb")}
        self.assertSame(v["dots"][0].xreplace(bar).xreplace({sp.Symbol("δF"): sp.Symbol("dF"), sp.Symbol("δp_m"): sp.Symbol("dp"), sp.Symbol("δq_k"): sp.Symbol("dq")}),
                        "dF - b*dp/m - (k + 3*k3*qb**2)*dq")
        self.assertEqual([str(f) for f in v["operating"]], ["F‾ - b*p_m‾/m - k*q_k‾ - k3*q_k‾**3", "p_m‾/m"])
        self.assertEqual([name for name, _ in v["outputs"]], ["δx", "δv"])
        # In power variables: dδe_k/dt = (k + 3 k3 q‾^2) δf_m
        Q = sp.Symbol("q_k‾")
        v = bg.view(m, c, "power", linearize=True)
        self.assertSame(v["dots"][1].xreplace({sp.Symbol("δf_m"): sp.Symbol("df"), Q: sp.Symbol("qb")}), "(k + 3*k3*qb**2)*df")
        out = bg.latex_report(m, c, linearize=True)
        self.assertIn(r"\delta \dot{p}_{m} &= ", out["separate"])
        self.assertIn(r"0 &= \bar{F} - ", out["operating"])

    def test_views_of_a_linear_model(self):
        m = bg.parse(str(Path(__file__).parents[1] / "examples" / "dc_motor.bond"))
        c = bg.causality(m)
        v = bg.view(m, c, linearize=True)
        self.assertFalse(v["nonlinear"] or v["linearized"])  # a linear model is its own approximation
        self.assertEqual(bg.latex_report(m, c), bg.latex_report(m, c, "energy", True))

    def test_power_variables_of_a_nonlinear_law_take_its_inverse(self):
        m = bg.parse_text(self.HARD_SPRING)
        v = bg.view(m, bg.causality(m), "power")
        self.assertEqual([str(x) for x in v["states"]], ["f_m", "e_k"])
        self.assertIn("linear approximation is simpler", v["hint"])

    def test_power_variables_that_depend_on_an_input(self):
        # An electrostatic actuator driven by a voltage source: port 1 of the C-field
        # is in derivative causality, so its charge q1 = Ae*u/(l - q2) follows the
        # source, and port 2's force e2 = -q1^2/(2*Ae) depends on the input too.
        m = bg.parse_text("junction J1 1: CF1.2+, Se1+, I1-\njunction bond5 1: F-, CF1.1+\nelement F Se value=F0*Heaviside(t)\n"
                          "element CF1 Cfield e1=q1*(l-q2)/(Ae) e2=-q1^2/(2*Ae)\nelement I1 I I=I_I1\nelement Se1 Se value=1\n")
        c = bg.causality(m)
        v = bg.view(m, c, "power")
        self.assertEqual((v["note"], [str(x) for x in v["states"]]), ("", ["e_CF1_2", "f_I1"]))
        self.assertIn("rate of change", v["hint"])
        # de2/dt = de2/dq2 * dq2/dt + de2/du * du/dt, written with the inverse that keeps
        # the gap l - q2 positive. Check it at a point, from q2 and u.
        l, Ae, q, u, du, f = 0.8, 0.6, 0.1, 0.3, 0.7, 0.2
        e2 = -Ae * u**2 / (2 * (l - q)**2)
        expected = -Ae * u**2 / (l - q)**3 * f - Ae * u / (l - q)**2 * du
        slot, rate = (next(col["slot"] for col in v["columns"] if col["name"] == name) for name in ("u_F", "u_F˙"))
        point = {sp.Symbol("l"): l, sp.Symbol("Ae"): Ae, slot: u, rate: du, sp.Symbol("e_CF1_2"): e2, sp.Symbol("f_I1"): f}
        self.assertAlmostEqual(complex(v["dots"][0].xreplace(point).evalf()).real, expected)
        out = bg.latex_report(m, c, "power")
        self.assertIn(r"\delta\left(t\right)", out["separate"])  # the step input's rate, an impulse
        # Linearized, the input's rate is one more input: dδe2/dt = ... - Ae*u‾/(l - q2‾)^2 * δu_F˙
        v = bg.view(m, c, "power", linearize=True)
        self.assertEqual(v["note"], "")
        self.assertIn(sp.Symbol("δu_F˙"), v["dots"][0].free_symbols)
        self.assertEqual(bg.tex_name("δu_F˙"), r"\delta \dot{u}_{F}")

    def test_advanced_examples_solve(self):
        for path in sorted((Path(__file__).parents[1] / "examples").glob("*.bond")):
            if "#@title" not in path.read_text():
                continue
            with self.subTest(path.name):
                m = bg.parse(str(path))
                c = bg.causality(m)
                outputs = {}
                _, _, _, solution = bg.derive(m, c, outputs=outputs)
                self.assertIsNotNone(solution)
                self.assertEqual(bg.unsolved(m, bg.derive(m, c)[0], outputs), [])
                self.assertTrue(bg.view(m, c, "lagrange", linearize=True)["note"] == "")

    def test_names_of_deviations_and_operating_points(self):
        self.assertEqual(bg.tex_name("δq_C1", dot=True), r"\delta \dot{q}_{\mathrm{C1}}")
        self.assertEqual(bg.tex_name("omega‾"), r"\bar{\omega}")

    def test_field_mistakes_are_rejected(self):
        base = self.TWO_MASSES
        for law, why in (("e1=q1 f2=q2", r"Cfield takes a law per port, like e1=, q1="),
                         ("e1=q1", r"port 2 needs a law \(e2= or q2=\)"),
                         ("e1=q1 e2=q2 e3=q1", "has a law for port 3"),
                         ("e1=q1 q1=e1 e2=q2", "port 1 has two laws"),
                         ("e1=q1 e2=p2", "uses p2, which isn't one of its variables")):
            with self.assertRaisesRegex(ValueError, why):
                bg.parse_text(base + f"element S Cfield {law}\n")
        with self.assertRaisesRegex(ValueError, "port 2 isn't attached"):
            bg.parse_text("junction a 1: F-, S.1+, S.3+\nelement F Se value=1\nelement S Rfield e1=f1 e2=f2 e3=f3\n")
        with self.assertRaisesRegex(ValueError, "needs a port number, like S.1"):
            bg.parse_text("junction a 1: F-, S+\nelement F Se value=1\nelement S Rfield e1=f1\n")


if __name__ == "__main__":
    unittest.main()
