#!/usr/bin/env python3
"""Assign causality and derive state equations for the project's .bond format."""
from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import sympy as sp
except ImportError:  # helpful instead of a traceback for a fresh checkout
    sys.exit("This program needs SymPy. Run: python3 -m pip install -r requirements.txt")


@dataclass
class Element:
    name: str
    kind: str
    params: dict[str, str] = field(default_factory=dict)
    ports: list[tuple[str, int]] = field(default_factory=list)  # junction, sign


@dataclass
class Model:
    junctions: dict[str, str] = field(default_factory=dict)
    attachments: dict[str, list[tuple[str, int, int]]] = field(default_factory=dict)
    elements: dict[str, Element] = field(default_factory=dict)


def parse(path: str) -> Model:
    m = Model()
    pending: list[tuple[str, list[tuple[str, int]]]] = []
    connects: list[tuple[str, str, str]] = []
    for number, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        jm = re.fullmatch(r"junction\s+(\w+)\s+([01])\s*:\s*(.*)", line, re.I)
        em = re.fullmatch(r"element\s+(\w+)\s+(I|C|R|Se|Sf|TF|GY)(?:\s+(.*))?", line, re.I)
        cm = re.fullmatch(r"connect\s+(\w+)\s+(\w+)\s+(\w+)", line, re.I)
        if jm:
            name, kind, members = jm.groups()
            if name in m.junctions:
                raise ValueError(f"line {number}: duplicate junction {name}")
            m.junctions[name] = kind
            parsed = []
            for token in filter(None, map(str.strip, members.split(","))):
                mm = re.fullmatch(r"(\w+)(?:\.(1|2))?\s*([+-])?", token)
                if not mm:
                    raise ValueError(f"line {number}: bad attachment {token!r}")
                ename, port, sign = mm.groups()
                parsed.append((ename + ("." + port if port else ""), 1 if sign != "-" else -1))
            pending.append((name, parsed))
        elif em:
            name, kind, rest = em.groups()
            if name in m.elements:
                raise ValueError(f"line {number}: duplicate element {name}")
            params = {}
            for word in (rest or "").split():
                if "=" not in word:
                    raise ValueError(f"line {number}: expected key=value, got {word!r}")
                key, value = word.split("=", 1)
                params[key] = value
            m.elements[name] = Element(name, kind, params)
        elif cm:
            connects.append(cm.groups())
        else:
            raise ValueError(f"line {number}: expected junction, element, or connect")
    # create one-port connections from junction member lists
    for jname, members in pending:
        for token, sign in members:
            bits = token.split(".")
            ename = bits[0]
            if ename not in m.elements:
                raise ValueError(f"junction {jname}: unknown element {ename}")
            elem = m.elements[ename]
            port = int(bits[1]) if len(bits) == 2 else 1
            if elem.kind.upper() in {"TF", "GY"}:
                if len(bits) != 2:
                    raise ValueError(f"junction {jname}: two-port {ename} needs .1 or .2")
            elif len(bits) != 1:
                raise ValueError(f"junction {jname}: one-port {ename} cannot have a port suffix")
            while len(elem.ports) < port:
                elem.ports.append(("", 1))
            if elem.ports[port - 1][0]:
                raise ValueError(f"element {ename}.{port} attached twice")
            elem.ports[port - 1] = (jname, sign)
    for ename, j1, j2 in connects:
        if ename not in m.elements or m.elements[ename].kind.upper() not in {"TF", "GY"}:
            raise ValueError(f"connect {ename}: must name a TF or GY")
        elem = m.elements[ename]
        if elem.ports:
            raise ValueError(f"connect {ename}: use either connect or junction attachments, not both")
        if j1 not in m.junctions or j2 not in m.junctions:
            raise ValueError(f"connect {ename}: unknown junction")
        elem.ports = [(j1, 1), (j2, -1)]  # power flows j1 -> element -> j2
    for elem in m.elements.values():
        required = 2 if elem.kind.upper() in {"TF", "GY"} else 1
        if len(elem.ports) != required or any(not p[0] for p in elem.ports):
            raise ValueError(f"element {elem.name}: needs {required} attached port(s)")
    for elem in m.elements.values():
        for p, (j, sign) in enumerate(elem.ports, 1):
            m.attachments.setdefault(j, []).append((elem.name, p, sign))
    return m


def port_name(name: str, port: int) -> str:
    return f"{name}_{port}" if port != 1 else name


def causality(m: Model):
    """Return effort/flow produced onto each junction, maximizing integral storage."""
    ports = [(e.name, p) for e in m.elements.values() for p in range(1, len(e.ports) + 1)]
    fixed: dict[tuple[str, int], str] = {}
    preferred: dict[tuple[str, int], str] = {}
    pair_rules = []
    for e in m.elements.values():
        k = e.kind.upper()
        key = (e.name, 1)
        if k == "SE": fixed[key] = "e"
        elif k == "SF": fixed[key] = "f"
        elif k == "I": preferred[key] = "f"
        elif k == "C": preferred[key] = "e"
        elif k == "R" and e.params.get("causality"):
            mode = e.params["causality"].lower()
            if mode not in {"effort", "flow"}:
                raise ValueError(f"element {e.name}: R causality must be effort or flow")
            fixed[key] = "e" if mode == "effort" else "f"
        elif k in {"TF", "GY"}:
            pair_rules.append(((e.name, 1), (e.name, 2), k == "GY"))  # GY same, TF opposite
    choices = [[fixed[p]] if p in fixed else ["e", "f"] for p in ports]
    best, best_score = None, -1
    # Kept deliberately bounded: an oversized model needs a proper constraint solver.
    if len(ports) > 22:
        raise ValueError("model has over 22 ports; split it or add a constraint-solver backend")
    for vals in itertools.product(*choices):
        result = dict(zip(ports, vals))
        if any((result[a] == result[b]) != same for a, b, same in pair_rules):
            continue
        valid = True
        for j, kind in m.junctions.items():
            required = "e" if kind == "0" else "f"
            members = m.attachments.get(j, [])
            if not members or sum(result[(n, p)] == required for n, p, _ in members) != 1:
                valid = False; break
        if not valid:
            continue
        score = sum(result[p] == wanted for p, wanted in preferred.items())
        if score > best_score:
            best, best_score = result, score
    if best is None:
        raise ValueError("no valid causality: check source/non-invertible constraints and junction connections")
    return best


def sym(value: str):
    # Give names used like u(t) function semantics; all other undeclared names
    # remain symbolic parameters, which is convenient for R, C, I, and n.
    # SymPy functions (Piecewise, Abs, Max, sin, sqrt, ...) keep their meaning.
    functions = {
        name: sp.Function(name)
        for name in re.findall(r"\b([A-Za-z_]\w*)\s*\(", value)
        if name not in {"t", "sqrt"} and not isinstance(getattr(sp, name, None), sp.FunctionClass)
    }
    return sp.sympify(value.replace("^", "**"), locals={"t": sp.Symbol("t"), **functions})


def derive(m: Model, c):
    e, f, fin = {}, {}, {}
    equations, states, dots = [], [], []
    for el in m.elements.values():
        for p, (_, sign) in enumerate(el.ports, 1):
            suffix = port_name(el.name, p)
            e[el.name, p], f[el.name, p] = sp.symbols(f"e_{suffix} f_{suffix}")
            # A bond's sign only sets the reference direction of its flow; element
            # laws use the flow *into* the element, so the physics never depends on it.
            fin[el.name, p] = sign * f[el.name, p]
    # Junction equations
    for j, kind in m.junctions.items():
        a = m.attachments[j]
        n0, p0, _ = a[0]
        base = e[n0, p0] if kind == "0" else f[n0, p0]
        for n, p, _ in a[1:]:
            equations.append((e[n, p] if kind == "0" else f[n, p]) - base)
        summed = sum(sign * (f[n, p] if kind == "0" else e[n, p]) for n, p, sign in a)
        equations.append(summed)
    # Element laws and storage state equations
    for el in m.elements.values():
        k = el.kind.upper(); key = (el.name, 1)
        if k == "I":
            I = sym(el.params.get("I", "I_" + el.name)); p = sp.Symbol("p_" + el.name)
            equations.append(fin[key] - p / I)
            if c[key] == "f":
                states.append(p); dots.append(e[key])
        elif k == "C":
            C = sym(el.params.get("C", "C_" + el.name)); q = sp.Symbol("q_" + el.name)
            equations.append(e[key] - q / C)
            if c[key] == "e":
                states.append(q); dots.append(fin[key])
        elif k == "R":
            R = sym(el.params.get("R", "R_" + el.name)); equations.append(e[key] - R * fin[key])
        elif k == "SE": equations.append(e[key] - sym(el.params.get("value", "u_" + el.name)))
        elif k == "SF": equations.append(-fin[key] - sym(el.params.get("value", "u_" + el.name)))  # flow delivered
        # Two-ports: power enters port 1 and leaves port 2 (f2 out = -fin[2]).
        elif k == "TF":
            n = sym(el.params.get("n", "n_" + el.name)); equations += [e[el.name,1] - n*e[el.name,2], -fin[el.name,2] - n*fin[el.name,1]]
        elif k == "GY":
            n = sym(el.params.get("n", "n_" + el.name)); equations += [e[el.name,1] + n*fin[el.name,2], e[el.name,2] - n*fin[el.name,1]]
    variables = list(e.values()) + list(f.values())
    # Piecewise element laws (check valves, piecewise resistors) are still linear in
    # the unknowns, but sp.solve evaluates their conditions during elimination and
    # can collapse to 0/0 -> nan. Swap each for a placeholder, solve, then put it back.
    placeholders: dict[sp.Expr, sp.Expr] = {}

    def linearize(eq):
        for pw in eq.atoms(sp.Piecewise):
            placeholders.setdefault(pw, sp.Dummy(f"Piecewise{len(placeholders)}"))
        return eq.xreplace(placeholders)

    solution = sp.solve([linearize(eq) for eq in equations], variables, dict=True, simplify=True)
    if not solution:
        return states, dots, equations, None
    sol = solution[0]
    if placeholders:
        restore = {dummy: pw for pw, dummy in placeholders.items()}
        sol = {var: sp.simplify(val.xreplace(restore)) for var, val in sol.items()}
    return states, [sp.simplify(x.subs(sol)) for x in dots], equations, sol


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file", help="input .bond file")
    ap.add_argument("--show-equations", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    try:
        m = parse(args.file); c = causality(m); states, dots, eqs, sol = derive(m, c)
    except (OSError, ValueError, sp.SympifyError) as exc:
        sys.exit(f"Error: {exc}")
    report = {"causality": {f"{n}.{p}": v for (n,p),v in c.items()}, "states": [str(x) for x in states], "state_equations": [f"d({x})/dt = {y}" for x,y in zip(states,dots)]}
    if args.json:
        print(json.dumps(report, indent=2)); return
    print("Causality (element port -> quantity caused onto junction):")
    for key, value in report["causality"].items(): print(f"  {key}: {value}ffort" if value == "e" else f"  {key}: flow")
    print("\nStates:", ", ".join(report["states"]) or "none (all storage has derivative causality)")
    print("\nState equations:")
    for line in report["state_equations"]: print("  " + line)
    if sol is None:
        print("  Unable to eliminate algebraic variables; full equations follow.")
    if args.show_equations or sol is None:
        print("\nFull equations (each equals zero):")
        for equation in eqs: print("  " + str(equation))


if __name__ == "__main__":
    main()
