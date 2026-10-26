#!/usr/bin/env python3
"""Assign causality and derive state equations for the project's .bond format."""
from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import sympy as sp
    from sympy.core.function import AppliedUndef
except ImportError:  # helpful instead of a traceback for a fresh checkout
    sys.exit("This program needs SymPy. Run: python3 -m pip install -r requirements.txt")


# Law keys per element kind. The first is the default: a coefficient (R, C, I, n)
# or a source value. R, C and I can instead give the law as a function, keyed by
# the variable it sets: R e=phi(f) or f=phi(e); C e=phi(q) or q=phi(e); I f=phi(p)
# or p=phi(f). Inside any R, C or I law, e, f, q and p are the element's own
# effort, inflow and state.
LAW_FORMS = {"R": ("R", "e", "f"), "C": ("C", "e", "q"), "I": ("I", "f", "p"),
             "SE": ("value",), "SF": ("value",), "TF": ("n",), "GY": ("n",)}

# Multiport fields have any number of ports and one law per port, keyed by the
# variable it sets and the port: e1=phi(...), q2=phi(...). These are the letters
# each kind takes. In a C, I or IC field the letter also sets the port's type: e or
# q makes a C port (state q), f or p an I port (state p). In a field's laws, e1,
# f1, q1, p1, ... are its own effort, inflow and state at each port.
FIELDS = {"CFIELD": "eq", "IFIELD": "fp", "RFIELD": "ef", "ICFIELD": "eqfp"}


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
    outputs: dict[str, str] = field(default_factory=dict)  # name -> expression


def parse(path: str) -> Model:
    return parse_text(Path(path).read_text())


def parse_text(text: str) -> Model:
    m = Model()
    pending: list[tuple[str, list[tuple[str, int]]]] = []
    connects: list[tuple[str, str, str]] = []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        jm = re.fullmatch(r"junction\s+(\w+)\s+([01])\s*:\s*(.*)", line, re.I)
        em = re.fullmatch(r"element\s+(\w+)\s+(ICfield|Cfield|Ifield|Rfield|I|C|R|Se|Sf|TF|GY)(?:\s+(.*))?", line, re.I)
        cm = re.fullmatch(r"connect\s+(\w+)\s+(\w+)\s+(\w+)", line, re.I)
        om = re.fullmatch(r"output\s+(\w+)\s*=\s*(.+)", line, re.I)
        if jm:
            name, kind, members = jm.groups()
            if name in m.junctions:
                raise ValueError(f"line {number}: duplicate junction {name}")
            m.junctions[name] = kind
            parsed = []
            for token in filter(None, map(str.strip, members.split(","))):
                mm = re.fullmatch(r"(\w+)(?:\.(\d+))?\s*([+-])?", token)
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
            if kind.upper() in FIELDS:
                letters = FIELDS[kind.upper()]
                for key in params:
                    if not re.fullmatch(rf"[{letters}][1-9]\d*", key):
                        laws = ", ".join(f"{letter}1=" for letter in letters)
                        raise ValueError(f"line {number}: {kind} takes a law per port, like {laws} (the variable it sets and the port), not {key}=")
                m.elements[name] = Element(name, kind, params)
                continue
            forms = LAW_FORMS[kind.upper()]
            allowed = forms + (("causality",) if kind.upper() == "R" else ())
            unknown = [key for key in params if key not in allowed]
            if unknown:
                raise ValueError(f"line {number}: {kind} takes {', '.join(k + '=' for k in allowed)} not {unknown[0]}=")
            given = [key for key in forms if key in params]
            if len(given) > 1:
                raise ValueError(f"line {number}: {name} has two laws ({given[0]}= and {given[1]}=); give one")
            m.elements[name] = Element(name, kind, params)
        elif cm:
            connects.append(cm.groups())
        elif om:
            name, expr = om.groups()
            if name in m.outputs:
                raise ValueError(f"line {number}: duplicate output {name}")
            m.outputs[name] = expr
        else:
            raise ValueError(f"line {number}: expected junction, element, connect, or output")
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
                if len(bits) != 2 or port not in {1, 2}:
                    raise ValueError(f"junction {jname}: two-port {ename} needs .1 or .2")
            elif elem.kind.upper() in FIELDS:
                if len(bits) != 2 or port < 1:
                    raise ValueError(f"junction {jname}: field {ename} needs a port number, like {ename}.1")
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
        if elem.kind.upper() in FIELDS:
            gap = next((p for p, (j, _) in enumerate(elem.ports, 1) if not j), None)
            if not elem.ports or gap:
                raise ValueError(f"element {elem.name}: " + (f"port {gap} isn't attached; number a field's ports 1, 2, 3, ..." if gap else "needs at least one attached port"))
            field_ports(elem)  # each port has one law
            continue
        required = 2 if elem.kind.upper() in {"TF", "GY"} else 1
        if len(elem.ports) != required or any(not p[0] for p in elem.ports):
            raise ValueError(f"element {elem.name}: needs {required} attached port(s)")
    return link(m)


def link(m: Model) -> Model:
    """Rebuild m.attachments from the elements' ports."""
    m.attachments = {}
    for elem in m.elements.values():
        for p, (j, sign) in enumerate(elem.ports, 1):
            m.attachments.setdefault(j, []).append((elem.name, p, sign))
    return m


def port_name(name: str, port: int) -> str:
    return f"{name}_{port}" if port != 1 else name


def field_ports(el: Element) -> list[tuple[int, str, str, str]]:
    """A field's (port, letter, law, type) for each port in order; type is C, I or R.

    Checks that each port has one law and that the laws use only the field's own
    variables among names like e1, f2, q1 and p2.
    """
    n = len(el.ports)
    laws: dict[int, tuple[str, str]] = {}
    for key, value in el.params.items():
        port = int(key[1:])
        if port > n:
            raise ValueError(f"element {el.name}: has a law for port {port} ({key}=) but {n} port{'s' if n != 1 else ''}")
        if port in laws:
            raise ValueError(f"element {el.name}: port {port} has two laws ({laws[port][0]}{port}= and {key}=); give one")
        laws[port] = (key[0], value)
    letters, out = FIELDS[el.kind.upper()], []
    for p in range(1, n + 1):
        if p not in laws:
            raise ValueError(f"element {el.name}: port {p} needs a law ({' or '.join(f'{letter}{p}=' for letter in letters)})")
        letter, law = laws[p]
        out.append((p, letter, law, "R" if el.kind.upper() == "RFIELD" else "C" if letter in "eq" else "I"))
    own = {f"{v}{p}" for p in range(1, n + 1) for v in "ef"} | {f"{'q' if t == 'C' else 'p'}{p}" for p, _, _, t in out if t != "R"}
    for p, _, law, _ in out:
        stray = sorted(set(re.findall(r"\b[efqp]\d+\b", law)) - own)
        if stray:
            raise ValueError(f"element {el.name}: the law for port {p} uses {stray[0]}, which isn't one of its variables")
    return out


def storage_ports(m: Model) -> list[tuple[Element, int, str]]:
    """(element, port, C or I) for every storage port: C and I elements, and the C and I ports of fields."""
    out = []
    for el in m.elements.values():
        k = el.kind.upper()
        if k in {"C", "I"}:
            out.append((el, 1, k))
        elif k in FIELDS:
            out += [(el, p, t) for p, _, _, t in field_ports(el) if t != "R"]
    return out


def port_state(el: Element, port: int, kind: str) -> sp.Symbol:
    """State of a C (q) or I (p) storage port: q_X, or q_X_2 for a field's port 2."""
    return sp.Symbol(("q_" if kind == "C" else "p_") + port_name(el.name, port))


def integral_port(c, el: Element, port: int, kind: str) -> bool:
    return c.get((el.name, port)) == ("e" if kind == "C" else "f")


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
        elif k == "R" and e.params.get("causality"):
            mode = e.params["causality"].lower()
            if mode not in {"effort", "flow"}:
                raise ValueError(f"element {e.name}: R causality must be effort or flow")
            fixed[key] = "e" if mode == "effort" else "f"
        elif k in {"TF", "GY"}:
            pair_rules.append(((e.name, 1), (e.name, 2), k == "GY"))  # GY same, TF opposite
    for el, p, kind in storage_ports(m):
        preferred[el.name, p] = "e" if kind == "C" else "f"
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


def sym(value: str, own: dict | None = None):
    # Give names used like u(t) function semantics; all other undeclared names
    # remain symbolic parameters, which is convenient for R, C, I, and n.
    # SymPy functions (Piecewise, Abs, Max, sin, sqrt, ...) keep their meaning.
    # `own` binds an element's own e, f, q, p inside its law.
    functions = {
        name: sp.Function(name)
        for name in re.findall(r"\b([A-Za-z_]\w*)\s*\(", value)
        if name not in {"t", "sqrt"} and not isinstance(getattr(sp, name, None), sp.FunctionClass)
    }
    return sp.sympify(value.replace("^", "**"), locals={"t": sp.Symbol("t"), **(own or {}), **functions})


def integral_arg(expr: str) -> str | None:
    """The inside of integral(...) when that is the whole expression, else None."""
    mm = re.match(r"\s*integral\s*\(", expr)
    if not mm:
        return None
    depth = 0
    for i in range(mm.end() - 1, len(expr)):
        depth += {"(": 1, ")": -1}.get(expr[i], 0)
        if depth == 0:
            return None if expr[i + 1:].strip() else expr[mm.end():i]
    return None


def derive(m: Model, c, inputs: dict | None = None, reduction: dict | None = None, outputs: dict | None = None):
    """derive_states(), plus the outputs that need a state of their own.

    An output written integral(expr) is the integral of expr from 0, so it becomes
    an extra state (named after the output, with that derivative) and isn't put in
    `outputs`. Outputs may use each other's names: a plain one is substituted, an
    integral one is that state.
    """
    if outputs is None:
        return derive_states(m, c, inputs, reduction)
    plain, integrals = {}, {}
    for name, expr in m.outputs.items():
        arg = integral_arg(expr)
        (plain if arg is None else integrals)[name] = sym(expr if arg is None else arg)
    for y in [*plain.values(), *integrals.values()]:
        if any(f.func.__name__ == "integral" for f in y.atoms(AppliedUndef)):
            raise ValueError("integral(...) must be a whole output on its own; give it an output of its own and use that output's name")
    work = {**plain, **{"integral " + n: a for n, a in integrals.items()}}
    states, dots, equations, sol = derive_states(m, c, inputs, reduction, work)
    done: dict = {}

    def resolve(key, path=()):
        if key not in done:
            if key in path:
                raise ValueError(f"outputs {' and '.join(sorted(set(path)))} depend on each other")
            done[key] = work[key].xreplace({s: resolve(s.name, path + (key,)) for s in work[key].free_symbols if s.name in plain})
        return done[key]

    result = {name: resolve(name) for name in plain}
    extra, extra_dots = [], []
    for name in integrals:
        arg = resolve("integral " + name)
        if unsolved(m, states, {name: arg}):
            result[name] = arg  # reported as unsolved, like any output
        else:
            extra.append(sp.Symbol(name)); extra_dots.append(arg)
    outputs.update(result)
    return list(states) + extra, list(dots) + extra_dots, equations, sol


def derive_states(m: Model, c, inputs: dict | None = None, reduction: dict | None = None, outputs: dict | None = None):
    """Return states, their derivatives, the full equations and the solution (None if unsolved).

    Source values that don't depend on the system are solved as input symbols.
    Pass a dict as `inputs` to get derivatives in terms of those symbols (it's
    filled with symbol -> (source name, value)); otherwise the values are put back.

    Storage in derivative causality is reduced away first where it can be (see
    reduce_graph); the states and derivatives are still the original elements'.
    Pass a dict as `reduction` to get the steps taken ("steps") and each storage
    element that isn't a state as (its state, value in terms of the states) ("dependent").

    Pass a dict of name -> expression as `outputs`: it's updated in place with
    each one in terms of the states and inputs, as far as it solves (see unsolved).
    """
    ys = dict(outputs or {})
    red = reduce_graph(m, c)
    if red is None:
        # Storage left in derivative causality is solved for alongside the outputs.
        follow = [port_state(el, p, kind) for el, p, kind in storage_ports(m) if not integral_port(c, el, p, kind)]
        work = {**ys, **{("follows", x): x for x in follow}}
        result = derive_graph(m, c, inputs, work)
        if reduction is not None:
            dependent = [(x, work["follows", x]) for x in follow]
            reduction.update(steps=[], dependent=[(x, y) for x, y in dependent if not unsolved(m, result[0], {x: y})])
        if outputs is not None:
            outputs.update({name: work[name] for name in ys})
        return result
    states, dots, equations, sol = derive_graph(red.model, red.causality, inputs)
    stuck = set(derivative_storage(m, c))
    integral = {port_state(el, p, kind) for el, p, kind in storage_ports(m) if integral_port(c, el, p, kind)}
    order = {port_state(el, p, kind): i for i, (el, p, kind) in enumerate(storage_ports(m))}
    # Each reduced state is a multiple of one original state: keep an integral one.
    # A field's states aren't reduced, so they're their own.
    picks, dependent = [], []
    for x in states:
        origin = red.origin.get(x.name[2:], {x: sp.Integer(1)})
        pick = next((y for y in origin if y in integral), next(iter(origin)))
        picks.append((pick, origin[pick]))
        dependent += [(y, sp.simplify(k / origin[pick] * pick)) for y, k in origin.items() if y != pick]
    swap = {x: pick / k for x, (pick, k) in zip(states, picks)}
    result = sorted(((pick, sp.simplify(k * d.xreplace(swap))) for (pick, k), d in zip(picks, dots)), key=lambda r: order[r[0]])
    if reduction is not None:
        reduction.update(steps=red.steps, dependent=sorted(dependent, key=lambda r: order[r[0]]))
    if ys:
        ys = {name: y.xreplace(dict(dependent)) for name, y in ys.items()}
        if sol is not None:
            # The reduced graph's bond variables aren't the original's, so solve the
            # original graph again. Every storage element's rate is known now, and each
            # one in derivative causality imposes its rate in place of its law.
            rate = dict(result)
            rate.update((y, sum(sp.diff(value, x) * dx for x, dx in result)) for y, value in dependent)
            derive_graph(m, c, inputs, ys, {x: r for x, r in rate.items() if x.name[2:] in stuck})
    if outputs is not None:
        outputs.update(ys)
    return [x for x, _ in result], [d for _, d in result], equations, sol


def derive_graph(m: Model, c, inputs: dict | None = None, outputs: dict | None = None, rates: dict | None = None):
    """derive() without the reduction.

    A source already in `inputs` keeps its symbol. Each expression in `outputs`
    (name -> expression) is replaced by its value in terms of the states. `rates`
    maps the state of a storage element in derivative causality to its known rate,
    which the element imposes in place of its law.
    """
    rates = rates or {}
    e, f, fin = {}, {}, {}
    equations, states, dots = [], [], []
    given: dict = {}
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
    system = set(e.values()) | set(f.values()) | {port_state(el, p, kind) for el, p, kind in storage_ports(m)}
    for el in m.elements.values():
        k = el.kind.upper(); key = (el.name, 1)
        if k in {"R", "C", "I"}:
            form = next((f for f in LAW_FORMS[k] if f in el.params), k)
            own = {"e": e[key], "f": fin[key]}
            if k != "R":
                state = sp.Symbol(("q_" if k == "C" else "p_") + el.name)
                own["q" if k == "C" else "p"] = state
            law = sym(el.params.get(form, f"{k}_{el.name}"), own)
            if k != "R" and state in rates: equations.append((fin[key] if k == "C" else e[key]) - rates[state])
            elif form != k: equations.append(own[form] - law)  # e=, f=, q= or p=: the law itself
            elif k == "R": equations.append(e[key] - law * fin[key])
            elif k == "C": equations.append(e[key] - state / law)
            else: equations.append(fin[key] - state / law)
            if k != "R" and (k, c[key]) in {("C", "e"), ("I", "f")}:
                states.append(state); dots.append(fin[key] if k == "C" else e[key])
        elif k in FIELDS:
            ports = field_ports(el)
            own = {f"{v}{p}": x for p in range(1, len(el.ports) + 1) for v, x in (("e", e[el.name, p]), ("f", fin[el.name, p]))}
            own.update({f"{'q' if t == 'C' else 'p'}{p}": port_state(el, p, t) for p, _, _, t in ports if t != "R"})
            for p, letter, law, t in ports:
                equations.append(own[f"{letter}{p}"] - sym(law, own))
                if t != "R" and integral_port(c, el, p, t):
                    states.append(port_state(el, p, t)); dots.append(fin[el.name, p] if t == "C" else e[el.name, p])
        elif k in {"SE", "SF"}:
            value = sym(el.params.get("value", "u_" + el.name))
            if not value.free_symbols & system:  # an input, not a modulated source
                same = [u for u, (src, _) in (inputs or {}).items() if src == el.name]
                u = same[0] if same else sp.Dummy("u_" + el.name)
                given[u] = (el.name, value)
                value = u
            equations.append((e[key] if k == "SE" else -fin[key]) - value)  # Sf: flow delivered
        # Two-ports: power enters port 1 and leaves port 2 (f2 out = -fin[2]).
        elif k == "TF":
            n = sym(el.params.get("n", "n_" + el.name)); equations += [e[el.name,1] - n*e[el.name,2], -fin[el.name,2] - n*fin[el.name,1]]
        elif k == "GY":
            n = sym(el.params.get("n", "n_" + el.name)); equations += [e[el.name,1] + n*fin[el.name,2], e[el.name,2] - n*fin[el.name,1]]
    # Storage in derivative causality isn't a state: its law fixes its q or p, so
    # that's solved for too (a flow source driving a mass fixes the mass's momentum).
    stuck = [port_state(el, p, kind) for el, p, kind in storage_ports(m) if not integral_port(c, el, p, kind)]
    variables = list(e.values()) + list(f.values()) + [x for x in stuck if x not in rates]
    # Piecewise element laws (check valves, piecewise resistors) are still linear in
    # the unknowns, but sp.solve evaluates their conditions during elimination and
    # can collapse to 0/0 -> nan, and it can't solve through Abs or sign at all.
    # Swap each for a placeholder, solve, then put it back.
    placeholders: dict[sp.Expr, sp.Expr] = {}

    def linearize(eq):
        for pw in eq.atoms(sp.Piecewise, sp.Abs, sp.sign, sp.Heaviside, sp.Max, sp.Min):
            placeholders.setdefault(pw, sp.Dummy(f"Switch{len(placeholders)}"))
        return eq.xreplace(placeholders)

    try:
        solutions = sp.solve([linearize(eq) for eq in equations], variables, dict=True, simplify=True)
    except NotImplementedError:  # e.g. inverting e=a*f*Abs(f); report it as unsolved
        solutions = []
    values = {u: value for u, (_, value) in given.items()}
    back = {} if inputs is not None else values
    if inputs is not None:
        inputs.update(given)
    equations = [eq.xreplace(values) for eq in equations]
    if not solutions:
        return states, dots, equations, None
    restore = {dummy: pw for pw, dummy in placeholders.items()}
    unknowns = set(variables)
    candidates = []
    for sol in solutions:
        sol = settle({var: val.xreplace(restore) for var, val in sol.items()}, unknowns)
        result = [x.xreplace(sol) for x in dots]
        solved = [r for r in result if not r.free_symbols & unknowns]
        candidates.append((len(solved), solved, result, sol))
    # Prefer the solution that eliminates the most unknowns. A law with several
    # inverses (p=Psi*tanh(f/i0) solved for f) gives several solutions; then
    # keep the first that's real wherever the others are.
    count, _, result, sol = max(candidates, key=lambda cand: (cand[0], real_samples(cand[1])))
    if outputs:
        outputs.update({name: sp.simplify(y.xreplace(sol).xreplace(back)) for name, y in outputs.items()})
    return states, [sp.simplify(r.xreplace(back)) for r in result], equations, sol if count == len(dots) else None


def settle(sol: dict, unknowns: set) -> dict:
    """Substitute solved unknowns into each other, as far as they resolve.

    A restored placeholder can still mention unknowns: a diode
    f=Piecewise((e/R_on,e>0),...) depends on its own e. Values that depend on an
    unsolved unknown, or on themselves through a loop, are left as they are.
    """
    done: dict = {}  # var -> resolved value, or None when it can't be resolved

    def resolve(var, path):
        if var in path or var not in sol:
            return None
        if var not in done:
            reps = {}
            for u in sol[var].free_symbols & unknowns:
                reps[u] = resolve(u, path | {var})
                if reps[u] is None:  # every var on `path` depends on var, so this is a loop through var
                    break
            done[var] = None if None in reps.values() else sol[var].xreplace(reps)
        return done[var]

    resolved = {var: resolve(var, frozenset()) for var in sol}
    return {var: val if resolved[var] is None else resolved[var] for var, val in sol.items()}


def real_samples(exprs, count: int = 8) -> int:
    """Count random points (every symbol and input in 0.1..1) where all exprs are real."""
    rng = random.Random(0)
    atoms = set().union(*(x.free_symbols | x.atoms(AppliedUndef) for x in exprs))
    hits = 0
    for _ in range(count):
        point = {a: sp.Float(rng.uniform(0.1, 1)) for a in atoms}
        try:
            values = [complex(x.xreplace(point).evalf()) for x in exprs]
        except (TypeError, ValueError, ZeroDivisionError):
            continue
        hits += all(math.isfinite(v.real) and abs(v.imag) < 1e-9 for v in values)
    return hits


def storage_causality(m: Model, integral: dict[str, bool]):
    """Causality for derive() from which storage is in integral causality.

    `integral` is keyed by element name, or "name.port" for a field's port.
    """
    c = {}
    for el, p, kind in storage_ports(m):
        wanted, other = ("e", "f") if kind == "C" else ("f", "e")
        c[el.name, p] = wanted if integral.get(el.name if el.kind.upper() in {"C", "I"} else f"{el.name}.{p}", True) else other
    return c


def state_symbol(el: Element) -> sp.Symbol:
    return sp.Symbol(("q_" if el.kind.upper() == "C" else "p_") + el.name)


def derivative_storage(m: Model, c) -> list[str]:
    return [el.name for el in m.elements.values() if (el.kind.upper(), c.get((el.name, 1))) in {("C", "f"), ("I", "e")}]


def unsolved(m: Model, states, outputs: dict) -> list[str]:
    """Outputs that still mention a bond variable, or a storage element that isn't a state."""
    left = {port_state(el, p, kind) for el, p, kind in storage_ports(m)} - set(states)
    left |= {sp.Symbol(f"{v}_{port_name(el.name, p)}") for el in m.elements.values() for p in range(1, len(el.ports) + 1) for v in "ef"}
    return [name for name, y in outputs.items() if y.free_symbols & left]


# ------------------------------------------------------------------ reduction
# Storage in derivative causality is often two I's (or two C's) competing for one
# flow (or effort), directly or across a TF or GY. Reflecting the elements on one
# side of the two-port onto the other side and combining the pair leaves a graph
# with fewer storage elements in derivative causality, solved like any other.

OWN = {v: sp.Symbol(v) for v in "efqp"}  # an element's own variables inside its law
KAPPA = {1: -1, 2: 1}  # GY: e at port P = KAPPA[P]*n*(inflow at the other port)


@dataclass
class Reduced:
    model: Model
    causality: dict | None
    steps: list  # (what was done, [(name, value)])
    origin: dict  # storage element -> {original state: k}, where original state = k * this element's state

    def copy(self) -> Reduced:
        els = {n: Element(el.name, el.kind, dict(el.params), list(el.ports)) for n, el in self.model.elements.items()}
        return Reduced(Model(dict(self.model.junctions), {}, els), None, list(self.steps),
                       {n: dict(o) for n, o in self.origin.items()})


def law_of(el: Element) -> tuple[str, sp.Expr]:
    """An R, C or I law as (key, expression in OWN); a coefficient becomes e=R*f, e=q/C or f=p/I."""
    k = el.kind.upper()
    form = next((f for f in LAW_FORMS[k] if f in el.params), k)
    if form != k:
        return form, sym(el.params[form], OWN)
    coef = sym(el.params.get(k, f"{k}_{el.name}"))
    return {"R": ("e", coef * OWN["f"]), "C": ("e", OWN["q"] / coef), "I": ("f", OWN["p"] / coef)}[k]


def coefficient(el: Element):
    """R, C or I of a linear law (e = R*f, q = C*e, p = I*f), else None."""
    key, law = law_of(el)
    x, y = {"R": "ef", "C": "qe", "I": "pf"}[el.kind.upper()]
    if law == 0:
        return None
    coef = sp.simplify(law / OWN[y] if key == x else OWN[x] / law)
    return None if coef.free_symbols & set(OWN.values()) or coef == 0 else coef


def plain(el: Element) -> bool:
    """el's laws don't read t or any element's variables, and its causality isn't fixed."""
    return "causality" not in el.params and not any(re.search(r"\bt\b|\b[efqp]_\w", v) for v in el.params.values())


def read_by_laws(m: Model) -> set[str]:
    """Elements whose e, f, q or p some law or modulus reads."""
    words = {w for el in m.elements.values() for v in el.params.values() for w in re.findall(r"\b[efqp]_(\w+)", v)}
    return {w if w in m.elements else re.sub(r"_\d+$", "", w) for w in words}


def join_names(els) -> str:
    names = [el.name for el in els]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def combine(r: Reduced, j: str) -> Reduced | None:
    """Merge the I's on a 1-junction (one flow) or the C's on a 0-junction (one effort) into one."""
    m = r.model
    kind = "I" if m.junctions[j] == "1" else "C"
    group = [m.elements[n] for n, _, _ in m.attachments.get(j, []) if m.elements[n].kind.upper() == kind]
    if len(group) < 2 or not all(map(plain, group)) or read_by_laws(m) & {el.name for el in group}:
        return None
    coefs = [coefficient(el) for el in group]
    if None in coefs:
        return None
    total = sp.simplify(sum(coefs))  # inertias sharing a flow add, and so do capacitances sharing an effort
    out = r.copy()
    name = "_".join(el.name for el in group)
    while name in m.elements:
        name += "_"
    origin = {}
    for el, coef in zip(group, coefs):
        # p_k = sign_k*I_k*f on a 1-junction; q_k = C_k*e on a 0-junction
        scale = coef / total * (el.ports[0][1] if kind == "I" else 1)
        origin.update({x: sp.simplify(k * scale) for x, k in out.origin.pop(el.name).items()})
        del out.model.elements[el.name]
    out.model.elements[name] = Element(name, kind, {kind: str(total)}, [(j, 1)])
    out.origin[name] = origin
    out.steps.append((f"combined {join_names(group)} on {j} into one {'inertia' if kind == 'I' else 'capacitor'}",
                      [(kind, total)]))
    link(out.model)
    return out


def reflect(r: Reduced, T: Element, P: int) -> Reduced | None:
    """Move the R, C and I on T's far junction (port 3-P) across T onto its near one, merging the two junctions."""
    m, Q = r.model, 3 - P
    (near, a), (far, b) = T.ports[P - 1], T.ports[Q - 1]
    gy = T.kind.upper() == "GY"
    if near == far or (m.junctions[near] == m.junctions[far]) == gy:  # a TF joins like junctions, a GY unlike
        return None
    moved = [m.elements[n] for n, _, _ in m.attachments[far] if n != T.name]
    if (not moved or not plain(T) or not all(map(plain, moved)) or any(el.kind.upper() not in {"R", "C", "I"} for el in moved)
            or read_by_laws(m) & {T.name, *(el.name for el in moved)}):
        return None
    n = sym(T.params.get("n", "n_" + T.name))
    out = r.copy()
    del out.model.elements[T.name], out.model.junctions[far]
    laws = []
    for el in moved:
        s = el.ports[0][1]
        # sub: far variable -> (near variable, k) with far = k * near
        if gy:  # efforts and flows swap, so I <-> C
            sigma = KAPPA[P] * b * s if m.junctions[far] == "1" else 1
            sign = KAPPA[Q] * a if m.junctions[far] == "0" else 1
            sub = {"e": ("f", sigma * n), "f": ("e", sigma / n), "p": ("q", sigma * n), "q": ("p", sigma / n)}
            kind = {"R": "R", "C": "I", "I": "C"}[el.kind.upper()]
        else:  # e1 = n*e2, f2 = n*f1
            g = n if Q == 2 else 1 / n
            sign = -s * a * b
            sub = {"e": ("e", 1 / g), "f": ("f", g), "q": ("q", g), "p": ("p", 1 / g)}
            kind = el.kind.upper()
        key, law = law_of(el)
        new_key, k = sub[key]
        law = sp.simplify(law.xreplace({OWN[v]: kv * OWN[w] for v, (w, kv) in sub.items()}) / k)
        new = out.model.elements[el.name] = Element(el.name, kind, {new_key: str(law)}, [(near, sign)])
        if el.kind.upper() != "R":
            k = sub["q" if el.kind.upper() == "C" else "p"][1]
            out.origin[el.name] = {x: sp.simplify(kx * k) for x, kx in out.origin[el.name].items()}
        coef = coefficient(new)
        laws.append((f"{kind}_{el.name}", coef) if coef is not None else
                    (f"{new_key}_{el.name}", law.xreplace({OWN[v]: sp.Symbol(f"{v}_{el.name}") for v in OWN})))
    out.steps.append((f"reflected {join_names(moved)} through {T.name} onto {near}", laws))
    link(out.model)
    return out


def reductions(r: Reduced):
    for j in r.model.junctions:
        yield combine(r, j)
    for T in r.model.elements.values():
        if T.kind.upper() in {"TF", "GY"}:
            for P in (1, 2):
                moved = reflect(r, T, P)
                yield moved and combine(moved, T.ports[P - 1][0])


def reduce_graph(m: Model, c) -> Reduced | None:
    """Reflect and combine storage while that leaves fewer elements in derivative causality.

    None when nothing is in derivative causality or no reduction helps.
    """
    stuck = len(derivative_storage(m, c))
    best = Reduced(m, c, [], {el.name: {state_symbol(el): sp.Integer(1)}
                              for el in m.elements.values() if el.kind.upper() in {"C", "I"}})
    while stuck:
        for trial in filter(None, reductions(best)):
            try:
                trial.causality = causality(trial.model)
            except ValueError:
                continue
            left = len(derivative_storage(trial.model, trial.causality))
            if left < stuck:
                best, stuck = trial, left
                break
        else:
            break
    return best if best.steps else None


def describe(step) -> str:
    text, laws = step
    return f"{text}: " + ", ".join(f"{name} = {value}" for name, value in laws)


def state_space(states, dots, inputs):
    """A, B, d with dots = A*x + B*u + d, or None if dots aren't linear in x and u."""
    x, u, f = sp.Matrix(states), sp.Matrix(len(inputs), 1, list(inputs)), sp.Matrix(dots)
    A = f.jacobian(x).applyfunc(sp.simplify)
    B = f.jacobian(u).applyfunc(sp.simplify) if inputs else sp.zeros(len(dots), 0)
    d = (f - A * x - B * u).applyfunc(sp.simplify)
    if any(M.free_symbols & (set(states) | set(inputs)) for M in (A, B, d)):
        return None
    return A, B, d


GREEK = set("alpha beta gamma delta epsilon varepsilon zeta eta theta vartheta iota kappa lambda mu nu xi "
            "pi rho sigma tau upsilon phi varphi chi psi omega Gamma Delta Theta Lambda Xi Pi Sigma "
            "Upsilon Phi Psi Omega".split())


def tex_name(name: str, dot: bool = False, bar: bool = False) -> str:
    """LaTeX for a name: q_C1 -> q_{\\mathrm{C1}}, R_a -> R_{a}, F0 -> F_{0}, rho -> \\rho.

    SymPy's own printer mangles some names (p_I_Electrical ends in "cal", which it
    reads as \\mathcal), so every symbol goes through this instead. A deviation
    δq_C1 is \\delta q_{\\mathrm{C1}} and an operating point q_C1‾ is \\bar{q}_{\\mathrm{C1}}.
    """
    if name.startswith(DELTA):
        return r"\delta " + tex_name(name[1:], dot)
    if name.endswith(BAR):
        return tex_name(name[:-1], dot, bar=True)
    if name.endswith(DOT):
        return tex_name(name[:-1], True, bar)
    base, _, sub = name.partition("_")
    digits = re.fullmatch(r"([A-Za-z]+?)(\d+)", base)
    if not sub and digits:
        base, sub = digits.groups()
    head = "\\" + base if base in GREEK else rf"\mathit{{{base}}}" if base.isalpha() and len(base) > 1 else base
    if dot:
        head = rf"\dot{{{head}}}"
    if bar:
        head = rf"\bar{{{head}}}"
    if not sub:
        return head
    sub = sub if len(sub) == 1 or sub.isdigit() else r"\mathrm{%s}" % sub.replace("_", r"\_")
    return f"{head}_{{{sub}}}"


def tex(expr) -> str:
    return sp.latex(expr, symbol_names={s: tex_name(s.name) for s in expr.free_symbols if isinstance(s, sp.Symbol)})


def tex_column(items) -> str:
    return r"\left[\begin{matrix}" + r" \\ ".join(items) + r"\end{matrix}\right]"


# ------------------------------------------------------------------ views
# The state and output equations can be shown in energy variables (the states p
# and q, as derived), in power variables (each storage port's effort or flow in
# their place) or in Lagrange variables (displacements q and flows f: each I
# port's flow in place of its p); and a nonlinear model can be replaced by its
# linear approximation about an operating point first.

DELTA, BAR, DOT = "δ", "‾", "˙"  # δx: a deviation from an operating point; x‾: the point's value; u˙: a rate
TIME = sp.Symbol("t")
INPUT_RATES = ("Some of these variables depend on an input, because a source drives storage in derivative causality, "
               "so their equations include that input's rate of change.")
VARIABLES = ("energy", "power", "lagrange")


def derive_system(m: Model, c, power: bool = False) -> dict:
    """One solve for any view: derive()'s states, derivatives, outputs and inputs.

    With `power`, each storage port's power variable (the effort of a C port, the
    flow of an I port) is solved for too, in terms of the states: "power" maps each
    storage state to (its power variable, that value).
    """
    hidden = {}
    if power:
        hidden = {f"\0{port_state(el, p, kind)}": f"{'e' if kind == 'C' else 'f'}_{port_name(el.name, p)}"
                  for el, p, kind in storage_ports(m)}
    given: dict = {}
    reduction: dict = {}
    ys: dict = {}
    states, dots, eqs, sol = derive(Model(m.junctions, m.attachments, m.elements, {**m.outputs, **hidden}), c, given, reduction, ys)
    powers = {sp.Symbol(k[1:]): (sp.Symbol(v), ys.pop(k)) for k, v in hidden.items()}
    return {"states": states, "dots": dots, "equations": eqs, "solved": sol is not None, "inputs": given,
            "reduction": reduction, "outputs": ys, "power": powers if power else None}


def input_columns(given: dict) -> list[dict]:
    """The inputs as the entries of u, in order.

    Each source is an entry: its value when that's a lone name or input like v(t),
    else u_<source>, defined by a "where" line. Sources fed by one v(t) share one.
    An entry has slot (the symbol for it in the equations), shown (how u shows it),
    value (what it stands for), name (shown, as a string) and sources (their input
    symbols in derive's results).
    """
    cols: list[dict] = []
    for u, (src, value) in given.items():
        shown = value if isinstance(value, (sp.Symbol, AppliedUndef)) else sp.Symbol("u_" + src)
        col = next((col for col in cols if col["shown"] == shown), None)
        if col is None:
            name = shown.func.__name__ if isinstance(shown, AppliedUndef) else shown.name
            col = {"slot": sp.Dummy(name), "shown": shown, "value": value, "name": name, "sources": []}
            cols.append(col)
        col["sources"].append(u)
    return cols


def jacobian(exprs, wrt, at: dict) -> sp.Matrix:
    rows = sp.Matrix(len(exprs), 1, list(exprs))
    if not wrt:
        return sp.zeros(len(exprs), 0)
    return rows.jacobian(list(wrt)).xreplace(at).applyfunc(sp.simplify)


def view(m: Model, c, variables: str = "energy", linearize: bool = False, system: dict | None = None) -> dict:
    """The state and output equations to show, in the variables asked for.

    variables: energy (the states p and q, as derived), power (each storage port's
    effort or flow in place of its state) or lagrange (each I port's flow in place
    of its p; displacements q stay). linearize: replace a nonlinear model by its linear approximation
    about an operating point (x‾, u‾), in deviations δx, δu and δy from it, with
    "operating" the state derivatives there (each 0 at an equilibrium).

    Inputs appear as the slots of "columns" (see input_columns). "note" says why
    the view asked for can't be shown, and is "" when it can.
    """
    sysd = system or derive_system(m, c, power=variables != "energy")
    cols = input_columns(sysd["inputs"])
    to_slot = {u: col["slot"] for col in cols for u in col["sources"]}
    slots = [col["slot"] for col in cols]
    X = list(sysd["states"])
    F = [d.xreplace(to_slot) for d in sysd["dots"]]
    stuck = unsolved(m, X, sysd["outputs"])
    Y = [(name, y.xreplace(to_slot)) for name, y in sysd["outputs"].items()]
    dep = [(x, y.xreplace(to_slot)) for x, y in sysd["reduction"].get("dependent", [])]
    out = {"states": X, "dots": F, "outputs": Y, "dependent": dep, "columns": cols, "unsolved_outputs": stuck,
           "solved": sysd["solved"], "equations": sysd["equations"], "steps": sysd["reduction"].get("steps", []),
           "variables": variables, "linearized": False, "operating": [], "note": "", "hint": ""}
    ready = bool(sysd["solved"] and X)
    linear_ys = lambda: not Y or bool(stuck) or state_space(X, [y for _, y in Y], slots) is not None  # noqa: E731
    out["nonlinear"] = ready and not (state_space(X, F, slots) is not None and linear_ys())
    if variables == "energy" and not (linearize and out["nonlinear"]):
        return out
    if not ready:
        out["note"] = ("There are no states." if sysd["solved"] else
                       "This needs every equation solved, and SymPy couldn't eliminate every algebraic variable. "
                       "Energy variables show how far it got.")
        return out
    # The variables that take the states' places: every storage port's for power
    # variables, the I ports' flows for Lagrange variables (a C port keeps its q).
    powers = {x: (w, value.xreplace(to_slot)) for x, (w, value) in (sysd["power"] or {}).items()
              if variables == "power" or w.name.startswith("f_")}

    if linearize and out["nonlinear"]:
        # Jacobians at the operating point, as deviations from it.
        point = {**{x: sp.Symbol(x.name + BAR) for x in X}, **{col["slot"]: sp.Symbol(col["name"] + BAR) for col in cols}}
        Ys = [] if stuck else Y
        A, B = jacobian(F, X, point), jacobian(F, slots, point)
        C, D = jacobian([y for _, y in Ys], X, point), jacobian([y for _, y in Ys], slots, point)
        Jd, Ed = jacobian([y for _, y in dep], X, point), jacobian([y for _, y in dep], slots, point)
        names = [x.name for x in X]
        dep_names = [x.name for x, _ in dep]
        rates = None  # with power variables, the coefficients of the inputs' rates
        if powers:
            # δw = T δx + S δu with T and S the power variables' Jacobians at the
            # operating point, so dδw/dt = T dδx/dt + dT/dt δx + S dδu/dt + dS/dt δu.
            Wall = [powers[x][1] if x in powers else x for x in X]
            T, S = jacobian(Wall, X, point), jacobian(Wall, slots, point)
            if T.det() == 0:
                out["note"] = "The power variables don't determine the states at this operating point."
                return out
            Ti = T.inv().applyfunc(sp.simplify)
            A2 = ((T * A + T.diff(TIME)) * Ti).applyfunc(sp.simplify)
            A, B, rates = A2, (T * B + S.diff(TIME) - A2 * S).applyfunc(sp.simplify), S
            C, D = (C * Ti).applyfunc(sp.simplify), (D - C * Ti * S).applyfunc(sp.simplify)
            Wd = [powers[x][1] if x in powers else y for x, y in dep]
            Jw = jacobian(Wd, X, point)
            Jd, Ed = (Jw * Ti).applyfunc(sp.simplify), (jacobian(Wd, slots, point) - Jw * Ti * S).applyfunc(sp.simplify)
            names = [powers[x][0].name if x in powers else x.name for x in X]
            dep_names = [powers[x][0].name if x in powers else x.name for x, _ in dep]
        dX = sp.Matrix([sp.Symbol(DELTA + n) for n in names])
        dU = sp.Matrix([sp.Symbol(DELTA + col["name"]) for col in cols])
        lin = lambda M, N: list(M * dX + N * dU) if M.rows else []  # noqa: E731
        columns = [{**col, "slot": u, "shown": u, "value": u, "name": u.name, "defines": col} for col, u in zip(cols, dU)]
        dots = lin(A, B)
        used = [j for j in range(len(cols)) if rates is not None and any(rates[:, j])]
        if used:  # the rates of the inputs a power variable depends on
            dR = [sp.Symbol(DELTA + cols[j]["name"] + DOT) for j in used]
            dots = [d + sum(rates[i, j] * r for j, r in zip(used, dR)) for i, d in enumerate(dots)]
            columns += [{"slot": r, "shown": r, "value": r, "name": r.name, "sources": [], "rate": True} for r in dR]
            out["hint"] = INPUT_RATES
        out.update(states=list(dX), dots=dots, linearized=True, columns=columns,
                   outputs=[(DELTA + name, y) for (name, _), y in zip(Ys, lin(C, D))],
                   dependent=[(sp.Symbol(DELTA + n), y) for n, y in zip(dep_names, lin(Jd, Ed))],
                   operating=[sp.simplify(f.xreplace(point)) for f in F])
    elif powers:
        # Exact: invert the new variables for the states, w = W(x, u) -> x = X(w, u),
        # and differentiate: dw/dt = dW/dx dx/dt + dW/du du/dt + dW/dt. W depends on an
        # input when a source drives storage in derivative causality that W reads,
        # and on time when a storage law is modulated by it.
        storage = [x for x in X if x in powers]
        W = [powers[x][1] for x in storage]
        # SymPy can take forever on laws that couple states nonlinearly (a cubic system
        # for a capacitive transducer), so it's only asked for an inverse that's linear,
        # or where each new variable depends on one state.
        coupled = any(len(w.free_symbols & set(storage)) > 1 for w in W)
        if coupled and sp.Matrix(W).jacobian(storage).free_symbols & set(storage):
            out["note"] = ("These storage laws couple the states nonlinearly, so they can't be inverted to write the states "
                           "in these variables. The linear approximation can.")
            return out
        try:
            inverses = sp.solve([powers[x][0] - w for x, w in zip(storage, W)], storage, dict=True)
        except NotImplementedError:
            inverses = []
        inverses = [inv for inv in inverses if set(inv) == set(storage)]
        if not inverses:
            out["note"] = "SymPy couldn't invert a storage law to write the states in these variables. Try the linear approximation."
            return out
        inv = max(inverses, key=lambda i: (round_trips(i, storage, [powers[x][0] for x in storage], W, slots), real_samples(list(i.values()))))
        hints = []
        if any(p.exp.is_Rational and not p.exp.is_Integer for value in inv.values() for p in value.atoms(sp.Pow)):
            hints.append("Writing the states in these variables took the inverse of a nonlinear storage law, which has roots in it. "
                         "Where it has several, the one kept gives back states that are small next to the parameters. "
                         "The linear approximation is simpler.")
        used = [col for col in cols if any(w.has(col["slot"]) for w in W)]
        rate_cols = [{"slot": sp.Dummy(col["name"] + DOT), "shown": sp.Symbol(col["name"] + DOT), "value": sp.diff(col["value"], TIME),
                      "name": col["name"] + DOT, "sources": [], "rate": True} for col in used]
        Wm = sp.Matrix(W) if W else sp.zeros(0, 1)
        rate = Wm.jacobian(storage) * sp.Matrix([F[X.index(x)] for x in storage]) + Wm.diff(TIME) if storage else sp.Matrix([])
        if used:
            rate += Wm.jacobian([col["slot"] for col in used]) * sp.Matrix([r["slot"] for r in rate_cols])
            hints.append(INPUT_RATES)
        out.update(states=[powers[x][0] if x in powers else x for x in X], columns=cols + rate_cols, hint=" ".join(hints),
                   dots=[sp.simplify((rate[storage.index(x)] if x in powers else f).xreplace(inv)) for x, f in zip(X, F)],
                   outputs=[(name, sp.simplify(y.xreplace(inv))) for name, y in Y],
                   dependent=[(powers[x][0], sp.simplify(powers[x][1].xreplace(inv))) if x in powers else (x, sp.simplify(y.xreplace(inv)))
                              for x, y in dep])

    return out


def round_trips(inverse: dict, states: list, ws: list, W: list, others) -> int:
    """At how many random points x = inverse(w) gives back the x that w = W(x) came from.

    Several inverses (x = l - r and x = l + r) are told apart this way. The states and
    inputs are sampled small next to the parameters, as in most models (a
    displacement smaller than the gap it closes), so the physical one wins.
    """
    rng = random.Random(1)
    small = set(states) | set(others)
    atoms = (set().union(*(w.free_symbols | w.atoms(AppliedUndef) for w in W), *(v.free_symbols for v in inverse.values())) | small) - set(ws)
    hits = 0
    for _ in range(8):
        point = {a: sp.Float(rng.uniform(0.01, 0.1) if a in small else rng.uniform(0.5, 1)) for a in atoms}
        try:
            w = {wi: Wi.xreplace(point) for wi, Wi in zip(ws, W)}
            back = [complex(inverse[x].xreplace(point).xreplace(w).evalf()) for x in states]
        except (TypeError, ValueError, ZeroDivisionError):
            continue
        hits += all(abs(b - complex(point[x])) < 1e-6 * (1 + abs(b)) for b, x in zip(back, states))
    return hits


def latex_report(m: Model, c, variables: str = "energy", linearize: bool = False, system: dict | None = None) -> dict:
    """The state and output equations as LaTeX, one per line and in matrix form, for --latex and the editor.

    `variables` and `linearize` pick the view (see view); `system` is a
    derive_system() result, so that several views share one solve.
    """
    v = view(m, c, variables, linearize, system)
    X, F, Y, cols = v["states"], v["dots"], v["outputs"], v["columns"]
    inline = {col["slot"]: col["value"] for col in cols}
    pairs = [(tex_name(x.name, dot=True), tex(d.xreplace(inline))) for x, d in zip(X, F)]
    outs = [(tex_name(name), tex(y.xreplace(inline))) for name, y in Y]
    # A deviation δu_F stands for u_F's, so its definition still applies; a rate's follows from it.
    defines = [col.get("defines", col) for col in cols if not col.get("rate")]
    where = [f"{tex(d['shown'])} = {tex(d['value'])}" for d in defines if d["shown"] != d["value"]]
    out = {
        "variables": v["variables"], "linearized": v["linearized"], "nonlinear": v["nonlinear"], "note": v["note"], "hint": v["hint"],
        "states": [tex_name(x.name) for x in X],
        "equations": pairs,
        "outputs": outs,
        "unsolved_outputs": v["unsolved_outputs"],
        "separate": aligned([pairs, outs], where if v["linearized"] else []),
        "solved": v["solved"],
        "full": [tex(eq) + " = 0" for eq in v["equations"]],
        "matrix": None, "linear": False, "outputs_linear": False, "inputs": [],
        "reductions": [{"text": text, "math": r",\quad ".join(f"{tex_name(name)} = {tex(value)}" for name, value in laws)}
                       for text, laws in v["steps"]],
        "dependent": [f"{tex_name(str(x))} = {tex(value.xreplace(inline))}" for x, value in v["dependent"]],
        "operating": aligned([[("0", tex(f)) for f in v["operating"]]]) if v["operating"] else None,
    }
    if v["note"] or not X or not v["solved"]:
        return out
    slots = [col["slot"] for col in cols]
    uses_u = False

    def form(lhs, exprs):
        """lhs = M x + N u (+ the rest), or lhs = the vector of exprs when they aren't linear in x and u."""
        nonlocal uses_u
        split = state_space(X, exprs, slots)
        if split is None:
            return f"{lhs} = {tex(sp.Matrix([y.xreplace(inline) for y in exprs]))}", False
        M, N, rest = split
        rhs = [f"{tex(M)} {tex_column(tex_name(x.name) for x in X)}"]
        if any(N):
            rhs.append(f"{tex(N)} {tex_column(tex(col['shown']) for col in cols)}")
            uses_u = True
        if any(rest):
            rhs.append(tex(rest))
        return f"{lhs} = " + " + ".join(rhs), True

    rows = []
    row, out["linear"] = form(tex_column(tex_name(x.name, dot=True) for x in X), F)
    rows.append(row)
    if Y and not v["unsolved_outputs"]:
        row, out["outputs_linear"] = form(tex_column(tex_name(name) for name, _ in Y), [y for _, y in Y])
        rows.append(row)
    if uses_u:
        out["inputs"] = where
    if out["inputs"]:
        rows.append("\\text{where } " + ",\\quad ".join(out["inputs"]))
    out["matrix"] = rows[0] if len(rows) == 1 else "\\begin{gathered}\n" + " \\\\\n".join(rows) + "\n\\end{gathered}"
    return out


def aligned(groups, where=()) -> str:
    """Groups of (lhs, rhs) as one aligned block, a little apart, then a where line."""
    blocks = [" \\\\\n".join(f"{l} &= {r}" for l, r in group) for group in groups if group]
    if where:
        blocks.append(r"&\text{where } " + ",\\quad ".join(where))
    return "\\begin{aligned}\n" + " \\\\[1ex]\n".join(blocks) + "\n\\end{aligned}"


def view_text(v: dict) -> dict:
    """A view's sections as plain text, for the CLI: heading -> lines."""
    cols = v["columns"]
    inline = {col["slot"]: col["value"] for col in cols}
    show = lambda e: str(e.xreplace(inline))  # noqa: E731
    defines = [col.get("defines", col) for col in cols if not col.get("rate")]
    where = [f"{d['shown']} = {d['value']}" for d in defines if d["shown"] != d["value"]]
    out = {"Note": [v["hint"]]} if v["hint"] else {}
    if v["linearized"]:
        out["Linear approximation about an operating point (x‾, u‾); the state derivatives there, each 0 at an equilibrium"] = \
            [f"{show(f)} = 0" for f in v["operating"]]
    kind = {"power": "power variables", "lagrange": "Lagrange variables"}.get(v["variables"], "deviations" if v["linearized"] else "")
    out[f"States ({kind})" if kind else "States"] = [", ".join(map(str, v["states"]))]
    out["State equations"] = [f"d({x})/dt = {show(f)}" for x, f in zip(v["states"], v["dots"])]
    if v["dependent"]:
        out["Storage in derivative causality (follows the states)"] = [f"{x} = {show(y)}" for x, y in v["dependent"]]
    if v["outputs"]:
        out["Outputs"] = [f"{name} = {show(y)}" for name, y in v["outputs"]]
    if where and v["linearized"]:
        out["Inputs"] = where
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file", help="input .bond file")
    ap.add_argument("--show-equations", action="store_true")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--latex", nargs="?", const="equations", choices=["equations", "matrix"],
                    help="print the state equations as LaTeX, one per state or in matrix form")
    ap.add_argument("--variables", choices=VARIABLES, default="energy",
                    help="energy: states p and q (the default); power: each storage element's effort or flow "
                         "in their place; lagrange: displacements q and flows f, each I element's flow in place of its p")
    ap.add_argument("--linearize", action="store_true",
                    help="show a nonlinear model's linear approximation about an operating point")
    args = ap.parse_args()
    try:
        m = parse(args.file); c = causality(m)
        if args.latex:
            latex = latex_report(m, c, args.variables, args.linearize)
            if args.json:
                print(json.dumps(latex, indent=2)); return
            if latex["note"]:
                sys.exit(latex["note"])
            print(latex["separate"] if args.latex == "equations" else latex["matrix"] or latex["separate"])
            for line in [f"{step['text']}: {step['math']}" for step in latex["reductions"]] + latex["dependent"]:
                print("% " + line)
            if not latex["solved"]:
                sys.exit("Unable to eliminate algebraic variables; run without --latex for the full equations.")
            if latex["operating"]:
                print("% Linear approximation about an operating point; there, each of these is 0 at an equilibrium:")
                print("% " + latex["operating"].replace("\n", " "))
            if latex["hint"]:
                print("% " + latex["hint"])
            if latex["unsolved_outputs"]:
                print(f"% Unable to write {', '.join(latex['unsolved_outputs'])} in terms of the states and inputs.")
            return
        if args.variables != "energy" or args.linearize:
            v = view(m, c, args.variables, args.linearize)
            if v["note"]:
                sys.exit(v["note"])
            sections = view_text(v)
            if args.json:
                print(json.dumps(sections, indent=2)); return
            print("\n\n".join(heading + ":\n" + "\n".join("  " + line for line in lines) for heading, lines in sections.items()))
            if not v["nonlinear"] and args.linearize:
                print("\n(The model is linear, so it's its own linear approximation.)")
            return
        reduction: dict = {}
        outputs: dict = {}
        states, dots, eqs, sol = derive(m, c, reduction=reduction, outputs=outputs)
    except (OSError, ValueError, sp.SympifyError) as exc:
        sys.exit(f"Error: {exc}")
    report = {"causality": {f"{n}.{p}": v for (n,p),v in c.items()}, "states": [str(x) for x in states], "state_equations": [f"d({x})/dt = {y}" for x,y in zip(states,dots)],
              "reductions": [describe(step) for step in reduction.get("steps", [])],
              "dependent_states": [f"{x} = {y}" for x, y in reduction.get("dependent", [])],
              "outputs": [f"{name} = {y}" for name, y in outputs.items()]}
    if args.json:
        print(json.dumps(report, indent=2)); return
    print("Causality (element port -> quantity caused onto junction):")
    for key, value in report["causality"].items(): print(f"  {key}: {value}ffort" if value == "e" else f"  {key}: flow")
    if report["reductions"]:
        print("\nReduced the graph to remove derivative causality:")
        for line in report["reductions"]: print("  " + line)
    print("\nStates:", ", ".join(report["states"]) or "none (all storage has derivative causality)")
    print("\nState equations:")
    for line in report["state_equations"]: print("  " + line)
    if report["dependent_states"]:
        print("\nStorage in derivative causality (follows the states):")
        for line in report["dependent_states"]: print("  " + line)
    if sol is None:
        print("  Unable to eliminate algebraic variables; full equations follow.")
    if report["outputs"]:
        print("\nOutputs:")
        for line in report["outputs"]: print("  " + line)
        stuck = unsolved(m, states, outputs)
        if stuck and sol is not None:
            print(f"  Unable to write {', '.join(stuck)} in terms of the states and inputs.")
    if args.show_equations or sol is None:
        print("\nFull equations" + (" of the reduced graph" if report["reductions"] else "") + " (each equals zero):")
        for equation in eqs: print("  " + str(equation))


if __name__ == "__main__":
    main()
