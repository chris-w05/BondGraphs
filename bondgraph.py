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


def derive(m: Model, c, inputs: dict | None = None, reduction: dict | None = None):
    """Return states, their derivatives, the full equations and the solution (None if unsolved).

    Source values that don't depend on the system are solved as input symbols.
    Pass a dict as `inputs` to get derivatives in terms of those symbols (it's
    filled with symbol -> (source name, value)); otherwise the values are put back.

    Storage in derivative causality is reduced away first where it can be (see
    reduce_graph); the states and derivatives are still the original elements'.
    Pass a dict as `reduction` to get the steps taken ("steps") and each storage
    element that isn't a state as (its state, value in terms of the states) ("dependent").
    """
    red = reduce_graph(m, c)
    if red is None:
        return derive_graph(m, c, inputs)
    states, dots, equations, sol = derive_graph(red.model, red.causality, inputs)
    stuck = set(derivative_storage(m, c))
    integral = {state_symbol(el) for el in m.elements.values() if el.kind.upper() in {"C", "I"} and el.name not in stuck}
    order = {state_symbol(el): i for i, el in enumerate(m.elements.values()) if el.kind.upper() in {"C", "I"}}
    # Each reduced state is a multiple of one original state: keep an integral one.
    picks, dependent = [], []
    for x in states:
        origin = red.origin[x.name[2:]]
        pick = next((y for y in origin if y in integral), next(iter(origin)))
        picks.append((pick, origin[pick]))
        dependent += [(y, sp.simplify(k / origin[pick] * pick)) for y, k in origin.items() if y != pick]
    swap = {x: pick / k for x, (pick, k) in zip(states, picks)}
    result = sorted(((pick, sp.simplify(k * d.xreplace(swap))) for (pick, k), d in zip(picks, dots)), key=lambda r: order[r[0]])
    if reduction is not None:
        reduction.update(steps=red.steps, dependent=sorted(dependent, key=lambda r: order[r[0]]))
    return [x for x, _ in result], [d for _, d in result], equations, sol


def derive_graph(m: Model, c, inputs: dict | None = None):
    """derive() without the reduction."""
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
    system = set(e.values()) | set(f.values()) | {
        sp.Symbol(("q_" if el.kind.upper() == "C" else "p_") + el.name)
        for el in m.elements.values() if el.kind.upper() in {"C", "I"}}
    for el in m.elements.values():
        k = el.kind.upper(); key = (el.name, 1)
        if k in {"R", "C", "I"}:
            form = next((f for f in LAW_FORMS[k] if f in el.params), k)
            own = {"e": e[key], "f": fin[key]}
            if k != "R":
                state = sp.Symbol(("q_" if k == "C" else "p_") + el.name)
                own["q" if k == "C" else "p"] = state
            law = sym(el.params.get(form, f"{k}_{el.name}"), own)
            if form != k: equations.append(own[form] - law)  # e=, f=, q= or p=: the law itself
            elif k == "R": equations.append(e[key] - law * fin[key])
            elif k == "C": equations.append(e[key] - state / law)
            else: equations.append(fin[key] - state / law)
            if k != "R" and (k, c[key]) in {("C", "e"), ("I", "f")}:
                states.append(state); dots.append(fin[key] if k == "C" else e[key])
        elif k in {"SE", "SF"}:
            value = sym(el.params.get("value", "u_" + el.name))
            if not value.free_symbols & system:  # an input, not a modulated source
                given[sp.Dummy("u_" + el.name)] = (el.name, value)
                value = list(given)[-1]
            equations.append((e[key] if k == "SE" else -fin[key]) - value)  # Sf: flow delivered
        # Two-ports: power enters port 1 and leaves port 2 (f2 out = -fin[2]).
        elif k == "TF":
            n = sym(el.params.get("n", "n_" + el.name)); equations += [e[el.name,1] - n*e[el.name,2], -fin[el.name,2] - n*fin[el.name,1]]
        elif k == "GY":
            n = sym(el.params.get("n", "n_" + el.name)); equations += [e[el.name,1] + n*fin[el.name,2], e[el.name,2] - n*fin[el.name,1]]
    variables = list(e.values()) + list(f.values())
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
    """Causality for derive() from which storage elements are in integral causality."""
    c = {}
    for el in m.elements.values():
        if el.kind.upper() in {"C", "I"}:
            wanted, other = ("e", "f") if el.kind.upper() == "C" else ("f", "e")
            c[el.name, 1] = wanted if integral.get(el.name, True) else other
    return c


def state_symbol(el: Element) -> sp.Symbol:
    return sp.Symbol(("q_" if el.kind.upper() == "C" else "p_") + el.name)


def derivative_storage(m: Model, c) -> list[str]:
    return [el.name for el in m.elements.values() if (el.kind.upper(), c.get((el.name, 1))) in {("C", "f"), ("I", "e")}]


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
    return {w if w in m.elements else re.sub(r"_[12]$", "", w) for w in words}


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


def tex_name(name: str, dot: bool = False) -> str:
    """LaTeX for a name: q_C1 -> q_{\\mathrm{C1}}, R_a -> R_{a}, F0 -> F_{0}, rho -> \\rho.

    SymPy's own printer mangles some names (p_I_Electrical ends in "cal", which it
    reads as \\mathcal), so every symbol goes through this instead.
    """
    base, _, sub = name.partition("_")
    digits = re.fullmatch(r"([A-Za-z]+?)(\d+)", base)
    if not sub and digits:
        base, sub = digits.groups()
    head = "\\" + base if base in GREEK else rf"\mathit{{{base}}}" if base.isalpha() and len(base) > 1 else base
    if dot:
        head = rf"\dot{{{head}}}"
    if not sub:
        return head
    sub = sub if len(sub) == 1 or sub.isdigit() else r"\mathrm{%s}" % sub.replace("_", r"\_")
    return f"{head}_{{{sub}}}"


def tex(expr) -> str:
    return sp.latex(expr, symbol_names={s: tex_name(s.name) for s in expr.free_symbols if isinstance(s, sp.Symbol)})


def tex_column(items) -> str:
    return r"\left[\begin{matrix}" + r" \\ ".join(items) + r"\end{matrix}\right]"


def latex_report(m: Model, c) -> dict:
    """The state equations as LaTeX, one per state and in matrix form, for --latex and the editor."""
    given: dict = {}
    reduction: dict = {}
    states, dots, eqs, sol = derive(m, c, given, reduction)
    values = {u: value for u, (_, value) in given.items()}
    pairs = [(tex_name(str(x), dot=True), tex(d.xreplace(values))) for x, d in zip(states, dots)]
    out = {
        "states": [tex_name(str(x)) for x in states],
        "equations": pairs,
        "separate": "\\begin{aligned}\n" + " \\\\\n".join(f"{l} &= {r}" for l, r in pairs) + "\n\\end{aligned}",
        "solved": sol is not None,
        "full": [tex(eq) + " = 0" for eq in eqs],
        "matrix": None, "linear": False, "inputs": [],
        "reductions": [{"text": text, "math": r",\quad ".join(f"{tex_name(name)} = {tex(value)}" for name, value in laws)}
                       for text, laws in reduction.get("steps", [])],
        "dependent": [f"{tex_name(str(x))} = {tex(value)}" for x, value in reduction.get("dependent", [])],
    }
    if not states or sol is None:
        return out
    # Each source becomes an entry of u: its value when that's a lone name or
    # input like v(t), else u_<source> defined underneath.
    shown = {u: value if isinstance(value, (sp.Symbol, AppliedUndef)) else sp.Symbol("u_" + src)
             for u, (src, value) in given.items()}
    columns = list(dict.fromkeys(shown.values()))  # two sources fed by one v(t) share a column
    slot = {name: sp.Dummy() for name in columns}
    split = state_space(states, [d.xreplace({u: slot[shown[u]] for u in given}) for d in dots], list(slot.values()))
    lhs = tex_column(tex_name(str(x), dot=True) for x in states)
    if split is None:
        out["matrix"] = f"{lhs} = {tex(sp.Matrix([d.xreplace(values) for d in dots]))}"
        return out
    A, B, d = split
    rhs = [f"{tex(A)} {tex_column(tex_name(str(x)) for x in states)}"]
    if columns:
        rhs.append(f"{tex(B)} {tex_column(tex(name) for name in columns)}")
    if any(d):
        rhs.append(tex(d))
    inputs = [f"{tex(shown[u])} = {tex(value)}" for u, (_, value) in given.items() if shown[u] != value]
    matrix = f"{lhs} = " + " + ".join(rhs)
    if inputs:
        matrix = f"\\begin{{gathered}}\n{matrix} \\\\\n\\text{{where }} " + ",\\quad ".join(inputs) + "\n\\end{gathered}"
    out.update(matrix=matrix, linear=True, inputs=inputs)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file", help="input .bond file")
    ap.add_argument("--show-equations", action="store_true")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--latex", nargs="?", const="equations", choices=["equations", "matrix"],
                    help="print the state equations as LaTeX, one per state or in matrix form")
    args = ap.parse_args()
    try:
        m = parse(args.file); c = causality(m)
        if args.latex:
            latex = latex_report(m, c)
            if args.json:
                print(json.dumps(latex, indent=2)); return
            print(latex["separate"] if args.latex == "equations" else latex["matrix"] or latex["separate"])
            for line in [f"{step['text']}: {step['math']}" for step in latex["reductions"]] + latex["dependent"]:
                print("% " + line)
            if not latex["solved"]:
                sys.exit("Unable to eliminate algebraic variables; run without --latex for the full equations.")
            return
        reduction: dict = {}
        states, dots, eqs, sol = derive(m, c, reduction=reduction)
    except (OSError, ValueError, sp.SympifyError) as exc:
        sys.exit(f"Error: {exc}")
    report = {"causality": {f"{n}.{p}": v for (n,p),v in c.items()}, "states": [str(x) for x in states], "state_equations": [f"d({x})/dt = {y}" for x,y in zip(states,dots)],
              "reductions": [describe(step) for step in reduction.get("steps", [])],
              "dependent_states": [f"{x} = {y}" for x, y in reduction.get("dependent", [])]}
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
    if args.show_equations or sol is None:
        print("\nFull equations" + (" of the reduced graph" if report["reductions"] else "") + " (each equals zero):")
        for equation in eqs: print("  " + str(equation))


if __name__ == "__main__":
    main()
