/* Bond graph engine: expressions, causality and simulation.
 *
 * Shared by ui/index.html (browser) and the Node tests (ui/engine.test.js).
 * Conventions match bondgraph.py:
 *   - A bond's half-arrow is the reference direction of power. Element laws use
 *     the flow INTO the element, so flipping an arrow never changes the physics.
 *   - C: e = q/C, dq/dt = f_in      I: f_in = p/I, dp/dt = e      R: e = R*f_in
 *     or a law as a function: C e=φ(q) | q=φ(e), I f=φ(p) | p=φ(f), R e=φ(f) | f=φ(e).
 *     In an R, C or I law, e, f, q and p are the element's own effort, inflow and state.
 *   - Se: e = value                  Sf: flow delivered (out of the source) = value
 *   - TF: power in at port 1, out at port 2: e1 = n*e2, f2 = n*f1
 *   - GY: e1 = n*f2, e2 = n*f1
 */
(function (root) {
  "use strict";

  // ---------------------------------------------------------------- expressions

  const FUNCS = {
    sin: Math.sin, cos: Math.cos, tan: Math.tan, asin: Math.asin, acos: Math.acos,
    atan: Math.atan, atan2: Math.atan2, sinh: Math.sinh, cosh: Math.cosh, tanh: Math.tanh,
    exp: Math.exp, log: (x, b) => (b === undefined ? Math.log(x) : Math.log(x) / Math.log(b)),
    sqrt: Math.sqrt, Abs: Math.abs, abs: Math.abs, sign: (x) => (x > 0 ? 1 : x < 0 ? -1 : 0),
    floor: Math.floor, ceiling: Math.ceil, Max: Math.max, Min: Math.min, max: Math.max, min: Math.min,
    Heaviside: (x, h0) => (x > 0 ? 1 : x < 0 ? 0 : h0 === undefined ? 0.5 : h0),
    Mod: (a, b) => a - b * Math.floor(a / b),
  };
  // Bare names SymPy (and so the CLI) reads as something other than a parameter.
  const RESERVED = { I: "the imaginary unit", S: "a SymPy internal", N: "a SymPy function",
    Q: "SymPy's assumption system", O: "SymPy's order symbol", lambda: "a Python keyword" };
  const CONSTS = { pi: Math.PI, E: Math.E, True: 1, False: 0, oo: Infinity };

  function tokenize(src) {
    const out = [];
    let i = 0;
    while (i < src.length) {
      const ch = src[i];
      if (/\s/.test(ch)) { i++; continue; }
      const num = /^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/.exec(src.slice(i));
      if (num) { out.push({ t: "num", v: parseFloat(num[0]), at: i }); i += num[0].length; continue; }
      const id = /^[A-Za-z_]\w*/.exec(src.slice(i));
      if (id) { out.push({ t: "id", v: id[0], at: i }); i += id[0].length; continue; }
      const two = src.slice(i, i + 2);
      if (["**", "<=", ">=", "==", "!="].includes(two)) { out.push({ t: "op", v: two, at: i }); i += 2; continue; }
      if ("+-*/^()<>,&|~".includes(ch)) { out.push({ t: "op", v: ch, at: i }); i++; continue; }
      throw new SyntaxError(`unexpected character '${ch}' at position ${i + 1}`);
    }
    return out;
  }

  // Recursive-descent parser producing a small AST.
  function parse(src) {
    if (typeof src !== "string" || !src.trim()) throw new SyntaxError("empty expression");
    const toks = tokenize(src);
    let k = 0;
    const peek = () => toks[k];
    const isOp = (v) => toks[k] && toks[k].t === "op" && toks[k].v === v;
    const expect = (v) => {
      if (!isOp(v)) {
        const got = toks[k] ? `'${toks[k].v}'` : "end of expression";
        throw new SyntaxError(`expected '${v}' but found ${got}`);
      }
      k++;
    };
    function orExpr() {
      let a = andExpr();
      while (isOp("|")) { k++; a = { k: "or", a, b: andExpr() }; }
      return a;
    }
    function andExpr() {
      let a = notExpr();
      while (isOp("&")) { k++; a = { k: "and", a, b: notExpr() }; }
      return a;
    }
    function notExpr() {
      if (isOp("~")) { k++; return { k: "not", a: notExpr() }; }
      return cmp();
    }
    function cmp() {
      const a = sum();
      const t = peek();
      if (t && t.t === "op" && ["<", ">", "<=", ">=", "==", "!="].includes(t.v)) {
        k++;
        return { k: "cmp", op: t.v, a, b: sum() };
      }
      return a;
    }
    function sum() {
      let a = prod();
      while (isOp("+") || isOp("-")) { const op = toks[k++].v; a = { k: "bin", op, a, b: prod() }; }
      return a;
    }
    function prod() {
      let a = unary();
      while (isOp("*") || isOp("/")) { const op = toks[k++].v; a = { k: "bin", op, a, b: unary() }; }
      return a;
    }
    function unary() {
      if (isOp("-")) { k++; return { k: "neg", a: unary() }; }
      if (isOp("+")) { k++; return unary(); }
      return power();
    }
    function power() {
      const base = atom();
      if (isOp("**") || isOp("^")) { k++; return { k: "bin", op: "^", a: base, b: unary() }; }
      return base;
    }
    function atom() {
      const t = peek();
      if (!t) throw new SyntaxError("expression ends too early");
      if (t.t === "num") { k++; return { k: "num", v: t.v }; }
      if (t.t === "id") {
        k++;
        if (isOp("(")) {
          k++;
          const args = [];
          if (!isOp(")")) {
            args.push(orExpr());
            while (isOp(",")) { k++; args.push(orExpr()); }
          }
          expect(")");
          return { k: "call", name: t.v, args };
        }
        return { k: "id", name: t.v };
      }
      if (isOp("(")) {
        k++;
        const items = [orExpr()];
        while (isOp(",")) { k++; items.push(orExpr()); }
        expect(")");
        return items.length === 1 ? items[0] : { k: "tuple", items };
      }
      throw new SyntaxError(`unexpected '${t.v}'`);
    }
    const ast = orExpr();
    if (k < toks.length) throw new SyntaxError(`unexpected '${toks[k].v}' after the expression`);
    return ast;
  }

  // Classify every name an expression uses: variables (states, bond variables,
  // parameters, t) and called input functions such as u(t).
  function names(ast, acc) {
    acc = acc || { ids: new Set(), inputs: new Set() };
    (function walk(n) {
      switch (n.k) {
        case "id": if (!(n.name in CONSTS)) acc.ids.add(n.name); break;
        case "call":
          if (n.name !== "Piecewise" && !(n.name in FUNCS)) acc.inputs.add(n.name);
          n.args.forEach(walk); break;
        case "tuple": n.items.forEach(walk); break;
        case "num": break;
        default: if (n.a) walk(n.a); if (n.b) walk(n.b);
      }
    })(ast);
    return acc;
  }

  // Compile an AST to a closure fn(ctx). `resolve(name)` returns a closure for a
  // variable; `callInput(name, argFn)` returns a closure for an input call.
  function compile(ast, resolve, callInput) {
    function c(n) {
      switch (n.k) {
        case "num": { const v = n.v; return () => v; }
        case "id": {
          if (n.name in CONSTS) { const v = CONSTS[n.name]; return () => v; }
          return resolve(n.name);
        }
        case "neg": { const a = c(n.a); return (x) => -a(x); }
        case "not": { const a = c(n.a); return (x) => (a(x) ? 0 : 1); }
        case "and": { const a = c(n.a), b = c(n.b); return (x) => (a(x) && b(x) ? 1 : 0); }
        case "or": { const a = c(n.a), b = c(n.b); return (x) => (a(x) || b(x) ? 1 : 0); }
        case "cmp": {
          const a = c(n.a), b = c(n.b);
          switch (n.op) {
            case "<": return (x) => (a(x) < b(x) ? 1 : 0);
            case ">": return (x) => (a(x) > b(x) ? 1 : 0);
            case "<=": return (x) => (a(x) <= b(x) ? 1 : 0);
            case ">=": return (x) => (a(x) >= b(x) ? 1 : 0);
            case "==": return (x) => (a(x) === b(x) ? 1 : 0);
            default: return (x) => (a(x) !== b(x) ? 1 : 0);
          }
        }
        case "bin": {
          const a = c(n.a), b = c(n.b);
          switch (n.op) {
            case "+": return (x) => a(x) + b(x);
            case "-": return (x) => a(x) - b(x);
            case "*": return (x) => a(x) * b(x);
            case "/": return (x) => a(x) / b(x);
            default: return (x) => Math.pow(a(x), b(x));
          }
        }
        case "tuple": throw new SyntaxError("a bracketed list is only allowed inside Piecewise(...)");
        case "call": {
          if (n.name === "Piecewise") {
            if (!n.args.length) throw new SyntaxError("Piecewise needs at least one (value, condition) pair");
            const pairs = n.args.map((p) => {
              if (p.k !== "tuple" || p.items.length !== 2) throw new SyntaxError("each Piecewise argument must be (value, condition)");
              return [c(p.items[0]), c(p.items[1])];
            });
            // When x.fz is set (a step retried with frozen switches), branch
            // conditions read the start-of-step values so the law is smooth.
            return (x) => {
              const cx = x.fz || x;
              for (let i = 0; i < pairs.length; i++) if (pairs[i][1](cx)) return pairs[i][0](x);
              return NaN;
            };
          }
          const args = n.args.map(c);
          if (n.name === "Heaviside" || n.name === "sign") {
            const fn = FUNCS[n.name], a = args[0], b = args[1];
            return b ? (x) => fn(a(x.fz || x), b(x)) : (x) => fn(a(x.fz || x));
          }
          if (n.name in FUNCS) {
            const fn = FUNCS[n.name];
            if (args.length === 1) { const a = args[0]; return (x) => fn(a(x)); }
            if (args.length === 2) { const a = args[0], b = args[1]; return (x) => fn(a(x), b(x)); }
            return (x) => fn(...args.map((g) => g(x)));
          }
          if (args.length !== 1) throw new SyntaxError(`input ${n.name}(...) takes one argument, usually t`);
          return callInput(n.name, args[0]);
        }
      }
      throw new SyntaxError("unsupported expression");
    }
    return c(ast);
  }

  function checkReserved(ast) {
    const { ids } = names(ast);
    for (const id of ids) if (id in RESERVED) {
      throw new SyntaxError(`'${id}' is reserved (${RESERVED[id]} in SymPy); use another name`);
    }
  }

  // Evaluate a constant expression (parameter values, initial conditions).
  function evalConst(src, params) {
    const ast = parse(String(src));
    return compile(ast, (name) => {
      if (params && name in params) { const v = params[name]; return () => v; }
      throw new SyntaxError(`unknown name '${name}'`);
    }, (name) => { throw new SyntaxError(`input ${name}(t) can't be used here`); })({});
  }

  // ---------------------------------------------------------------- graph model

  const STORAGE = { C: true, I: true };
  const TWO_PORT = { TF: true, GY: true };
  const JUNCTION = { "0": true, "1": true };
  const isElement = (kind) => !JUNCTION[kind];
  const maxBonds = (kind) => (JUNCTION[kind] ? Infinity : TWO_PORT[kind] ? 2 : 1);

  // Law keys per element kind (same keys as the .bond format). The first is the
  // default: a coefficient (C, I, R, n) or a source value. The others give the
  // law as a function, keyed by the variable it sets.
  const LAW_FORMS = { C: ["C", "e", "q"], I: ["I", "f", "p"], R: ["R", "e", "f"], Se: ["value"], Sf: ["value"], TF: ["n"], GY: ["n"] };
  const LAW_KEY = Object.fromEntries(Object.entries(LAW_FORMS).map(([k, forms]) => [k, forms[0]]));
  // Names an R, C or I law uses for the element's own effort, inflow and state.
  const OWN_VARS = { C: ["e", "f", "q"], I: ["e", "f", "p"], R: ["e", "f"] };

  // Key of the law a node uses: the first form it has, else the default.
  function lawKeyOf(node) {
    const forms = LAW_FORMS[node.kind];
    if (!forms) return undefined;
    return forms.find((k) => node.law && node.law[k] !== undefined) || forms[0];
  }

  function stateName(node) {
    return node.kind === "C" ? "q_" + node.name : node.kind === "I" ? "p_" + node.name : null;
  }

  // Port of `node` used by bond b (1 for one-ports and junctions).
  function portOf(b, nodeId) {
    if (b.from === nodeId && b.fromPort) return b.fromPort;
    if (b.to === nodeId && b.toPort) return b.toPort;
    return 1;
  }

  // Names that refer to bond b: e_X/f_X for element X's port 1, e_X_2/f_X_2 for
  // port 2, and e_bondN/f_bondN for a bond between two junctions.
  function bondAliases(graph, b) {
    const out = [];
    for (const end of ["from", "to"]) {
      const n = graph.nodes.find((x) => x.id === b[end]);
      if (n && isElement(n.kind)) {
        const p = portOf(b, n.id);
        out.push(p === 2 ? n.name + "_2" : n.name);
      }
    }
    if (!out.length) out.push("bond" + b.id);
    return out;
  }

  function incident(graph, nodeId) {
    return graph.bonds.filter((b) => b.from === nodeId || b.to === nodeId);
  }

  // Structural problems that block simulation (warnings are shown but tolerated).
  function validate(graph) {
    const errors = [], warnings = [];
    const seen = new Map();
    for (const n of graph.nodes) {
      if (!/^[A-Za-z_]\w*$/.test(n.name)) errors.push({ node: n.id, msg: `“${n.name}” is not a valid name (letters, digits and _ only)` });
      if (seen.has(n.name)) errors.push({ node: n.id, msg: `Two nodes are named ${n.name}` });
      seen.set(n.name, n.id);
      const deg = incident(graph, n.id).length;
      if (TWO_PORT[n.kind] && deg === 1) errors.push({ node: n.id, msg: `${n.name} needs a bond on both ports` });
      else if (deg === 0) warnings.push({ node: n.id, msg: `${n.name} isn't bonded to anything, so it's left out` });
    }
    return { errors, warnings };
  }

  // ---------------------------------------------------------------- causality
  //
  // One boolean per bond: x = 1 when the causal stroke sits at the bond's `to`
  // end (that end receives effort and returns flow). Local rules become
  // counting constraints; storage elements prefer integral causality. A
  // branch-and-bound search minimises violated rules first, then maximises the
  // number of storage elements in integral causality.

  function assignCausality(graph, opts) {
    const limit = (opts && opts.limit) || 400000;
    const bonds = graph.bonds;
    const idx = new Map(bonds.map((b, i) => [b.id, i]));
    const nodeById = new Map(graph.nodes.map((n) => [n.id, n]));
    // literal: [bondIndex, valueThatMeansStrokeAtNode]
    const lit = (b, nodeId) => [idx.get(b.id), b.to === nodeId ? 1 : 0];
    const cons = [];   // { node, lits, allowed:Set }
    const prefs = [];  // { node, lit, want }  want = stroke-at-node wanted (1/0)

    for (const n of graph.nodes) {
      const inc = incident(graph, n.id);
      if (!inc.length) continue;
      const lits = inc.map((b) => lit(b, n.id));
      switch (n.kind) {
        case "Se": cons.push({ node: n.id, lits, allowed: new Set([0]), why: "an effort source must set the effort" }); break;
        case "Sf": cons.push({ node: n.id, lits, allowed: new Set([1]), why: "a flow source must set the flow" }); break;
        case "R": {
          const cz = n.law && n.law.causality;
          if (cz === "effort") cons.push({ node: n.id, lits, allowed: new Set([0]), why: "this resistor is fixed to set effort" });
          if (cz === "flow") cons.push({ node: n.id, lits, allowed: new Set([1]), why: "this resistor is fixed to set flow" });
          break;
        }
        case "C": prefs.push({ node: n.id, lit: lits[0], want: 0 }); break;
        case "I": prefs.push({ node: n.id, lit: lits[0], want: 1 }); break;
        case "0": cons.push({ node: n.id, lits, allowed: new Set([1]), why: "exactly one bond may set the effort on a 0-junction" }); break;
        case "1": cons.push({ node: n.id, lits, allowed: new Set([lits.length - 1]), why: "exactly one bond may set the flow on a 1-junction" }); break;
        case "TF": if (lits.length === 2) cons.push({ node: n.id, lits, allowed: new Set([1]), why: "a transformer passes effort straight through" }); break;
        case "GY": if (lits.length === 2) cons.push({ node: n.id, lits, allowed: new Set([0, 2]), why: "a gyrator turns effort on one side into flow on the other" }); break;
      }
    }

    const nb = bonds.length;
    const byBond = Array.from({ length: nb }, () => []);
    cons.forEach((c, ci) => c.lits.forEach(([bi]) => byBond[bi].push(ci)));
    const prefByBond = Array.from({ length: nb }, () => []);
    prefs.forEach((p, pi) => prefByBond[p.lit[0]].push(pi));

    // Visit bonds outward from sources and fixed elements so rules close early.
    const order = [], placed = new Array(nb).fill(false);
    const seeds = graph.nodes.filter((n) => n.kind === "Se" || n.kind === "Sf" || (n.kind === "R" && n.law && n.law.causality))
      .concat(graph.nodes.filter((n) => STORAGE[n.kind]), graph.nodes);
    for (const s of seeds) {
      const queue = [s.id];
      const seenNode = new Set();
      while (queue.length) {
        const id = queue.shift();
        if (seenNode.has(id)) continue;
        seenNode.add(id);
        for (const b of incident(graph, id)) {
          const bi = idx.get(b.id);
          if (!placed[bi]) { placed[bi] = true; order.push(bi); }
          queue.push(b.from === id ? b.to : b.from);
        }
      }
    }

    const x = new Int8Array(nb).fill(-1);
    const trues = new Int32Array(cons.length), open = new Int32Array(cons.length);
    cons.forEach((c, ci) => { open[ci] = c.lits.length; });
    const feasible = (ci) => {
      const lo = trues[ci], hi = trues[ci] + open[ci];
      for (const a of cons[ci].allowed) if (a >= lo && a <= hi) return true;
      return false;
    };
    let violations = 0, integral = 0, openStorage = prefs.length;
    let best = null, bestViol = Infinity, bestInt = -1, steps = 0, truncated = false;

    function setBond(bi, v, sign) {
      for (const ci of byBond[bi]) {
        const before = feasible(ci);
        for (const [lb, lv] of cons[ci].lits) if (lb === bi) { open[ci] -= sign; if (lv === v) trues[ci] += sign; }
        const after = feasible(ci);
        if (before !== after) violations += after ? -1 : 1;
      }
      for (const pi of prefByBond[bi]) {
        const p = prefs[pi];
        openStorage -= sign;
        if ((v === p.lit[1] ? 1 : 0) === p.want) integral += sign;
      }
    }

    function preferred(bi) {
      for (const pi of prefByBond[bi]) { const p = prefs[pi]; return p.want ? p.lit[1] : 1 - p.lit[1]; }
      for (const ci of byBond[bi]) {
        const c = cons[ci];
        if (c.lits.length === 1 && c.allowed.size === 1) {
          const want = [...c.allowed][0];
          return want ? c.lits[0][1] : 1 - c.lits[0][1];
        }
      }
      return 0;
    }

    function search(depth) {
      if (++steps > limit) { truncated = true; return; }
      if (violations > bestViol || (violations === bestViol && integral + openStorage <= bestInt)) return;
      if (depth === nb) { best = Array.from(x); bestViol = violations; bestInt = integral; return; }
      const bi = order[depth];
      const first = preferred(bi);
      for (const v of [first, 1 - first]) {
        x[bi] = v; setBond(bi, v, 1);
        search(depth + 1);
        setBond(bi, v, -1); x[bi] = -1;
        if (truncated) return;
      }
    }
    if (nb) search(0); else { best = []; bestViol = 0; bestInt = 0; }

    // Report: stroke end per bond, storage status, violated rules.
    const stroke = {};
    bonds.forEach((b, i) => { stroke[b.id] = best && best[i] === 1 ? b.to : b.from; });
    const storage = {};
    for (const p of prefs) {
      const atNode = best ? (best[p.lit[0]] === p.lit[1] ? 1 : 0) : 0;
      storage[p.node] = atNode === p.want ? "integral" : "derivative";
    }
    const conflicts = [];
    if (best) {
      cons.forEach((c) => {
        const t = c.lits.reduce((s, [bi, v]) => s + (best[bi] === v ? 1 : 0), 0);
        if (!c.allowed.has(t)) conflicts.push({ node: c.node, msg: `${nodeById.get(c.node).name}: ${c.why}` });
      });
    }
    return { stroke, storage, conflicts, truncated };
  }

  // ---------------------------------------------------------------- simulation

  function luFactor(A, n) {
    const piv = new Int32Array(n);
    for (let k = 0; k < n; k++) {
      let p = k, m = Math.abs(A[k * n + k]);
      for (let i = k + 1; i < n; i++) { const v = Math.abs(A[i * n + k]); if (v > m) { m = v; p = i; } }
      if (!(m > 1e-300)) return null;
      piv[k] = p;
      if (p !== k) for (let j = 0; j < n; j++) { const t = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = t; }
      const d = A[k * n + k];
      for (let i = k + 1; i < n; i++) {
        const f = (A[i * n + k] /= d);
        if (f) for (let j = k + 1; j < n; j++) A[i * n + j] -= f * A[k * n + j];
      }
    }
    return piv;
  }
  function luSolve(A, piv, n, b) {
    for (let k = 0; k < n; k++) { const p = piv[k]; if (p !== k) { const t = b[k]; b[k] = b[p]; b[p] = t; } }
    for (let i = 1; i < n; i++) { let s = b[i]; for (let j = 0; j < i; j++) s -= A[i * n + j] * b[j]; b[i] = s; }
    for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < n; j++) s -= A[i * n + j] * b[j]; b[i] = s / A[i * n + i]; }
    return b;
  }

  // Build the residual system for the bonded part of the graph.
  //   unknowns v = [e_b, f_b for every bond..., one state per C/I]
  //   residuals  = junction rules, element laws, and one integration rule per state
  function buildModel(graph, settings) {
    settings = settings || {};
    const nodes = graph.nodes.filter((n) => incident(graph, n.id).length > 0);
    const bonds = graph.bonds;
    const nB = bonds.length;
    const bi = new Map(bonds.map((b, i) => [b.id, i]));
    const E = (b) => 2 * bi.get(b.id), F = (b) => 2 * bi.get(b.id) + 1;

    const sym = new Map();          // name -> unknown index
    const symInfo = [];             // index -> { name, kind, node }
    bonds.forEach((b) => {
      const aliases = bondAliases(graph, b);
      aliases.forEach((a) => { sym.set("e_" + a, E(b)); sym.set("f_" + a, F(b)); });
      symInfo[E(b)] = { name: "e_" + aliases[0], kind: "effort", bond: b.id };
      symInfo[F(b)] = { name: "f_" + aliases[0], kind: "flow", bond: b.id };
    });
    const states = [];
    for (const n of nodes) if (STORAGE[n.kind]) {
      const i = 2 * nB + states.length;
      const name = stateName(n);
      states.push({ node: n.id, name, index: i, kind: n.kind });
      sym.set(name, i);
      symInfo[i] = { name, kind: "state", node: n.id };
    }
    const N = 2 * nB + states.length;

    // Parameters & inputs referenced by the laws.
    const paramNames = new Set(), inputNames = new Set();
    const lawAst = new Map(), errors = [];
    for (const n of nodes) {
      if (!LAW_KEY[n.kind]) continue;
      const src = (n.law && n.law[lawKeyOf(n)]) || "";
      const own = OWN_VARS[n.kind] || [];
      try {
        const ast = parse(src);
        checkReserved(ast);
        lawAst.set(n.id, ast);
        const nm = names(ast);
        nm.ids.forEach((id) => { if (id !== "t" && !sym.has(id) && !own.includes(id)) paramNames.add(id); });
        nm.inputs.forEach((id) => inputNames.add(id));
      } catch (err) {
        errors.push({ node: n.id, msg: `${n.name}: ${err.message}` });
      }
    }
    // Inputs may themselves use parameters.
    const inputAst = new Map();
    for (const name of inputNames) {
      const src = (settings.inputs && settings.inputs[name]) || "";
      try {
        const ast = parse(src || "1"); // undefined inputs default to a unit step
        checkReserved(ast);
        inputAst.set(name, ast);
        names(ast).ids.forEach((id) => { if (id !== "t" && !sym.has(id)) paramNames.add(id); });
        if (names(ast).inputs.size) throw new SyntaxError("an input can't call another input");
      } catch (err) {
        errors.push({ input: name, msg: `${name}(t): ${err.message}` });
      }
    }
    const params = [...paramNames].sort();
    const P = new Float64Array(params.length);
    const pIndex = new Map(params.map((p, i) => [p, i]));
    params.forEach((p, i) => {
      const raw = settings.params && settings.params[p];
      try { P[i] = raw === undefined || raw === "" ? 1 : evalConst(raw); }
      catch (err) { errors.push({ param: p, msg: `${p}: ${err.message}` }); P[i] = NaN; }
      if (!Number.isFinite(P[i]) && !errors.some((e) => e.param === p)) errors.push({ param: p, msg: `${p} must be a finite number` });
    });

    const resolve = (name) => {
      if (name === "t") return (c) => c.t;
      if (sym.has(name)) { const i = sym.get(name); return (c) => c.v[i]; }
      const i = pIndex.get(name); return (c) => c.P[i];
    };
    const inputFns = new Map();
    for (const [name, ast] of inputAst) {
      try {
        const scratch = { v: null, P: null, t: 0 };
        const body = compile(ast, resolve, () => () => NaN);
        inputFns.set(name, (c, tt) => { scratch.v = c.v; scratch.P = c.P; scratch.t = tt; return body(scratch); });
      } catch (err) { errors.push({ input: name, msg: `${name}(t): ${err.message}` }); }
    }
    const callInput = (name, argFn) => (c) => { const f = inputFns.get(name); return f ? f(c, argFn(c)) : NaN; };
    // In an R, C or I law, e and f are the element's own effort and inflow (so
    // the law doesn't depend on its bond's arrow) and q or p its own state.
    function ownVars(n) {
      if (!OWN_VARS[n.kind]) return null;
      const b = incident(graph, n.id)[0], e = E(b), f = F(b), s = b.to === n.id ? 1 : -1;
      const own = new Map([["e", (c) => c.v[e]], ["f", (c) => s * c.v[f]]]);
      if (STORAGE[n.kind]) { const i = sym.get(stateName(n)); own.set(n.kind === "C" ? "q" : "p", (c) => c.v[i]); }
      return own;
    }
    const law = new Map();
    for (const [id, ast] of lawAst) {
      const n = graph.nodes.find((x) => x.id === id), own = ownVars(n);
      try { law.set(id, compile(ast, own ? (name) => own.get(name) || resolve(name) : resolve, callInput)); }
      catch (err) { errors.push({ node: id, msg: `${n.name}: ${err.message}` }); }
    }

    // Residuals, each r(c) with c = { v, P, t }.
    const res = [];
    for (const n of nodes) {
      const inc = incident(graph, n.id);
      // +1 when the bond points into n
      const into = (b) => (b.to === n.id ? 1 : -1);
      if (JUNCTION[n.kind]) {
        const common = n.kind === "0" ? E : F, summed = n.kind === "0" ? F : E;
        const b0 = common(inc[0]);
        for (let k = 1; k < inc.length; k++) { const j = common(inc[k]); res.push((c) => c.v[j] - c.v[b0]); }
        // signed sum of power-carrying variable, + for bonds leaving the junction
        const terms = inc.map((b) => [summed(b), -into(b)]);
        res.push((c) => { let s = 0; for (const [j, sg] of terms) s += sg * c.v[j]; return s; });
        continue;
      }
      const L = law.get(n.id) || (() => NaN);
      if (TWO_PORT[n.kind]) {
        const p1 = inc.find((b) => portOf(b, n.id) === 1), p2 = inc.find((b) => portOf(b, n.id) === 2);
        if (!p1 || !p2) continue; // reported by validate()
        const e1 = E(p1), f1 = F(p1), s1 = into(p1), e2 = E(p2), f2 = F(p2), s2 = into(p2);
        // flow into port 1 and out of port 2
        const fi1 = (c) => s1 * c.v[f1], fo2 = (c) => -s2 * c.v[f2];
        if (n.kind === "TF") {
          res.push((c) => c.v[e1] - L(c) * c.v[e2]);
          res.push((c) => fo2(c) - L(c) * fi1(c));
        } else {
          res.push((c) => c.v[e1] - L(c) * fo2(c));
          res.push((c) => c.v[e2] - L(c) * fi1(c));
        }
        continue;
      }
      const b = inc[0], e = E(b), f = F(b), s = into(b), form = lawKeyOf(n);
      switch (n.kind) {
        case "C": {
          const q = sym.get(stateName(n));
          res.push(form === "e" ? (c) => c.v[e] - L(c) : form === "q" ? (c) => c.v[q] - L(c) : (c) => c.v[e] - c.v[q] / L(c));
          break;
        }
        case "I": {
          const p = sym.get(stateName(n));
          res.push(form === "f" ? (c) => s * c.v[f] - L(c) : form === "p" ? (c) => c.v[p] - L(c) : (c) => s * c.v[f] - c.v[p] / L(c));
          break;
        }
        case "R":
          res.push(form === "e" ? (c) => c.v[e] - L(c) : form === "f" ? (c) => s * c.v[f] - L(c) : (c) => c.v[e] - L(c) * s * c.v[f]);
          break;
        case "Se": res.push((c) => c.v[e] - L(c)); break;
        case "Sf": res.push((c) => -s * c.v[f] - L(c)); break;
      }
    }
    // State derivatives: dq/dt = f_in for C, dp/dt = e for I.
    const deriv = states.map((st) => {
      const n = nodes.find((x) => x.id === st.node);
      const b = incident(graph, n.id)[0];
      const s = b.to === n.id ? 1 : -1;
      return n.kind === "C" ? ((j) => (c) => s * c.v[j])(F(b)) : ((j) => (c) => c.v[j])(E(b));
    });

    return { N, nB, states, deriv, res, P, params, inputs: [...inputNames].sort(), symInfo, sym, errors };
  }

  // Implicit simulation: BDF2 with a backward-Euler start, Newton on the full
  // set of efforts, flows and states. Storage in derivative causality and
  // algebraic loops are handled by the same solve.
  function simulate(graph, settings, causality) {
    const t0 = Date.now();
    const model = buildModel(graph, settings);
    if (model.errors.length) return { ok: false, errors: model.errors, model };
    const { N, nB, states, deriv, res, P } = model;
    if (!states.length && !nB) return { ok: false, errors: [{ msg: "Add elements and bonds to simulate." }], model };
    const tEnd = Number(settings.tEnd), steps = Math.max(10, Math.min(20000, Math.round(Number(settings.steps) || 1000)));
    if (!(tEnd > 0)) return { ok: false, errors: [{ msg: "End time must be a positive number." }], model };
    const h = tEnd / steps;
    const M = res.length + states.length;
    if (M !== N) return { ok: false, errors: [{ msg: `The equations don't balance (${M} equations, ${N} unknowns). Check that every element is fully bonded.` }], model };

    const integral = new Set(states.filter((s) => !causality || causality.storage[s.node] !== "derivative").map((s) => s.index));
    const x0 = new Float64Array(N);
    for (const st of states) {
      if (!integral.has(st.index)) continue;
      const raw = settings.init && settings.init[st.name];
      try { x0[st.index] = raw === undefined || raw === "" ? 0 : evalConst(raw, Object.fromEntries(model.params.map((p, i) => [p, P[i]]))); }
      catch (err) { return { ok: false, errors: [{ state: st.name, msg: `Initial ${st.name}: ${err.message}` }], model }; }
    }

    const ctx = { v: null, P, t: 0, fz: null };
    const fz = { v: new Float64Array(N), P, t: 0, fz: null }; // frozen-switch context
    const r = new Float64Array(N), rp = new Float64Array(N);
    const J = new Float64Array(N * N);
    let piv = null, jacFor = null;

    // Residual of the step equations. mode: {init:true} or {alpha, hist} where
    // the integration rule is x - hist - alpha*xdot = 0; mode.frozen holds
    // switch conditions at the start-of-step values in fz.v.
    function residual(v, t, mode, out) {
      ctx.v = v; ctx.t = t; fz.t = t;
      ctx.fz = mode.frozen ? fz : null;
      for (let i = 0; i < res.length; i++) out[i] = res[i](ctx);
      for (let k = 0; k < states.length; k++) {
        const st = states[k];
        if (mode.init) out[res.length + k] = integral.has(st.index) ? v[st.index] - x0[st.index] : deriv[k](ctx);
        else out[res.length + k] = v[st.index] - mode.hist[k] - mode.alpha * deriv[k](ctx);
      }
      return out;
    }
    function jacobian(v, t, mode) {
      residual(v, t, mode, r);
      const w = Float64Array.from(v);
      for (let j = 0; j < N; j++) {
        const d = 1e-7 * Math.max(1, Math.abs(v[j]));
        w[j] = v[j] + d;
        residual(w, t, mode, rp);
        for (let i = 0; i < N; i++) J[i * N + j] = (rp[i] - r[i]) / d;
        w[j] = v[j];
      }
      return luFactor(J, N);
    }
    const norm = (a) => { let m = 0; for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i])); return m; };

    // Modified Newton: reuse the factored Jacobian while it keeps converging.
    function newton(v, t, mode, key) {
      const start = Float64Array.from(v);
      if (mode.frozen) key += "/frozen";
      for (let attempt = 0; attempt < 2; attempt++) {
        if (attempt === 1) v.set(start);
        if (!piv || jacFor !== key || attempt === 1) {
          piv = jacobian(v, t, mode); jacFor = key;
          if (!piv) return { ok: false, singular: true };
        }
        let prev = Infinity;
        for (let it = 0; it < 30; it++) {
          residual(v, t, mode, r);
          if (!Number.isFinite(norm(r))) break;
          const dx = luSolve(J, piv, N, Float64Array.from(r));
          for (let i = 0; i < N; i++) v[i] -= dx[i];
          const step = norm(dx);
          if (step <= 1e-10 * norm(v) + 1e-300) return { ok: true };
          if (it > 3 && step > 0.5 * prev) break; // slow: refresh the Jacobian
          prev = step;
        }
      }
      return { ok: false };
    }

    // Consistent start: integral states fixed; dependent storage starts at rest.
    const v = Float64Array.from(x0);
    let out = newton(v, 0, { init: true }, "init");
    if (!out.ok && !out.singular) {
      v.set(x0); fz.v.set(x0);
      out = newton(v, 0, { init: true, frozen: true }, "init");
    }
    if (!out.ok) {
      return { ok: false, model, errors: [{ msg: out.singular
        ? "The equations are singular at t = 0. Usually this means a causal conflict (for example two sources fighting over one junction) or a law that's zero where it divides."
        : "Couldn't find consistent starting values at t = 0. Check the laws for division by zero or undefined Piecewise branches." }] };
    }
    const T = new Float64Array(steps + 1);
    const X = [Float64Array.from(v)];
    let xPrev = null, xCur = Float64Array.from(v);
    const hist = new Float64Array(states.length);
    for (let n = 1; n <= steps; n++) {
      const t = n * h;
      T[n] = t;
      const guess = Float64Array.from(xCur);
      if (xPrev) for (let i = 0; i < N; i++) guess[i] = 2 * xCur[i] - xPrev[i];
      let mode;
      if (xPrev) {
        states.forEach((st, k) => { hist[k] = (4 * xCur[st.index] - xPrev[st.index]) / 3; });
        mode = { alpha: (2 / 3) * h, hist };
      } else {
        states.forEach((st, k) => { hist[k] = xCur[st.index]; });
        mode = { alpha: h, hist };
      }
      const key = xPrev ? "bdf2" : "be";
      const start = Float64Array.from(guess);
      out = newton(guess, t, mode, key);
      if (!out.ok) {
        // A switch inside the step can leave the implicit equation without a
        // solution. Hold the switches at their start-of-step state and retry.
        guess.set(start); fz.v.set(xCur);
        out = newton(guess, t, Object.assign({}, mode, { frozen: true }), key);
      }
      if (!out.ok) {
        // Retry the interval with smaller backward-Euler substeps.
        let ok = true;
        const sub = Float64Array.from(xCur);
        const m = 8;
        for (let s = 1; s <= m && ok; s++) {
          states.forEach((st, k) => { hist[k] = sub[st.index]; });
          fz.v.set(sub);
          ok = newton(sub, t - h + (s * h) / m, { alpha: h / m, hist, frozen: true }, "sub").ok;
        }
        jacFor = null;
        if (!ok) {
          return { ok: true, partial: true, model, T: T.slice(0, n), X,
            errors: [{ msg: `The solver stopped at t = ${fmt(t - h)}: Newton's method didn't converge. A law may be discontinuous or undefined here, or the step is too large.` }],
            ms: Date.now() - t0 };
        }
        guess.set(sub);
        xPrev = null; // restart BDF2 from the substepped point
      } else {
        xPrev = xCur;
      }
      xCur = guess;
      X.push(Float64Array.from(xCur));
      if (!Number.isFinite(norm(xCur))) {
        return { ok: true, partial: true, model, T: T.slice(0, n + 1), X, ms: Date.now() - t0,
          errors: [{ msg: `The response became infinite near t = ${fmt(t)}.` }] };
      }
    }
    return { ok: true, model, T, X, errors: [], ms: Date.now() - t0 };
  }

  function fmt(x) {
    if (!Number.isFinite(x)) return String(x);
    const a = Math.abs(x);
    if (a !== 0 && (a < 1e-3 || a >= 1e5)) return x.toExponential(2);
    return String(+x.toPrecision(4));
  }

  // ---------------------------------------------------------------- .bond files

  // Export: plain bonds between two junctions become unit transformers named
  // bondN; bonds between two elements go through a two-bond 1-junction named bondN.
  function toBond(graph, settings) {
    const byId = new Map(graph.nodes.map((n) => [n.id, n]));
    const lines = ["# Bond graph exported from the bond graph editor (ui/index.html).", ""];
    const members = new Map();     // junction name -> [tokens]
    const extraEls = [], connects = [];
    const addMember = (j, tok) => { if (!members.has(j)) members.set(j, []); members.get(j).push(tok); };
    for (const n of graph.nodes) if (JUNCTION[n.kind]) members.set(n.name, []);
    for (const b of graph.bonds) {
      const A = byId.get(b.from), B = byId.get(b.to);
      const tok = (n) => (TWO_PORT[n.kind] ? `${n.name}.${portOf(b, n.id)}` : n.name);
      if (JUNCTION[A.kind] && JUNCTION[B.kind]) {
        extraEls.push(`element bond${b.id} TF n=1  # plain bond ${A.name} -> ${B.name}`);
        connects.push(`connect bond${b.id} ${A.name} ${B.name}`);
      } else if (JUNCTION[A.kind]) addMember(A.name, tok(B) + "+");
      else if (JUNCTION[B.kind]) addMember(B.name, tok(A) + "-");
      else {
        const j = "bond" + b.id;
        members.set(j, []);
        addMember(j, tok(A) + "-"); addMember(j, tok(B) + "+");
        lines.push(`# bond${b.id}: direct bond ${A.name} -> ${B.name}`);
      }
    }
    const linked = new Set(connects.flatMap((c) => c.split(" ").slice(2)));
    for (const [name, toks] of members) {
      const n = graph.nodes.find((x) => x.name === name && JUNCTION[x.kind]);
      const kind = n ? n.kind : "1";
      if (toks.length || linked.has(name)) lines.push(`junction ${name} ${kind}: ${toks.join(", ")}`);
    }
    lines.push("");
    for (const n of graph.nodes) {
      if (JUNCTION[n.kind] || !incident(graph, n.id).length) continue;
      const key = lawKeyOf(n);
      let line = `element ${n.name} ${n.kind} ${key}=${String((n.law && n.law[key]) || "").replace(/\s+/g, "")}`;
      if (n.kind === "R" && n.law && n.law.causality) line += ` causality=${n.law.causality}`;
      lines.push(line);
    }
    extraEls.forEach((l) => lines.push(l));
    connects.forEach((l) => lines.push(l));
    // Editor-only data the CLI ignores: positions, parameter values, inputs, settings.
    lines.push("", "# --- editor data (comments; bondgraph.py ignores these) ---");
    for (const n of graph.nodes) {
      if (incident(graph, n.id).length) continue;
      const key = lawKeyOf(n);
      lines.push(`#@node ${n.name} ${n.kind}` + (key ? ` ${key}=${String((n.law && n.law[key]) || "").replace(/\s+/g, "")}` : ""));
    }
    for (const n of graph.nodes) lines.push(`#@pos ${n.name} ${Math.round(n.x)} ${Math.round(n.y)}`);
    if (settings) {
      for (const [k, v] of Object.entries(settings.params || {})) lines.push(`#@param ${k}=${String(v).replace(/\s+/g, "")}`);
      for (const [k, v] of Object.entries(settings.inputs || {})) lines.push(`#@input ${k}=${String(v).replace(/\s+/g, "")}`);
      for (const [k, v] of Object.entries(settings.init || {})) lines.push(`#@init ${k}=${String(v).replace(/\s+/g, "")}`);
      lines.push(`#@sim tEnd=${settings.tEnd} steps=${settings.steps}`);
    }
    return lines.join("\n") + "\n";
  }

  function fromBond(text) {
    const nodes = [], bonds = [];
    const pos = {}, settings = { params: {}, inputs: {}, init: {}, tEnd: 10, steps: 1000 };
    const junctionLines = [], connects = [], unbonded = [];
    let nextId = 1;
    const els = new Map(), juncs = new Map();
    text.split(/\r?\n/).forEach((raw, i) => {
      const meta = /^#@(\w+)\s+(.*)$/.exec(raw.trim());
      if (meta) {
        const [, tag, rest] = meta;
        const into = { param: "params", input: "inputs", init: "init" }[tag];
        if (tag === "pos") { const [nm, x, y] = rest.split(/\s+/); pos[nm] = [+x, +y]; }
        else if (tag === "node") unbonded.push(rest.trim().split(/\s+/));
        else if (into) {
          const eq = rest.indexOf("=");
          if (eq > 0) settings[into][rest.slice(0, eq).trim()] = rest.slice(eq + 1).trim();
        } else if (tag === "sim") {
          rest.split(/\s+/).forEach((kv) => { const [k, v] = kv.split("="); if (k === "tEnd" || k === "steps") settings[k] = Number(v); });
        }
        return;
      }
      const line = raw.split("#")[0].trim();
      if (!line) return;
      let m;
      if ((m = /^junction\s+(\w+)\s+([01])\s*:\s*(.*)$/i.exec(line))) {
        if (juncs.has(m[1])) throw new SyntaxError(`line ${i + 1}: duplicate junction ${m[1]}`);
        const n = { id: nextId++, kind: m[2], name: m[1], x: 0, y: 0, law: {} };
        juncs.set(m[1], n); junctionLines.push([n, m[3], i + 1]);
      } else if ((m = /^element\s+(\w+)\s+(I|C|R|Se|Sf|TF|GY)(?:\s+(.*))?$/i.exec(line))) {
        const kind = { i: "I", c: "C", r: "R", se: "Se", sf: "Sf", tf: "TF", gy: "GY" }[m[2].toLowerCase()];
        if (els.has(m[1])) throw new SyntaxError(`line ${i + 1}: duplicate element ${m[1]}`);
        const law = {};
        (m[3] || "").split(/\s+/).filter(Boolean).forEach((w) => {
          const eq = w.indexOf("=");
          if (eq < 0) throw new SyntaxError(`line ${i + 1}: expected key=value, got '${w}'`);
          law[w.slice(0, eq)] = w.slice(eq + 1);
        });
        const forms = LAW_FORMS[kind], allowed = kind === "R" ? forms.concat("causality") : forms;
        const unknown = Object.keys(law).find((k) => !allowed.includes(k));
        if (unknown) throw new SyntaxError(`line ${i + 1}: ${m[2]} takes ${allowed.map((k) => k + "=").join(", ")} not ${unknown}=`);
        const given = forms.filter((k) => k in law);
        if (given.length > 1) throw new SyntaxError(`line ${i + 1}: ${m[1]} has two laws (${given[0]}= and ${given[1]}=); give one`);
        if (!given.some((k) => law[k])) {
          given.forEach((k) => delete law[k]);
          law[LAW_KEY[kind]] = (kind === "Se" || kind === "Sf" ? "u_" : kind === "TF" || kind === "GY" ? "n_" : kind + "_") + m[1];
        }
        els.set(m[1], { id: nextId++, kind, name: m[1], x: 0, y: 0, law });
      } else if ((m = /^connect\s+(\w+)\s+(\w+)\s+(\w+)$/i.exec(line))) {
        connects.push([m[1], m[2], m[3], i + 1]);
      } else {
        throw new SyntaxError(`line ${i + 1}: expected junction, element, or connect`);
      }
    });
    // Unit transformers named bondN with connect are plain junction-junction bonds.
    const plain = new Set();
    for (const [name, a, b] of connects) {
      const el = els.get(name);
      if (el && el.kind === "TF" && /^bond\d+$/.test(name) && String(el.law.n) === "1") plain.add(name);
    }
    for (const [n, list, ln] of junctionLines) {
      list.split(",").map((s) => s.trim()).filter(Boolean).forEach((tok) => {
        const mm = /^(\w+)(?:\.(1|2))?\s*([+-])?$/.exec(tok);
        if (!mm) throw new SyntaxError(`line ${ln}: bad attachment '${tok}'`);
        const el = els.get(mm[1]);
        if (!el) throw new SyntaxError(`junction ${n.name}: unknown element ${mm[1]}`);
        const port = TWO_PORT[el.kind] ? Number(mm[2]) : undefined;
        if (TWO_PORT[el.kind] && !port) throw new SyntaxError(`junction ${n.name}: two-port ${el.name} needs .1 or .2`);
        const intoEl = mm[3] !== "-";
        bonds.push(intoEl ? { id: 0, from: n.id, to: el.id, toPort: port } : { id: 0, from: el.id, to: n.id, fromPort: port });
      });
    }
    for (const [name, a, b, ln] of connects) {
      const el = els.get(name), A = juncs.get(a), B = juncs.get(b);
      if (!el || !TWO_PORT[el.kind]) throw new SyntaxError(`line ${ln}: connect ${name} must name a TF or GY`);
      if (!A || !B) throw new SyntaxError(`line ${ln}: connect ${name}: unknown junction`);
      if (plain.has(name)) { bonds.push({ id: Number(name.slice(4)) || 0, from: A.id, to: B.id }); continue; }
      bonds.push({ id: 0, from: A.id, to: el.id, toPort: 1 }, { id: 0, from: el.id, to: B.id, fromPort: 2 });
    }
    for (const name of plain) els.delete(name);
    // A two-bond 1-junction named bondN that only links two elements is a direct bond.
    for (const [name, j] of [...juncs]) {
      if (!/^bond\d+$/.test(name)) continue;
      const inc = bonds.filter((b) => b.from === j.id || b.to === j.id);
      if (inc.length !== 2) continue;
      const [b1, b2] = inc;
      const inBond = b1.to === j.id ? b1 : b2.to === j.id ? b2 : null;
      const outBond = inBond === b1 ? b2 : b1;
      if (!inBond || outBond.from !== j.id) continue;
      bonds.splice(bonds.indexOf(b1), 1); bonds.splice(bonds.indexOf(b2), 1);
      bonds.push({ id: Number(name.slice(4)) || 0, from: inBond.from, fromPort: inBond.fromPort, to: outBond.to, toPort: outBond.toPort });
      juncs.delete(name);
    }
    // Nodes the editor saved without bonds (the CLI would reject them).
    for (const [name, kind, ...kv] of unbonded) {
      if (!name || !kind || juncs.has(name) || els.has(name)) continue;
      const law = {};
      kv.forEach((w) => { const eq = w.indexOf("="); if (eq > 0) law[w.slice(0, eq)] = w.slice(eq + 1); });
      (JUNCTION[kind] ? juncs : els).set(name, { id: nextId++, kind, name, x: 0, y: 0, law });
    }
    nodes.push(...juncs.values(), ...els.values());
    // Give every bond a unique id, keeping bondN ids where possible.
    const used = new Set();
    for (const b of bonds) if (b.id && !used.has(b.id)) used.add(b.id); else b.id = 0;
    let nid = 1;
    for (const b of bonds) if (!b.id) { while (used.has(nid)) nid++; b.id = nid; used.add(nid); }
    const missing = nodes.filter((n) => !pos[n.name]);
    nodes.forEach((n) => { if (pos[n.name]) { n.x = pos[n.name][0]; n.y = pos[n.name][1]; } });
    if (missing.length) autoLayout({ nodes, bonds }, new Set(missing.map((n) => n.id)));
    return { graph: { nodes, bonds }, settings };
  }

  // Force-directed placement for nodes without stored positions.
  function autoLayout(graph, movable) {
    const ns = graph.nodes;
    const free = (n) => !movable || movable.has(n.id);
    ns.forEach((n, i) => {
      if (!free(n)) return;
      const a = (2 * Math.PI * i) / ns.length;
      n.x = 420 + 220 * Math.cos(a); n.y = 300 + 180 * Math.sin(a);
    });
    const byId = new Map(ns.map((n) => [n.id, n]));
    for (let it = 0; it < 400; it++) {
      const temp = 30 * (1 - it / 400) + 1;
      const disp = new Map(ns.map((n) => [n.id, [0, 0]]));
      for (let i = 0; i < ns.length; i++) for (let j = i + 1; j < ns.length; j++) {
        const a = ns[i], b = ns[j];
        let dx = a.x - b.x, dy = a.y - b.y; const d = Math.hypot(dx, dy) || 0.01;
        const f = 9000 / (d * d);
        dx /= d; dy /= d;
        disp.get(a.id)[0] += dx * f; disp.get(a.id)[1] += dy * f;
        disp.get(b.id)[0] -= dx * f; disp.get(b.id)[1] -= dy * f;
      }
      for (const b of graph.bonds) {
        const A = byId.get(b.from), B = byId.get(b.to);
        let dx = A.x - B.x, dy = A.y - B.y; const d = Math.hypot(dx, dy) || 0.01;
        const f = (d - 110) * 0.12;
        dx /= d; dy /= d;
        disp.get(A.id)[0] -= dx * f; disp.get(A.id)[1] -= dy * f;
        disp.get(B.id)[0] += dx * f; disp.get(B.id)[1] += dy * f;
      }
      for (const n of ns) {
        if (!free(n)) continue;
        const [dx, dy] = disp.get(n.id); const d = Math.hypot(dx, dy) || 1;
        n.x += (dx / d) * Math.min(d, temp); n.y += (dy / d) * Math.min(d, temp);
      }
    }
    // Shift into view.
    const minX = Math.min(...ns.map((n) => n.x)), minY = Math.min(...ns.map((n) => n.y));
    ns.forEach((n) => { n.x += 90 - minX; n.y += 80 - minY; });
  }

  const api = { parse, names, compile, evalConst, checkReserved, FUNCS, RESERVED, assignCausality, buildModel,
    simulate, validate, toBond, fromBond, autoLayout, bondAliases, portOf, stateName, incident, maxBonds,
    lawKeyOf, STORAGE, TWO_PORT, JUNCTION, LAW_KEY, LAW_FORMS, OWN_VARS, fmt };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.BondEngine = api;
})(typeof self !== "undefined" ? self : this);
