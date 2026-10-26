// Run with: node --test ui/engine.test.js
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const E = require("./engine.js");

const near = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ""} expected ${b}, got ${a}`);

// Tiny graph builder: nodes as [name, kind, law], bonds as "A>B" (power A -> B).
function G(nodes, bonds) {
  const ns = nodes.map(([name, kind, law], i) => ({ id: i + 1, name, kind, x: 0, y: 0, law: law || {} }));
  const id = (nm) => ns.find((n) => n.name === nm.split(".")[0]).id;
  const port = (nm) => (nm.includes(".") ? Number(nm.split(".")[1]) : undefined);
  const bs = bonds.map((s, i) => {
    const [a, b] = s.split(">");
    return { id: i + 1, from: id(a), to: id(b), fromPort: port(a), toPort: port(b) };
  });
  return { nodes: ns, bonds: bs };
}
function run(graph, settings) {
  const c = E.assignCausality(graph);
  const r = E.simulate(graph, Object.assign({ tEnd: 5, steps: 1000, params: {}, inputs: {}, init: {} }, settings), c);
  assert.ok(r.ok, JSON.stringify(r.errors));
  const last = r.X[r.X.length - 1];
  const val = (name) => last[r.model.sym.get(name)];
  const at = (name, t) => r.X[Math.round((t / settings.tEnd) * (r.X.length - 1))][r.model.sym.get(name)];
  return { c, r, val, at };
}

test("expressions follow SymPy syntax", () => {
  const ev = (s, p) => E.evalConst(s, p);
  assert.equal(ev("2+3*4"), 14);
  assert.equal(ev("-2**2"), -4);
  assert.equal(ev("2^3^2"), 512);
  assert.equal(ev("Piecewise((1,x>0),(2,True))", { x: 3 }), 1);
  assert.equal(ev("Piecewise((1,x>0),(2,True))", { x: -3 }), 2);
  assert.equal(ev("Piecewise((1,(x>0)&(x<1)),(0,True))", { x: 0.5 }), 1);
  assert.equal(ev("Heaviside(0)"), 0.5);
  assert.equal(ev("Max(-1,Min(1,5))"), 1);
  near(ev("sqrt(Abs(-4))+sin(pi/2)"), 3, 1e-12);
  assert.throws(() => ev("1+"), /ends too early/);
  assert.throws(() => E.checkReserved(E.parse("I*2")), /reserved/);
  const nm = E.names(E.parse("R0*(1+a*q_C1)+u(t)"));
  assert.deepEqual([...nm.ids].sort(), ["R0", "a", "q_C1", "t"]);
  assert.deepEqual([...nm.inputs], ["u"]);
});

test("mass-spring-damper is all integral and matches the analytic response", () => {
  // Undamped: x'' + x = 0, x(0) = 1  =>  q = cos t
  const g = G([["F", "Se", { value: "0" }], ["v", "1"], ["m", "I", { I: "1" }], ["k", "C", { C: "1" }]],
    ["F>v", "v>m", "v>k"]);
  const { c, at } = run(g, { tEnd: 10, steps: 4000, init: { q_k: "1" } });
  assert.deepEqual(Object.values(c.storage), ["integral", "integral"]);
  assert.equal(c.conflicts.length, 0);
  near(at("q_k", 10), Math.cos(10), 2e-3, "q_k(10)");
  near(at("p_m", 10), -Math.sin(10), 2e-3, "p_m(10)");
});

test("RC circuit matches q = RC(1 - exp(-t/RC)) and ignores arrow direction", () => {
  const exact = 0.5 * 2 * (1 - Math.exp(-5 / 1)); // R = 0.5, C = 2, u = 2  => q = u*R*C*(1-e^{-t/RC})
  for (const arrows of [["Is>n", "n>C1", "n>R1"], ["n>Is", "C1>n", "R1>n"]]) {
    const g = G([["Is", "Sf", { value: "u(t)" }], ["n", "0"], ["C1", "C", { C: "C" }], ["R1", "R", { R: "R" }]], arrows);
    const { val } = run(g, { tEnd: 5, steps: 2000, params: { R: "0.5", C: "2" }, inputs: { u: "2" } });
    near(val("q_C1"), 2 * exact, 2e-4, `q_C1 with arrows ${arrows}`);
  }
});

test("R written as e = φ(f) or f = φ(e) matches the linear law whichever way the arrows point", () => {
  const exact = 2 * 0.5 * 2 * (1 - Math.exp(-5 / 1));
  for (const law of [{ e: "R*f" }, { f: "e/R" }]) {
    for (const arrows of [["Is>n", "n>C1", "n>R1"], ["n>Is", "C1>n", "R1>n"]]) {
      const g = G([["Is", "Sf", { value: "u(t)" }], ["n", "0"], ["C1", "C", { C: "C" }], ["R1", "R", law]], arrows);
      const { r, val } = run(g, { tEnd: 5, steps: 2000, params: { R: "0.5", C: "2" }, inputs: { u: "2" } });
      assert.deepEqual(r.model.params, ["C", "R"], "e and f are the resistor's own variables, not parameters");
      near(val("q_C1"), exact, 2e-4, `q_C1 with ${JSON.stringify(law)} and arrows ${arrows}`);
    }
  }
});

test("C and I laws can be written as functions in either direction", () => {
  for (const [k, m] of [[{ e: "q" }, { f: "p" }], [{ q: "e" }, { p: "f" }]]) {
    const g = G([["F", "Se", { value: "0" }], ["v", "1"], ["m", "I", m], ["k", "C", k]], ["F>v", "v>m", "v>k"]);
    const { c, at } = run(g, { tEnd: 10, steps: 4000, init: { q_k: "1" } });
    assert.deepEqual(Object.values(c.storage), ["integral", "integral"]);
    near(at("q_k", 10), Math.cos(10), 2e-3, `q_k(10) with C ${JSON.stringify(k)}, I ${JSON.stringify(m)}`);
  }
});

test("nonlinear laws are inverted when the causality needs it", () => {
  // Cubic spring e = q^3 behind a damper: settles where q^3 = F.
  const spring = G([["F", "Se", { value: "8" }], ["v", "1"], ["k", "C", { e: "q**3" }], ["b", "R", { R: "1" }]], ["F>v", "v>k", "v>b"]);
  near(run(spring, { tEnd: 20, steps: 2000 }).val("q_k"), 2, 1e-4, "cubic spring deflection");
  // Saturating inductor p = tanh(f) in integral causality: the solver inverts it for f.
  const coil = G([["V", "Se", { value: "1" }], ["a", "1"], ["L", "I", { p: "tanh(f)" }], ["R1", "R", { R: "1" }]], ["V>a", "a>L", "a>R1"]);
  const { val } = run(coil, { tEnd: 20, steps: 2000 });
  near(val("f_L"), 1, 1e-4, "steady current");
  near(val("p_L"), Math.tanh(1), 1e-4, "flux linkage");
  // Diode written f = φ(e) in flow causality, the way it's naturally given.
  const diode = G([["V", "Se", { value: "sin(t)" }], ["l", "1"], ["D", "R", { f: "Piecewise((e/0.01,e>0),(e/100,True))" }], ["C1", "C", { C: "1" }]],
    ["V>l", "l>D", "l>C1"]);
  const d = run(diode, { tEnd: 6, steps: 1200 });
  assert.ok(!d.r.partial, JSON.stringify(d.r.errors));
  near(d.at("q_C1", Math.PI / 2), 1, 0.05, "capacitor charges to the peak through the diode");
});

test("function-form laws survive .bond export and import", () => {
  const g = G([["Is", "Sf", { value: "1" }], ["n", "0"], ["C1", "C", { q: "C*e" }], ["R1", "R", { f: "e/R", causality: "flow" }],
    ["m", "I", { p: "m*f" }], ["lonely", "R", { e: "b*f" }]], ["Is>n", "n>C1", "n>R1", "n>m"]);
  const text = E.toBond(g);
  assert.match(text, /^element C1 C q=C\*e$/m);
  assert.match(text, /^element R1 R f=e\/R causality=flow$/m);
  assert.match(text, /^#@node lonely R e=b\*f$/m);
  const back = E.fromBond(text).graph.nodes;
  const law = (nm) => back.find((n) => n.name === nm).law;
  assert.deepEqual(law("C1"), { q: "C*e" });
  assert.deepEqual(law("R1"), { f: "e/R", causality: "flow" });
  assert.deepEqual(law("m"), { p: "m*f" });
  assert.deepEqual(law("lonely"), { e: "b*f" });
  assert.throws(() => E.fromBond("junction n 0: R1+\nelement R1 R R=1 f=e\n"), /two laws \(R= and f=\)/);
  assert.throws(() => E.fromBond("junction n 0: C1+\nelement C1 C f=e\n"), /C takes C=, e=, q= not f=/);
});

test("a junction joined only to other junctions survives export", () => {
  const g = G([["F", "Se", { value: "1" }], ["v", "1"], ["m", "I", { I: "1" }], ["J1", "1"], ["J2", "1"], ["I1", "I", { I: "1" }]],
    ["F>v", "v>m", "v>J1", "J1>J2", "J2>I1"]);
  const text = E.toBond(g);
  assert.match(text, /^junction J1 1: *$/m);
  const back = E.fromBond(text).graph;
  assert.deepEqual(back.nodes.map((n) => n.name).sort(), ["F", "I1", "J1", "J2", "m", "v"]);
  assert.equal(back.bonds.length, 5);
});

test("DC motor (gyrator) reaches K v / (R b + K^2)", () => {
  const g = G([["V", "Se", { value: "1" }], ["arm", "1"], ["Ra", "R", { R: "1" }], ["La", "I", { I: "0.5" }],
    ["G", "GY", { n: "0.5" }], ["shaft", "1"], ["J", "I", { I: "0.1" }], ["b", "R", { R: "0.05" }]],
  ["V>arm", "arm>Ra", "arm>La", "arm>G.1", "G.2>shaft", "shaft>J", "shaft>b"]);
  const { c, val } = run(g, { tEnd: 20, steps: 4000 });
  assert.equal(c.conflicts.length, 0);
  near(val("p_J") / 0.1, 0.5 / (1 * 0.05 + 0.25), 1e-4, "omega");
  near(val("p_La") / 0.5, (1 - 0.5 * (0.5 / 0.3)) / 1, 1e-4, "current"); // v = R i + K w
});

test("transformer passes effort and flow by n", () => {
  // Lever: Se force F -> TF n -> spring k. Static: e1 = n e2 => spring force = F/n; q = F/(n k)
  const g = G([["F", "Se", { value: "3" }], ["a", "1"], ["d", "R", { R: "1" }], ["T", "TF", { n: "2" }], ["k", "C", { C: "0.25" }]],
    ["F>a", "a>d", "a>T.1", "T.2>k"]);
  const { val } = run(g, { tEnd: 30, steps: 3000 });
  near(val("q_k"), (3 / 2) * 0.25, 1e-5, "spring deflection");
});

test("parallel capacitors put one in derivative causality and still simulate", () => {
  const g = G([["S", "Sf", { value: "1" }], ["n", "0"], ["C1", "C", { C: "1" }], ["C2", "C", { C: "3" }]],
    ["S>n", "n>C1", "n>C2"]);
  const { c, val } = run(g, { tEnd: 2, steps: 400 });
  assert.deepEqual(Object.values(c.storage).sort(), ["derivative", "integral"]);
  near(val("q_C1"), 2 * 0.25, 1e-6, "q_C1"); // charge splits 1:3
  near(val("q_C2"), 2 * 0.75, 1e-6, "q_C2");
});

test("two effort sources on one 0-junction are a causal conflict", () => {
  const g = G([["A", "Se", { value: "1" }], ["B", "Se", { value: "2" }], ["n", "0"], ["R1", "R", { R: "1" }]],
    ["A>n", "B>n", "n>R1"]);
  const c = E.assignCausality(g);
  assert.ok(c.conflicts.length >= 1);
});

test("a fixed-flow resistor pushes a storage element into derivative causality", () => {
  const g = G([["D", "R", { R: "1", causality: "flow" }], ["n", "1"], ["m", "I", { I: "1" }], ["F", "Se", { value: "1" }]],
    ["F>n", "n>D", "n>m"]);
  const c = E.assignCausality(g);
  assert.equal(c.conflicts.length, 0);
  assert.equal(c.storage[3], "derivative");
});

test("piecewise law on a state switches the dynamics", () => {
  // Tank with outlet resistance that drops above V_max.
  const g = G([["Q", "Sf", { value: "1" }], ["t", "0"], ["Tank", "C", { C: "1" }],
    ["Out", "R", { R: "Piecewise((0.5,q_Tank>1),(100,True))" }]], ["Q>t", "t>Tank", "t>Out"]);
  const { val } = run(g, { tEnd: 20, steps: 4000 });
  near(val("q_Tank"), 1, 0.02, "level settles at the switch point region");
});

test("a state-switched outlet passes through its switch level", () => {
  // Same numbers as the editor's overflow example: settles at C*R_low*Q above V_max.
  const g = G([["Q", "Sf", { value: "0.01" }], ["t", "0"], ["Tank", "C", { C: "1/9810" }],
    ["Out", "R", { R: "Piecewise((6e5,q_Tank>0.5),(2e6,True))" }]], ["Q>t", "t>Tank", "t>Out"]);
  const { r, val } = run(g, { tEnd: 600, steps: 3000 });
  assert.ok(!r.partial, JSON.stringify(r.errors));
  near(val("q_Tank"), 6e5 * 0.01 / 9810, 2e-3, "overflow level");
});

test("a switch that chatters (sliding mode) still completes", () => {
  // Below 0.5 the outlet can't keep up; above it, it drains faster than the inflow.
  const g = G([["Q", "Sf", { value: "0.01" }], ["t", "0"], ["Tank", "C", { C: "1/9810" }],
    ["Out", "R", { R: "Piecewise((2e4,q_Tank>0.5),(2e6,True))" }]], ["Q>t", "t>Tank", "t>Out"]);
  const { r, val } = run(g, { tEnd: 600, steps: 3000 });
  assert.ok(!r.partial, JSON.stringify(r.errors));
  // It chatters in a band about one step of drain wide: h * |dq/dt above| = 0.2 * 0.235.
  const q = val("q_Tank");
  assert.ok(q > 0.5 - 0.05 && q < 0.5 + 0.01, `held at the switch level, got ${q}`);
});

test("nonlinear resistor in flow causality (diode-like) converges", () => {
  const g = G([["V", "Se", { value: "sin(t)" }], ["l", "1"], ["D", "R", { R: "Piecewise((0.01,f_D>0),(100,True))" }], ["C1", "C", { C: "1" }]],
    ["V>l", "l>D", "l>C1"]);
  const { r } = run(g, { tEnd: 6, steps: 1200 });
  assert.ok(!r.partial, JSON.stringify(r.errors));
});

test(".bond export runs through the CLI format and imports back", () => {
  const g = G([["V", "Se", { value: "v(t)" }], ["arm", "1"], ["Ra", "R", { R: "R_a" }], ["La", "I", { I: "L_a" }],
    ["G", "GY", { n: "K_t" }], ["shaft", "1"], ["J", "I", { I: "J" }], ["b", "R", { R: "b" }], ["j2", "0"], ["Cx", "C", { C: "1" }],
    ["lonely", "R", { R: "5" }]],
  ["V>arm", "arm>Ra", "arm>La", "arm>G.1", "G.2>shaft", "shaft>J", "shaft>j2", "j2>b", "j2>Cx"]);
  const text = E.toBond(g, { params: { R_a: "1" }, inputs: { v: "1" }, init: {}, tEnd: 10, steps: 500 });
  assert.match(text, /element bond7 TF n=1/);
  assert.match(text, /connect bond7 shaft j2/);
  assert.match(text, /#@node lonely R R=5/);
  const back = E.fromBond(text);
  assert.equal(back.graph.nodes.length, g.nodes.length);
  assert.equal(back.graph.bonds.length, g.bonds.length);
  assert.equal(back.settings.params.R_a, "1");
  assert.equal(back.settings.inputs.v, "1");
  fs.writeFileSync(path.join(require("node:os").tmpdir(), "ui_export.bond"), text);
});

test("outputs follow the response, and a broken one doesn't stop it", () => {
  // Series RC from a 2 V source, RC = 1: v_C = 2(1 - e^-t), i = (2 - v_C)/R
  const g = G([["V", "Se", { value: "2" }], ["a", "1"], ["R1", "R", { R: "R" }], ["C1", "C", { C: "C" }]], ["V>a", "a>R1", "a>C1"]);
  g.outputs = [{ name: "i", expr: "f_R1" }, { name: "v_C", expr: "q_C1/C" }, { name: "P", expr: "k*e_R1*f_R1" }, { name: "bad", expr: "f_R1+" }];
  const { r } = run(g, { tEnd: 5, steps: 2000, params: { R: "0.5", C: "2", k: "3" } });
  assert.deepEqual(r.model.params, ["C", "R", "k"], "a parameter only an output uses is listed");
  assert.match(r.model.outputs[3].error, /ends too early/);
  const end = (k) => r.Y[k][r.Y[k].length - 1];
  near(end(1), 2 * (1 - Math.exp(-5)), 2e-4, "v_C(5)");
  near(end(0), 4 * Math.exp(-5), 2e-4, "i(5)");
  near(end(2), 3 * 0.5 * end(0) ** 2, 1e-9, "P = k R i^2");
  assert.ok(r.Y[3].every(Number.isNaN));
});

test("an integral output is a state with its own initial value", () => {
  // Series RL from 2 V, L/R = 1: i = 2/R (1 - e^-t), and the charge through it is its integral.
  const g = G([["V", "Se", { value: "2" }], ["a", "1"], ["R1", "R", { R: "1" }], ["L1", "I", { I: "1" }]], ["V>a", "a>R1", "a>L1"]);
  g.outputs = [{ name: "pw", expr: "k*e_R1*f_R1" }, { name: "q", expr: "integral(f_L1)" }, { name: "wasted", expr: "integral( pw )" },
    { name: "loopA", expr: "loopB" }, { name: "loopB", expr: "loopA" }, { name: "nested", expr: "2*integral(f_L1)" }];
  const { r, val } = run(g, { tEnd: 5, steps: 4000, params: { k: "3" }, init: { q: "10" } });
  assert.deepEqual(r.model.states.map((s) => s.name), ["p_L1", "q", "wasted"]);
  const y = (name) => r.Y[r.model.outputs.findIndex((o) => o.name === name)];
  const end = (name) => y(name)[y(name).length - 1];
  near(end("q"), 10 + 2 * (5 - (1 - Math.exp(-5))), 2e-3, "q(5) = 10 + integral of i");
  near(y("q")[0], 10, 1e-9, "q starts at its initial value");
  near(end("wasted"), 3 * 4 * (5 - 2 * (1 - Math.exp(-5)) + (1 - Math.exp(-10)) / 2), 5e-2, "an integral of another output");
  assert.match(r.model.outputs[3].error, /loop/);
  assert.match(r.model.outputs[5].error, /whole output/);
  assert.ok(r.Y[5].every(Number.isNaN));
  assert.deepEqual(r.model.params, ["k"]);
});

test("outputs survive .bond export and import", () => {
  const g = G([["Is", "Sf", { value: "1" }], ["n", "0"], ["C1", "C", { C: "C" }]], ["Is>n", "n>C1"]);
  g.outputs = [{ name: "v", expr: " q_C1 / C " }, { name: "unfinished", expr: "" }];
  const text = E.toBond(g);
  assert.match(text, /^output v = q_C1 \/ C$/m);
  assert.doesNotMatch(text, /unfinished/);
  assert.deepEqual(E.fromBond(text).graph.outputs, [{ name: "v", expr: "q_C1 / C" }]);
  assert.throws(() => E.fromBond("output y = 1\noutput y = 2\n"), /line 2: duplicate output y/);
});

// Two runs agree on the named variables at a few times.
function same(a, b, names, tol) {
  for (const k of [0.25, 0.5, 1]) {
    const i = Math.round(k * (a.r.X.length - 1));
    for (const nm of names) near(a.r.X[i][a.r.model.sym.get(nm[0])], b.r.X[i][b.r.model.sym.get(nm[1] || nm[0])], tol, `${nm.join(" vs ")} at ${k}`);
  }
}

test("a C-field couples its ports like separate springs would", () => {
  const params = { m1: "1", m2: "2", k1: "3", k2: "1", k12: "2" }, init = { p_m1: "1" };
  const springs = G([["v1", "1"], ["v2", "1"], ["m1", "I", { I: "m1" }], ["m2", "I", { I: "m2" }], ["k1", "C", { C: "1/k1" }],
    ["k2", "C", { C: "1/k2" }], ["J0", "0"], ["k12", "C", { C: "1/k12" }]],
  ["v1>m1", "v1>k1", "v2>m2", "v2>k2", "v1>J0", "J0>v2", "J0>k12"]);
  const field = G([["v1", "1"], ["v2", "1"], ["m1", "I", { I: "m1" }], ["m2", "I", { I: "m2" }],
    ["S", "Cfield", { e1: "k1*q1+k12*(q1-q2)", e2: "k2*q2+k12*(q2-q1)" }]], ["v1>m1", "v2>m2", "v1>S.1", "v2>S.2"]);
  const a = run(springs, { tEnd: 6, steps: 1200, params, init }), b = run(field, { tEnd: 6, steps: 1200, params, init });
  assert.deepEqual(b.r.model.states.map((s) => s.name), ["p_m1", "p_m2", "q_S", "q_S_2"]);
  assert.deepEqual(b.c.ports[5], { 1: "integral", 2: "integral" });
  same(a, b, [["p_m1"], ["p_m2"], ["q_k1", "q_S"], ["q_k2", "q_S_2"]], 1e-9);
});

test("an I-field given as p = M f moves like the rigid bar it describes", () => {
  // Bounce and pitch of a bar on two springs: v_front = vc + a w, v_rear = vc - b w.
  const params = { m: "2", J: "1", a: "0.6", b: "0.9", k_f: "30", k_r: "40", c_f: "1.5", c_r: "2", F0: "5" };
  const ends = [["F", "Se", { value: "F0*Heaviside(t-0.2)" }], ["vf", "1"], ["vr", "1"], ["kf", "C", { C: "1/k_f" }],
    ["cf", "R", { R: "c_f" }], ["kr", "C", { C: "1/k_r" }], ["cr", "R", { R: "c_r" }]];
  const endBonds = ["F>vf", "vf>kf", "vf>cf", "vr>kr", "vr>cr"];
  const bar = G(ends.concat([["vc", "1"], ["w", "1"], ["m", "I", { I: "m" }], ["Jb", "I", { I: "J" }], ["nf", "0"], ["nr", "0"],
    ["TA", "TF", { n: "a" }], ["TB", "TF", { n: "1/b" }]]),
  endBonds.concat(["vc>m", "w>Jb", "vc>nf", "w>TA.1", "TA.2>nf", "nf>vf", "vc>nr", "nr>vr", "nr>TB.1", "TB.2>w"]));
  const field = G(ends.concat([["Body", "Ifield", { p1: "((m*b^2+J)*f1+(m*a*b-J)*f2)/(a+b)^2", p2: "((m*a*b-J)*f1+(m*a^2+J)*f2)/(a+b)^2" }]]),
    endBonds.concat(["vf>Body.1", "vr>Body.2"]));
  const a = run(bar, { tEnd: 3, steps: 3000, params }), b = run(field, { tEnd: 3, steps: 3000, params });
  assert.equal(Object.values(b.c.storage).filter((s) => s === "derivative").length, 0);
  same(a, b, [["q_kf"], ["q_kr"], ["f_kf"], ["e_cr"]], 1e-8);
});

test("an R-field couples its ports like a resistor network", () => {
  const params = { C1: "1", C2: "0.5", R1: "2", R2: "3", Rm: "0.5" }, init = { q_C1: "1" };
  const net = G([["n1", "0"], ["n2", "0"], ["C1", "C", { C: "C1" }], ["C2", "C", { C: "C2" }], ["R1", "R", { R: "R1" }],
    ["R2", "R", { R: "R2" }], ["J1", "1"], ["Rm", "R", { R: "Rm" }]], ["n1>C1", "n1>R1", "n1>J1", "J1>Rm", "J1>n2", "n2>C2", "n2>R2"]);
  const field = G([["n1", "0"], ["n2", "0"], ["C1", "C", { C: "C1" }], ["C2", "C", { C: "C2" }],
    ["G", "Rfield", { f1: "(1/R1+1/Rm)*e1-e2/Rm", f2: "-e1/Rm+(1/R2+1/Rm)*e2" }]], ["n1>C1", "n1>G.1", "n2>C2", "n2>G.2"]);
  same(run(net, { tEnd: 3, steps: 600, params, init }), run(field, { tEnd: 3, steps: 600, params, init }), [["q_C1"], ["q_C2"]], 1e-9);
});

test("a field port in derivative causality follows its source", () => {
  // A flow source drives port 1, so that port's momentum is m1 times the source's flow.
  const g = G([["S", "Sf", { value: "sin(t)" }], ["a", "1"], ["V", "Se", { value: "1" }], ["b", "1"],
    ["X", "Ifield", { p1: "m1*f1", f2: "p2/m2" }]], ["S>a", "a>X.1", "V>b", "b>X.2"]);
  const { c, r, val } = run(g, { tEnd: 2, steps: 2000, params: { m1: "3", m2: "2" } });
  assert.deepEqual(c.ports[5], { 1: "derivative", 2: "integral" });
  assert.equal(c.storage[5], "derivative");
  near(val("p_X"), 3 * Math.sin(2), 1e-6, "p1 = m1 f1");
  near(val("p_X_2"), 2, 1e-6, "p2 integrates the unit effort");
  assert.ok(!r.partial);
});

test("fields survive .bond export and import, and bad ones are reported", () => {
  const g = G([["v1", "1"], ["v2", "1"], ["m1", "I", { I: "m1" }], ["m2", "I", { I: "m2" }],
    ["S", "ICfield", { e1: "k*q1 + a*p2", f2: "p2/m+a*q1" }], ["lonely", "Rfield", { e1: "R*f1" }]], ["v1>m1", "v2>m2", "v1>S.1", "S.2>v2"]);
  const text = E.toBond(g);
  assert.match(text, /^element S ICfield e1=k\*q1\+a\*p2 f2=p2\/m\+a\*q1$/m);
  assert.match(text, /^junction v1 1: m1\+, S\.1\+$/m);
  assert.match(text, /^junction v2 1: m2\+, S\.2-$/m);
  assert.match(text, /^#@node lonely Rfield e1=R\*f1$/m);
  const back = E.fromBond(text).graph;
  const S = back.nodes.find((n) => n.name === "S");
  assert.deepEqual(S.law, { e1: "k*q1+a*p2", f2: "p2/m+a*q1" });
  assert.deepEqual(E.fieldPorts(back, S).map((fp) => [fp.port, fp.type]), [[1, "C"], [2, "I"]]);
  const bad = (law) => () => E.fromBond(`junction a 1: F-, S.1+, S.2+\nelement F Se value=1\nelement S Cfield ${law}\n`);
  assert.throws(bad("e1=q1 f2=q2"), /Cfield takes a law per port, like e1=, q1=/);
  assert.throws(bad("e1=q1"), /port 2 needs a law \(e2= or q2=\)/);
  assert.throws(bad("e1=q1 q1=e1 e2=q2"), /port 1 has two laws/);
  assert.throws(() => E.fromBond("junction a 1: F-, S.1+, S.3+\nelement F Se value=1\nelement S Rfield e1=f1 e2=f2 e3=f3\n"), /port 2 isn't attached/);
  // A law that reaches for another port's variable is caught before simulating.
  const stray = G([["v", "1"], ["m", "I", { I: "1" }], ["S", "Cfield", { e1: "q2" }]], ["v>m", "v>S.1"]);
  assert.match(E.buildModel(stray, {}).errors[0].msg, /S port 1: q2 isn't one of its variables/);
});

test("the equations window's copy of bondgraph.py is up to date", () => {
  const embed = require("./embed-python.js");
  assert.ok(fs.existsSync(embed.OUT) && fs.readFileSync(embed.OUT, "utf8") === embed.render(),
    "ui/bondgraph-py.js is stale. Run: node ui/embed-python.js");
});

test("the editor's copy of the examples is up to date", () => {
  const embed = require("./embed-examples.js");
  assert.ok(fs.existsSync(embed.OUT) && fs.readFileSync(embed.OUT, "utf8") === embed.render(),
    "ui/examples.js is stale. Run: node ui/embed-examples.js");
});

test("every advanced example imports and simulates to the end", () => {
  const dir = path.join(__dirname, "..", "examples");
  const titled = fs.readdirSync(dir).filter((f) => f.endsWith(".bond")).map((f) => [f, fs.readFileSync(path.join(dir, f), "utf8")])
    .filter(([, text]) => /^#@title /m.test(text));
  assert.ok(titled.length >= 6);
  for (const [file, text] of titled) {
    const { graph, settings } = E.fromBond(text);
    const c = E.assignCausality(graph);
    assert.equal(c.conflicts.length, 0, file);
    const r = E.simulate(graph, settings, c);
    assert.ok(r.ok && !r.partial, `${file}: ${JSON.stringify(r.errors)}`);
    assert.ok(r.model.outputs.every((o) => !o.error), file);
  }
});

test("examples/dc_motor.bond imports and simulates", () => {
  const text = fs.readFileSync(path.join(__dirname, "..", "examples", "dc_motor.bond"), "utf8");
  const { graph } = E.fromBond(text);
  const params = { R_a: "1", L_a: "0.5", K_t: "0.5", J: "0.1", b: "0.05" };
  const { val } = run(graph, { tEnd: 20, steps: 4000, params, inputs: { v: "1" } });
  near(val("p_J") / 0.1, 0.5 / 0.3, 1e-4, "omega from the CLI example");
});
