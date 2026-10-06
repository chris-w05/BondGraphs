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

test("examples/dc_motor.bond imports and simulates", () => {
  const text = fs.readFileSync(path.join(__dirname, "..", "examples", "dc_motor.bond"), "utf8");
  const { graph } = E.fromBond(text);
  const params = { R_a: "1", L_a: "0.5", K_t: "0.5", J: "0.1", b: "0.05" };
  const { val } = run(graph, { tEnd: 20, steps: 4000, params, inputs: { v: "1" } });
  near(val("p_J") / 0.1, 0.5 / 0.3, 1e-4, "omega from the CLI example");
});
