# Bond graph state-equation tool

`bondgraph.py` reads a small `.bond` file, assigns causalities, and produces
symbolic state equations.  It is intentionally a command-line program so a
model can live alongside a homework/project and be version controlled.

## Install and run

```sh
python3 -m pip install -r requirements.txt
python3 bondgraph.py examples/rc.bond
```

Use `--json` for machine-readable output.  `--show-equations` also prints the
complete constitutive/junction system used to obtain the result.  `--latex`
prints the state equations as LaTeX, one per state; `--latex=matrix` prints
them as `ẋ = A x + B u`, with each source's value as an entry of `u`.

## Editor

Open [`ui/index.html`](ui/index.html) in a browser (no install or server
needed) to draw a bond graph instead of typing it:

- **Draw:** pick an element or junction from the toolbar (or press its key:
  `E` Se, `F` Sf, `C`, `I`, `R`, `T` TF, `G` GY, `0`, `1`) and click to place it.
  Switch to **Bond** (`B`), or hold Shift, and drag from one node to another.
  Power flows the way you drag.
- **Causality:** strokes update as you draw.  Storage elements in integral
  causality are green with an `∫` badge, derivative causality is amber
  (`d/dt`), and causal conflicts are red.
- **Laws:** click any node to set its constitutive law.  Laws use the same
  SymPy-style expressions as `.bond` files, including `Piecewise` and state
  references like `q_C1`.  For `R`, `C` and `I`, **Law form** switches between
  the linear coefficient and the law as a function either way round (`e = φ(f)`
  or `f = φ(e)` for `R`).  Common laws (diode, orifice, hardening spring...)
  are one click away.
- **Simulation:** the right pane plots every state over time and re-runs on each
  change.  It also lists parameters, inputs and initial conditions to edit.
  It uses an implicit BDF2 solver, so stiff models and derivative causality
  work too.
- **Equations:** the **Equations** button opens a window with the solved
  state equations, typeset, either one per state or in matrix form
  (`ẋ = A x + B u`; a nonlinear model is shown as one vector).  **Copy LaTeX**
  copies what it shows.  It runs this CLI's solver with SymPy inside the
  browser (via [Pyodide](https://pyodide.org)), so it shows exactly what
  `bondgraph.py --latex` prints, using the causality drawn in the editor.  The
  first time it opens it downloads about 11 MB, so that needs an internet
  connection; after that the browser's cache serves it.  Drag the window by
  its title and resize it from its corner.
- **Files:** **Export .bond** produces a file this CLI reads; **Import .bond**
  opens one (layout is kept in `#@` comment lines the CLI ignores).

The engine's tests run with `node --test ui/engine.test.js`.  The equations
window runs an embedded copy of `bondgraph.py` (`ui/bondgraph-py.js`, since a
page opened from disk can't read the `.py` file).  After changing
`bondgraph.py`, run `node ui/embed-python.js` to refresh the copy; the engine
tests fail until you do.

## `.bond` format

The full reference covers syntax, every element kind, gyrators, modulated
elements, piecewise laws, causality, output and errors:
[`docs/bond-format.md`](docs/bond-format.md).  A summary follows.

Blank lines and text after `#` are ignored.  A junction lists attached
one-port elements; each attachment creates one bond automatically.  `+`
(the default) points the bond's half-arrow from the junction into the element;
`-` points it from the element into the junction, as for sources.  Signs only
set reference directions: element laws use the flow into the element, so
flipping a sign never changes the physics.

```text
junction electrical 0: Is-, C1+, R1+

element Is Sf value=u(t)
element C1 C C=0.01
element R1 R R=100
```

Supported one-port elements are `I`, `C`, `R`, `Se`, and `Sf`.  A two-port
`TF` or `GY` is connected with a `connect` line; this also creates its two
bonds.  Its ports are named `name.1` and `name.2` in diagnostics.

```text
junction left 1: I1+, T.1+
junction right 0: T.2-, R1+
element I1 I I=2
element R1 R R=4
element T TF n=3
```

Alternatively, omit `T.1`/`T.2` from the junction lists and write `connect T
left right`; the two bonds are then attached automatically, with power flowing
`left -> T -> right`.

Parameters and source values are SymPy expressions.  For linear `R`, use
`R=...`; for `I` and `C`, use `I=...` and `C=...`.  To give the law itself,
name the variable it sets: `R` takes `e=φ(f)` or `f=φ(e)`, `C` takes `e=φ(q)`
or `q=φ(e)`, and `I` takes `f=φ(p)` or `p=φ(f)`.  Inside these laws `e`, `f`,
`q` and `p` are the element's own effort, inflow and state, e.g.
`element D R f=Piecewise((e/R_on,e>0),(0,True))`.  Arbitrary symbolic input
names such as `u(t)` are accepted, and SymPy functions such as `Piecewise`,
`Abs`, `Max`, `sin` and `sqrt` keep their usual meaning.  Parameters may
reference states (`q_C1`, `p_I1`) for modulated or piecewise laws, e.g.
`R=Piecewise((R_low,q_Tank>V_max),(R_high,True))` (no spaces).
For a non-invertible resistance, add `causality=effort` (it imposes effort) or
`causality=flow` (it imposes flow).  The latter is appropriate for the diode
example described in the notes, written with `f=`.

### Conventions and scope

At a 0-junction efforts are equal and signed flows sum to zero.  At a
1-junction flows are equal and signed efforts sum to zero.  Bond variables are
referenced as `e_Element_port` and `f_Element_port` in the full equations.
With power entering port 1 and leaving port 2, the transformer convention is
`e1=n*e2`, `f2=n*f1` and the gyrator convention is `e1=n*f2`, `e2=n*f1`.
Both conserve power whatever signs the two bonds carry.

Integral-causality `I` elements use momentum `p_Name` as their state;
integral-causality `C` elements use displacement/charge `q_Name`.  If a
storage element must take derivative causality it is reported and is not
included in the state vector.  The solver is intended for ordinary, solvable
symbolic models; unsupported algebraic loops are reported with the unresolved
equations rather than silently guessed.
