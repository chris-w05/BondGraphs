# Writing `.bond` files

A `.bond` file describes a bond graph as a list of **junctions** and
**elements**. `bondgraph.py` reads it, assigns causality (preferring integral
causality on storage elements), and derives symbolic state equations.

```sh
python3 bondgraph.py model.bond                  # causality, states, state equations
python3 bondgraph.py model.bond --show-equations # also print every junction/element equation
python3 bondgraph.py model.bond --json           # machine-readable output
python3 bondgraph.py model.bond --latex          # state equations as LaTeX, one per state
python3 bondgraph.py model.bond --latex=matrix   # ... or as x' = A x + B u
```

Contents:

1. [File layout](#1-file-layout)
2. [Junctions](#2-junctions)
3. [Elements](#3-elements)
4. [Two-ports: TF and GY](#4-two-ports-tf-and-gy)
5. [Expressions](#5-expressions)
6. [Modulated elements](#6-modulated-elements)
7. [Piecewise and nonlinear laws](#7-piecewise-and-nonlinear-laws)
8. [Causality](#8-causality)
9. [Reading the output](#9-reading-the-output)
10. [Errors and troubleshooting](#10-errors-and-troubleshooting)
11. [Quick reference](#11-quick-reference)

---

## 1. File layout

- One statement per line. There are three statements: `junction`, `element`
  and `connect`.
- `#` starts a comment that runs to the end of the line. Blank lines are
  ignored.
- Keywords and element kinds are case-insensitive (`GY`, `gy`). Names are
  case-sensitive.
- Names are letters, digits and `_` only (`R1`, `L_a`, `tank`).
- Statements can appear in any order. An element can be used in a junction
  line before its `element` line appears.

```text
# Current source driving an RC network.
junction node 0: Is-, C1+, R1+

element Is Sf value=u(t)
element C1 C  C=C
element R1 R  R=R
```

## 2. Junctions

```text
junction <name> <0|1>: <attachment>, <attachment>, ...
```

| Type | Common variable | Summed variable |
|------|-----------------|-----------------|
| `0`  | effort (all efforts equal) | signed flows sum to zero |
| `1`  | flow (all flows equal)     | signed efforts sum to zero |

Each attachment creates **one bond** between the junction and an element
port:

- **One-port element:** write its name: `R1`.
- **Two-port element (`TF`/`GY`):** write the port: `T.1` or `T.2`.
- **Sign:** an optional `+` or `-` after the name sets the bond's half-arrow
  (the direction of positive power). `+` (the default) points from the
  junction **into the element**, the usual choice for `I`, `C` and `R`. `-`
  points from the element into the junction, the usual choice for sources and
  for the output port of a two-port.

```text
junction armature 1: V-, Ra+, La+, G.1+
```

This means `-e_V + e_Ra + e_La + e_G = 0`, i.e. `e_V = e_Ra + e_La + e_G`.

The sign only sets the **reference direction** of that bond's variables.
Element laws are written in terms of the flow *into* the element (or, for a
source, the flow it delivers), so flipping a sign never changes the physics.
It only flips the sign of that bond's `f_X` in the full equations.

Each element port must be attached exactly once in the whole file.

### Junction-to-junction bonds

There is no direct junction-to-junction bond. To join two junctions with a
plain bond, place a unit transformer between them:

```text
element link TF n=1
connect link j_left j_right
```

The [editor](../ui/index.html) writes these for you, named `bondN`. A direct
bond between two elements is exported through a two-bond 1-junction, also named
`bondN`.

## 3. Elements

```text
element <name> <kind> [key=value ...]
```

`key=value` pairs are separated by whitespace, so **an expression must not
contain spaces** (`R=R0*(1+a*q_C1)`, not `R = R0 * (1 + a*q_C1)`).

| Kind | Law used by the tool | Parameter | Default if omitted | State |
|------|----------------------|-----------|--------------------|-------|
| `I`  | `f_in = p / I`, `dp/dt = e` | `I=`, or `f=` / `p=` | `I_<name>` | `p_<name>` (momentum) |
| `C`  | `e = q / C`, `dq/dt = f_in` | `C=`, or `e=` / `q=` | `C_<name>` | `q_<name>` (displacement) |
| `R`  | `e = R * f_in` | `R=`, or `e=` / `f=`; optional `causality=` | `R_<name>` | none |
| `Se` | `e = value` | `value=` | `u_<name>` | none |
| `Sf` | `f_out = value` (flow delivered) | `value=` | `u_<name>` | none |
| `TF` | `e1 = n*e2`, `f2 = n*f1` | `n=` | `n_<name>` | none |
| `GY` | `e1 = n*f2`, `e2 = n*f1` | `n=` | `n_<name>` | none |

`f_in` is the flow into the element and `f_out` the flow out of it. With the
usual signs (`+` on `I`, `C`, `R`; `-` on sources) these are just the bond's
`f_X`. For two-ports, see [Sign conventions](#sign-conventions).

Examples:

```text
element m   I  I=m
element k   C  C=1/k            # spring: compliance is 1/stiffness
element b   R  R=b
element F   Se value=F0*sin(w*t)
element Q   Sf value=Q(t)
element D   R  R=R_d causality=flow
```

`causality=effort` or `causality=flow` on an `R` forces it to impose that
variable on its junction (see [Causality](#8-causality)).

### Laws as functions

`R=`, `C=` and `I=` give a coefficient, so the law is linear in it. To give
any law instead, name the variable it sets. Use one key per element.

| Kind | Key | Law | Natural causality |
|------|-----|-----|-------------------|
| `R` | `e=` | `e = φ(f)` | sets effort (it is given the flow) |
| `R` | `f=` | `f_in = φ(e)` | sets flow (it is given the effort) |
| `C` | `e=` | `e = φ(q)` | integral |
| `C` | `q=` | `q = φ(e)` | derivative |
| `I` | `f=` | `f_in = φ(p)` | integral |
| `I` | `p=` | `p = φ(f_in)` | derivative |

Inside any `R`, `C` or `I` law, these names mean the element's own variables:

| Name | Meaning |
|------|---------|
| `e` | its effort |
| `f` | the flow into it (`f_in`), whichever way its bond points |
| `q` | its displacement (`C` only), the same as `q_<name>` |
| `p` | its momentum (`I` only), the same as `p_<name>` |

```text
element b  R  e=c*f*Abs(f)                           # quadratic damper / orifice
element D  R  f=Piecewise((e/R_on,e>0),(0,True))     # ideal-ish diode
element D2 R  f=I_s*(exp(e/V_T)-1)                   # Shockley diode
element k  C  e=k*q+k3*q**3                          # hardening spring
element L  I  p=Psi*tanh(f/i0)                       # saturating inductor
```

The law form doesn't constrain causality. When the graph gives an element the
other causality, the tool solves the law for the variable it needs (see
[what works well](#what-works-well-and-what-doesnt)). Because of that, prefer
`f` over `f_<name>` inside an element's own law: `f_<name>` follows the bond's
half-arrow, so flipping the arrow would flip its sign.

## 4. Two-ports: TF and GY

A transformer or gyrator has ports `.1` and `.2` and needs **both** attached.
You can attach them in either of two ways, but not both for the same element.

**Option 1: in the junction lists**, which lets you choose each sign:

```text
junction armature 1: V-, Ra+, La+, G.1+
junction shaft    1: G.2-, J+, b+
element G GY n=K_t
```

**Option 2: with `connect`.** Power flows from the first junction, through the
element, to the second (port 1 signed `+`, port 2 signed `-`):

```text
junction armature 1: V-, Ra+, La+
junction shaft    1: J+, b+
element G GY n=K_t
connect G armature shaft
```

### Sign conventions

Power enters port 1 and leaves port 2:

- **Transformer:** `e1 = n*e2`, `f2 = n*f1`.
- **Gyrator:** `e1 = n*f2`, `e2 = n*f1`.

Here `f1` is the flow into the element at port 1 and `f2` is the flow out at
port 2. These are the textbook conventions (Karnopp, Margolis & Rosenberg).
Because the laws use those flows rather than the raw bond variables, they
conserve power (`e1*f1 = e2*f2`) whatever signs you give the two bonds. The
usual drawing signs port 1 `+` and port 2 `-`, as in the DC motor example.

For the motor, `G GY n=K_t` gives back-EMF `e1 = K_t*ω` and torque
`e2 = K_t*i`.

In the variables, port 1 of element `G` is `e_G`/`f_G` and port 2 is
`e_G_2`/`f_G_2`.

## 5. Expressions

Every parameter value and source value is parsed as a
[SymPy](https://www.sympy.org) expression.

- **Operators:** `+ - * / **`. `^` is also accepted as a power.
- **Symbols:** any undeclared name becomes a symbolic parameter: `R=R_a`,
  `C=A/(rho*g)`.
- **Time:** `t` is time: `value=F0*sin(w*t)`.
- **Inputs:** an unknown name followed by `(` becomes an undefined function,
  which is how inputs are written: `u(t)`, `v(t)`, `Q(t)`.
- **SymPy functions** keep their meaning: `sin`, `cos`, `exp`, `log`, `sqrt`,
  `tanh`, `Abs`, `sign`, `Max`, `Min`, `Heaviside`, `Piecewise`, and so on.
  (Consequence: an input called `beta(t)` or `gamma(t)` means the SymPy
  function, not an unknown input. Pick another name.)
- **Numbers** are fine anywhere: `R=100`, `C=1e-3`.

### Reserved names

These bare names mean something to SymPy and **must not be used as your own
parameter names**:

| Name | SymPy meaning |
|------|---------------|
| `I` | imaginary unit (so write `I=m`, never `I=I`) |
| `E` | Euler's number e |
| `pi` | π (fine to use *as* π) |
| `oo` | infinity |
| `S`, `N`, `Q`, `O` | SymPy internals; using them causes an error or nonsense |
| `lambda` | Python keyword; causes a parse error |

`Q(t)` *with* parentheses is fine: it becomes an input function. Only the bare
names are a problem.

Inside an `R`, `C` or `I` law, `e`, `f`, `q` and `p` are that element's own
variables ([Laws as functions](#laws-as-functions)), so they can't be used as
parameter names there.

## 6. Modulated elements

Any parameter can depend on **states**, **inputs**, **time**, or **bond
variables**, so modulated (`MTF`, `MGY`) and modulated-parameter elements
need no special syntax.

| To refer to | Write |
|-------------|-------|
| State of `C` element `X` | `q_X` |
| State of `I` element `X` | `p_X` |
| Effort / flow on one-port `X` | `e_X` / `f_X` (`f_X` along the bond's half-arrow) |
| The element's own effort / inflow / state, inside an `R`, `C` or `I` law | `e` / `f` / `q` or `p` |
| Effort / flow on port 2 of two-port `X` | `e_X_2` / `f_X_2` |
| Input | any `name(t)` |

```text
element G  GY n=k*q_C1             # gyrator modulated by a capacitor's state
element T  TF n=r(t)               # transformer with a time-varying ratio
element Rv R  R=R0*(1+alpha*q_C1)  # state-dependent resistance
```

Caveats:

- **Integral causality:** `q_X`/`p_X` is only a state if `X` ends up in
  integral causality. If `X` is forced into derivative causality, the symbol
  is treated as an unknown and the solver usually can't eliminate it.
- **Power:** a modulated TF or GY still conserves power. Don't use modulation
  to add or remove energy; model that with sources or resistors.

## 7. Piecewise and nonlinear laws

Write piecewise values with SymPy's `Piecewise`. It takes `(value, condition)`
pairs, checked in order, and the last condition is usually `True`:

```text
Piecewise((value_1,condition_1),(value_2,condition_2),...,(value_else,True))
```

Remember: **no spaces**.

```text
# Outlet resistance switches when the tank volume passes V_max
element Outlet R R=Piecewise((R_low,q_Tank>V_max),(R_high,True))

# A spring that stiffens beyond a gap
element k C C=Piecewise((1/k1,Abs(q_k)<gap),(1/k2,True))

# Saturated input
element F Se value=Max(-F_max,Min(F_max,u(t)))
```

Conditions can use `<`, `>`, `<=` and `>=`, and can be combined with `&`
(and) and `|` (or): `(q_C1>0)&(q_C1<1)`. Use `Heaviside(x)` or `sign(x)` for
step or sign behavior.

### What works well and what doesn't

A piecewise or nonlinear expression can be a **coefficient** (`R`, `C`, `I`,
`n`), a **source value**, or, for `R`, `C` and `I`, the **whole law**
([Laws as functions](#laws-as-functions)).

- **Works well:** conditions on **states, inputs and time** (`q_*`, `p_*`,
  `u(t)`, `t`). The solver treats them as known, and the result is a
  piecewise state equation:

  ```text
  d(q_Tank)/dt = Piecewise((Q(t) - g*q_Tank*rho/(A*R_low), V_max < q_Tank),
                           (Q(t) - g*q_Tank*rho/(A*R_high), True))
  ```

- **Works well:** a law on the element's **own variables**, written for the
  causality it ends up in. A diode the graph gives an effort is written
  `f=Piecewise((e/R_on,e>0),(0,True))`, and its state equation comes out
  explicit.

- **Usually fails:** a `Piecewise`, `Abs` or `sign` law that has to be
  **inverted** for the causality it ends up in (the same diode where the graph
  gives it a flow, or `R=Piecewise((R_on,f>0),(R_off,True))`, which hides the
  flow inside the coefficient). SymPy can't invert these, so the tool prints
  the equations it could reduce, then "Unable to eliminate algebraic
  variables" and the full equations. Those are still correct, so you can
  finish the reduction by hand. Smooth laws such as `exp` or `tanh` are often
  inverted fine. If one has several inverses, the tool keeps the one that's
  real. Alternatives:
  - **Write the law the other way round:** `e=` instead of `f=`, or
    the reverse.
  - **Condition on a state:** rewrite the condition in terms of a state that
    determines the same switch.
  - **One mode per file:** analyze each mode in its own file with a constant
    `R`.

The [editor](../ui/index.html) solves every law numerically, so it simulates
laws in either direction. The exception is a law whose slope is infinite at
the starting point, such as `f=sign(e)*sqrt(Abs(e)/a)` at `e=0`. Its solver
can't start there; write it as `e=a*f*Abs(f)` instead.

## 8. Causality

The tool assigns causality automatically.

- **Sources:** `Se` imposes effort and `Sf` imposes flow.
- **Junctions:** a `0` junction has exactly one port that sets the effort, and
  a `1` junction has exactly one port that sets the flow.
- **Two-ports:** a `TF` passes causality straight through. A `GY` flips it, so
  both of its ports impose the same kind of variable.
- **Storage:** among all valid assignments, the one with the most `I` and `C`
  elements in **integral causality** wins. `I` imposes flow and `C` imposes
  effort.
- **Resistors:** an `R` takes whatever the rest of the graph leaves it, unless
  you fix it with `causality=effort` or `causality=flow`. Use this when its
  law can only be solved one way. Writing the law as `e=` or `f=` doesn't fix
  the causality on its own.

A storage element left in **derivative causality** is reported in the
causality list and left out of the state vector.

### Reducing derivative causality away

Derivative causality usually means two storage elements compete for one
variable: two `I`s on a 1-junction (one flow), two `C`s on a 0-junction (one
effort), or the same pair on either side of a `TF` or `GY`. Before solving, the
tool reduces the graph the way you would by hand, keeping each step only if it
leaves fewer elements in derivative causality:

- **Reflect** the `R`, `C` and `I` elements on one side of a two-port onto the
  other side, then merge the two junctions. Through a `TF` with
  `g = n` (port 2 reflected onto port 1) or `g = 1/n` (port 1 onto port 2), an
  `R` becomes `g²R`, an `I` becomes `g²I` and a `C` becomes `C/g²`. Through a
  `GY`, an `I` becomes a `C` of `I/n²`, a `C` becomes an `I` of `n²C`, and an
  `R` becomes `n²/R`. Laws written as functions are reflected too.
- **Combine** the competing pair: inertias on a 1-junction add, and so do
  capacitances on a 0-junction.

The equations are still in the original states: the combined element's state
is rewritten in terms of the member that was in integral causality, and the
other member's state is reported as a function of it. For
[`fluid_piston.bond`](../examples/fluid_piston.bond) (fluid inertia `If`
driving a mass `m` through `TF1 n=A`):

```text
Reduced the graph to remove derivative causality:
  reflected b, C and m through TF1 onto J1: R_b = A**2*b, C_C = 1/(A**2*k), I_m = A**2*m
  combined m and If on J1 into one inertia: I = A**2*m + If

States: q_C, p_If

State equations:
  d(q_C)/dt = A*p_If/If
  d(p_If)/dt = (-A**2*b*p_If - A*If*k*q_C + If*P(t) - Rf*p_If)/(A**2*m + If)

Storage in derivative causality (follows the states):
  p_m = A*m*p_If/If
```

A side is reflected only if it holds nothing but `R`, `C` and `I` elements, and
only if no law reads `t` or another element's variables (no modulation, no
`causality=`). Derivative causality that a source forces, such as an `Sf` on a
1-junction with an `I`, can't be reduced away. That storage just follows the
source.

The search is exhaustive, so models are limited to **22 ports** in total
(each one-port counts 1, each TF/GY counts 2). For larger systems, split the
model.

## 9. Reading the output

```text
Causality (element port -> quantity caused onto junction):
  V.1: effort
  La.1: flow
  ...

States: p_La, p_J

State equations:
  d(p_La)/dt = v(t) - R_a*p_La/L_a - K_t*p_J/J
  d(p_J)/dt = K_t*p_La/L_a - b*p_J/J
```

- **Causality:** what each port imposes on its junction. For `I` and `C`
  ports, `flow` and `effort` respectively mean integral causality.
- **Reduced the graph:** only when storage was in derivative causality; the
  reflections and combinations made first (see
  [Reducing derivative causality away](#reducing-derivative-causality-away)).
- **States:** the state vector, built from integral-causality storage only.
- **State equations:** derivatives written in terms of states, inputs and
  parameters.
- **Storage in derivative causality (follows the states):** after a
  reduction, each storage element that isn't a state, written in terms of the
  states.
- **Unable to eliminate algebraic variables:** SymPy couldn't reduce the
  system. The full equations (each `= 0`) follow so you can see what's left.
  `--show-equations` always prints them.

### LaTeX and matrix form

`--latex` prints the same state equations as a LaTeX `aligned` block.
`--latex=matrix` prints them as `ẋ = A x + B u`:

```text
\left[\begin{matrix}\dot{p}_{\mathrm{La}} \\ \dot{p}_{J}\end{matrix}\right] = \left[\begin{matrix}- \frac{R_{a}}{L_{a}} & - \frac{K_{t}}{J}\\\frac{K_{t}}{L_{a}} & - \frac{b}{J}\end{matrix}\right] \left[\begin{matrix}p_{\mathrm{La}} \\ p_{J}\end{matrix}\right] + \left[\begin{matrix}1\\0\end{matrix}\right] \left[\begin{matrix}v{\left(t \right)}\end{matrix}\right]
```

- **Inputs:** `u` has one entry per `Se`/`Sf`. When a source's value is a
  single name or input, such as `v(t)` or `Q_in`, that's the entry; otherwise
  it's `u_<source>`, defined under the matrices ("where `u_F = F_0 θ(t − 0.5)`").
  A source whose value depends on a state or bond variable is part of the
  dynamics, not an input.
- **Nonlinear models:** if the equations aren't linear in the states and inputs
  (a state-switched `Piecewise`, `q**3`, ...), there's no `A` and `B`, and the
  matrix form is the vector `ẋ = f(x, u)`.
- **Names:** a one-letter subscript stays italic (`R_a` → R<sub>a</sub>) and a
  longer one is upright (`q_C1` → q<sub>C1</sub>, `R_on` → R<sub>on</sub>).
  Greek names (`rho`, `omega_n`) become Greek letters.

The [editor](../ui/index.html)'s **Equations** window shows the same output.

## 10. Errors and troubleshooting

| Message | Cause / fix |
|---------|-------------|
| `line N: expected junction, element, or connect` | Typo in a keyword, or a missing `:` after the junction type. |
| `line N: expected key=value, got '...'` | A space inside an expression. Remove it. |
| `line N: R takes R=, e=, f=, causality= not X=` | A key that this kind of element doesn't have; see the [Elements](#3-elements) tables. |
| `line N: X has two laws (R= and f=); give one` | An element was given its law twice. Keep one form. |
| `line N: bad attachment '...'` | Junction member isn't `Name`, `Name.1`/`Name.2`, optionally followed by `+`/`-`. |
| `junction J: unknown element X` | `X` is attached but never declared with `element`. |
| `two-port X needs .1 or .2` | A TF/GY is attached in a junction list without a port. |
| `one-port X cannot have a port suffix` | `.1`/`.2` used on a one-port element. |
| `element X.N attached twice` | The same port appears in two junctions (or twice in one). |
| `element X: needs N attached port(s)` | Element declared but not connected, or a TF/GY with only one port attached. |
| `connect X: use either connect or junction attachments, not both` | Choose one way of attaching the two-port. |
| `no valid causality` | Conflicting sources (e.g. two `Se` on one `0` junction), or a junction with nothing to set its common variable. |
| `model has over 22 ports` | Split the model. |
| SymPy error / traceback while parsing | A [reserved name](#reserved-names) was used as a parameter, or the expression is malformed. |
| "Unable to eliminate algebraic variables" | See [Section 7](#what-works-well-and-what-doesnt). Also check for storage in derivative causality and algebraic loops. |

## 11. Quick reference

```text
# comment
junction <name> 0: A+, B-, T.1+       # 0: common effort, signed flows sum to 0
junction <name> 1: C+, T.2-           # 1: common flow, signed efforts sum to 0

element <name> I  I=<expr>            # state p_<name>; or f=φ(p) / p=φ(f)
element <name> C  C=<expr>            # state q_<name>; or e=φ(q) / q=φ(e)
element <name> R  R=<expr> [causality=effort|flow]   # or e=φ(f) / f=φ(e)
                                      # in R, C, I laws: e, f, q, p are its own variables
element <name> Se value=<expr>
element <name> Sf value=<expr>
element <name> TF n=<expr>            # e1 = n e2,  f2 = n f1  (power in at 1, out at 2)
element <name> GY n=<expr>            # e1 = n f2,  e2 = n f1

connect <TF|GY name> <junction for port 1> <junction for port 2>   # power flows 1 -> 2
```

Working examples are in [`examples/`](../examples):

- `rc.bond`: RC circuit driven by a current source.
- `dc_motor.bond`: armature-controlled DC motor, using a gyrator.
- `piecewise_leak.bond`: tank with a piecewise, state-dependent outlet
  resistance.
