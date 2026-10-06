# Writing `.bond` files

A `.bond` file describes a bond graph as a list of **junctions** and
**elements**. `bondgraph.py` reads it, assigns causality (preferring integral
causality on storage elements), and derives symbolic state equations.

```sh
python3 bondgraph.py model.bond                  # causality, states, state equations
python3 bondgraph.py model.bond --show-equations # also print every junction/element equation
python3 bondgraph.py model.bond --json           # machine-readable output
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
| `I`  | `f_in = p / I`, `dp/dt = e` | `I=` | `I_<name>` | `p_<name>` (momentum) |
| `C`  | `e = q / C`, `dq/dt = f_in` | `C=` | `C_<name>` | `q_<name>` (displacement) |
| `R`  | `e = R * f_in` | `R=`, optional `causality=` | `R_<name>` | none |
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

## 6. Modulated elements

Any parameter can depend on **states**, **inputs**, **time**, or **bond
variables**, so modulated (`MTF`, `MGY`) and modulated-parameter elements
need no special syntax.

| To refer to | Write |
|-------------|-------|
| State of `C` element `X` | `q_X` |
| State of `I` element `X` | `p_X` |
| Effort / flow on one-port `X` | `e_X` / `f_X` (`f_X` along the bond's half-arrow) |
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

The element laws stay in the forms in the [Elements](#3-elements) table. A
piecewise or nonlinear expression is the **coefficient** (`R`, `C`, `I`, `n`)
or the **source value**. It can't replace the law itself.

- **Works well:** conditions on **states, inputs and time** (`q_*`, `p_*`,
  `u(t)`, `t`). The solver treats them as known, and the result is a
  piecewise state equation:

  ```text
  d(q_Tank)/dt = Piecewise((Q(t) - g*q_Tank*rho/(A*R_low), V_max < q_Tank),
                           (Q(t) - g*q_Tank*rho/(A*R_high), True))
  ```

- **Usually fails:** conditions on the element's **own bond variables** (for
  example a diode written as `R=Piecewise((R_on,f_D>0),(R_off,True))`). SymPy
  generally can't invert such a law, so the tool prints "Unable to eliminate
  algebraic variables" and the full equations. Those equations are still
  correct, so you can finish the reduction by hand. Alternatives:
  - **Smooth approximation:** replace the switch with a smooth function, e.g.
    `R=R_off+(R_on-R_off)*(1+tanh(f_D/eps))/2`. SymPy may still not solve it
    in closed form.
  - **Condition on a state:** rewrite the condition in terms of a state that
    determines the same switch.
  - **One mode per file:** analyze each mode in its own file with a constant
    `R`.

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
  law can only be solved one way.

A storage element left in **derivative causality** is reported in the
causality list and left out of the state vector.

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
  d(p_J)/dt = -K_t*p_La/L_a - b*p_J/J
```

- **Causality:** what each port imposes on its junction. For `I` and `C`
  ports, `flow` and `effort` respectively mean integral causality.
- **States:** the state vector, built from integral-causality storage only.
- **State equations:** derivatives written in terms of states, inputs and
  parameters.
- **Unable to eliminate algebraic variables:** SymPy couldn't reduce the
  system. The full equations (each `= 0`) follow so you can see what's left.
  `--show-equations` always prints them.

## 10. Errors and troubleshooting

| Message | Cause / fix |
|---------|-------------|
| `line N: expected junction, element, or connect` | Typo in a keyword, or a missing `:` after the junction type. |
| `line N: expected key=value, got '...'` | A space inside an expression. Remove it. |
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

element <name> I  I=<expr>            # state p_<name>
element <name> C  C=<expr>            # state q_<name>
element <name> R  R=<expr> [causality=effort|flow]
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
