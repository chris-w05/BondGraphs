# Summary

Bond graphs are build around 3 objects: Bonds, Junctions, and elements. A bond carries power flow between a junction and element, a junction distributes power flow, and an element does 'something' with the power. 

## Bonds

A bond connects a junction to an element. It has 4 traits: Direction, causality, effort, and flow. The direction for a bond is the notional direction of power flow, used when deriving state equations later. Causality is how to determine where the source of the effort and flow are for the bond, which allows traceability of efforts/flows when deriving state equations. An effort is the 'force' on a bond, i.e. Voltage, Pressure, Force, Torque, temperature, etc.. A flow is the rate that produces a power when multiplied with an effort, this is in the form of mass flow rate, volumentric flow rate, velocity, current, rotational velocity, entropy change rate, etc. When a bond has causality, it means it causes effort in one direction, which means it ALWAYS causes flow in the opposite direction. 

## Junctions

There are two types of junctions: 0 and 1 junctions:

### 0 Junctions:

A zero junction can have an unlimited number of elements connected to it. All of the bonds connected to a 0 junction have the same effort. The sum of all bonds connected to a 0 junction is always 0. Lastly, a zero junction can only have 1 bond that causes effort on the 1 junction, and all other bonds on the junction must cause a flow. 

### 1 Junctions:

A one junction is similar to a 0 junction, but has the exact same flow for all bonds connects. Similarly, the sum of all efforts going into a 1 junction is always zero. A one junction can only have one bond that causes flow on it, and all other bonds must apply effort. 

## Elements:

### I: Inertance: 

An inductance relates a flow to momentum. This is defined in the generic constitutive law P = If where I is the inductance, P is the momentum of the element and f is the flow of the element. Examples are things such as inductors, masses, and long pipes in hydraulics. For an inertance to have a unique state, it must be in integral causality, meaning that it causes flow on whatever bond it is connected to. 

### C: Capacitance:

Capacitances relate forces to displacements. This is defined by the equation: q = Ce where 1 is the charge in the element, C is the capacitance (think electrical capacitor C, or mechanical spring 1/k) and e is the effort carried by the bond. Capacitive elements must cause effort on the junction they are connected to in order to have a unique state. 

### R: Resistance:

Resistances relate efforts and flows, generally through the form e = rf. When a resistance has a non-invertible function it creates a causality constraint. For example, a diode is a resistive element that is exclusively in conductive causality because it creates a flow for any given effort, but does not create an effort with a flow as a input. This means a diode must always cause flow on whatever junction it is connected to. Generalized, if a resistive element has a non-reversible constitutive law, it specifies the causality of the resistor depending on whether its function is defined across all values for effort/flow

### TF: Transformer

A transformer relates two efforts or two flows. The only rule for a transformer is that it must be power-conservative, meaning e1f1 = e2f2 where e1, f1 are the effort and flwo into the transformer, and e2, f2 are the effort and flows leaving the trasnformer. Tranformers generally have a constant modulus that determines the scale, meaning e1 = n * e2, f1 = (1/n) * f2. Transformers can also be modulated, where the scale of multiplication is controlled via a function of another variable( usually a state of an I or C somewhere else on the graph). this would look like n = function(P). Notably, this function does not have to be invertible. Some real worl examples include gears, levers, and electric transformers. When assigning causality, it simply transmits the causality, so if an element causes flow on the transformer, the transformer  causes flow on its other side, and vice versa.

### GY: Gyrator:

A gyrator behaves very similarly to a Transformer, but instead of mapping effort->effort it maps effort->flow (meaning e1 = n * f2, and f1 = (1/n) * e2). The same function requirements apply as the transformer. Gyrators do the opposite of transformers when it comes to causality; when an effort is caused onto a gyrator it causes a flow on its other side, and when a flow is caused on one side, a gyrator causes an effort on the other side.

### Se, Sf: Sources/sinks

Se (effort sources) and Sf (flow sources) are simple. They impose efforts or flows into junctions. As such an effort source will always cause effort on a junction, and a flow source will always cause flow on a junction. 

# Solving Causality:

Once a bond graph exists, causality can be assigned. This includes the following steps:

1: Assign causality to elements where their causality is already defined (Non-invertible R, Se, Sf)

2: Apply 'preferred' (Integral) causality: Inertances cause flow, and capacitances cause effort

2: Propogate causality through bond graph, using rules for 0 and 1 junctions (see junction section)

3: If the propogation disagrees with the preferred causality of I's and C's, then switch the causality of the I or C elements as necessary this will put them into derivative causality, and they will no longer be states. If the causality for the graph is indeterminate, then arbirarily choose a solution, with the note that there will be algebraic loops in the state equation solutions. If an I or C is in conflict with a Non-invertible R, then the graph will need to be modified as no state equations can be derived. 

The number of state equations for a system is the same as the number of elements in integral causality (I causing flow, C causing effort). After solving state equations, there will be a matrix of n first order differential equations. With source and sink functions defined, the simulation can be simulated. 