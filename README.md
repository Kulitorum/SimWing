# SimWing

[![CI](https://github.com/Kulitorum/SimWing/actions/workflows/ci.yml/badge.svg)](https://github.com/Kulitorum/SimWing/actions/workflows/ci.yml)

SimWing is a separate GPL-3.0 engineering project for coupled XPBD-CFD
simulation of ram-air paragliders, including the long-term goals of dynamic
inflation, collapse, and reinflation. It is bootstrapped from
[LEparagliding Studio](https://github.com/Kulitorum/LeParagliding) and preserves
that project's exact wing-design engine, Qt/OCCT application, and soft-wing
structural baseline. SimWing starts from an independent root commit; upstream
provenance is recorded in [`UPSTREAM.md`](UPSTREAM.md), not copied Git history.

The imported application and executables still use their LEparagliding names
while the solver-independent scene format, structural adapter, CFD worker, and
Flight Lab are developed. The implementation roadmap and validation boundary
are recorded in
[`docs/coupled-fsi-architecture.md`](docs/coupled-fsi-architecture.md).
SimWing reuses the isolated XPBD structural primitives, but its fluid and flight
simulation is a remake: inherited Playground aerodynamics, pressure, cell-air,
flight, and scenario behavior are not compatibility targets.

Numerical development includes a standalone Qt/OpenGL diagnostic viewer from
the first implementation slice. Interactive runs open it by default so the
live structure, loads, contact, coupling residuals, and progressively added CFD
fields can be inspected; automated runs remain available through `--no-viewer`.
The separate Qt-free `simwing-mimetic-conductance-audit` tool runs deliberately
opt-in, graph-free pressure-cell conductance samples at one resolution and one
canonical grid phase (or all eight phases explicitly). It is an offline
convergence-evidence runner and has no path to worker pressure or loads.
For example, `simwing-mimetic-conductance-audit --resolution 64 --phase 0`
runs one expensive selected trajectory; `--all-phases` is intentionally
required for a complete ensemble and streams one validated phase product at a
time to keep retained audit storage bounded.
The current complete `64^3` developer run accepts all eight placements but
still fails the established trend screen: mean drift barely contracts and
three phase trajectories reverse direction. It remains evidence against a
live pressure-owner switch, not authorization for one.

The first remake foundations are now present: scene-v2 validation and bounded
serialization, direct export from captured analytical wing geometry,
deterministic scene-to-XPBD membrane/cable/bending assembly, production
suspension/payload/contact checkpoints, rigid-pilot/harness assembly, immutable
structural diagnostic frames, a
standalone Qt/OpenGL trace viewer, and a canonical Qt-free worker that opens
the viewer by default and streams accepted XPBD steps through a replayable
growing trace. The first Qt-free fluid kernel adds a periodic staggered-grid
pressure projection with deterministic rollback behavior, exact discrete
gradient/divergence pairing, Taylor-Green preservation, and manufactured
second-order pressure convergence. Its first velocity-evolution operator adds
explicit laminar viscosity directly on the periodic MAC components, preserves
momentum and solenoidal Fourier modes, dissipates kinetic energy, and rejects
steps above the sharp diffusion-number limit without mutating the field. It is
a stage oracle for a two-stage SSPRK2 viscosity path. The SSPRK2 result is a
convex average of the old field and two accepted Euler candidates, retains the
same per-stage stability/conservation contract, and exhibits second-order
temporal convergence against the exact discrete Fourier decay.
A companion uniform-flow transport oracle applies an unsplit conservative
donor-cell update to all periodic MAC components. At total absolute CFL one it
performs an exact one-cell periodic shift; below that limit it preserves
component momentum, the maximum principle, non-increasing kinetic energy, and
discrete solenoidal modes. It is intentionally first order and does not yet
claim higher order. A variable-flow companion reconstructs one shared upwind
flux on every periodic face of each translated component control volume. It
requires a divergence-free MAC advector, delegates uniform fields bit-exactly
to the oracle, and supports safe nonlinear self-advection. Under its local
outgoing-CFL limit it preserves component momentum and bounds without adding
kinetic energy; a periodic shear fixture observes first-order convergence.
A selectable monotonized-central reconstruction limits a second-order slope on
those same shared faces. Its forward-Euler stages may add the expected
second-order-in-time energy defect only while private inside SSPRK2; the
committed convex update rechecks original component bounds, momentum, and
non-increasing energy. A discontinuous pulse remains bounded and smooth
uniform-flow L1 refinement approaches second order.
A pressure-projected nonlinear SSPRK2 operator now uses two self-advection
stages without feeding a divergent intermediate field back as its own
advector. It projects stage one, advances stage two, convexly averages that
prediction with the old field, and projects the result. All four candidates
remain private until pressure, velocity, momentum, energy, and continuity pass;
a fixed-grid vortical refinement observes second-order temporal behavior. A
Galilean-translated Taylor-Green vortex, run with `dt` proportional to `h^2`,
separates that temporal error and observes first-order donor-cell versus
near-second-order limited-MC L1 spatial refinement on 16/32/64 grids. This is a
smooth periodic nonlinear canonical, not yet a cut-cell or moving-interface
accuracy claim. A symmetric full flow path now reconciles those pressure stages
with viscosity: SSPRK2 diffusion
advances half a step on each side of the full projected nonlinear transport
step. The resulting Strang composition is transactional across all three
sub-integrators, closes their energy-loss sum, and observes second-order
fixed-grid temporal refinement. An independent continuously viscous,
Galilean-translated Taylor-Green solution also observes near-second-order L1
refinement for the complete limited-MC path on 12/24/48 grids with `dt`
proportional to `h^2`. This remains a smooth periodic oracle, not a cut-cell or
moving-interface accuracy result.
A bounded wrapper advances a requested outer interval through equal Strang
substeps. It pre-sizes the schedule from the explicit viscous limit, then
restarts the whole private interval with more substeps only when transport
reports a CFL or limited-reconstruction bound rejection. Projection and ledger
failures are never retried, the configured substep count is capped, and caller
pressure/velocity commit only after the complete outer conservation ledger
passes.
The `simwing-fsi --case periodic-flow` worker now exercises that complete path
on an `18 x 18 x 2` Galilean-shifted Taylor-Green field. Every accepted state
is copied into an immutable trace frame with one cell-centre point per cell,
exact cell pressure, averaged MAC velocity, speed, the projection operator's
finite-volume divergence, diagnostic centred-curl vorticity, and conservation
ledgers.
The standalone viewer renders those otherwise-unconnected points with selectable
scalar colouring. Its independent vector selector adds normalized velocity or
vorticity arrows while retaining the chosen scalar colour field.
Arrow construction is
Qt-free, owning, deterministically bounded for large fields, and has no path
back into solver data. This makes the verification flow inspectable without
making it a cut-cell, canopy, or aerodynamic-truth case.
The periodic worker advances candidates through frame validation before
committing them. Its immutable in-memory checkpoint binds the exact grid and
case definition to pressure, velocity, last accepted diagnostics, step, and
time; the initial state and later accepted states resume bit-for-bit in the
same or an equivalent rebuilt worker. A bounded, versioned, checksummed
little-endian file codec preserves that complete payload. The worker accepts it
through `--checkpoint-in` and atomically replaces `--checkpoint-out`; requested
steps are additional after restore. For long runs, `--checkpoint-every N` also
saves after accepted absolute step multiples and at the final state, so a
resumed run keeps the original cadence and never writes a rejected candidate.
Moving-porous-flow checkpoint validation regenerates bounded canonical history;
batch checkpoint output is therefore limited to 10,000 total accepted steps.
Runs that would exceed that bound fail before opening the trace or checkpoint,
including after adding the restored step count.
A first transport-neutral worker-control protocol now defines versioned,
checksummed, byte-bounded `advance`, `checkpoint`, and `stop` commands with
correlated `ready`, `advanced`, `checkpointed`, `stopped`, and coded-error
responses. Every response carries the absolute accepted step and simulated
time. A shared Qt-free control session executes those decoded commands on the
worker owner thread, with typed periodic-flow, moving-porous-flow, open-piston,
strong-piston, and porous-sheet adapters. Each
publishes immutable accepted frames through an injected sink, delegates its
own complete checkpoint codec, and makes stop terminal. `simwing-fsi
--control-stdio` exposes that boundary for all five cases as self-framed binary
stdin/stdout messages. It launches no viewer,
puts no human-readable text on protocol stdout, completes the trace before the
stopped response, and uses `--checkpoint-out` as the checkpoint-command target.
Supplying `--checkpoint-in` makes the initial Ready response expose the restored
absolute step/time; the new trace starts with the next accepted frame rather
than replaying earlier frames. The moving-porous-flow regression advances
through its second periodic wrap, persists step 101, proves exact step-102
continuation in the original trace, and resumes with that one exact next frame.
The porous-sheet end-to-end control regression
checkpoints its rebased step 330 and verifies exact step-331 continuation in
both the original and resumed traces, then drives the original worker to its
pump-collision rejection and proves that only preceding accepted frames were
published. The terminal accepted endpoint is itself persistable even when it
has moved beyond the face segment used by that step's constitutive midpoint.
A separate process-level regression checkpoints that terminal state, restores
it into a fresh controlled worker, and verifies the same rejection with an
empty completed continuation trace.
Ordinary step-driven CLI behavior is unchanged.
Both the projected transport and Strang flow paths select donor-cell or limited
monotonized-central reconstruction explicitly; donor-cell remains the default.
The original single-stage transport and Euler/SSPRK2 viscosity modes compose
with the zero-mean
pressure solve in one transactional periodic fluid step. Every stage runs on
candidate fields; velocity and pressure commit together only when transport
bounds, viscous stability, projection convergence, and the final
momentum/energy ledger all pass. Focused rollback cases fail each stage
independently without changing either caller field. This remains a periodic
verification kernel; general cut-cell evolution and whole-wing CFD are not yet
implemented.
A validated ordered sharp-interface field now preserves prescribed two-sided
static pressure jumps without smearing or spurious flow, including across the
periodic domain boundary. Multiple crossings on one face-normal cell segment
retain their stable surfaces and region sequence while the paired gradient and
Poisson stencils use their deterministic signed sum. A split-region slab is
bit-identical to its compact jump and a balanced folded subcell pocket creates
no spurious pressure or flow.
The axis-aligned face epoch used by moving porous sheets is now a reusable
planar topology primitive, with the established porous API delegating to the
same arithmetic. A thin multi-layer pressure oracle builds complete X/Y/Z
closed region chains and rigidly translates them through adjacent faces and
periodic wraps. Its `+70/-70 Pa` pocket deliberately records the current
limitation: both same-face layers and fractions remain visible, but their dense
stencil sum is zero. A separate static regional profile now partitions the
exact one-period interval volumes and reconstructs a consistent pressure
potential, recovering the physical thin-pocket volume and 70 Pa pressure
difference under a declared volume-mean gauge. It creates no regional velocity
degree of freedom and does not change the projection, so this remains
moving-topology groundwork rather than a general folded-fabric or leakage
solve.
A bounded two-epoch regional sweep ledger now binds those profiles by stable
layer identity. It derives each layer velocity, independently integrates the
lower/upper boundary sweep of every interval, and closes that value against
geometric volume change per interval, per region, and globally. Rigid and
breathing pockets close on X/Y/Z through positive and negative periodic
rebases. The ledger still contains no Eulerian regional flux, so it proves the
moving geometry/GCL prerequisite without claiming no-leakage dynamics.
A companion regional flux compatibility screen now fits one uniform axial
fluid velocity to each interval's two moving boundaries. Rigid motion matches
both exactly and is sealed. A full-plane breathing pocket instead reports the
least unavoidable one-sided slip and material-relative flow while independently
closing `volume change + integrated outward relative flow = 0`. This keeps
impermeability distinct from continuity and provides evidence for the next
regional velocity discretization; it still does not advance fluid state.
The screen is published as a nonzero semantic-fingerprinted artifact. Its
validator reconstructs every interval result, stable-ID chain, sorted region
summary, aggregate, tolerance flag, and storage ledger before accepting a
copy, so diagnostics cannot be silently reused after mutation.
A bounded planar opening-flow oracle now provides the compatible counterpart.
It treats authored openings as an oriented region graph and solves the minimum
area-weighted relative flow needed by the regional volume changes. The sealed
breathing pocket remains locally infeasible despite global cancellation; one
`0.5 m²` exterior-to-pocket opening carries `1.6 m³/s` at `3.2 m/s`, parallel
openings split by area, and a serial three-region graph balances its
intermediate region. This is an offline feasibility result, not an embedded
opening discretization, pressure solve, or fluid-state update.
A companion midpoint pressure-power audit prevents “feasible” from being read
as “passively driven.” The inflating 70 Pa pocket closes continuity but moves
`1.6 m³/s` uphill and therefore needs `112 W` from an external aerodynamic or
kinetic source; deflation releases the same power. Opening pressure power plus
regional `p*dV/dt` closes to roundoff. Local uphill deficits and the net graph
deficit are reported separately, without inventing an opening constitutive law
or energy source.
A one-epoch planar regional fragment builder now performs the missing spatial
split. It cuts every Cartesian cell at the current pressure layers and emits
one full-dimensional control per transverse tile, with stable boundary
identity, physical volume/centroid, and regional pressure. The canonical
16-cell grid becomes 24 controls: four cells each contain separate
exterior/pocket/exterior fragments, and the thin pocket remains four `0.6 m³`
layer-to-layer controls. Per-cell volume and first moment, per-region volume,
and the periodic domain all close through X/Y/Z and signed rebases. No
connectivity, velocity basis, or pressure solve uses these controls yet.
A companion fragment topology now pairs every one of those control faces. The
canonical 24 controls close as 72 two-sided links: 64 periodic same-region grid
links and eight pressure-layer walls, producing separate 13.6/2.4 `m³`
exterior/pocket components. Each wall tile retains its authored surface,
orientation, area, centroid, center distance, and signed 70 Pa jump, but has
exactly zero conductance. All six face incidences and boundary area close per
fragment through X/Y/Z, rigid in-segment motion, and signed periodic rebases.
This is connectivity evidence only; it still owns no regional velocity or
pressure solve.
A bounded symmetric pressure operator now acts on that graph without crossing
fabric. Its 24 rows contain 128 directed entries from the 64 same-region links;
the eight layer walls contribute no entry. The canonical regional pressure is
therefore an exact two-component null mode while arbitrary fields satisfy
symmetry, nonnegative link energy, and zero integrated row sum per component.
The unique `area / distance` weight is `160/3 m` and the two-sided diagonal sum
is `320/3 m`. Stable row/entry/gauge identity survives motion within one
topology segment. This remains an ungauged offline operator with no RHS, solve,
velocity update, or production pressure ownership.
A transactional conjugate-gradient oracle now solves compatible corrections on
that operator. It removes only roundoff-sized component RHS imbalance, keeps
the static regional potential separate, and commits each correction with
roundoff-zero volume mean per component. A manufactured two-component field is
recovered to `3e-11 Pa` with a recomputed residual below `2e-12 Pa·m`; a zero
RHS removes correction gauges while preserving every authored 70 Pa wall jump.
Incompatible or iteration-truncated attempts leave the warm start bit-for-bit
unchanged. No velocity divergence, pressure RHS, or production step invokes
this solve yet.
A first static regional face projection now supplies that physical RHS without
changing production. One oriented normal velocity is owned per topology link;
the 64 same-region Cartesian links contribute area-weighted outward flow and
receive the matching `dt/rho` pressure-gradient correction, while all eight
pressure-layer walls must remain exactly zero-flow. Manufactured divergence is
cancelled below `3e-14 m3/s` on X/Y/Z, uniform wall-tangential flow is
preserved, and a failed pressure solve publishes neither corrected velocity
nor its pressure warm start. Projection work has an explicit storage ceiling.
The base overload rejects moving layers, and this opt-in oracle still owns no
face momentum mass, opening conductance, kinetic-energy claim, or production
worker state.
A separate topology-stable fragment volume-rate ledger now supplies that local
moving geometry. It reconstructs each
current control's previous volume from its two exact layer-boundary
displacements and makes `change = current - previous` bit-exact. Canonical
breathing assigns `+1.6/-1.6 m3/s` to the pocket/exterior components while
every fixed Cartesian cell and the global periodic domain remain closed;
rigid motion instead exchanges volume locally with zero component change.
The result closes independently through fragments, cells, regions, and graph
components on X/Y/Z. It is bounded, fingerprinted, and rejects a layer that
crosses a grid segment because appearance/retirement ownership is not yet
defined. The static overload remains unchanged.
The topology-stable moving projection overload now does consume it. Its
continuity residual is exactly `dV/dt + net outward grid flow`, and its
physical RHS is `-rho/dt` times that residual. A rigid `0.1 m` translation
starts with `0.1 m3/s` local geometry demand and pressure-projects the required
same-region redistribution on X/Y/Z while all layer-wall relative flows remain
zero. The corrected residual closes below `1e-11 m3/s`. A sealed breathing
pocket instead exposes its `1.6 m3/s` component deficit, maps it to `3.84
Pa*m`, and rolls back as incompatible rather than fabricating a wall or opening
flux. Duration and source identities are exact, and static projection behavior
remains a separate overload.
A first diagonal face-inertia metric now closes the missing geometric mass
ownership without yet creating a velocity state. Each same-region Cartesian
link has one shared normal-velocity degree of freedom with dual volume `area *
center distance`. Each pressure-layer wall instead has two independent
one-sided trace degrees of freedom, each owning only its adjacent half-volume;
velocity and momentum are never averaged across fabric. The canonical case has
64 shared and 16 trace degrees of freedom, whose `44.8 + 3.2 = 48 m3` total
closes the `16 m3` periodic domain independently on X/Y/Z. Every fragment and
the separate `13.6/2.4 m3` exterior/pocket components close on each axis,
including rigid motion and breathing geometry. This immutable bounded metric
contains no velocity values, density, kinetic-energy acceptance, advection, or
production state.
A source-bound regional velocity state now applies density to that metric and
publishes scalar normal momentum plus kinetic energy for every degree of
freedom. Half-volume contributions reconstruct fragment and component mass on
each axis, while component/global momentum and energy are summed without
collapsing the two wall traces. At `1.25 kg/m3`, the canonical uniform
`[2, -0.5, 0.25] m/s` field has `20 kg` of diagonal mass per axis,
`[40, -10, 5] kg*m/s` momentum, and `43.125 J` kinetic energy. Assigning
different `+1/-2 m/s` values to the two sides of one fabric surface keeps their
exterior/pocket ledgers separate and closes `-2.5 kg*m/s` and `3.25 J`
globally. Motion updates the source fingerprint and component mass while a
uniform field preserves global momentum and energy. This immutable state does
not prescribe wall velocity, prove projection acceptance, transport momentum,
apply pressure work, or enter production.
A projection-energy audit now certifies both static and topology-stable moving
corrections without changing the projector. For each of the 64 shared velocity
degrees it reconstructs `dt/rho * (p_minus-p_plus)/distance`, proves pressure
impulse equals diagonal momentum change, and closes midpoint pressure work to
kinetic-energy change. Static geometry keeps all 16 one-sided fabric traces at
exact zero, closes continuity below `3e-14 m3/s` through X/Y/Z, strips a pure
gradient below `1e-26 J`, and preserves the tangential null field bit-exactly.
The moving overload binds the local volume-rate artifact and sets each trace to
its material-wall velocity. Its affine identity is `delta-K = geometry-pressure
work - correction kinetic energy`: a rigid `0.1 m/s` translation from zero grid
flow receives twice the correction energy from geometry and retains half as new
kinetic energy, with work residuals below `4e-13 J` on every axis. Sealed
breathing remains incompatible. The bounded fingerprinted certificate still
excludes the authored 70 Pa jump, topology rebase, momentum transport, and
production ownership.
A separate authored-jump energy audit now resolves the absolute pressure force
on each fluid side of every fabric tile and proves that their sum is the signed
`area * pressure-jump * normal` load. It publishes the opposite sheet force and
impulse without applying either one. Static walls do exactly zero work; a rigid
`0.1 m/s` X translation transfers `+28 J` through one layer and `-28 J` through
the other, closing globally. The deliberately projection-incompatible breathing
pocket is still useful as a work oracle: its expanding `70 Pa` pocket gives
`-56 J` to the fluid and `+56 J` to the sheet, while exterior/pocket geometry
work closes at `-11.2/-44.8 J`. The same bounded source-bound checks cover
X/Y/Z, but remain diagnostic only: no impulse is applied and no transport,
topology rebase, or production state is introduced.
An immutable regional pressure state now composes the authored and correction
owners only after both independent audits validate against the same after-state,
metric, topology, epoch, and time step. Every fragment retains all three scalar
pressures; every fabric tile retains authored, correction, and total jumps plus
the corresponding sheet-force split. Static manufactured correction exposes a
real transient fabric jump while doing exact zero wall work. Under rigid X/Y/Z
motion, total material-wall work closes per component and globally to authored
jump work plus correction geometry work within `5e-13 J`. This is still a
diagnostic endpoint: it neither applies the published sheet force nor advances
fluid momentum, transport, rebase, or production state.
That endpoint can now be captured as a minimal authored-surface load ledger.
All eight canonical wall tiles retain stable link/surface identity, wrapped
centroid, area, normal, authored/correction/total pressure traction, force,
impulse, and work. Deterministic surface summaries recover the two `4 m2`
planes at `x=-0.8/-0.2 m` and their authored `-280/+280 N` sheet forces; moving
X/Y/Z ledgers close exactly back to the composed pressure state. This defines
the future structural handoff data, but deliberately performs no nodal load
distribution or XPBD/worker mutation.
The projected nonlinear SSPRK2 operator now applies one immutable jump field at
both internal pressure stages, and the first-order, Strang, and retrying
subcycled full-flow paths carry that same topology through every private step.
The balanced slab stays sharp without spurious velocity throughout; empty-jump
overloads remain bit-exact to the original evolution APIs. Powered or moving
jump work still belongs in an explicit coupled energy ledger.
The disconnected moving-interface projection also accepts the same immutable
sharp field on unconstrained faces. Its Poisson source is gauged independently
inside each sealed fluid component, while impermeable moving faces retain sole
ownership of their MAC velocity and reaction; overlapping ownership is rejected
before either field can change. A translating sealed-slab canonical carries a
second static pressure slab inside one component with analytic pressure and no
spurious flow. This is coexistence of two fixed-grid boundary types, not yet a
general moving porous-sheet coupling.
The immutable diagnostic adapter keeps the owning cell pressure and velocity
samples and adds a separate oriented quad for every authored crossing, so
same-face layers remain visible with their region pairing, signed jump, normal,
and subcell fraction. The `pressure-jump` worker exercises that path with 48
ordered transitions in a periodic split slab, preserving the analytic
`-125/+125 Pa` state without spurious velocity. It is a static visibility
oracle, not moving folded topology or cut-cell evolution.
A calibrated Darcy-Forchheimer adapter now converts resolved normal fluid slip
relative to a fabric tile into the same signed sharp-jump representation. Its
exact monotone inverse, canonical X/Y/Z MAC sampling, tile volume flux, and
nonnegative pressure dissipation are tested. A periodic prescribed-flux plane
with an explicit balancing pressure source retains its analytic `195 Pa` loss
without spurious velocity. A transactional Picard solve now closes that
constitutive law against the full periodic sharp projection at either the
endpoint or temporal midpoint, including separately prescribed pressure-source
jumps, moving-sheet-relative velocity, and heterogeneous tile resistance. Its
uniform endpoint Darcy canonical accelerates to the analytic `0.4 m/s` under
`20 Pa` and closes `9.6 W` of porous loss. Midpoint mode exactly closes the
grid pressure-work/porous-loss/kinetic-energy identity and matches the scalar
nonlinear plug oracle. Every accepted result separately exposes total oriented
jump force/impulse on the fluid, jump power/work at the selected constitutive
time, and porous dissipation. Exhausting the nonlinear solve changes neither
field. The first complete periodic macro-step now uses those ledgers to admit
interface-driven momentum and energy without weakening its conservation gate:
the uniform midpoint drive and unforced decay both close analytically, failed
nonlinear projection rolls back the whole step, and empty porous topology is
bit-exact to the original evolution path.
A second-order fixed-grid path now symmetrically brackets the complete
viscosity/transport/viscosity Strang step with two implicit-midpoint porous
half-steps. It closes the summed interface impulse/work ledger, retains porous
dissipation separately, matches the exact documented three-operator
composition, rolls back all prior operators on a later failure, and observes
second-order temporal refinement for the uniform driven canonical. Empty
topology delegates to the original bulk Strang path without changing its
fields or nested diagnostics.
A stage-resolved planar overload now samples one immutable sheet definition at
each porous half-step midpoint. The two stages must keep stable identity,
oriented regions, and resistance, while the pure topology selector permits the
second sample to retain or move by one adjacent segment. Their physical
displacement must also match trapezoidal integration of the two normal-velocity
samples across the `dt/2` interval within explicit absolute/relative
tolerances. A wrapped `3 -> 0`, image `0 -> 1` regression retains both complete
unwrapped epochs and the closed kinematic residual; inconsistent motion,
identity, placement, skipped topology, or a later numerical failure commits
neither field. This advances moving planar porous sources through the complete
symmetric flow step; it does not claim a cut-cell remap.
The `MovingPorousFlowCase` turns that operator into an immutable-frame
canonical: a prescribed `0.4 m/s` sheet begins just before the positive domain
wrap, crosses `face=3,image=0` to `face=0,image=1` inside its first macro-step,
and continues through a second periodic wrap without resetting the fluid. Each
frame owns the final sheet and pump planes plus both stage epochs, kinematic
residual, porous loss, jump work, and flow conservation residual. The Qt-free
`--case moving-porous-flow` worker publishes accepted frames to a completed trace
headlessly. Its immutable in-memory checkpoint owns the fluid fields, accepted
jumps and diagnostics, sheet kinematics, and unwrapped topology epoch, with
bit-identical replay at the initial state and through the second wrap. A
bounded/checksummed `SWMF` codec stores fluid fields and sharp crossings and
regenerates the accepted diagnostics only after bounded canonical replay
matches every stored value. The headless checkpoint flags provide additional-
step resume, absolute-step autosave cadence, final-write deduplication, and
atomic same-path replacement. Coupled structure remains open.
The nonlinear porous iteration can also use the disconnected moving-interface
projector as its inner solve. A translating two-region slab retains exact wall
velocities while endpoint or midpoint porous slip closes on unconstrained
interior faces; the final moving-wall reaction and porous jump ledgers remain
separate. Empty topology is exact to the moving-only path, and a face claimed
by both an impermeable constraint and a porous jump is rejected transactionally.
This is fixed-grid coexistence with moving boundaries, not yet moving porous
cut-cell topology.
A separate pure topology selector now tracks a planar porous sheet by axis,
wrapped MAC face, and signed periodic image. It deterministically retains or
rebases one segment in either direction on X/Y/Z and rejects exact face
placement or skipped segments before any fluid state changes.
Its complete-plane assembler expands that unwrapped physical epoch into one
deterministically ordered porous crossing per transverse MAC tile on X, Y, or
Z, validating stable IDs, two-sided regions, material resistance, kinematics,
and strict segment placement first. The static porous-flow and moving
porous-sheet workers now share this boundary instead of owning X-only loops.
The accepted wrapper can cross the planar open-piston GCL and physical
cut-surface boundaries without discarding its nonlinear status. In the porous
opening canonical, the summed tile flux exactly fills the swept chamber volume,
matches the independently sampled opening transport, and retains the final
moving-surface reaction; rejected or internally inconsistent wrappers cannot
expose either downstream ledger.
A pressure-driven companion advances a uniform fluid plug with an
implicit-midpoint nonlinear solve; every step independently closes pressure
impulse, driving work, porous dissipation, and kinetic energy.
Its viewable worker converges to the analytic `1.74165739 m/s` speed and `250 Pa`
loss while carrying its endpoint interfaces through both pressure stages of a
complete Strang/SSPRK2 step. The grid iteration still offers first-order
endpoint coupling for comparison and a verified midpoint pressure impulse;
moving cut-cell topology remains future work.
Face-aligned moving membranes can now partition stable fluid regions, impose
an exact normal MAC velocity, and project each region transactionally while
retaining its prior pressure gauge. A translating sealed-slab canonical closes
the same `264 N*s` and `66 J` per-wall impulse/work values as the structural
piston case; incompatible sealed fixed-grid volume change remains explicitly
rejected. Equal-sided interface labels can also describe a nonseparating sheet
whose fluid remains connected around a resolved grid path. The first open
planar control volume follows that sheet through one partial cell and
independently closes geometric volume change, surface sweep, and projected
opening transport on X, Y, and Z grids. At an exact cell crossing it can build
a candidate on the next periodic MAC plane, require old/new chamber volumes to
match, and reject skipped planes, changed identities, or collision with the
opening before commit. A topology-bound conservative transfer layer integrates
uniform triangle
traction or explicit barycentric quadrature patches into structural loads while
independently closing force, moment, and rigid-motion power ledgers. A versioned
macro-step coupling layer now integrates nonuniform temporal samples into nodal
impulse, angular impulse,
and work, and applies the result transactionally across XPBD substeps. Its
prescribed moving-piston canonical closes analytic pressure-volume work and
delivers the same total impulse to structural momentum. The same boundary now
owns a definition-bound vector Aitken delta-squared relaxer: its first factor
and dynamic factors are explicitly bounded, degenerate residual changes fall
back deterministically, and an in-memory checkpoint restores the exact
iteration history without accepting foreign vector definitions. This is the
rollback primitive for future repeated fluid/structure iterations, not yet a
strong-coupled worker. A companion convergence decision requires displacement,
velocity, and traction to each pass both an absolute tolerance and a
floor-stabilized relative tolerance, enforces minimum/maximum iteration counts,
and reports exhaustion separately. A macro-step iteration owner now composes
that decision with the Aitken state and current relaxed vector. Its checkpoint
replays the exact next update, failed advances/restores are transactional, and
converged or exhausted states are terminal; solver checkpoints remain beside
it. A separate composite rollback owner now binds that algorithm checkpoint to
one real `Structure` checkpoint and one accepted moving-interface fluid epoch.
It prebuilds and validates all three replacements, then commits them through
no-throw moves, so a foreign or corrupt member changes none of the owners and a
valid restore replays the exact next mutation. A bounded macro-step retry
policy now distinguishes the active attempt, an explicit restore-before-retry
handshake, acceptance, and terminal failure. Exhaustion reduces the step by a
configured factor down to a hard minimum and retry-count budget; its
definition-bound checkpoint reproduces the exact pending or terminal decision.
A one-macro-step state owner composes that policy with the real three-owner
rollback boundary. It accepts only a fresh iteration baseline and exposes one
restore-and-begin operation, ensuring a reduced attempt cannot start from the
discarded Structure, fluid, or Aitken state. Within an active attempt it also
owns a solver-only checkpoint: Structure and fluid rewind before the next
fixed-point solve while the advanced Aitken iterate remains intact. A generic
Qt-free loop driver now invokes a typed physical-solver callback, validates its
accepted epoch, advances Aitken/convergence, performs those iteration rewinds,
and carries exhaustion through bounded full-state retries. Terminal numerical
failure or a callback exception restores the accepted macro-step baseline. Its
bounded result retains one terminal decision, iteration result, time step, and
solver-run count per attempt, so an exhausted large step remains observable
after a smaller step succeeds. The strong-piston canonical supplies the first
real fluid/Structure callback. The
macro-step surface supplies the residual norms from
stable-ID-bound baseline/previous/current kinematics and consecutive nodal
traction results. It uses maximum physical nodal updates, references motion to
the saved macro-step baseline so coordinate and bulk-velocity shifts cancel,
and rejects foreign topology before reduction. Alongside the strict
uniform fluid-to-structure subset, a planar face-resolved bridge clips canonical
MAC tiles against structural triangles and conservatively maps nonuniform face
traction while closing per-face and aggregate area, force, moment, and power
ledgers. The moving-interface overload can explicitly select the complete
fluid constraint reaction—adjacent pressure plus direct constrained-face
impulse—without relabelling it as pressure. Fixed same-region panels are
supported for resolved flow-around paths. Its rigid-normal mode keeps those
material patches while the physical
plane moves and the Eulerian grid plane rebases, provided transverse geometry
remains fixed and fluid/structural normal velocities agree.
`simwing-fsi --case flag` uses the fixed mode for the first flexible
CFD-to-XPBD slice: a sinusoidally accelerating periodic gust is projected around
a finite one-metre reference panel, and its complete face reaction is mapped
conservatively into a 5-by-5 membrane with an edge clamp. The accepted trace
shows pressure, direct constraint traction, complete reaction, mapped nodal
force, divergence, and normal deformation. The displaced fabric is not fed
back to the fixed CFD reference, so this is one-way coupling—not a moving
cut-cell, classical tangential-flow flag, or two-way energy result. The
`simwing-fsi --case ram-cell` canonical extends that same honest one-way
boundary to an open five-panel cavity. Back, side, top, and bottom reference
walls each receive their complete pressure-plus-direct MAC constraint reaction,
and five independent conservative bridges add those loads to one shared-node
89-node XPBD shell. Two perimeter rows clamp the open mouth; the remaining
panels share seam nodes and flex under the projected cross-flow. Frames expose
panel identity, pressure and complete-reaction traction, mapped nodal force,
deformation, net and RMS mouth flow, and divergence. Incompressibility closes
the dead-ended cavity's net mouth flux while permitting local exchange. This
is still fixed-reference one-way deformation: displaced panels do not yet
rebuild CFD cut cells or close a two-way energy ledger. The
`simwing-fsi --case piston` harness crosses that
face-resolved bridge, temporal coupling, XPBD acceptance, and replayable viewer
frames with visible deterministic motion. The same production target now also
contains a headless strong-coupling canonical: a `6 kg` tributary-mass piston
against `28.8 kg` of projected fluid. Its scalar interface speed is solved by
the generic Aitken loop, while every iteration reruns the actual moving-interface
projection, face-resolved pressure transfer, temporal integration, and XPBD
step from the accepted baseline. Only the converged Structure/fluid epoch
persists. It is selectable as `simwing-fsi --case strong-piston`; accepted
immutable frames expose the real coupling count, residuals, retry count,
interface closure, pressure traction, and transferred impulse/work. The
projection's endpoint force slope recovers its analytic `10.8 kg` discrete
added mass. Because force is trapezoidally integrated from the unchanged start
load, the accepted speed matches `F0*dt/(6 kg + 10.8 kg/2)` rather than the
unstable weak one-pass acceleration; both mass and speed residual are published.
Its in-memory accepted-state checkpoint composes the real Structure and
moving-interface fluid payloads, rejects foreign public identities/topology
transactionally, and reproduces the exact next strong-coupled result and frame.
No in-progress iteration or rejected attempt is checkpointable, and this
accepted boundary now has a deterministic bounded `SWSPCKP1` envelope that
nests the existing Structure and moving-interface fluid codecs. It validates a
fresh canonical owner before encode/decode, checks an outer checksum and nested
length limits, rejects corrupt/truncated/trailing data without changing output,
and resumes the next coupled step exactly. The headless `--checkpoint-in`,
`--checkpoint-out`, and absolute `--checkpoint-every` workflow now routes this
accepted boundary through atomic same-file replacement. Its typed safe-point
control adapter also exposes only whole accepted macro-steps; the process-level
stdio regression checkpoints step two, verifies the exact third frame, and
resumes that continuation in a fresh worker.
The
`simwing-fsi --case porous-sheet` harness drives fluid through a translating
linear-resistance sheet, transfers only the sheet reaction to XPBD, and closes
the pump impulse/work, porous dissipation, and fluid/structure kinetic-energy
ledgers before publishing a frame. It uses that reusable selector to commit six
consecutive MAC-face rebases, including the wrapped `7 -> 0` epoch with signed
periodic image `1`, then stops before colliding with the next periodic image of
the prescribed pump surface. Its in-memory composite checkpoint restores the
complete accepted state with exact next-frame replay;
the checksummed persistent form reuses Structure's codec, stores the complete
topology epoch (version, axis, wrapped face, and signed periodic image) plus all
MAC fields, and validates them against bounded deterministic replay before
decode. Every pre-pump rebased epoch is persistable and continues
bit-identically. Checkpoint topology is validated at the accepted constitutive
midpoint, so the last accepted endpoint before pump collision can also be
persisted and deterministically repeats the same transactional rejection after
restore. A process-level restart at step 900 additionally resumes the ordinary
wrapped `face=0, image=1` epoch and advances only new accepted frames.
A focused reverse-direction instance of the same worker flips the pump and
material motion together, crosses the `0 -> 7` wrap into signed image `-1`,
persists and resumes that epoch, and rejects its next negative pump image. The
immutable direction has distinct provenance and checkpoint identity;
cross-direction restore is rejected. The CLI remains the positive-direction
canonical.
The
`simwing-fsi --case open-piston` harness drives a `6000 kg` plate at
`0.05 m/s`, exposes the resisting CFD pressure separately from its actuator,
and publishes accepted partial-cell and geometric-conservation ledgers. Before
that load reaches XPBD, a fluid-side planar cut-surface operator places each
face-resolved complete constraint reaction on the congruent physical plane and
closes area, force, moment, power, and periodic-image ledgers. The complete load
adds the direct MAC-velocity enforcement reaction to adjacent pressure traction.
Because projection produces a macro-step-average reaction, the same accepted
force is sampled at both endpoint kinematics before temporal integration. This
closes the worker's structure, fluid, actuator, and total momentum and kinetic-
energy ledgers; it remains a bounded reaction-geometry operator, not cell-
pressure interpolation. At
step 1200 it crosses its first `0.5 m` cell, rebases without a chamber-volume
jump, and continues in the new topology. At step 2400 it crosses the periodic
boundary while retaining an unwrapped `4 m` physical position. Its accepted
cut-surface pressure loads use the moving face-resolved bridge rather than the
uniform-only subset. A versioned in-memory fluid checkpoint binds accepted
pressure, velocity, projection diagnostics, grid geometry, and interface
topology behind an immutable payload. The open-piston worker composes it with
the full XPBD checkpoint, partial-cell epoch, and conservation ledgers; an
equivalent rebuilt worker resumes bit-for-bit, including immediately after
both ordinary and periodic topology rebases. The contained accepted
moving-interface fluid epoch now has a deterministic bounded/checksummed file
codec that preserves fields, interface topology/kinematics, and complete
projection diagnostics. Both structural-core halves now have deterministic,
bounded/checksummed codecs: the opaque SoftBody payload overlays committed
mutable state onto an equivalent rebuilt topology without trusting wire
topology, while the suspension payload preserves rigid-payload, line, control,
ledger, identity, and diagnostic state behind its complete-state fingerprint.
The enclosing Structure codec now combines those payloads with accepted
step/time and pending/applied nodal-load state, validating both directions
through an equivalent rebuilt Structure and rejecting public/opaque node-state
disagreement. The full composite open-piston checkpoint now also has a
deterministic bounded/checksummed codec. It nests those Structure and fluid
payloads with the partial-cell/rebase epoch and every committed transfer and
conservation ledger, validates through an equivalent rebuilt worker, and
recomputes the canonical diagnostic identities, geometry, acceptance state,
and cross-ledger force/work relationships before restore. It resumes ordinary
or freshly rebased arithmetic bit-for-bit. The open-piston
worker accepts the same `--checkpoint-in`, `--checkpoint-out`, and absolute
`--checkpoint-every` file workflow as periodic flow, including atomic same-file
resume. Its typed checkpoint adapter also supports the binary standard-I/O
control protocol without launching the viewer. A
real 3.28 fixture
also crosses export, structural assembly, one coupled pilot/suspension step,
checkpoint replay, and completed trace using synthetic physical settings.
An authoritative Qt-free surface bridge now compacts the scene's used fluid
regions and porous materials with stable-ID-sorted vertices, oriented
two-sided triangles, and ordered openings, then captures the corresponding
motion from one accepted Structure checkpoint. It deliberately does not infer
grid crossings or cut cells, and it rejects opening-only vertices until they
have a trusted structural motion rule. A companion binding carries that exact
oriented surface and accepted kinematics into the existing conservative
uniform or barycentric-quadrature traction transfer, retaining both scene and
Structure identities through validated XPBD load application. The first
grid-facing geometry stage now bins padded current-triangle AABBs into bounded,
cell-major candidates tied to the complete accepted surface-state fingerprint.
It covers both adjacent cells on internal grid planes and rejects geometry
outside the selected domain. A normalized separating-axis narrow phase now
removes triangle/AABB false positives, retains exact or explicitly tolerated
contact, and revalidates every expected pair. Exact zero-tolerance pairs are
now clipped against all six cell planes with original-triangle barycentric
coordinates, analytic area/centroid data, and explicit point/segment contact.
Every material and virtual-cap cell bound uses the canonical
`gridLower + faceIndex * spacing` coordinate, so independently clipped
adjacent cells cannot acquire different shared planes through non-associative
floating-point addition.
Coincident grid-plane area remains duplicated and flagged on both adjacent
cells in the raw clipping output. A separate ownership stage pairs those exact
duplicates into one canonical MAC-face area with authored side-region IDs and
explicit winding sign; ordinary areas remain cell-owned and zero-area contacts
remain diagnostic. Those unique owners now become stable barycentric
quadrature points which retain the authored material, sheet, role, and both
fluid regions and feed ordered traction through the conservative Structure
transfer without integrating shared-plane area twice. The reciprocal adapter
samples accepted Structure position and velocity at those identical stable
barycentric points. This boundary does not
invent traction. Explicit one-sided pressure samples can now produce exactly
one `(p_negative-p_positive)*normal` load through that same path; equal
pressures cancel and no polar or pressure stamp is added. Ordinary clipped polygons also expose their exact boundary
segments: independently clipped copies from adjacent cells must agree inside a
fixed machine-roundoff envelope and become one stable lower-cell transverse
face crossing with an in-face authored-side direction, while unpaired edges
remain contact and coplanar area remains separately face-owned. Unpaired
periodic-domain boundary area and grid-edge crossing ambiguity are rejected.
A bounded sparse index now groups every crossing and coplanar owner by stable
grid-bound MAC face while keeping multiple sheets separate; its summed length
and area are diagnostics, not union coverage. A face-local graph now stitches
adjacent triangle segments through stable authored vertex/edge provenance,
recomputes shared-edge intersections canonically, and keeps opening-boundary
and grid-edge endpoints explicit. Segments are now extracted as
winding-directed open chains or closed loops independently per authored region
pair. A valid multi-region junction may terminate several pair-specific chains
at one physical higher-degree node; branching within one region pair still
rejects. Open chains are
not falsely treated as closed regions. Simple closed loops now carry signed
area, centroid, and winding-derived enclosed/exterior region identity;
self-intersections and degenerate loops reject, while nested loops remain
separate until a bounded containment stage proves they do not touch, requires
authored region continuity through every parent, and closes exact per-region
area over an eligible MAC face. Simple directed open-chain arrangements whose
leaves all reach the rectangular face boundary also close exact per-region
areas while retaining every source chain. A bounded half-edge traversal
enumerates left faces from authored winding and rejects unstitched crossings or
conflicting labels. Opening-ended chains, coplanar sheets, and
boundary-touching loops remain explicitly unresolved. The versioned
`simwing_scene_fluid_grid_epoch` boundary now composes this entire chain plus
unique conservative quadrature into one immutable accepted-Structure remap.
Its fingerprint prevents candidates, patches, topology, partitions, or loads
from being mixed across structural steps; every stage keeps its own bounds and
the aggregate owned payload has an additional byte ceiling. A focused moving
regression crosses a fabric triangle through a MAC plane, preserves authored
region/material/sheet identity, and retains force/moment closure. A first
`simwing_scene_fluid_cell_volume` subset now reconstructs deterministic sparse
per-cell region volumes for closed, consistently wound region cycles.
Each oriented interface triangle forms a signed tetrahedron against the grid
origin; deterministic convex clipping distributes that chain across exact
cells, including cells wholly inside a region and cases whose face-local chains
remain open at tile boundaries. Every tetrahedron closes across its clipped
cells, every cell closes across its regions, and the global sum is compared with
an independent whole-surface divergence volume. Decomposition closure and
sparse-region publication have separate tolerances. Their defaults retain the
established worker arithmetic; offline refinement audits can retain smaller
resolved slivers with cell-local first moments and bounded centroid repair for
cancellation-scale complements. Scene-v2.2 may now attach an
oriented boundary-vertex disk to an intake or crossport. A separate immutable
opening-cap owner closes those explicit planar or nonplanar facets for volume
accounting only. When that disk is absent, planar convex loops retain their
exact fan and planar concave loops use deterministic reference-geometry ear
clipping, so triangle identity stays fixed under accepted motion. Winding comes
from the adjacent fabric boundary; the individual authored facet normals
survive, and no cap becomes fabric or traction. Every finite-area edge must
form one closed oriented region cycle, including valid three-region
sheet/cap junctions. A cap-free directed material loop is accepted only while
both its reference and live geometry remain collapsed. The analytic open-tetrahedron
regression therefore has a finite cell volume. Stable one-point triangle
samples now exactly integrate the accepted piecewise-linear cap velocity into
per-opening surface-sweep rates; a rigid material-plus-cap surface closes that
geometric ledger. The exact barycentric
triangle/box clipper is shared with material geometry and partitions each cap
into bounded positive-area grid patches. Off-face pieces have unique cell
owners; paired grid-plane copies become one canonical non-periodic face owner,
and their polygon, area, accepted velocity, and sweep remain explicit.
Cell-owned polygons now also expose exact positive-length cap intersections
with internal Cartesian faces. Independently clipped segments from both
adjacent cells must agree before one stable, winding-directed crossing is
published; face-owned aperture area remains area, and grid-edge ownership is
not guessed. A bounded planar half-edge arrangement combines these virtual-cap
crossings with the material chains on each touched face. It supports
disconnected signed cycles, publishes exact same-region areas when they close,
their global first moments and centroids, and otherwise retains an explicit
unresolved face record. The coarse real wing
closes all nine cap-crossed faces through this path. Exact coordinates already
owned by an earlier clip plane survive later-axis clipping, so a mathematical
face segment cannot disappear through one-ulp interpolation drift. Face-owned
aperture area remains with its existing owner. Every rejected face carries a
stable reason and, when one source is uniquely responsible, its material-chain
or opening-crossing stable ID. The pressure epoch retains the complete
product without turning caps into fabric. A
read-only flux epoch binds the complete MAC field, reads exact resolved face
flow or deterministic off-face staggered quadrature, and reports signed fluid
flow, cap sweep, and relative volume flow. It also maps every intake or
crossport into equal-and-opposite outward balances for its two authored
regions. Uniform co-moving air and mouth give zero relative flow, and all
region balances cancel globally. A bounded two-epoch audit now pairs those
balances with sparse region volumes and evaluates `delta volume + outward
relative flow` using endpoint trapezoidal integration. The expanding-cell
regression closes only with its matching inlet transport and reports a local
mismatch when that flow is removed. Authored openings now also define
deterministic pressure-connectivity components: an intake joins its cell to
Outside, while a crossport joins only its adjacent cells. Every component owns
one stable gauge ID and a component-wise continuity/source check. Thus two
sealed cells with equal-and-opposite volume changes are rejected locally even
though their global volume residual is zero. Applying these components to the
grid now starts with an immutable pressure-control-volume topology: every
positive sparse cell/region volume owns one pressure unknown, stable fixed-grid
cell/region ID, volume weight, exact centroid, component, and deterministic
gauge owner. The clipped-tetrahedron ledger also closes each cell's first
moment and recovers analytic region moments under translated grids. Mixed cut
cells and full interior cells close independently and exactly recover all
region and domain volumes. A complementary immutable face topology now links
those unknowns only through exact same-region Cartesian areas and their exact
subface centroids. Resolved nested
interfaces retain one link per region partition. Untouched faces use their one
unambiguous common region when available; when sparse supports overlap, the
accepted material-plus-cap closed surface may instead prove one region at the
face centroid by oriented solid angle, and that region must still own controls
in both adjacent cells. On a transversely cap-crossed face, the exact capped
material-plus-opening partition supersedes the material-only result; an
unsupported arrangement has its own unresolved status and no fabricated link.
The analytic open tetrahedron carries its `0.105 m²` Cell section and
`0.895 m²` Outside complement into pressure links, while the coarse real wing
now has ten resolved partition faces rather than one. A refined 4-by-4-by-4
real-wing audit closes every pressure-active Cartesian face: same-region
internal-sheet chains remain audited material topology but do not split a
pressure region, and stitched interior multi-region chains raise the result to
58 resolved partitions (15 material-only and 43 cap-touched). Its remaining
294 embedded-opening rejections are explicitly non-admissible centroid
stencils, not unresolved material faces. Every one has correctly sided
same-region Cartesian support within one ring on both authored sides after
periodic-image unwrapping. That support is retained as bounded immutable
provenance; it is evidence for a future multipoint reconstruction, not a
replacement conductance. A face-aligned authored opening now contributes its
exact oriented cross-region aperture area plus an unambiguous same-region link
over the rest of that Cartesian face. The analytic intake therefore retains
`0.18 m²` of opening and `0.82 m²` of Outside flow area instead of smearing the
whole face into either one. A cell-owned opening patch contributes an embedded
cross-region link inside its cut cell. Its authored normal and exact area are
retained, while conductance distance is the positive projected separation of
the two cell-region pressure centroids. When that separation is not positive,
the patch count and exact area remain explicitly unresolved and no conductance
link is fabricated. The rejection additionally owns deterministic ranges of
all same-region Cartesian neighbors reached from each side control, their
source link/control identities, periodic-image offsets from the aperture
centroid, signed normal projections, and correctly-sided flags.
Material/interface ambiguity remains unlinked, and the
pressure operator rejects either incomplete topology. Physical pressure
coefficients enter only in the downstream fixed-epoch projection. For fully
resolved geometry,
those links assemble into a bounded symmetric integrated graph Laplacian. Its exact
constant nullspace, positive link energy, component conservation, and retained
gauge owners are tested. Assembly requires one connected link graph per
authored pressure component. Face-aligned and tilted off-face intakes now join
Cell and Outside with exact aperture graph energy instead of acting like
sealed fabric. A deterministic component-wise conjugate-gradient layer now solves
manufactured integrated systems on both the nested closed regions and that open
intake. It rejects incompatible component sources, removes only admitted
roundoff, explicitly recomputes the final residual, sets every canonical gauge
to exactly zero, and leaves warm pressure unchanged on failure. A first physical fixed-epoch
adapter now maps the predicted MAC field into one oriented volume flow per
resolved link, substitutes accepted fluid-minus-cap-sweep flow on each exact
face-aligned or embedded aperture patch, assembles `-rho/dt` times the
integrated outward
flow as the pressure RHS, and applies the solved pressure difference back to
each link. It publishes pressure and corrected flow only after explicit
control-volume continuity closes; non-convergence exposes diagnostics without
a partial corrected state. The analytic intake retains separate aperture and
complement flows rather than collapsing its cut face back to one MAC value. A
consecutive-epoch owner now derives exact per-control-volume `dV/dt` from
stable cell/region identities. The moving projection overload targets
`dV/dt + net relative flow = 0`; with zero predicted air velocity, the
expanding analytic cell develops pressure and draws corrected flow inward
through its intake. A newly positive cell/region row is now marked explicitly
and receives a zero-volume previous endpoint. A disappeared row can retire its
complete previous volume to one unique retained same-region neighbour;
ambiguous or donorless retirement rejects. One bounded, versioned
pressure-topology transition now pairs every retained row and owns both
appearance donors and retirement recipients. Geometry-volume rates,
transported momentum, and pressure warm state consume that same fingerprinted
decision. General swept-volume remap, a direct Cartesian representation of
embedded opening flow, and general moving-boundary projection remain future
work. Embedded opening flow now continues through explicit-normal region
momentum; the derived bulk MAC field deliberately contains Cartesian faces
only.
Accepted projection pressure now also returns through the authoritative
surface path. Quadrature-v2 retains the exact cell owner of each authored
surface side; a bounded sampler resolves those cell/region pressure unknowns,
requires both sides to share one gauge component, and passes their difference
to the existing conservative pressure-traction transfer. The expanding open
tetra therefore produces a nonzero, force/moment-conservative XPBD load from
its scene-derived pressure. Uniform pressure-gauge warm starts give identical
samples and loads, while a sealed inside/outside pair with independent gauges
rejects rather than inventing an arbitrary pressure jump.
A versioned pressure epoch now composes the entire accepted geometry side of
that solve atomically: grid remap, opening caps, patches, transverse face
crossings, capped face partitions, sparse region
volumes, pressure controls, conservative links, and the ungauged operator all
share one Structure-state identity and an aggregate storage bound. Downstream
macro-step code can no longer accidentally mix those products across accepted
surface states.
The first strong feedback owner now closes that pressure path
against XPBD. Each Aitken iteration restores the same structural baseline,
applies the trapezoidal average of the accepted start pressure load and a
relaxed end-load guess, rebuilds the pressure epoch, projects moving-volume
flow, and returns the conservative scene pressure load. Only a converged
iterate commits. The expanding open-tetra regression proves pressure opposes
the imposed expansion, continues deterministically into the next macro-step,
and restores both Structure and pressure ownership on iteration exhaustion or
projection failure. The MAC predictor is held fixed inside each nonlinear
feedback loop. Accepted corrected link flows can now be conservatively
area-collapsed back onto one absolute bulk MAC velocity per Cartesian face,
including oriented intake-cap sweep. Embedded opening links reject that
Cartesian collapse until their oriented flow is carried by region momentum.
The visible pressure-cell canonical uses
that field as its next predictor. Region-resolved transport and a local
material-wall exchange now continue the accepted state after bootstrap;
general immersed-boundary advection/boundary layers and general topology
swept-volume remapping remain later work.
A separate topology-bound link-flow continuation primitive now proves exactly
what the collapse loses. It restores accepted opening-cap sweep, carries each
link's absolute velocity deviation from its previous face mean, recentres those
deviations under changed subface areas, and supplies an explicit per-link
pressure predictor while preserving the new bulk MAC face total to roundoff.
It is deliberately not enabled in the worker: a static carry has no
region-resolved convective update and strongly drains the repeated pressure
load. The regression records that failure mode so this bookkeeping adapter
cannot be mistaken for the missing cut-region momentum equation.
Accepted corrected links can now also be reconstructed into an immutable
cell/region momentum state. Each vector component is the area-weighted mean of
incident absolute link velocities; components without an incident face retain
the cell-centred value from the exact MAC predictor bound to the projection.
The state reports momentum, kinetic energy, fallback coverage, and the loss in
mapping staggered link normals to collocated region vectors. It is the explicit
input boundary for conservative region transport and does not itself advance
or apply a wall model.
A region transport kernel now advances that state with first-order
donor-cell vector-momentum flux over the corrected pressure links, followed by
equal-and-opposite graph-viscous impulses. It deterministically subcycles from
the smallest local outgoing-volume and viscous stability limits, conserves
internal momentum, requires non-increasing post-forcing stage energy, and
publishes no candidate on a failed ledger. Moving pressure volumes advance
linearly through their accepted geometry rate, which gives uniform flow an
exact discrete GCL update. An optional cell-centred delta between bound bulk
MAC predictors carries pump/advection/viscosity impulse and work into the
region state without erasing cut-region differences. A topology-stable
moving-epoch adapter then retains transported cell/region velocities while
remapping momentum to current physical volumes, averages them onto current
pressure links, subtracts exact opening-cap sweep, and supplies the
fingerprinted predictor to pressure projection alongside dV/dt.
A bounded first crossing adapter now handles one-ring appearance and
retirement through that shared topology transition. Retained controls keep
transported velocity; each newly positive control receives the transition's
area-weighted retained same-region velocity. Momentum is recomputed from
current volume with an explicit geometric-change ledger. A
disappeared control transfers its complete source volume and momentum to one
unique previous same-region neighbour, and that source mapping closes before
geometric correction. The pressure warm start seeds appearances, preserves
retained rows, and drops retired values. The strong
coupling owner composes this rebase with wall exchange and pressure projection
transactionally. Cross-material donation and missing or ambiguous one-ring
ownership still reject; this is not a general swept-volume remap.
A separate material-wall operator remaps that fixed transport to each current
strong-coupling geometry, applies two-sided tangential viscous exchange at the
authoritative material quadrature, and returns both adjusted region momentum
and equal-and-opposite Structure traction. A local `0.5*volume/wall-area`
distance closure and deterministic explicit subcycling bound the exchange;
fluid/Structure impulse, wall work, and nonnegative dissipation are checked
before publication. Pressure remains the sole owner of normal traction. The
wall-adjusted current-link predictor feeds pressure projection, while the
combined pressure-plus-shear load follows the existing conservative XPBD
transfer. This is a local coarse cut-region closure, not a resolved boundary
layer.
Its in-memory composite checkpoint retains Structure plus the accepted sparse
pressure projection; restore rebuilds and validates the complete pressure
epoch and conservative load before committing either owner. Initial and
post-step checkpoints reproduce the exact next coupled result in the same or
an equivalent owner, while foreign settings and corrupt or missing accepted
pressure reject transactionally.
The selectable `pressure-cell` worker now owns this path after its first
bootstrap solve. A soft
three-panel tetrahedral cell has one triangular intake and three fixed mouth
vertices. A prescribed `-0.85 m/s` periodic mean wind replaces the former
mechanical apex actuator: a uniform mean-flow correction pumps the accepted
bulk MAC field, the existing symmetric SSPRK2 viscosity/projected-nonlinear-
advection/viscosity operator advances it, and the scene pressure solve supplies
the only load on the free apex. The bulk velocity is advanced on a private
candidate and commits only after scene feedback accepts. Its
accepted region momentum receives the bound bulk-MAC delta, advances through
moving-volume GCL transport, and predicts every current pressure link during
each strong geometry iterate. Its
immutable frames publish deformation, area-averaged triangle pressure jump,
separate pressure/wall/total-fluid nodal and global forces, mean-flow pump
force, bulk-flow change and divergence, viscous energy loss, region transport
loss/GCL change/momentum residual, wall loss/momentum residual,
strong-iteration count, bulk MAC speed, and
the maximum
mixed-subface velocity spread discarded by the area collapse, plus the count
of embedded opening links retained outside the Cartesian bulk field. The final
scene projection remains cut-region-aware, but the two intermediate bulk-advection
projections see only one velocity per Cartesian face; this is a bootstrap
diagnostic, not a general immersed-boundary CFD or wing-aerodynamics claim. A
600-step headless run remains topology-stable for 10 simulated seconds and
reports a relaxed `0.00309 Pa`, `0.0260 mm` deformation, and `0.908 m/s`
maximum MAC speed; the two-second transient reaches about `2.12 Pa` and
`17.8 mm`. The molecular-viscosity wall reaction is conservative but tiny at
this coarse resolution (about `1.4e-6 N` at two seconds), so sustaining a
bluff-body pressure wake still requires general immersed-boundary advection
and a validated boundary-layer closure rather than hidden predictor reset.
An explicit `--mimetic-pressure-audit` mode now runs the mixed-hybrid pressure
path only after the normal graph iteration converges. It bootstraps from the
same MAC/opening predictor, then uses transported wall-adjusted flow and the
consecutive pressure warm remap. The accepted shadow retains its complete
topology, physical source, solve diagnostics, pressure state, and material
samples but never applies those samples to Structure; default frames and graph
checkpoints remain byte-identical. The same endpoint crosses the coarse real
wing's 138 controls and 42,927 shared traces without requiring the old graph
operator to accept inadmissible embedded-opening coefficients. The bounded
`SWPCELL10` checkpoint now composes the compact accepted `SWMP` rows with the
graph restart. Restore rebuilds trusted control/full/condensed topology from
the Structure payload before decoding those rows, then resumes the exact
consecutive wall-predicted endpoint without persisting transient topology.
Before any load switch, the same opt-in path now compares every graph and
mimetic material pressure jump and evaluates both through the unchanged
conservative Structure transfer. It retains bounded per-sample and per-node
deltas and reports pressure, force, moment, and power disagreement without
applying the shadow load. The diagnostic is decisive on the current pressure
cell: relative pressure-jump and net-force differences are about `0.60` both
at step 4 and after 600 steps. Mimetic load selection therefore remains
disabled rather than treating solver acceptance as physical agreement.
That evidence now feeds a separate immutable pressure-owner decision. Its
explicit default policy requires the graph and mimetic source vectors to agree,
then bounds pressure magnitude and shape, nodal-force magnitude and shape, net
force, moment, power, and conservative-transfer closure. Every failed check has
a typed rejection bit and selects the graph owner; only a zero-rejection result
names the mimetic candidate. The real-wing exact self-comparison proves the
positive path when source evidence is explicitly waived, while the live cell
fails pressure magnitude/scale, nodal-force scale, and net-force checks. The
decision is diagnostic only: no worker uses it to apply mimetic loads, and
agreement with the current graph path is not presented as a continuum-truth
criterion. The opt-in coupling now publishes this decision atomically beside
the accepted comparison. Its CLI line reports `owner=graph` with rejection
mask `0xeb00` at step 4 (pressure magnitude/scale, nodal-force scale, net
force, moment, and power) and `0x6b00` after 600 steps, when the power-delta
check has cleared but the other five blockers remain. Even a deliberately
permissive policy that selects the mimetic candidate leaves the production
frame byte-identical; the selection is not yet a load-application command.
The mismatch has now been narrowed further. The graph and mimetic source
vectors agree to relative roundoff (`5.3e-16` at step 4 and `1.7e-16` after
600 steps), while the shadow pressure and every transferred nodal load are
almost exact scalar multiples of the graph result: gains `2.53035` and
`2.50693`, with post-fit relative shape residuals near `2e-16`. That locates
the open question in the cut-cell graph versus mixed-hybrid operator response,
not pressure-source units, gauge-safe surface sampling, or load transfer. It
does not by itself establish which spatial operator is the better physical
reference.
An offline inverse-response audit now probes that operator question with the
accepted source and six independent manufactured modes on the same 65-control
cut-cell topology. The accepted direction remains exceptional: its
gauge-aligned control pressures have about `2.54` shadow gain and `2.4%`
post-fit shape residual. The coordinate/mixed/high-frequency probes instead
have gains from `0.999` to `1.007`; their shape residuals range from `1.2%` to
`16.5%`. A piecewise-constant inside/outside mode has no source on same-region
graph links—only on the authored intake—and it reproduces the outlier with
`2.562` gain and `1.98%` residual. A gauge-invariant source-work audit makes
that comparison dimensional: the graph gives the intake a `0.18 m`
two-terminal conductance, while the mixed-hybrid cut-cell response gives
`0.0700821 m`, a `2.5684` ratio. This rules out a global pressure-unit or
normalization repair and localizes the transition question to authored-opening
cut-cell closure rather than bulk response. Neither value is yet treated as
the physical oracle.
A separate rest-geometry refinement audit now removes the visible case's
face-aligned special position. It uses a fixed skew tetrahedral intake and
rebuilds the complete graph and mixed-hybrid topologies at `2^3`, `4^3`, and
`8^3`. The graph/shadow conductance ratios are `28.024`, `33.672`, and `8.830`.
After multiplying conductance by nominal cell width and dividing by intake
area, the shadow sequence rises smoothly from `0.102` through `0.358` to
`0.503`; the graph sequence jumps from `2.857` to `12.063` and back to `4.446`.
This demonstrates strong grid sensitivity in the embedded two-point graph
coefficient. It is not a convergence proof for the shadow method, so live
loads remain on the established graph path.
A fixed-`4^3` grid-placement audit reaches the same conclusion independently.
It moves the grid through all eight combinations of zero and negative
half-cell phase while keeping the skew tetrahedron fixed. Six phases assemble
complete operators; two are rejected with typed diagnostics because one and
two embedded intake patches, respectively, have no complete pressure-face
ownership. Over the six accepted placements, normalized graph conductance
ranges from `2.5575` to `13.9854` with population coefficient of variation
`0.77175`; the shadow response ranges from `0.13695` to `0.40878` with
coefficient `0.34030`. The shadow result is less placement-sensitive in this
small coarse ensemble, but it remains an offline comparator rather than the
live pressure owner or a demonstrated continuum oracle.
Repeating that same eight-phase set at `2^3`, `4^3`, and `8^3` exposes a more
fundamental limit: only `4/8`, `6/8`, and `2/8` placements, respectively,
assemble complete paired products. On those surviving subsets the normalized
graph means are `3.9025`, `6.2525`, and `4.5187`, while shadow means rise from
`0.10168` through `0.26909` to `0.45753`. The shadow sequence is suggestive,
but not a convergence claim—the graph ownership gate censors six of eight
fine-grid phases before a paired response can be recorded. The next offline
step is therefore a graph-independent shadow phase spectrum, not a live solver
switch or a repair that silently seals unresolved intake patches. Its generic
terminal primitive is now in place: it finds permeable cross-region traces in
the trusted mixed-hybrid topology, represents a face-aligned authored opening
as a Cartesian trace and an embedded opening as an authored-opening trace,
excludes material walls, and distributes a balanced fixed transfer uniformly
by aperture area. The direct condensed solve reports conductance from
gauge-invariant source work. It gives `0.0700820848335194 m` on the
face-aligned case, within `3.6e-14 m` of the earlier graph-manufactured shadow
measurement, and `0.0608388978079475 m` on the embedded two-point fixture.
The resulting graph-independent matrix now retains and solves all 40 requested
phase/refinement attempts. The normalized mean/CV values at `2^3`, `4^3`, and
`8^3` are `0.100660/0.01043`, `0.240930/0.39050`, and
`0.521248/0.38828`; at `16^3` they are `0.902570/0.20104`, with range
`0.624544`-`1.141234`; at `32^3` they are `1.091977/0.07436`, with range
`0.956232`-`1.207995`. The five former `8^3`
fine-grid failures were not ill-conditioned inversions. Tiny Cartesian
subfaces were measured with absolute-coordinate shoelace moments, and
triangle/grid-edge vertices inherited whichever sequential clipping path
created them. Face-local moment charts plus canonical authored-triangle/two-
grid-plane nodes remove that precision loss without fitting closure or
relaxing the unchanged `1e-10` algebraic tolerance. At `16^3`, tiny
full-minus-positive boundary complements additionally take their centroid
from an independently integrated reverse polygon when the subtraction error
exceeds a fixed coordinate-ULP envelope. The resulting four-face/two-wall
sliver closes its divergence moment to `2.14e-21 m^3`; if its normal `7 x 7`
wall auxiliary core is numerically singular, a bounded direct principal-block
fallback retains the exact Schur metric without persistent dense storage. All
eight phases then solve at `16^3`. Extending to `32^3` exposes a real
`2.83e-15 m^3` complementary region that the former shared `1e-12 m^3`
closure/publication tolerance erased. The audit retains it with an independent
publication envelope and cell-local moment path. Two valid four-face local
operators report algebraic errors `5.52e-10` and `1.31e-10`, so this offline
matrix explicitly fingerprints a `1e-9` algebraic-consistency tolerance;
production retains the default `1e-10`. Phase CV contracts by about 48% from
`8^3` to `16^3` and another 63% from `16^3` to `32^3`, while mean drift also
contracts. These two intervals are useful continuum evidence, not a convergence
claim. A separate immutable three-level assessment now makes that caveat
machine-checkable. The aggregate `8^3`/`16^3`/`32^3` means have apparent order
`1.0095`, extrapolate to `1.27891`, and leave a `14.6%` fine-grid gap, all
inside its screening policy. The per-phase trajectories are not asymptotic:
four reverse increment direction, four fail the `0.75` contraction bound, and
only one of eight passes both. The strict outcome is therefore
`InsufficientEvidence` with rejection mask `0x300`; even an explicitly
permissive zero-rejection result names only a read-only trend candidate and
cannot select or apply a pressure field. The live 600-step trace
and audited checkpoint remain byte-identical. The complete phase matrix still
does not establish shadow convergence or authorize a live solver switch; the
uniform area-weighted multi-opening source also differs intentionally from the
earlier graph-manufactured source.
A separate origin-invariance oracle moves that same `8^3`, phase
`[0,-0.5,0]` sample and its grid together by `[256,-512,1024] m`. All 524
controls, full/reduced traces, and four opening traces retain their topology;
the intake area changes by only `5.93e-14 m^2` and normalized conductance by
`5.59e-12`. Adjacent-cell crossings now compare geometric position plus
barycentric zero/nonzero provenance, canonical face nodes admit only a fixed
coordinate-ULP envelope, and capped arrangements switch to a local chart when
that envelope exceeds their declared minimum tolerance. The ordinary-origin
path and both production hashes remain byte-identical.
Its checkpoint also stores the trusted Structure state, complete
accepted sparse pressure projection, accepted wall-traction endpoint, and
accepted region momentum. Initial and
accepted files
round-trip deterministically, reject foreign/corrupt input transactionally,
resume through the normal worker checkpoint flags, and reconstruct the exact
derived MAC continuation without duplicating either it or transient bulk
pressure on the wire. Version 10 preserves transported-region and wall-exchange
projection provenance, including explicit-normal embedded-opening momentum,
and bounds/revalidates every momentum control volume, material quadrature
traction, and optional nested `SWMP` pressure row before publication.
Nonplanar openings without authored cap triangles, self-intersecting or folded
caps, opening-only interior construction vertices, branching or inconsistently
wound surface junctions, periodic-boundary ambiguity, and general
moving-boundary fluid equations still
reject or remain open; the grid epoch itself continues to own geometry and
transfer only. The analytical model exporter now supplies deterministic
boundary-vertex cap disks for its nonplanar intakes. It reuses the actual skin
lip and rib-mesh boundary vertices, so every exported opening vertex reaches
the live Structure-to-fluid surface state. The cap owner now validates the
full real wing's three-region skin/rib cycles and consumes its complete opening
set. Pairwise signed cell-volume accounting also closes the complete real-wing
region ledger on a centered coarse grid that crosses the canopy. One simple
boundary-to-boundary interface now resolves to exact oriented face areas;
junctions survive as pair-specific open chains and explicit unresolved
partitions. Opening quadrature/grid patches preserve the full cap area, and
authored connectivity plus sparse pressure controls assemble. Pressure-link
construction now reaches that real geometry and retains every coarse embedded
opening without admissible projected cell-region centroid separation as a
typed rejection with patch/opening identity, signed centroid-to-patch
distances, exact count and area, and no conductance link. The current coarse
fixture has 24 non-positive patches across its two mirrored intake openings;
both region centroids lie on the positive side of each local cap plane, so an
absolute-value coefficient would hide a non-admissible two-point stencil.
All 24 publish correctly sided one-ring support on both authored sides,
including validated periodic-image geometry. Direct donor substitution would
place flux in the wrong control-volume row. A new isolated mixed-hybrid
mimetic local-cell kernel now supplies that formulation's first numerical
building block without changing the production graph solve. From exact volume,
cell centroid, and oriented half-face area/centroid/normal geometry it builds a
fingerprinted SPD inverse flux inner product satisfying linear consistency.
The exact diagonal-plus-two-rank-three factorization uses seven stored doubles
per half-face and applies matrix-free; it no longer materializes an `n x n`
cell matrix. It eliminates the cell scalar with roundoff-only integrated
conservation, reproduces linear normal fluxes on a skew tetrahedron, and
reduces exactly to the existing `area / centre-distance` stencil after
Cartesian trace condensation. A companion exact local Schur kernel now removes
prescribed-zero-flux material-wall trace equations without a dense wall
matrix. It exploits the cell trace operator's diagonal-plus-seven-low-rank
form, retains one equilibrated `7 x 7` Woodbury inverse and linear face data,
and matches independent dense action, diagonal, right-hand-side, and wall-trace
reconstruction oracles while preserving the active constant null mode. An
immutable scene-side audit now assembles those periodic-image-unwrapped
half-face shells from Cartesian traces, two-sided material-wall quadrature, and
cell-owned opening patches. Nested, face-aligned-opening, and deliberately
two-point-rejected embedded-opening fixtures close every area vector and
`N^T R = volume I`; every completed shell builds the generic SPD kernel.
The same audit path now assembles one bounded, fingerprinted global
mixed-hybrid trace system. Shared Cartesian and authored-opening half-faces
pair into two-incidence traces, every impermeable material side owns a
one-incidence zero-flux wall trace, and one deterministic gauge is retained
per pressure component after rejecting disconnected component topology. Its
matrix-free condensed action is symmetric
positive semidefinite, preserves one exact constant null mode per component,
and builds a source-compatible right-hand side without materializing either a
global matrix or dense local matrices.
An immutable global wall-condensed adapter now composes the exact local Schur
kernels into a separate shared-trace system while retaining the full topology
for reconstruction. Its matrix-free action is symmetric positive
semidefinite, preserves the same component constant modes, publishes exact
assembled diagonals and deterministic shared gauges, condenses a full trace
right-hand side, and reconstructs wall values that close every original row.
An isolated gauge-fixed solver now pins the retained trace in every component
and applies deterministic Jacobi-preconditioned conjugate gradients directly
to that matrix-free action. It admits and removes only bounded component-sum
roundoff, uses the exact condensed row diagonals, freshly recomputes the full
residual before publication, and leaves its warm start unchanged on
incompatibility, non-finite arithmetic, or iteration exhaustion. Manufactured
multi-component fields recover their gauge-normalized traces, and
source-driven cases close both local conservation and every shared/wall trace.

The coarse real-wing audit now closes all 138 sparse cell/region shells. Ten
untouched Cartesian faces whose adjacent cells share many authored region IDs
are proven Outside by the accepted material-plus-cap closed surface; no
dominant-volume or closure-derived label is used. All controls are topology
complete and build-ready. All 42,826 cell-owned opening patches still appear
as 85,652 exactly paired half-faces, so the 24 non-admissible two-point
rejections are geometrically complete rather than missing aperture geometry.
A manual `4 x 4 x 4` audit likewise resolves its six formerly ambiguous faces
and closes all 358 controls. It safely omits 240 material quadrature sides
whose corresponding cell/region volume is exactly absent, misses no opening
side, retains 95,984 paired opening half-faces, and reaches 2,947 total
half-faces in its largest control. Every coarse real-wing shell now also builds
that compact local factorization within the linear storage bound. The complete
coarse system has 191,579 trace unknowns and retains 13,132,336 bytes of compact
local factors; its full matrix-free component-constant action is roundoff-null.
The real-wing regression also executes one finite, residual-reducing PCG step
over all 191,579 rows and proves truncated non-publication. All 138 local
controls now accept exact wall condensation: 148,652 one-sided wall traces
reduce to a separately assembled 42,927-row shared system using 3,986,602 bytes
of linear Schur data, and its assembled diagonals are positive. Solving that
reduced system, any further local/multilevel preconditioning, physical
right-hand-side assembly, and production integration remain open. The
production graph operator and worker are unchanged.
Export still requires explicit physical material/pilot settings;
manufacturing-pattern UVs, exact authored attachment vertices, structural seam
assembly, curved or transversely deforming grid-to-surface correspondence,
general cut-cell pressure metrics, nonplanar topology events, multiple crossings
in moving/cut-cell topology, and AMR CFD remain open. These worker cases
validate the pipeline;
they are not yet wing CFD or aerodynamic truth.
The inherited Playground is not used by these targets.

## Inherited LEparagliding Studio baseline

This project is a C++ port of Pere Casellas' LEparagliding 3.28 “Jardins”
engineering program with a Qt 6 desktop interface. It is free software under
the GNU GPL 3.0, like the
[original program](https://www.laboratoridenvol.com/leparagliding/lep.en.html)
it derives from.

The application accepts the same design text file as the original program.
Airfoil file references remain relative to the selected design file. Generated
files can be written to a separate folder:

- `leparagliding.dxf` — 2D manufacturing plans
- `lep-3d.step` — exact OCCT NURBS wing model
- `lep-3d.dxf` — legacy 3D wireframe retained as a reference
- `lep-out.txt` — calculated design data
- `lines.txt` — suspension line data
- `run-log.txt` — calculation progress and diagnostics

Version 3.28 can also create `stl` and `xflr5` subdirectories when the
corresponding design options are enabled.

The desktop application is a complete design studio:

- every numbered block in the selected design has its own syntax-highlighted
  editor and independent Undo/Redo history;
- every Save embeds the wing's complete version history in the design file, so
  Undo/Redo can continue across restarts and the `Versions...` window
  can restore the whole wing to any saved state;
- the `?` button on a section opens format guidance and a link to the full
  manual;
- opening a design automatically calculates a fresh OCCT NURBS model in a
  temporary folder, preventing an older exported model from being shown;
- **Build paraglider** validates the current editors and refreshes that
  temporary preview without saving the design or writing user output files;
- **Export files...** writes the manufacturing plans, 3D geometry, reports,
  and line data to the selected Output folder;
- the viewport reads `lep-3d.step` with OCCT, triangulates its exact NURBS
  surfaces with OCCT, and renders them with the native OCCT OpenGL viewer;
  isometric, front, back, left, right, top, and bottom views are available in
  perspective or orthographic projection;
- **Preferences…** adjusts the viewport triangulation resolution from very
  coarse to ultra fine; the setting applies immediately, is remembered between
  sessions, and never affects exported files.

Viewport navigation follows the slicer convention: drag with the left mouse
button to orbit, drag with the right or middle mouse button to pan, and use the
wheel to zoom. `Shift` + left drag also pans; double-click or press `F` to fit
the model. Number keys `0`–`6` select the preset views and `P` toggles the
projection. Double-click a completed output file to open it in its associated
viewer.

While editing a section, press `Enter` to build and refresh the temporary 3D
preview. Use `Shift+Enter` when you intentionally need to insert another input
record. Preview builds and exports use the current editors but do not save them;
use **Save** when the change should become a new embedded version.
`Ctrl+Z`/`Ctrl+Y` and the section's Undo/Redo buttons affect only the currently
visible section; switching sections does not merge or clear their histories.
Once the live editor history is exhausted, Undo continues through that
section's saved versions. Restoring an older whole-wing version does not delete
newer versions; the restored state becomes a new latest version when it is next
saved.

Version history is stored as a marked comment trailer at the end of the same
design file. Each entry is a compressed full-wing snapshot with a UTC
timestamp, changed-section list, parent identifier, and SHA-256 identifier.
The calculation engine removes this trailer in a temporary input copy before
calling the strict Fortran-compatible parser; the editable design file and its
history remain intact. A design without embedded history is treated as version
1, preserving the wing exactly as it was first opened.

The original
[LEparagliding user manual](https://www.laboratoridenvol.com/leparagliding/manual.en.html)
documents every input section. In particular, its record order is strict and
blank lines are not valid records.

## Build on Windows

The CMake preset mirrors the compiler, Qt, and Open CASCADE setup used by
`C:\CODE\cobod-slicer`: Visual Studio 2022, the newest compatible Qt 6 MSVC kit
under `C:\Qt`, and OCCT under `C:\OpenCASCADE-8.0\build2`. CMake auto-detects
that OCCT installation. A different compatible build can be selected with
`-DLEP_OCCT_ROOT=C:\path\to\occt`.

```powershell
cmake --preset windows
cmake --build --preset release --parallel
ctest --preset release
```

Run:

```powershell
.\build\bin\Release\LEparagliding.exe
.\build\bin\Release\simwing-fsi.exe --case hemisphere --steps 600
.\build\bin\Release\simwing-fsi.exe --case flag --steps 600
.\build\bin\Release\simwing-fsi.exe --case ram-cell --steps 600
.\build\bin\Release\simwing-fsi.exe --case strong-piston --steps 120 --no-viewer
.\build\bin\Release\simwing-fsi.exe --case strong-piston --steps 120 --no-viewer --checkpoint-out strong-piston.swsp --checkpoint-every 60
.\build\bin\Release\simwing-fsi.exe --case strong-piston --checkpoint-in strong-piston.swsp --steps 60 --no-viewer --checkpoint-out strong-piston.swsp --checkpoint-every 60
.\build\bin\Release\simwing-fsi.exe --case periodic-flow --steps 600
.\build\bin\Release\simwing-fsi.exe --case periodic-flow --steps 600 --no-viewer
.\build\bin\Release\simwing-fsi.exe --case periodic-flow --steps 600 --no-viewer --checkpoint-out periodic-flow.swpc --checkpoint-every 60
.\build\bin\Release\simwing-fsi.exe --case periodic-flow --checkpoint-in periodic-flow.swpc --steps 600 --checkpoint-out periodic-flow.swpc --checkpoint-every 60
.\build\bin\Release\simwing-fsi.exe --case open-piston --steps 1200 --no-viewer --checkpoint-out open-piston.swop --checkpoint-every 600
.\build\bin\Release\simwing-fsi.exe --case open-piston --checkpoint-in open-piston.swop --steps 600 --checkpoint-out open-piston.swop --checkpoint-every 600
.\build\bin\Release\simwing-fsi.exe --case pressure-jump --steps 4 --no-viewer
.\build\bin\Release\simwing-fsi.exe --case pressure-cell --steps 600
.\build\bin\Release\simwing-fsi.exe --case pressure-cell --mimetic-pressure-audit --steps 600 --no-viewer --checkpoint-out pressure-cell-audit.swpcell --checkpoint-every 600
.\build\bin\Release\simwing-fsi.exe --case pressure-cell --mimetic-pressure-audit --steps 300 --no-viewer --checkpoint-in pressure-cell-audit.swpcell --checkpoint-out pressure-cell-audit.swpcell --checkpoint-every 300
.\build\bin\Release\simwing-fsi.exe --case pressure-cell --steps 600 --no-viewer --checkpoint-out pressure-cell.swpcell --checkpoint-every 300
.\build\bin\Release\simwing-fsi.exe --case pressure-cell --checkpoint-in pressure-cell.swpcell --steps 600
.\build\bin\Release\simwing-fsi.exe --case porous-flow --steps 120 --no-viewer
.\build\bin\Release\simwing-fsi.exe --case moving-porous-flow --steps 101 --no-viewer
.\build\bin\Release\simwing-fsi.exe --case moving-porous-flow --steps 101 --no-viewer --checkpoint-out moving-porous-flow.swmf --checkpoint-every 50
.\build\bin\Release\simwing-fsi.exe --case moving-porous-flow --checkpoint-in moving-porous-flow.swmf --steps 4 --no-viewer --checkpoint-out moving-porous-flow.swmf --checkpoint-every 2
.\build\bin\Release\simwing-fsi.exe --case porous-sheet --steps 120 --no-viewer
.\build\bin\Release\simwing-fsi.exe --case porous-sheet --steps 330 --no-viewer --checkpoint-out porous-sheet.swps --checkpoint-every 165
.\build\bin\Release\simwing-fsi.exe --case porous-sheet --checkpoint-in porous-sheet.swps --steps 30 --checkpoint-out porous-sheet.swps --checkpoint-every 30
```

The hemisphere command launches the standalone trace viewer by default and
shows a soft fabric dome held at three equatorial points with a compliant rim.
A spatially alternating analytic pressure mode deforms four lobes without
exciting a free rigid-body rotation. It is a structural canonical, not a CFD
result. The following commands show
the other canonical cases. The flag command is the first real CFD-load-driven
fabric animation while retaining a stationary reference surface on the fluid
side. Commands with `--no-viewer` run unthrottled for
tests and scripted verification.
The ram-cell command shows the same conservative complete-reaction path acting
on a multi-panel open fabric shell; it is a deformation/inflation precursor,
not a resolved moving-cavity or aerodynamic-wing result.
The checkpoint commands save and resume exact accepted worker state. Resumed
steps are additional, autosave cadence uses absolute accepted-step multiples,
and input/output may name the same atomically replaced file.
For periodic flow, moving porous flow, open piston, or porous sheet,
programmatic controllers can
replace
`--steps` with `--control-stdio`; that mode requires binary protocol input and
an explicit stop command, so it is not an interactive terminal interface.

`windeployqt` and the OCCT runtime deployment run after the build, so the build
output is directly runnable. The install target includes the required OCCT
toolkit and third-party DLLs as well.

The calculation engine can also be used without the GUI:

```powershell
.\build\bin\Release\leparagliding-engine.exe <design-file> <output-directory>
```

The main Qt executable exposes the same operation in headless mode:

```powershell
.\build\bin\Release\LEparagliding.exe --headless <design-file> <output-directory>
```

Both commands return the engine's exit code and generate `leparagliding.dxf`,
`lep-3d.step`, the reference `lep-3d.dxf`, `lep-out.txt`, `lines.txt`, and
`run-log.txt` in the selected output directory. Relative airfoil paths are
resolved from the design file's directory.

The 3.28 input format adds sections 33–37 for detailed risers, line
characteristics, equilibrium calculations, XFLR5 export, and special
parameters. Older 3.17 designs that end at section 32 remain usable: the
command-line boundary appends disabled defaults for the five new sections to a
temporary input file. It never rewrites the selected design.

## Build on Linux and macOS

The same CMake project builds with GCC or Clang. Two dependencies are
required:

- **Qt 6.5+** (Widgets) — distro packages, the Qt online installer, or
  [aqtinstall](https://github.com/miurahr/aqtinstall).
- **Open CASCADE 8.0** — newer than any distro or Homebrew package today, so
  build it from source once and install it to a prefix:

```sh
curl -fsSL https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V8_0_0_p1.tar.gz | tar xz
cmake -S OCCT-* -B occt-build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DINSTALL_DIR="$HOME/occt-install" \
    -DCMAKE_INSTALL_RPATH='$ORIGIN;@loader_path' \
    -DBUILD_MODULE_Draw=OFF -DBUILD_MODULE_DETools=OFF \
    -DUSE_TK=OFF -DUSE_TCL=OFF -DBUILD_DOC_Overview=OFF
cmake --build occt-build --target install
```

On Debian/Ubuntu the OCCT build needs
`ninja-build libfreetype-dev libfontconfig1-dev libx11-dev libxext-dev
libxi-dev libgl1-mesa-dev libglu1-mesa-dev libxkbcommon-dev`; on macOS
`brew install ninja freetype fontconfig`.

Then configure the project with both prefixes visible:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$HOME/occt-install;/path/to/Qt/6.x.y/gcc_64"
cmake --build build
ctest --test-dir build --output-on-failure
```

The release builds in `.github/workflows/ci.yml` run exactly this on
Ubuntu 22.04 (so AppImages run on any distro with glibc 2.35+) and macOS
15, caching the OCCT install tree between runs. A
third job builds with Ninja/MSVC on Windows against the official OCCT
prebuilt vc14-64 binaries (with their third-party runtime DLLs merged in)
and runs the Fortran byte-compare reference test.

## Releasing

The application version has a single source of truth: the `project(...
VERSION x.y.z)` line at the top of `CMakeLists.txt` (the ported calculation
core remains LEparagliding 3.28 regardless). CI reads it and stamps every
artifact name with it. To cut a release:

1. Bump the version in `CMakeLists.txt`, commit, and push.
2. `git tag vX.Y.Z && git push origin vX.Y.Z`

The tag build packages all three platforms — Windows zip (self-contained
Release folder), Linux AppImage, macOS DMG — creates the GitHub release,
and attaches the downloads. CI runs only on `v*` tags; for a build
without a release, use the Actions tab ("Run workflow", optionally with
"Upload artifacts to the GitHub release").

The macOS app is signed with the Developer ID certificate and notarized
when the repository's signing secrets are configured
(`MACOS_DEVELOPER_ID_P12_BASE64` and friends); without them CI produces an
unsigned DMG (right-click → Open on first launch). The signed Windows
installer remains the local `installer/build_installer.ps1` flow described in
`docs/legacy/leparagliding/CLAUDE.md`.

## Port architecture

- `src/legacy/leparagliding_core.cpp` is the mechanically translated numerical
  and drawing core. It is built as C++, with Fortran indexing and I/O behavior
  retained for compatibility.
- `src/engine` supplies a small typed C++ boundary, input/output path handling,
  validation, and a command-line entry point.
- `src/model` interprets the fully transformed airfoil stations and analytical
  circular ballooning law directly. It converts each span arc to an exact
  rational B-spline, lofts semantic upper/vent/lower panel faces with OCCT,
  mirrors the calculated half wing, sews matching faces into shared shell
  topology while preserving designed intake openings, and writes an AP242
  STEP model in millimetres. The old tessellation is used only as a numerical
  regression oracle and remains available to the legacy DXF/STL exporters.
- `src/gui` is the Qt Widgets application. It runs the engine in a child
  process so the interface stays responsive and legacy input failures are
  isolated. Its viewport uses OCCT for STEP import, triangulation, and OpenGL
  presentation; it contains no application-side polygon model builder.
- `third_party/libf2c` is the portable runtime required by the translated I/O
  statements. Its original notice is included in that directory.

The translation boundary handles the Fortran features that `f2c` cannot
translate directly:

1. the Fortran 90 array-based word count was expressed as a character scan;
2. two whole-array negations were expanded to loops;
3. whole-array assignments used by the equilibrium solver were expanded to
   loops;
4. dynamic XFLR5/STL paths and directory creation were routed through the C++
   output boundary;
5. `kini=1` in `datair` was made explicit as `kini(i)=1`, matching how the
   array is subsequently consumed.

The C++ compatibility boundary also supplies reliable Windows `BACKSPACE`
record handling for LF and CRLF inputs and preserves GNU Fortran's formatted
output conventions. The regression fixture was generated by compiling
`leparagliding3.28.f` with native Windows `gfortran`: both DXFs, `lines.txt`,
and `run-log.txt` match byte-for-byte. The calculation report is compared
field-by-field with a 0.00015 display tolerance to accommodate signed zero and
two last-decimal rounding differences.

The active reference implementation is `leparagliding3.28.f` at the repository
root; the previous source is retained as `leparagliding3.17.f`. The original
source identifies itself as GNU GPL 3.0 software; the translated core is a
derivative under the same terms.
