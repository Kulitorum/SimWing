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
and grid-edge endpoints explicit. Degree-two components are now extracted as
winding-directed open chains or closed loops only when their authored region
pair stays consistent; branches and conflicts are rejected. Open chains are
not falsely treated as closed regions. Simple closed loops now carry signed
area, centroid, and winding-derived enclosed/exterior region identity;
self-intersections and degenerate loops reject, while nested loops remain
separate until a bounded containment stage proves they do not touch, requires
authored region continuity through every parent, and closes exact per-region
area over an eligible MAC face. Faces with open chains, coplanar sheets, or
boundary-touching loops remain explicitly unresolved. The versioned
`simwing_scene_fluid_grid_epoch` boundary now composes this entire chain plus
unique conservative quadrature into one immutable accepted-Structure remap.
Its fingerprint prevents candidates, patches, topology, partitions, or loads
from being mixed across structural steps; every stage keeps its own bounds and
the aggregate owned payload has an additional byte ceiling. A focused moving
regression crosses a fabric triangle through a MAC plane, preserves authored
region/material/sheet identity, and retains force/moment closure. A first
`simwing_scene_fluid_cell_volume` subset now reconstructs deterministic sparse
per-cell region volumes for closed, consistently wound two-sided manifolds.
Each oriented interface triangle forms a signed tetrahedron against the grid
origin; deterministic convex clipping distributes that chain across exact
cells, including cells wholly inside a region and cases whose face-local chains
remain open at tile boundaries. Every tetrahedron closes across its clipped
cells, every cell closes across its regions, and the global sum is compared with
an independent whole-surface divergence volume. A separate immutable opening-cap
owner now closes an explicitly authored planar, strictly convex intake or
crossport loop for volume accounting only. Its winding comes from the adjacent
fabric boundary, follows accepted Structure motion, and never becomes fabric or
traction. The analytic open-tetrahedron regression therefore has a finite cell
volume. Stable one-point triangle samples now exactly integrate the accepted
piecewise-linear cap velocity into per-opening surface-sweep rates; a rigid
material-plus-cap surface closes that geometric ledger. The exact barycentric
triangle/box clipper is shared with material geometry and partitions each cap
into bounded positive-area grid patches. Off-face pieces have unique cell
owners; paired grid-plane copies become one canonical non-periodic face owner,
and their polygon, area, accepted velocity, and sweep remain explicit. A
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
cell/region ID, volume weight, component, and deterministic gauge owner. Mixed
cut cells and full interior cells close independently and exactly recover all
region and domain volumes. A complementary immutable face topology now links
those unknowns only through exact same-region Cartesian areas. Resolved nested
interfaces retain one link per region partition; untouched faces require one
unambiguous common region. A face-aligned authored opening now contributes its
exact oriented cross-region aperture area plus an unambiguous same-region link
over the rest of that Cartesian face. The analytic intake therefore retains
`0.18 m²` of opening and `0.82 m²` of Outside flow area instead of smearing the
whole face into either one. Material/interface ambiguity remains unlinked, and
physical pressure coefficients enter only in the downstream fixed-epoch
projection. For fully resolved geometry,
those links assemble into a bounded symmetric integrated graph Laplacian. Its exact
constant nullspace, positive link energy, component conservation, and retained
gauge owners are tested. Assembly requires one connected link graph per
authored pressure component. The face-aligned intake now joins Cell and Outside
with exact aperture graph energy, while a complete Cartesian grid with a
missing off-face intake link still rejects instead of acting like sealed
fabric. A deterministic component-wise conjugate-gradient layer now solves
manufactured integrated systems on both the nested closed regions and that open
intake. It rejects incompatible component sources, removes only admitted
roundoff, explicitly recomputes the final residual, sets every canonical gauge
to exactly zero, and leaves warm pressure unchanged on failure. Off-face
opening transmissibility remains future work. A first physical fixed-epoch
adapter now maps the predicted MAC field into one oriented volume flow per
resolved link, substitutes accepted fluid-minus-cap-sweep flow on each exact
face-aligned aperture patch, assembles `-rho/dt` times the integrated outward
flow as the pressure RHS, and applies the solved pressure difference back to
each link. It publishes pressure and corrected flow only after explicit
control-volume continuity closes; non-convergence exposes diagnostics without
a partial corrected state. The analytic intake retains separate aperture and
complement flows rather than collapsing its cut face back to one MAC value. A
consecutive-epoch owner now derives exact per-control-volume `dV/dt` from
stable cell/region identities. The moving projection overload targets
`dV/dt + net relative flow = 0`; with zero predicted air velocity, the
expanding analytic cell develops pressure and draws corrected flow inward
through its intake. Sparse cell/region appearance or disappearance rejects as
a topology rebase. Conservative rebase/remap, corrected-MAC reconstruction
across partitioned faces, off-face opening transmissibility, and general
moving-boundary projection remain future work.
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
that solve atomically: grid remap, opening caps and patches, sparse region
volumes, pressure controls, conservative links, and the ungauged operator all
share one Structure-state identity and an aggregate storage bound. Downstream
macro-step code can no longer accidentally mix those products across accepted
surface states.
The first topology-stable strong feedback owner now closes that pressure path
against XPBD. Each Aitken iteration restores the same structural baseline,
applies the trapezoidal average of the accepted start pressure load and a
relaxed end-load guess, rebuilds the pressure epoch, projects moving-volume
flow, and returns the conservative scene pressure load. Only a converged
iterate commits. The expanding open-tetra regression proves pressure opposes
the imposed expansion, continues deterministically into the next macro-step,
and restores both Structure and pressure ownership on iteration exhaustion or
projection failure. The MAC predictor is still held fixed during this first
feedback loop; momentum evolution and topology rebasing remain later work.
Its in-memory composite checkpoint retains Structure plus the accepted sparse
pressure projection; restore rebuilds and validates the complete pressure
epoch and conservative load before committing either owner. Initial and
post-step checkpoints reproduce the exact next coupled result in the same or
an equivalent owner, while foreign settings and corrupt or missing accepted
pressure reject transactionally.
The selectable `pressure-cell` worker makes this path visible. A soft
three-panel tetrahedral cell has one triangular intake, three fixed mouth
vertices, and a sinusoidally driven apex. Its immutable frames publish
deformation, area-averaged triangle pressure jump, nodal/total pressure force,
actuator force, and strong-iteration count. A 600-step headless run remains
topology-stable for 10 simulated seconds and reaches centimetre-scale motion;
this is a diagnostic of the new feedback path, not wing aerodynamics.
Its bounded `SWPCELL1` checkpoint stores the trusted Structure state and the
complete accepted sparse pressure projection. Initial and accepted files
round-trip deterministically, reject foreign/corrupt input transactionally,
and resume through the normal worker checkpoint flags.
Nonplanar or concave openings, surface
junctions, periodic-boundary ambiguity, and general moving-boundary fluid
equations still reject or remain open; the grid epoch itself continues to own
geometry and transfer only.
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
