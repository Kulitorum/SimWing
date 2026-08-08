# Coupled XPBD-CFD architecture and delivery plan

## Purpose and claim boundary

The goal is a two-way fluid-structure simulation of a ram-air paraglider that
can eventually reproduce launch inflation, trimmed flight, asymmetric collapse,
deep collapse, reopening, and reinflation. The wing structure remains a
large-deformation XPBD system. A new CFD system owns both external and internal
air and exchanges conservative surface traction and motion with XPBD.

This is an engineering simulation program, not a promise of a "perfect"
simulation. Absolute predictions require measured material properties and
validation data. Until that evidence exists, results must be labelled by their
validation level and include numerical and parameter uncertainty.

This document is the target architecture. The imported Playground simulation
is not a baseline for SimWing physics. Its archived records may explain the
old code, but its aerodynamic loads, pressure model, cell-air model, flight
dynamics, metrics, scenarios, and numerical results are not regression oracles
for the remake.

## Imported baseline

The initial SimWing tree was copied from the working tree at
`C:\CODE\LeParagliding` on 2026-08-07:

- Git commit: `497b5632154fb09ef14c4b47943d3fa31a184309`.
- Imported tracked and non-ignored files: 1,182.
- Excluded: `.git`, generated build trees, and ignored generated outputs.
- SimWing begins with independent root commit
  `bec9843df172ffebe70cfb17f4d5871418e4e98e`; the original commit is recorded
  in `UPSTREAM.md` but is not part of SimWing history.
- Personal designs, the private planning transcript, local fixtures, and build
  output are not part of the public repository.
- Future upstream transfers must be deliberate source imports. There is no
  compatibility promise for inherited Playground behavior.

### Initial verification snapshot

The imported tree was configured on Windows with CMake 4.1.1, Visual Studio
2022/MSVC 19.44, Qt 6.11.1, and Open CASCADE 8.0. The focused Release targets
compiled and their tests passed. Only the XPBD material, membrane, cable,
suspension, volume, and contact results are relevant inputs to the remake. A
passing Playground simulation test means that the import is intact; it is not
a SimWing physics acceptance criterion.

A first run appeared to hang because the narrowly built executables could not
find `Qt6Core.dll`; Windows displayed a modal loader-error dialog outside
CTest's output. Prepending `C:\Qt\6.11.1\msvc2022_64\bin` to `PATH` corrected
the test environment. This was not a simulation failure. Only focused targets
were built, so the remaining configured tests were not executed in this
snapshot.

## Executive decision

Build one new full-FSI simulation path. It uses the proven XPBD structural
primitives but none of the inherited Playground simulation. An adaptive
Cartesian finite-volume flow solver with a projection step and a sharp,
two-sided immersed-interface treatment owns internal and external air. It uses
strong partitioned coupling to XPBD and is the sole path to flight, inflation,
collapse, and reinflation.

There is no preliminary port of the Playground pressure solver, polar forces,
cell-air network, flight dynamics, or SoftWingLab reduced-order aerodynamics.
If an interactive reduced model is useful later, derive and validate it against
the completed CFD system; do not make it a foundation of the remake.

The transcript's central intuition is sound: an Eulerian flow grid avoids
remeshing around extreme fabric motion, and the fabric must preserve distinct
pressures on its two sides. The recommended production discretization is not a
uniform-grid LBM, however. A low-Mach paraglider needs a pressure-accurate
solver, local refinement near thin fabric/intakes/wakes, and explicit handling
of folded interfaces and sub-grid gaps. Uniform LBM is valuable as a bounded
research prototype, but its domain-wide memory cost and common stair-step or
smeared membrane boundaries make it the wrong first production commitment.

The closest published precedents are sharp immersed-interface simulations of
thin compliant parachute shells and porous parachute fabric. AMReX is a useful
candidate substrate for block-structured AMR, GPU portability, and projection
solvers, but its standard embedded-boundary facility cannot directly represent
all multivalued cells created by folded, zero-thickness fabric. SimWing would
need its own discrete-surface interface/region treatment on top of the grid,
not merely turn on AMReX EB.

## What is reused, adapted, and replaced

| Existing component | Decision | Reason |
|---|---|---|
| Translated LEparagliding calculation core | Keep unchanged behind the engine process | It is the geometry/design oracle and compatibility boundary. |
| OCCT exact model and viewport | Keep | It remains the exact design model and production viewer. |
| `src/softwing` XPBD constraints, membranes, cables, contact, and suspension | Reuse behind a new structural adapter | These are the trusted large-deformation building blocks. Existing aerodynamic and pressure integration is outside this boundary. |
| `playground_sim`, pressure solve, cell air, metrics, analysis, and scenarios | Discard from the new path | Their behavior and results are not useful simulation truth. New diagnostics are specified from conservation and validation requirements. |
| Playground GUI and monolithic worker | Replace rather than migrate | The CFD worker and its UI are new components with no feature-parity gate. |
| `lep-sim.json` version 1 | Do not consume | Scene-v2 is exported directly from authoritative geometry; no compatibility adapter is required. |
| SoftWingLab reduced-order VLM/wake implementation | Do not port | Reduced aerodynamics must not delay or shape the full-CFD architecture. |
| XFLR5 | Keep as a polar and aerodynamic comparison tool | It is not a dynamic FSI solver. |

Only isolated XPBD primitives and their canonical unit tests cross the remake
boundary. No existing Playground tunnel, glide, collapse, trim, pressure, or
shape result is an oracle.

## Target process architecture

```text
LEparagliding design + airfoils + Studio-only material data
                         |
               child calculation engine
                         |
             exact OCCT model + scene-v2
                         |
                         |
              simwing-fsi worker process
                  checkpointable, no UI
                         |
             +-----------+-----------+
             |                       |
      structure adapter        fluid backend
      XPBD/contact/lines   adaptive FV + sharp IIM
             |                       |
             +-- strong coupling ----+
                         |
          immutable diagnostic frame stream
                         |
     +-------------------+-------------------+
     |                                       |
simwing-viewer                         Studio / Flight Lab
standalone Qt/OpenGL                  later control integration
```

The worker is a separate executable, following the successful calculation
engine boundary. A crash, out-of-memory condition, or long CFD step must not
take down the editor. The same executable must run locally, from a command line,
and later through WSL or a remote scheduler without changing simulation logic.

Use a small versioned control protocol. Commands include load scene, initialize,
advance, pause, checkpoint, resume, cancel, and request fields. Responses include
progress, convergence, metrics, diagnostics, preview frames, and fatal errors.
Do not stream full CFD fields to the GUI every step.

## Live diagnostic visualization

Visualization is part of the numerical-development workflow, not late UI
polish. Implement a separate `simwing-viewer` executable in the first vertical
slice. An interactive `simwing-fsi` run launches it by default; `--no-viewer`
keeps CI, benchmarks, remote jobs, and scripted runs headless. `--viewer` forces
it on when launch-mode detection is ambiguous.

The viewer subscribes to immutable sampled frames over a versioned local
protocol. The worker never waits for rendering: frames live in a bounded ring,
old frames may be dropped, and closing or crashing the viewer cannot change the
simulation state, time step, ordering, or result. Pause, single-step,
checkpoint, and stop are explicit control messages processed only at solver
safe points.

The first viewer renders:

- skin, ribs, seams, suspension lines, payload, face orientation, and region
  labels;
- structural strain, constraint error, contact pairs, sealing state, line
  tension, and applied traction;
- simulation time, time step, coupling iteration, residuals, and conservation
  ledgers;
- selectable clipping planes, camera presets, wireframe, normals, and stable
  entity IDs.

The first periodic CFD subset now renders scalar-coloured cell-centre points
and bounded selectable velocity glyphs. As later fields arrive, add AMR blocks,
streamlines, pressure and vorticity slices, divergence, interface pressure jump,
porous flux, and wake tracers. Field extraction is explicit and rate-limited;
the default stream is a surface/diagnostic subset rather than a copy of the
whole volume.

Support recording the same frame stream to a replayable trace. A visual anomaly
can then be inspected frame-by-frame, shared without rerunning an expensive
case, and correlated with the exact numerical log/checkpoint. The viewer must
always show the scene checksum, solver commit, step, and simulated time so a
screenshot is meaningful evidence.

## Scene-v2 contract

Create one authoritative, solver-independent scene representation. JSON is
appropriate for metadata and inspection; dense arrays should live in a compact
binary payload to avoid parsing and size costs. Every quantity has an explicit
unit. The simulation boundary converts the application's millimetres to SI.

Required scene data:

- right-handed, Z-up coordinates in metres;
- stable vertex, edge, face, rib, cell, seam, and attachment identifiers;
- oriented skin triangles with material coordinates and surface roles;
- a region on each face side, such as `outside/cell-12` or
  `cell-12/cell-13`;
- explicit intake and vent boundary loops represented as openings, not faces;
- triangulated ribs with their actual crossport holes;
- diagonal and mini-rib membranes, including attachments to the skin;
- seam/attachment graph and intended coincident-node relationships;
- warp, weft, shear, bending, areal-density, damping, permeability, porosity,
  and seam-property assignments;
- suspension graph, line type, diameter, density, elasticity, drag model,
  riser/brake roles, and attachment coordinates;
- pilot/harness mass, inertia, reference axes, and attachment model;
- design checksum, exporter version, source length unit, and coordinate-frame
  transform;
- optional initial shape, velocity, cell pressure, flow field reference, and
  restart/checkpoint identifiers.

The exporter must use authoritative geometry rather than reconstructing ribs
and diagonals from sampled curves in the simulator. The new worker never reads
the coarse v1 format.

Scene validation is a standalone library and command-line tool. It rejects:

- missing or inconsistent units and transforms;
- zero-area, duplicate, non-finite, or inconsistently oriented faces;
- invalid side-region assignments or a disconnected/unexpected region graph;
- accidental faces over intakes and missing rib-hole boundaries;
- missing material or mass assignments;
- dangling seams, suspension edges, or attachments;
- duplicate stable IDs and incompatible checkpoint/scene checksums.

## Structural model

The first structural backend wraps the trusted portions of `softwing_core`; it
does not rewrite XPBD during CFD development. Build a new Qt-free assembly path
from scene-v2 rather than extracting Playground body construction. It has these
responsibilities:

- build nodes, membrane/bending/seam constraints, ribs, diagonals, lines,
  payload, and contact from scene-v2;
- expose a contiguous coupling surface with stable IDs;
- accept time-integrated nodal forces or impulses from conservative traction
  transfer;
- predict and correct interface position and velocity;
- save and restore all state for implicit-coupling iterations;
- expose surface normals, region labels, contact state, cell volume, and
  topology through new SimWing interfaces;
- compute structural energy, constraint residuals, penetration, line loads,
  and material strain in every accepted step.

The accepted opaque `SoftBodyCheckpoint` now has its own versioned persistent
codec. It stores every committed mutable node, pressure, multiplier, contact
warm-start, record, diagnostic, and audit value in a deterministic,
checksummed, byte-bounded little-endian envelope. Decode is transactional and
requires a checkpoint from an equivalent rebuilt body as its topology template;
immutable masses, connectivity, constraint definitions, membrane materials,
and contact registration are copied from that trusted template rather than the
wire. This is the structural-core payload needed by an enclosing Structure
checkpoint format. The complementary `SuspensionCheckpoint` also has a bounded,
checksummed deterministic codec. It preserves rigid payload state, controls,
line and ground multipliers, pending work, stable identities, complete nested
diagnostics, and bounded UTF-8 text. Decode binds counts and identities to an
equivalent checkpoint template and verifies the canonical complete-state
fingerprint; the live owner's transactional `restore()` remains the final
semantic validator. The enclosing `StructureCheckpoint` codec composes both
solver-owned envelopes with public accepted-step, time, pending-load, and
last-applied-load state. It validates serialization and decoding through an
equivalent rebuilt Structure and explicitly compares the restored body nodes
with their public duplicate before publishing output. The enclosing
`OpenPistonCaseCheckpoint` codec then nests this Structure envelope with the
accepted moving-interface fluid envelope, partial-cell/rebase epoch, and all
control-volume, cut-surface, transfer, and conservation ledgers. It is likewise
deterministic, bounded, checksummed, and transactional. The rebuilt open-piston
worker is the final semantic validator: it reconstructs the canonical topology
epoch, stable identities, geometry, acceptance state, face aggregates, bridge
balances, and analytic conservation relationships before publishing the
decoded composition.

Material calibration is not optional. Replace prototype fallback constants with
identified parameter sets, while keeping an explicitly named synthetic material
for regression tests. Required measurements include biaxial fabric response,
bias shear, bending/folding response, hysteresis/damping, fabric and seam
permeability, seam stiffness, line load-extension, and line mass/drag.

Self-contact needs continuous or conservative broad-phase handling at the CFD
macro-step scale. Contact state also informs the fluid sealing model: two
surfaces in contact cannot leave a numerically open channel merely because the
Cartesian grid cannot resolve their separation.

## Fluid model

### Governing regime

Start with filtered incompressible Navier-Stokes in an inertial frame. Typical
paraglider speeds are deeply subsonic, so a compressible solver adds stiffness
without helping the dominant load mechanism. Use a finite-volume projection
scheme with second-order spatial and temporal targets, bounded convection, and
explicitly monitored kinetic-energy behavior.

Initial turbulence progression:

1. laminar and manufactured-solution verification;
2. implicit LES for early moving-interface tests;
3. wall-modelled LES or a documented hybrid RANS/LES closure for production
   Reynolds numbers;
4. model comparison and uncertainty reporting against measured data.

No turbulence choice is promoted to default solely because one free-flight run
looks plausible.

### Grid and interface

Use block-structured Cartesian refinement around:

- the canopy, ribs, intakes, vents, and suspension lines where modelled;
- high vorticity and wake gradients;
- close or contacting fabric sheets;
- cell openings and active collapse boundaries.

The fabric is a zero-thickness oriented discrete surface. The flow method must
retain separate one-sided pressure and shear values, enforce no-penetration or
prescribed interface velocity, and allow tangential wall modelling. A sharp
immersed-interface/jump formulation is preferred over a regularized force delta
for production because pressure jump is the primary canopy load.

Rib holes, intakes, and vents are genuine connectivity in the fluid region
graph. Fabric and seam leakage use a calibrated Darcy-Forchheimer pressure-jump
law. They are not represented by deleting random faces.

Folded fabric can put more than one interface crossing in one Cartesian cell.
That is a mandatory verification case. If the selected grid substrate cannot
represent such a cell directly, the discrete-surface layer needs regional
subcell reconstruction or a conservative unresolved-interface model. A
single-valued embedded-boundary assumption is insufficient.

### Far field and initialization

Support two physical setups from the beginning:

- **wind tunnel:** fixed or controlled wing, prescribed inflow, far-field
  outflow, deterministic force/shape sweeps;
- **free flight:** co-moving computational window with inertial bookkeeping,
  far-field air state, payload gravity, and domain recentering that preserves
  momentum and Galilean invariance.

Inflation cases need a controlled initial condition: packed or partially open
geometry, line state, cell pressure, ambient flow, pilot motion, and boundary
motion. Avoid an impulsive velocity field; ramped or physically reconstructed
initial conditions should be part of the scenario definition.

## Conservative transfer and strong coupling

Surface-to-grid and grid-to-surface operators must be paired and tested for
force, moment, and power consistency. Integrated one-sided pressure and shear
produce the only full-mode aerodynamic load. XPBD returns interface position
and velocity; no polar force, point force, pressure stamp, or hidden damping may
add the same physics a second time.

For each accepted fluid macro-step:

1. Save a rollback checkpoint of structure and fluid.
2. Predict the interface at the new time.
3. Rebuild or update interface geometry, side regions, refinement, and sealing.
4. Advance advection/diffusion and solve the pressure projection using the
   predicted moving boundary.
5. Integrate two-sided traction and conservatively transfer it to XPBD nodes.
6. Advance the structure through its smaller substeps, including suspension,
   payload, internal constraints, and contact.
7. Evaluate displacement/velocity and traction residuals at the interface.
8. Relax the coupled update with Aitken acceleration initially and IQN-ILS when
   enough history exists.
9. Repeat from the saved iteration state until absolute and relative residual
   tolerances pass; otherwise reduce the macro-step and retry.
10. Accept once, update metrics/energy ledgers, and emit a preview/checkpoint as
    requested.

Single-pass explicit coupling is allowed only as a named diagnostic. The
light, flexible canopy has precisely the added-mass conditions in which weak
partitioned coupling becomes unstable.

Multi-rate stepping is expected: multiple XPBD/contact substeps per fluid
macro-step, with all coupling exchanges defined in impulse/time-integrated form.
Time-step control considers fluid CFL, viscous limits, interface motion per
cell, contact events, and structural convergence.

## Inflation, collapse, and reinflation

These phenomena are architecture requirements, not late visual effects:

- Internal and external air are one Eulerian fluid domain separated by the
  moving skin and ribs.
- Intake mouths and rib crossports provide actual resolved mass and momentum
  exchange.
- Cell pressure is measured from the fluid solution. The inherited 0D gas model
  and pressure stamping are not used.
- Collapse is detected from region-resolved volume, pressure difference,
  contact topology, inlet effective area, fabric inversion, and suspension
  load—not one visual threshold.
- Contact and the flow solver share an unresolved-gap rule. Closed fabric must
  seal; an opening that grows beyond resolution must reconnect fluid regions
  conservatively.
- Reinflation requires correct intake momentum, cross-cell communication,
  fabric separation, and contact release. Artificial pressure floors are not
  permitted in a validation run.
- Topological fluid-region events are logged and checkpointed so a failed
  reinflation is diagnosable.

The first collapse experiments use controlled tunnel disturbances: asymmetric
inlet blockage, prescribed gust, brake/riser actuation, and angle-of-attack
ramps. Free-flight collapse follows only after those cases conserve mass,
momentum, and energy within defined tolerances.

## Proposed code boundaries

Names are provisional, but ownership should be explicit:

```text
src/fsi/
    scene.*                 schema, serialization, validation, stable IDs
    structure.*             adapter around softwing_core
    transfer.*              conservative traction/motion exchange
    coupling.*              rollback, Aitken, IQN-ILS, acceptance logic
    piston_case.*           sealed fixed-topology worker canonical
    open_piston_case.*      driven open-control-volume worker canonical
    open_piston_checkpoint_persistence.* bounded composite restart codec
    periodic_flow_case.*    visible periodic CFD verification worker
    pressure_jump_case.*    visible static layered-jump worker oracle
    periodic_flow_checkpoint.cpp bounded persistent worker restart codec
    worker_control_protocol.* transport-neutral safe-point messages
    worker_control_stream.*  bounded self-framing binary stream adapter
    worker_control_session.* case-neutral safe-point command execution
    periodic_flow_control.* typed periodic worker control adapter
    open_piston_control.*   typed open-piston worker control adapter
    diagnostics.*           conservation, residuals, events, profiling
    checkpoint.*            versioned restart state
    scenarios.*             tunnel, glide, inflation, collapse definitions
    fluid/
        grid.*              block hierarchy and field storage
        geometry.*          discrete surface and region reconstruction
        advection.*         bounded donor/limited-MC MAC transport
        projected_advection.* projected nonlinear SSPRK2 transport
        diffusion.*         bounded periodic MAC viscosity verification
        evolution.*         transactional flow steps and subcycling
        projection.*
        interface_jump.*
        porous_interface.* calibrated flux-driven porous jump and ledger
        moving_interface.*  grid-face constraints and region topology
        moving_control_volume.* open planar GCL and topology epochs
        checkpoint.*        accepted fluid/topology checkpoint + persistence
        turbulence.*
        boundaries.*

src/gui/
    flight_lab_page.*       controls and results, no solver arithmetic
    flight_lab_controller.* worker lifecycle and protocol

src/viewer/
    viewer_protocol.*       immutable frame/control schema and trace files
    fluid_frame.*           owning accepted MAC-to-cell diagnostic adapter
    pressure_jump_frame.*   owning cell and layered-interface diagnostics
    vector_glyphs.*         bounded deterministic vertex-vector arrows
    viewer_window.*         standalone Qt/OpenGL diagnostics
    viewer_layers.*         structure, interface, grid, field, and HUD layers

tools/
    simwing_fsi_main.cpp    headless worker executable
    simwing_viewer_main.cpp standalone live/replay viewer
    simwing_scene_main.cpp  inspect/validate/convert scenes
    simwing_fsi_bench.cpp   canonical deterministic/throughput cases
```

Likely CMake targets are `simwing_scene`, `simwing_structure`,
`simwing_transfer`, `simwing_fluid`, `simwing_coupling`,
`simwing_fluid_structure_bridge`, `simwing_piston_case`,
`simwing_open_piston_case`, `simwing_fluid_frame`,
`simwing_periodic_flow_case`, `simwing_pressure_jump_case`,
`simwing_worker_control_protocol`,
`simwing_worker_control_stream`,
`simwing_worker_control_session`, `simwing_periodic_flow_control`,
`simwing_open_piston_control`, `simwing_viewer_geometry`, `simwing-fsi`,
`simwing_viewer_protocol`,
`simwing-viewer`, and focused test executables. Keep Qt out of the numerical
targets; only the viewer/UI targets link it. Backend interfaces must not be so
abstract that they hide grid layout or force extra full-field copies.

## Delivery phases and gates

Current foundation status: the Qt-free scene-v2 core, XPBD structural adapter,
deterministic scene-to-structure assembly, immutable structural diagnostic
frames, replayable trace protocol, and standalone Qt/OpenGL trace viewer are
implemented with focused tests. A canonical Qt-free structural worker now
launches the viewer by default and publishes accepted steps through a bounded
growing-trace follower; it is a pipeline harness, not a fluid or flight model.
The same worker can now select a smooth periodic Taylor-Green CFD case. Its
owning adapter copies accepted cell pressure and averaged MAC velocity into
stable cell-centre points. It also publishes the projection operator's exact
finite-volume divergence and a diagnostic centred periodic curl of the sampled
velocity. The unchanged frame protocol carries those fields to the viewer for
scalar-coloured point rendering. A separate Qt-free geometry target
turns selected vertex vectors into normalized arrows with a deterministic glyph
budget, allowing scalar colouring and selected velocity or vorticity direction
to remain visible together. This is inspectable verification flow, not a
whole-wing CFD claim.
Periodic worker advance now keeps the solver candidates private through frame
validation and commits fields, diagnostics, step, and time together. Its
versioned in-memory checkpoint binds those committed values to exact grid
metadata and a case-definition fingerprint; initial-state, same-worker, and
equivalent rebuilt-worker continuations replay bit-for-bit. A deterministic
little-endian file codec now preserves every MAC component, pressure sample,
nested diagnostic, step, and time behind a versioned, byte-bounded, checksummed
envelope. Decode is transactional and reuses the worker restore validator.
`simwing-fsi --checkpoint-in/--checkpoint-out` exposes additional-step resume
for both periodic-flow and open-piston workers, with codec-specific validation
and a shared atomic output-replacement path. `--checkpoint-every N` runs at
accepted-step safe points using absolute step multiples, preserving cadence
across restarts; the final accepted state is always saved without a duplicate
write. The first transport-neutral control envelope is also present:
bounded versioned/checksummed `advance`, `checkpoint`, and `stop` commands use
nonzero request IDs, while ready/advanced/checkpointed/stopped/error responses
always expose absolute accepted step and time. It deliberately chooses no
named-pipe, socket, or scheduler transport. A shared Qt-free control session
executes decoded messages synchronously on a worker owner thread. Typed
periodic-flow and open-piston adapters bind numerical advance, absolute state,
and their distinct complete checkpoint payloads. Each publishes immutable
accepted frames, delegates checkpoint persistence, and makes stop terminal
without putting output or file policy in the numerical worker. A bounded
self-framing stream adapter now reads the
envelope payload length without a host-native prefix and flushes every response.
`simwing-fsi --control-stdio` binds that stream to binary stdin/stdout for
either adapter, suppresses viewer launch and textual stdout, and completes the
trace before acknowledging stop. A restored worker reports its checkpoint's
absolute step/time in Ready and writes only newly accepted frames to its trace.
The exact-model capture now exports validated scene-v2.1 skins, authored open
intakes, triangulated holed ribs, internal sheets, explicit suspension
junctions, and the uncollapsed line graph when supplied explicit physical
material and pilot settings. Scene assembly adds per-sheet bending and preserves
the junction graph. It now orients one pilot's line forest toward its harness
roots and assembles the rigid payload; contact remains an explicit worker policy
because scene-v2 has no authoritative contact material yet. `softwing_core`
provides complementary transactional checkpoints for complete SoftBody/contact
and suspension/rigid-payload state, and `simwing_structure` restores them as one
composite transaction. Both the opaque SoftBody payload and the public
suspension/rigid-payload checkpoint have bounded, checksummed deterministic
codecs whose decode is bound to equivalent rebuilt topology and stable
identities. The enclosing Structure format now composes both and validates via
an equivalent rebuilt adapter. The bounded/checksummed open-piston format nests
that Structure state with its accepted moving-interface fluid epoch and
semantically revalidates every committed diagnostic ledger. The
real 3.28 regression now reaches an accepted coupled
structural step and replayable diagnostic trace with synthetic physical export
settings. Manufacturing flat-pattern UVs, exact authored line-attachment
vertices, authored paired seams and stitch mechanics, live bidirectional
control, an authoritative settings source/engine CLI, and general cut-cell
moving interfaces, nonplanar topology events, curved or changing grid-side
correspondence, AMR, and full CFD evolution kernels remain open work.
Phase 2 has started with a dependency-free uniform periodic MAC-grid
verification kernel:
its finite-volume gradient and divergence are a tested adjoint pair, its
zero-mean conjugate-gradient pressure projection commits transactionally, and
the focused regression covers Taylor-Green invariance, a discretely
manufactured exact projection, observed second-order pressure convergence, and
periodic momentum/kinetic-energy budgets. A first explicit laminar diffusion
operator applies the centred seven-point Laplacian independently on the three
translated MAC-component lattices. The exact stability boundary is
`nu*dt*(1/dx^2+1/dy^2+1/dz^2) <= 0.5`; excessive steps preserve the input
bit-for-bit. Accepted steps preserve all three component momenta, do not add
kinetic energy, keep solenoidal Fourier modes divergence-free, and reproduce
the discrete Fourier eigenvalue. A two-stage SSPRK2 companion uses that exact
Euler path twice and convexly averages the twice-advanced candidate with the
old field. It retains the same sharp per-stage stability boundary, periodic
momentum and non-increasing energy contracts, returns the Nyquist mode exactly
after its two boundary sign flips, and shows second-order temporal convergence
against the discrete Fourier decay. Coupled operator splitting is not yet the
intended second-order production time integrator. A companion unsplit
donor-cell oracle transports every MAC component by one prescribed uniform
velocity. Its update is a conservative convex combination for
`sum(abs(U_i)*dt/h_i) <= 1`: it preserves periodic component momentum, cannot
create a new component extremum or increase kinetic energy, commutes with the
discrete divergence, and becomes an exact one-cell periodic translation at the
sharp CFL-one boundary. Its full-period sine regression converges at the
expected first order. A variable-flow companion averages a divergence-free MAC
advector onto every translated-component control-volume face and uses one
shared periodic upwind flux. Its default donor reconstruction delegates the
uniform subset bit-exactly to the oracle. General accepted donor steps preserve
all three component momenta, stay within old-time component bounds, and do not
add kinetic energy under a local outgoing-CFL limit; the aliased field path
supplies first-order nonlinear self-advection. A periodic shear fixture observes
the expected first-order refinement. The same conservative faces now support a
selectable monotonized-central reconstruction. It is enclosed in fixed-advector
SSPRK2 because an individual limited forward-Euler stage can add the expected
`O(dt^2)` energy even when the committed convex aggregate is dissipative. That
intermediate exception is explicit in diagnostics and does not relax the
committed step: the aggregate still enforces original component bounds, component
momentum, finite state, and non-increasing energy. A discontinuous pulse remains
bounded, while smooth uniform full-period transport shows near-second-order L1
refinement. Nonlinear spatial behavior requires the projected enclosure below.
A pressure-projected nonlinear SSPRK2 operator now keeps every intermediate
incompressible before it can become its own advector: stage one selected
transport is projected, stage two self-advects from that accepted field, the
twice-advanced prediction is convexly averaged with the old field, and the
result is projected again. Both pressure and velocity remain private until all
four stages and the aggregate momentum/energy/continuity ledger pass. Exact
manual-stage composition, repeated-step eligibility, and a fixed-grid vortical
refinement show deterministic second-order temporal behavior. Donor remains
the default; limited MC is selectable. A Galilean-translated Taylor-Green
vortex supplies the first analytic nonlinear spatial canonical. Scaling `dt`
with `h^2` suppresses the SSPRK2 time error; 16/32/64-grid L1 refinement is
first order for donor reconstruction and near second order for limited MC.
This smooth periodic result does not establish cut-cell or moving-interface
accuracy. A symmetric second-order temporal flow path now reconciles its
internal pressure stages with viscosity: half-step
SSPRK2 diffusion, full projected nonlinear SSPRK2 transport, and the matching
viscous half step.
Every sub-integrator works on the same private candidate, their individual
energy losses sum to the full-step loss, and pressure/velocity commit only
after the aggregate momentum, energy, and continuity ledger passes. Exact
manual composition, rollback at the viscous/transport/projection boundaries,
and fixed-grid refinement verify the Strang path. A continuous viscous,
Galilean-translated Taylor-Green solution independently verifies the complete
subcycled limited-MC path: with `dt <= 0.12 h^2`, successive 12/24/48-grid L1
error ratios must remain in `[3.0, 5.0]`. This is a smooth periodic full-flow
oracle and does not establish cut-cell or moving-interface accuracy. A bounded
outer wrapper now pre-sizes equal Strang substeps from the explicit
half-viscosity limit. An
explicit outgoing-CFL or limited-MC maximum-principle rejection restarts the
whole private interval with a finer equal schedule; projection and ledger
failures remain fatal. Diagnostics retain the final schedule and the last retry
trigger, the configured bound cannot exceed 4096 substeps, and caller pressure
and velocity commit only after an independent outer momentum/energy ledger.
Exact two-step viscous composition, CFL/limiter retry equivalence to the final
manual schedule, bounded-limit rollback, and projection-failure non-retry are
regressed. The original composed
periodic macro-step still selects a single-stage transport mode and either
Euler or SSPRK2 viscosity, then runs the zero-mean projection on
private candidates. Velocity and pressure commit together only after every
stage and an independent aggregate momentum/kinetic-energy ledger pass.
Regressions are bit-identical to their standalone uniform/nonlinear transport
and Euler/SSPRK2 diffusion stages, retain every stage diagnostic, and reject at
advection, diffusion, or projection without changing either input field. This
is the first complete nonlinear evolution path; its single-stage transport and
operator splitting remain first order. A canonical ordered grid-face field now
retains every stable surface ID, two-sided fluid-region transition, normal
crossing position, and signed pressure discontinuity. Multiple crossings on
one face-normal cell segment must form a continuous region chain; its paired
sharp gradient and Poisson source use their deterministic signed sum. A split
static slab is bit-identical to the compact jump, while a balanced folded
subcell pocket creates no pressure or flow. A separate physically grid-bound
face-aligned interface now
removes constrained faces from the pressure graph, partitions cells into stable
fluid regions, and holds prescribed normal MAC velocity exactly while solving
one zero-mean pressure correction per region. The prior mean pressure in each
region is retained because those means are independent incompressible null
modes. Nonzero fixed-grid region volume rate is rejected transactionally rather
than hidden by compatibility subtraction. A nonseparating plane instead uses
one region ID on both sides when the remaining grid graph supplies a resolved
fluid path around it. The first planar moving-control-volume operator binds one
such surface to one open MAC plane and follows its offset through one partial
cell. It computes geometric volume change, surface sweep, and projected opening
transport independently and requires all three to agree before acceptance.
At an explicitly terminal cell-boundary step it can construct a candidate on
the next positive-axis MAC plane. Stable identities and orientation must remain
unchanged, and the old terminal volume must equal the candidate reference
volume before the caller may commit the new epoch. Skipped planes and collision
with the resolved opening are rejected. This is a planar topology-epoch
verification subset, not a general cut-cell pressure or surface-reconstruction
scheme. Adjacent cell-centre pressure gives exact pressure force/impulse/work
for the compatible piecewise-constant slab canonical. When imposing a face
velocity also changes its predicted MAC degree of freedom, however, pressure
traction alone is not the coupled load. The complete fluid-on-interface
reaction on this uniform grid is
`A*(p_minus-p_plus)*n - rho*V_face*(U-u*)/dt`. Both terms remain explicit; the
second is exactly zero when the predicted velocity already meets the
constraint. Accepted diagnostics retain each canonical MAC tile's rectangle,
pressure traction, direct-enforcement force, complete reaction, constrained
velocity, impulse, and power ledger. This is deliberately not presented as
arbitrary surface reconstruction, general cut-cell pressure geometry,
folded/multiple-crossing motion, advection, turbulence, or a wing-flow solver.

The companion planar cut-surface operator now gives the accepted constraint
reaction an explicit physical application geometry within that partial cell.
It retains each MAC rectangle as the Eulerian source and translates a second
rectangle to an unwrapped physical plane only when grid-plane plus offset is a
matching periodic image. Independent face and surface area, force, physical
moment, and power ledgers must close. This remains the projection's Lagrange
reaction at a bounded planar cut, not an interpolation of new one-sided cell
pressures and not general cut-cell reconstruction. Projection produces that
constraint reaction as a macro-step average. A separate temporal adapter keeps
the accepted face force fixed while resampling congruent physical geometry and
rigid normal velocity at each coupling endpoint, so impulse and work use one
consistent step-average reaction.

The first structure-side conservative-transfer primitive is also present. A
canonical stable-ID surface is fingerprint-bound to one immutable structural
definition and rejects reversed or foreign topology. Uniform current-triangle
traction is integrated at the centroid and distributed in equal barycentric
shares. The same boundary also accepts explicit stable area/barycentric
quadrature patches, allowing nonuniform traction while exactly preserving
resultant, moment, and work against linear triangle kinematics. Surface and
nodal ledgers are accumulated independently, and immutable results are applied
additively through the real Structure nodal load path only after all bindings
validate. A separate versioned macro-step coupling layer trapezoidally
integrates nonuniform local-time samples into
immutable nodal impulses and independent surface/nodal impulse, angular
impulse, and work ledgers. Its acceptance boundary requires the exact exchange
duration, applies the equivalent average force through Structure's internal
XPBD substeps, and restores the checkpoint from before load application on any
failure. The first stable-ID fluid-to-structure bridge is present for the exact
uniform subset of those operators. It accepts one canonical face-aligned fluid
surface, requires its facewise pressure-traction deviation to meet an explicit
budget, reconstructs one uniform world-space traction, and then requires
independent fluid/structure area, force, and power ledgers to close before
returning the immutable structural transfer. A planar face-resolved extension
clips canonical nonoverlapping MAC tiles against fully covering, consistently
oriented reference triangles once. Every material overlap above the explicit
`1e-14 m^2` sliver threshold becomes a stable barycentric patch; current
nonuniform tile traction is mapped through those patches only after per-face
and aggregate area, force, moment, and power ledgers close. Its fixed-material
mode requires the original Eulerian face geometry exactly. Its
rigid-normal-translation mode instead keeps transverse tiling and material
patches immutable while explicitly separating the current Eulerian grid plane
from the unwrapped physical plane. It accepts only a matching rigid structural
translation and uniform normal velocity, with independent position and velocity
correspondence ledgers, and is verified for X, Y, and Z normals and periodic
grid-plane wrap. The `--case piston` worker now evaluates compatible
start/end face-resolved fluid samples, trapezoidally integrates them, accepts
the impulse through XPBD, and publishes the resulting moving two-triangle
surface and CFD ledger fields to the standalone viewer. The
`--case open-piston` worker adds the nonseparating connected-fluid projection,
partial-cell geometry, resolved-opening GCL ledger, an explicit plate actuator,
and a separately reported resisting CFD load. Its pressure reaction now crosses
the moving face-resolved bridge at the structural plate's physical plane rather
than through a surface-uniform bridge, but only after the fluid-side cut-surface
area, force, moment, power, and periodic-image ledgers are accepted. Both the
numerical state and published frame cross accepted boundaries only. The worker
also independently closes pressure/internal versus actuator/external impulse
and work against structure, fluid, and combined momentum and kinetic-energy
changes. On an exact cell crossing it also
verifies old/new chamber-volume continuity, remaps the constraint within a
written velocity budget, and commits the new fluid/control-volume epoch with
the completed viewer frame. This does not yet implement transverse or
non-rigid surface deformation, curved correspondence, general cut-cell pressure
metrics, nonplanar topology events, fluid/structure iteration, or a
strong-coupling convergence decision.

The first fluid checkpoint boundary is now implemented for an accepted
face-aligned moving-interface epoch. It captures pressure, velocity, exact
interface kinematics, and the complete projection diagnostics in an immutable
in-memory payload. Its deterministic bounded/checksummed little-endian codec
persists the complete state and reconstructs the grid/interface before reusing
the accepted-state validator transactionally. Restore rebinds exact grid
metadata and a deterministic topology fingerprint before returning a candidate
state. The codec regression requires byte-identical repeated and decode/re-encode
output for the complete nested state and for a rebased topology epoch. It
rejects bad magic/version/reserved bits, truncation, trailing bytes, checksum
damage, a recomputed-topology mismatch, and byte/sample/face/region/surface
limits without changing the caller's prior checkpoint.
The open-piston worker composes that payload with the complete Structure
checkpoint, partial-cell offset/topology epoch, and committed
coupling/conservation diagnostics. Its bounded, checksummed, deterministic
little-endian codec nests both solver-owned payloads and every committed
diagnostic ledger, then validates the decoded composition through an equivalent
rebuilt worker before publishing it. An
ordinary continuation and the first steps after both a normal and periodic
plane rebase replay bit-for-bit in the same or an equivalent rebuilt worker.
The worker CLI restores this composite before creating its trace and exposes
atomic same-file output, absolute accepted-step autosave, and final-save
deduplication. Its typed control adapter also delegates the same composite at
safe points through the shared transport-neutral command executor.

### Phase 0 — Establish the remake boundary

Work:

- define the exact `softwing_core` classes retained for membranes, constraints,
  cables, suspension, and contact;
- create new canonical XPBD tests for a membrane patch, bending strip,
  suspension graph, line load, self-contact, and closed volume under an
  analytic external load;
- define checkpoint/state and material-property interfaces needed by coupling;
- create the new SimWing CMake targets with no Playground dependencies;
- implement the versioned diagnostic-frame protocol and show the canonical
  XPBD cases in the standalone viewer;
- mark the inherited Playground targets for removal once the new structural
  target and worker skeleton build.

Gate: the new headless structural target exercises XPBD directly and does not
link any Playground simulation library. Its tests are based on analytic or
structural invariants, not inherited flight results.

### Phase 1 — Scene-v2 and structural extraction

Work:

- extend the exact model capture/export with authoritative simulation
  triangles, material coordinates, regions, rib holes, seams, and attachments;
- implement the v2 binary payload, reader, and validator with no v1 adapter;
- assemble the wing through a new structural adapter using `softwing_core`
  primitives directly;
- implement new conservation, strain, line-load, volume, and contact
  diagnostics from their definitions;
- add rollback/checkpoint support and structural conservation diagnostics;
- render stable scene IDs, materials, regions, strain, contact, and line loads
  live in `simwing-viewer`.

Gate: scene validation fixtures pass; a headless structural run has no Qt,
Playground, aerodynamic, pressure, or v1-format dependency; analytic structural
cases meet their declared tolerances.

### Phase 2 — CFD verification kernel

Status: in progress. `simwing_fluid` currently owns the uniform periodic
verification grid, bounded donor and limited monotonized-central uniform and
divergence-free variable-flow velocity transport, nonlinear self-advection,
pressure-projected second-order temporal transport, Euler and second-order
SSPRK2 laminar velocity diffusion, their
transactional selectable composed macro-step, the symmetric second-order
temporal Strang flow path and bounded transactional subcycling, pressure
projection, and fixed-topology face-aligned moving
constraints, plus the first open planar one-partial-cell control-volume
operator, its exact next-plane topology rebase, and the bounded physical
cut-surface reaction geometry described above. On solve failure it preserves
the input pressure and velocity bit-for-bit so future macro-step retry can be
transactional. Accepted moving-interface state now also has a versioned,
immutable checkpoint bound to exact grid geometry and a stable topology
fingerprint, plus deterministic bounded persistent serialization of its exact
fields, interface kinematics, topology, and nested projection diagnostics. It
also owns a canonical ordered sharp pressure-jump field with explicit
surface, two-sided region identity, and normal crossing position. Multiple
crossings on one grid face retain their individual metadata, require a
continuous region chain, and contribute through one deterministic aggregate
stencil jump. An owning frame adapter now retains every individual crossing as
an oriented quad at its subcell fraction with its two region IDs, signed jump,
normal, and reconstructed diagnostic pressure. A static headless worker makes
that layered field replayable without involving Qt. The regression budgets are:
gradient/divergence adjoint error at or below `2e-14` in the canonical
integral, Taylor-Green maximum divergence below `2e-14 1/s` with a
bit-identical zero correction, discretely manufactured post-projection L2
divergence below `2e-12 1/s`, pressure convergence ratios in `[3.9, 4.2]` for
two successive resolution doublings, viscous eigenvalue convergence ratios in
`[3.95, 4.05]` then `[3.98, 4.02]`, exact zero and uniform viscous modes, exact
acceptance of the non-amplifying diffusion-number `0.5` boundary, no
kinetic-energy increase for every committed transport aggregate, donor-cell
full-period error ratios in `[1.7, 2.2]` then `[1.8, 2.1]`, limited-MC L1
full-period error ratios in `[3.1, 4.8]` over two successive resolution
doublings, projected translating-Taylor-Green donor L1 ratios in `[1.7, 2.2]`
and limited-MC ratios in `[3.0, 5.0]` across 16/32/64 grids with `dt`
proportional to `h^2`, complete viscous translating-Taylor-Green limited-MC L1
ratios in `[3.0, 5.0]` across 12/24/48 grids with `dt <= 0.12 h^2`, exact
bounded CFL-one translation, and periodic
component-momentum sums preserved within `5e-14` in the projection mixed-mode
case; viscous physical-momentum residual stays below `2e-12 N*s`. Subcycling
must pre-size the viscous canonical to exactly two equal steps, reproduce that
manual schedule bit-for-bit, restart a rejected nonlinear CFL/limiter attempt
to another exact equal-step schedule, and roll back without retry on projection
failure or when the configured substep bound is exhausted. The
periodic-flow worker additionally replays 12 accepted steps bit-for-bit across
two independent instances. Each immutable frame owns exactly one stable point
per cell plus pressure, speed, velocity, exact MAC divergence, and diagnostic
vorticity fields and remains byte-identical after the solver advances again.
The initial Taylor-Green cell sampling has an analytic face-average and
discrete-curl oracle, with maximum divergence below `2e-14 1/s`; later accepted
frames must reproduce the finite-volume divergence bit-for-bit and retain exact
vorticity vector/magnitude agreement. A completed five-frame trace must replay
consecutive accepted steps, and the headless CLI must write a completed trace
without linking or launching Qt.
The static pressure-jump worker uses a `16 x 4 x 3` periodic slab with two
ordered transitions on every face at each of two boundaries: 48 crossings,
384 owning frame vertices, and 96 oriented triangles. Each fresh accepted
projection must reproduce the analytic `-125/+125 Pa` cell field, remain below
`2e-13 m/s` maximum velocity and `1e-12 1/s` frame divergence, and serialize
bit-identically across independent and repeated workers. This makes multiple
same-face layers inspectable end-to-end but does not claim moving folded or
cut-cell topology.
A calibrated Darcy-Forchheimer adapter now samples resolved X/Y/Z MAC normal
velocity relative to each authored sheet tile and emits the corresponding
signed sharp jump. It retains tile area, volume flow, and nonnegative pressure
dissipation, and its exact monotone inverse is stable across the pure-linear,
pure-quadratic, and combined fits. A periodic prescribed-flux plane with an
explicit balancing pressure source preserves a `195 Pa` porous loss without
changing its uniform velocity. This establishes the constitutive and ledger
boundary only; pressure-driven implicit porous coupling is still required for
the full porous-sheet canonical.
The case checkpoint must restore the initial state and a later accepted state,
then reproduce the next frame bit-for-bit in both the original and an equivalent
rebuilt worker. Version, case fingerprint, grid, sample count, step, time, and
payload presence are independently rejected before mutation. Its persistent
encoding must be byte-deterministic across repeats and decode/re-encode, retain
all nested diagnostics exactly, and continue bit-for-bit after either initial
or accepted-state decode. The default untrusted-input limits are `64 MiB`,
5,000,000 scalar samples, and 4096 substeps. Bad magic/version/reserved bits,
truncation, trailing data, checksum corruption, and every configured limit must
leave the destination checkpoint unchanged. The CLI integration writes four
steps, restores them, advances three additional steps, atomically replaces the
same restart file, decodes it again, and reports seven total accepted steps
without launching Qt. With interval two, the four-step run writes at steps two
and four, while the resumed run writes at absolute step six and final step
seven; each reports exactly two writes. Interval mode without an output path is
rejected before a worker or trace is created.
The open-piston CLI counterpart checkpoints at absolute steps 600 and 1200,
therefore persisting the first accepted topology rebase. It restores that file,
advances three additional steps, atomically replaces the same path at step 1202
and final step 1203, and decodes the result in a zero-step worker. A periodic
checkpoint presented to the open-piston worker is rejected before trace
creation. Persistence tests also recompute the outer checksum after changing
diagnostic stable identities, acceptance flags, stale rebase state, or physical
cut-plane geometry; semantic restore must reject each payload transactionally.
The separate control-protocol regression requires byte-identical repeat and
decode/re-encode results for all three command and five response kinds. It
rejects cross-decoding, zero request IDs, zero/oversized advances, misplaced
kind-specific fields, non-finite worker time, oversized error text/messages,
bad magic/version/reserved bits, corruption, truncation, and trailing bytes;
every failed decode leaves the caller's prior message unchanged.
The control-session regression then advances and publishes three consecutive
accepted frames, delegates a checkpoint that reproduces the exact next frame,
and requires checkpoint/protocol failure not to advance the solver. A sink
failure after acceptance reports that absolute committed step in a valid
bounded error response and continuation starts there. Stop prevents later
advance/checkpoint mutation, while a repeated stop is idempotent and receives
a newly correlated response.
The periodic and open-piston stdio integrations each send advance-two,
checkpoint, advance-one, and stop to a separate worker process. They require
an exact ready/advanced/checkpointed/
advanced/stopped response stream with no trailing stdout, a completed
three-frame trace, and a step-two checkpoint whose next frame is byte-identical
to trace step three. A second process restores that checkpoint, reports Ready at
step two, advances once, and produces a completed one-frame trace whose step
three is byte-identical to independent replay. The in-process open-piston
adapter additionally advances through the first topology rebase in one command,
delegates that composite checkpoint, and reproduces the next frame exactly.
Control mode is rejected for structural and sealed-piston cases.
The Qt-free viewer-geometry regression separately requires exact arrow shaft
direction and relative magnitude, one shaft plus two arrowhead segments per
retained nonzero vector, dimensional automatic scaling, owning deterministic
output, finite derived geometry, and bounded integer-stride sampling.
The
static 250 Pa slab must retain its two pressure levels within `2e-12 Pa` and
keep spurious velocity below `2e-13 m/s`; a jump placed on periodic face zero
has separate `1e-12 Pa` and `1e-13 m/s` budgets. These tolerances apply to the
serial CPU verification backend, not yet to a future parallel reduction
contract. The uniform-traction transfer canonical uses a `6 m^2` rectangle with
analytic force `[12, -18, 30] N` and moment `[12, -57, -39] N*m`; force,
moment, translation-power, and rotation-power residual budgets are `2e-14 N`,
`3e-14 N*m`, `1e-14 W`, and `2e-14 W`, respectively. It also verifies cyclic
topology canonicalization, bit-identical replay, foreign-result rejection
before mutation, and the resulting velocity through an accepted XPBD step. The
temporal canonical prescribes a `6 m^2` piston moving at `0.25 m/s` for `0.4 s`
under pressure rising linearly from `100 Pa` to `120 Pa`. Nonuniform samples at
`0`, `0.1`, and `0.4 s` must reproduce the analytic `0.6 m^3` swept volume,
`264 N*s` impulse, and `66 J` pressure-volume work. Impulse and angular-impulse
residuals are each below `3e-13` in their SI units, work residual is below
`1e-13 J`, and the accepted XPBD momentum gain is the same `264 kg*m/s`.
The fluid-side translating-slab canonical uses the same `6 m^2`, `0.25 m/s`,
`0.4 s`, and `110 Pa` values. Each wall must report `660 N`, `264 N*s`, and
`66 J` with opposite signs, while the closed slab totals remain zero. Interface
normal-velocity error is bit-exact zero on X, Y, and Z faces; disturbed regional
projection has L2 divergence below `2e-11 1/s` and compatibility volume-rate
roundoff below `2e-15 m^3/s`. A one-sided `1.5 m^3/s` fixed-grid volume request
is rejected without changing pressure or velocity.
The open-piston canonical uses one `6 m^2` nonseparating plane moving at
`0.25 m/s` for `0.4 s` in a connected periodic region. Projection routes a
uniform `0.25 m/s` flow through the explicit opening plane. The analytic
`0.6 m^3` chamber growth must agree independently with partial-cell geometry,
surface sweep within `3e-15 m^3`, and opening transport within `3e-12 m^3`;
the same operator is checked on X, Y, and Z axes. An explicitly terminal step
then rebases by one face, including across the periodic boundary, with exact
old/new volume continuity. Incomplete cells, skipped planes, changed region
identity, corrupted area/volume ledgers, and collision with the opening are
rejected.
The uniform bridge canonical maps the slab's stable right-wall ID to the same
`6 m^2` two-triangle structure, closes `660 N` and `165 W` at `0.25 m/s`, and
delivers `264 N*s`/`66 J` through the temporal exchange to accepted XPBD
momentum. It rejects nonuniform face traction, an absent surface ID, mismatched
area, failed fluid projection, and inconsistent interface power. The planar
face-resolved canonical then uses a disturbed projection with genuinely
nonuniform tile pressure. Deterministic tile/triangle clipping preserves its
force, moment, global power, and every face's local power within `2e-10` in the
corresponding SI units; gaps, tile or triangle overlap, orientation mismatch,
changed geometry, and locally inconsistent power are rejected. The canonical
additionally verifies rigid-normal translation on every grid axis. The
Eulerian grid plane may rebase, including a periodic index wrap, while the
physical plane remains unwrapped; stable material patches, transverse geometry,
matching node positions, and matching face/node normal velocities must remain
exact. Before this bridge transfer, the fluid-side cut-surface canonical checks
the same three axes and the terminal/rebased periodic images. Its complete
physical plane retains `16 m^2` area and closes source-versus-cut force and power;
its step-average reaction resampled from zero to terminal velocity closes the
analytic plug-flow kinetic-energy change. Invalid offsets, noncongruent physical
planes, failed projections, corrupted aggregates, and altered cut geometry are
rejected. The visible
piston worker repeats the face-resolved full chain at `120 Hz` with a synthetic
`6000 kg` tributary-mass plate so 600 default frames show about `1.38 m` of
deterministic translation without leaving the initial viewer scale. The open
piston worker uses the same mass and rate but drives at `0.05 m/s`: 600 frames
move `0.25 m`, grow the analytic chamber from `12` to `13.5 m^3`, expose the
actuator impulse separately from the moving face-resolved CFD pressure reaction,
and remain inside the first `0.5 m` cell. At frame 1200 the plate reaches
`3.5 m`, the chamber reaches `15 m^3`, and the worker transactionally advances
grid plane 6 to 7 with zero volume residual and less than `2e-12 m/s` velocity
remap; frame 1201 verifies continued partial-cell growth in the new epoch. At
frame 2400 the physical plate reaches the unwrapped `4.0 m` position while the
Eulerian plane wraps from index 7 to 0; frame 2401 verifies continued transfer
and growth after that periodic topology epoch. On the first acceleration step,
the `28.8 kg` fluid plug gains `1.44 N*s` and `0.036 J`; the complete reaction
delivers `-1.44 N*s` and `-0.036 J` to the structure, while the actuator supplies
`301.44 N*s` and `7.536 J`. Structure, fluid, and combined momentum residuals
must remain below `1e-8 N*s`, and their kinetic-energy residuals below `2e-9 J`,
including the two topology crossings.

Work:

- evaluate AMReX as grid/AMR/linear-solver infrastructure with a custom sharp
  discrete-surface layer;
- in parallel conceptually, use IBAMR or OpenFOAM/preCICE as an external
  reference for selected canonical cases, not as a required GUI dependency;
- implement arbitrary moving-interface/jump conditions, curved/changing
  paired grid-side correspondence, refinement, and broader case-specific
  control messages;
- extend the current cell-point and pressure-jump-layer diagnostics with
  rate-limited AMR blocks, slices, and traction as each field becomes available;
- verify CPU first, then add GPU kernels behind identical numerical tests.

Required canonical cases:

- Taylor-Green vortex and manufactured divergence/pressure solutions;
- channel/flat-plate flow and a static pressure jump across a membrane;
- flow through a calibrated porous sheet (the explicit prescribed-flux
  constitutive subset is complete; implicit pressure/flux coupling remains);
- moving piston/membrane with analytic volume and work balance;
- flexible flag and a thin shell with one interface per cell;
- two closely spaced and folded sheets with multiple crossings per cell;
- resolved opening that closes below grid scale and reopens;
- force, moment, and power transfer under rigid translation and rotation.

The structure side now verifies prescribed volume, interface work, temporal
impulse transfer, and transactional XPBD acceptance. The fluid side verifies a
volume-compatible translating sealed slab with exact face velocity and matching
per-wall pressure impulse/work. The planar face-resolved stable-ID subset now
also crosses those accepted fluid samples into XPBD and the viewer with explicit
interface residuals. It supports both exact fixed material correspondence and
rigid normal motion across grid rebases while preserving the transverse material
patches. The first open volume-changing piston now has one partial-cell geometry
update, an independently closed surface-sweep/opening-transport GCL ledger,
accepted physical cut-surface reaction geometry, and exact planar one-face
topology rebases including periodic wrap. Its analytic structure/fluid/actuator
momentum and kinetic-energy ledger is now complete for this driven planar case.
Its accepted fluid and structural state can also resume bit-identically from a
composite persistent checkpoint before or after a topology rebase.
General interpolated cut-cell pressure metrics, curved or transversely deforming
correspondence, nonplanar or opening topology events, sealed deforming chambers,
a complete coupled energy ledger for those general cases, and richer
case-specific control messages remain open, so the
full Phase 2 piston gate is not yet closed.

Gate: observed convergence is consistent with the intended order away from
interface singularities; mass, force, moment, and power errors satisfy written
budgets; folded-interface and unresolved-gap cases do not leak or change fluid
regions spuriously.

An LBM prototype proceeds beyond this phase only if it demonstrates adaptive
refinement, sharp two-sided pressure, folded-interface handling, acceptable
memory, a compatible license, and better measured cost/accuracy than the
projection solver.

### Phase 3 — One-way CFD on an inflated wing

Work:

- run frozen designed and XPBD-equilibrated shapes;
- perform mesh/refinement/domain/time-step studies;
- compare pressure maps, force/moment polars, wake, and intake mass flow with
  wind-tunnel data, XFLR5 where applicable, and an independent CFD reference;
- establish wall/turbulence model uncertainty.

Gate: integrated loads and pressure distributions are grid-converged within a
declared band and agree with validation data closely enough to justify two-way
coupling. A matching global lift coefficient alone is insufficient.

### Phase 4 — Two-way inflated-wing FSI

Work:

- implement rollback-capable Aitken coupling, then IQN-ILS;
- add adaptive macro-step control and coupled residual budgets;
- validate static inflation and steady deformed trim before maneuvers;
- compare structural strain, pressure, cell shape, line loads, and pilot motion.

Gate: results are insensitive to further coupling iterations, structural
substeps, and fluid time-step reduction within written tolerances; no hidden
force path is needed to keep the canopy inflated.

### Phase 5 — Resolved ram-air inflation

Work:

- enable internal/external region connectivity through intakes and crossports;
- start with restrained symmetric inflation, then moving-payload launch cases;
- calibrate permeability and inlet loss only against independent measurements;
- validate fill sequence, cell pressures, projected area, opening shock, line
  loads, and inflation time from synchronized experiment data.

Gate: symmetric inflation repeats across grid/time-step refinement and matches
the validation envelope without artificial pressure injection or geometry
stabilization.

### Phase 6 — Collapse and reinflation

Work:

- add fluid/contact sealing and reconnection events;
- validate controlled asymmetric disturbances in a tunnel configuration;
- progress through frontal collapse, asymmetric collapse, cravat-like contact,
  and recovery inputs;
- use ensembles for chaotic cases and compare distributions, not one trace.

Gate: mass does not cross sealed contact, topology events converge under
refinement, and measured collapse/reopening thresholds and time histories fall
within declared uncertainty bands.

### Phase 7 — Productization and performance

Work:

- GPU-port proven hotspots, add distributed block execution if required;
- mature progressive preview fields, trace replay, and robust
  cancellation/checkpoint resume;
- create presets for canonical verification, development resolution, and
  offline validation resolution;
- add result provenance: commit, scene checksum, backend, grid, time steps,
  tolerances, material set, validation level, and hardware;
- remove the inherited Playground from the SimWing build and UI; there is no
  behavioral or saved-scenario compatibility requirement.

Gate: production cases are restartable, reproducible within the declared
parallelism contract, bounded in memory, diagnosable, and covered by automated
and experimental regression suites.

## First implementation slice

Do not begin with a whole-wing CFD solver. The first mergeable vertical slice
should be:

1. `scene-v2` schema plus a deterministic exporter for one representative
   design;
2. scene validator and stable-ID round-trip tests;
3. a new Qt-free structural adapter that loads v2 and passes canonical XPBD
   membrane, suspension, contact, and analytic-load tests;
4. generic surface coupling view plus a synthetic, analytically known pressure
   jump;
5. conservative traction transfer tests for total force, moment, and work;
6. a minimal worker protocol that advances and checkpoints this structural
   case;
7. the standalone `simwing-viewer`, automatically opened for interactive runs,
   showing the live structure, loads, contact, residuals, and conservation
   state while also supporting trace recording/replay.

That slice attacks the largest integration risks—geometry truth, ownership,
load transfer, rollback, and process isolation—without pretending a CFD kernel
already exists. It also remains useful if the initial fluid-backend choice
changes.

## Verification strategy

Every accepted simulation records four independent ledgers:

- fluid mass and momentum by region and at domain boundaries;
- pressure, viscous, porous, gravity, line-drag, and contact impulses;
- kinetic, gravitational, structural, damping, and pressure-volume work;
- interface force, moment, and power transfer residuals.

Verification is layered:

1. unit and manufactured-solution tests;
2. canonical fluid, structure, contact, and coupling benchmarks;
3. whole-wing numerical convergence and conservation studies;
4. experimental validation with data not used in calibration;
5. scenario ensembles and sensitivity/uncertainty analysis.

Minimum physical validation data:

- biaxial and bias fabric tests, bending/folding, porosity, and seam leakage;
- line load-extension, line mass, and drag;
- exact inflated geometry or photogrammetry at known pressure;
- internal cell pressure taps and intake flow estimates;
- wind-tunnel pressure maps and six-component loads over angle and brake input;
- synchronized inflation/collapse video, pilot/wing trajectory, and line loads;
- instrumented flight airspeed, attitude, rates, acceleration, control input,
  and atmospheric conditions.

Calibration and validation datasets must be separated. Report parameter
sensitivity and confidence intervals rather than tuning one scenario until it
looks right.

## Principal risks and decision gates

| Risk | Early test or mitigation |
|---|---|
| Added-mass instability | Rollback and strong coupling on a light flexible-piston case before a wing. |
| Pressure leakage across thin/folded fabric | Two-sided pressure-jump and folded multi-interface cells are Phase 3 blockers. |
| Unresolved contact gaps | Couple contact state to a conservative seal/reopen model and refine locally. |
| Excessive full-domain memory | AMR from the first 3D kernel; record bytes per active cell and per block. |
| Nonphysical XPBD material response | Calibrate orthotropic, bending, seam, and damping parameters before absolute claims. |
| Double-counted aerodynamics | Exactly one fluid owner per mode; runtime assertions reject incompatible force backends. |
| GUI fragility or solver crash | Headless checkpointable worker with a versioned protocol. |
| A plausible but wrong whole-wing result | Canonical verification and conservation gates precede every fidelity increase. |
| Chaotic collapse comparisons | Use controlled disturbances, ensembles, and distributional metrics. |
| Framework lock-in | Make scene, structure, transfer tests, and protocol independent of the grid substrate. |

## External technical references

- Kolahdouz et al., [An immersed boundary method for fluid-structure
  interaction with thin, highly compliant shell structures](https://doi.org/10.1016/j.jcp.2021.110369),
  *Journal of Computational Physics*, 2021.
- Huang et al., [An immersed interface methodology for simulating supersonic
  spacecraft parachutes](https://doi.org/10.1016/j.jfluidstructs.2022.103742),
  including a Darcy-Forchheimer porous-fabric jump model, 2022.
- Qin, Kolahdouz, and Griffith, [An immersed interface method for incompressible
  flows and geometries represented by discrete surfaces](https://arxiv.org/abs/2003.12186),
  demonstrating sharp pressure-jump treatment in an interface/LBM setting.
- [AMReX documentation](https://amrex-codes.github.io/amrex/docs_html/) and its
  [embedded-boundary limitations](https://amrex-codes.github.io/amrex/docs_html/EB_Chapter.html).
- [IBAMR](https://ibamr.github.io/) as an established adaptive immersed-boundary
  research framework and comparison implementation.
- [preCICE](https://precice.org/) for mature partitioned-coupling algorithms and
  optional external-solver validation workflows.

These references support method selection; none is a turnkey paraglider solver.
The architecture remains conditional on the Phase 3 canonical gates.
