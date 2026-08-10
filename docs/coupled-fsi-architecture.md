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
    scene_structure.*       authoritative scene-to-XPBD assembly
    scene_fluid_surface.*   compact two-sided surface + accepted motion
    scene_fluid_surface_transfer.* authoritative conservative-load binding
    scene_fluid_quadrature.* unique owned area to stable transfer quadrature
    scene_fluid_pressure_traction.* one-sided pressure to normal traction
    scene_fluid_grid_epoch.* accepted full geometry/quadrature remap
    scene_fluid_opening_cap.* topology-only automatic/authored opening closure
    scene_fluid_opening_quadrature.* accepted cap motion/surface sweep
    scene_fluid_opening_patch.* exact grid-resolved opening polygons
    scene_fluid_opening_face_crossing.* paired transverse cap segments
    scene_fluid_capped_face_partition.* material-plus-cap face areas
    scene_fluid_opening_flux.* read-only relative MAC volume-flow ledger
    scene_fluid_cell_volume.* signed-chain sparse region volumes
    scene_fluid_region_continuity.* consecutive volume/flux compatibility
    scene_fluid_region_connectivity.* pressure components and gauges
    scene_fluid_pressure_control_volume.* sparse volume-weighted unknowns
    scene_fluid_pressure_face_link.* Cartesian and embedded pressure links
    scene_fluid_pressure_operator.* symmetric integrated graph Laplacian
    scene_fluid_pressure_epoch.* atomic accepted pressure geometry/operator
    scene_fluid_pressure_topology_transition.* shared crossing ownership
    scene_fluid_pressure_volume_rate.* consecutive sparse geometry rates
    scene_fluid_pressure_projection.* link-resolved pressure/flow correction
    scene_fluid_pressure_sampling.* gauge-safe surface pressure return path
    scene_fluid_region_momentum.* accepted collocated region momentum
    scene_fluid_region_transport.* conservative region-momentum advance
    scene_fluid_region_rebase.* bounded mapped current-topology rebase
    scene_fluid_region_wall.* two-sided material-wall momentum exchange
    scene_fluid_region_link_flow.* current-link region predictor
    scene_fluid_pressure_coupling.* strong pressure/shear feedback
    scene_pressure_cell_geometry.* shared visible/refinement tetrahedra
    scene_pressure_cell_case.* visible open-cell pressure-feedback canonical
    scene_pressure_cell_operator_phase_audit.h fixed-grid placement spectrum
    scene_pressure_cell_operator_phase_refinement_audit.* phase/refinement matrix
    scene_pressure_cell_operator_refinement_audit.* skew-intake grid spectrum
    scene_fluid_mimetic_region_conductance_audit.* graph-free terminal response
    scene_pressure_cell_mimetic_conductance_phase_refinement_audit.* uncensored shadow matrix
    scene_pressure_cell_mimetic_conductance_convergence_assessment.* typed three-level trend screen
    scene_pressure_cell_checkpoint_persistence.* bounded canonical restart
    transfer.*              conservative traction/motion exchange
    coupling.*              rollback, Aitken, IQN-ILS, acceptance logic
    coupled_state.*         composite rollback and bounded macro-step retries
    piston_case.*           sealed fixed-topology worker canonical
    projected_flag_case.*   fixed-reference projected-gust fabric canonical
    ram_air_cell_case.*     fixed-reference open five-panel cell canonical
    porous_sheet_case.*     coupled midpoint porous-sheet oracle
    porous_sheet_checkpoint_persistence.* bounded composite restart codec
    open_piston_case.*      driven open-control-volume worker canonical
    open_piston_checkpoint_persistence.* bounded composite restart codec
    periodic_flow_case.*    visible periodic CFD verification worker
    pressure_jump_case.*    visible static layered-jump worker oracle
    porous_flow_case.*      visible pressure-driven porous plug oracle
    periodic_flow_checkpoint.cpp bounded persistent worker restart codec
    worker_control_protocol.* transport-neutral safe-point messages
    worker_control_stream.*  bounded self-framing binary stream adapter
    worker_control_session.* case-neutral safe-point command execution
    periodic_flow_control.* typed periodic worker control adapter
    moving_porous_flow_control.* typed moving porous worker control adapter
    open_piston_control.*   typed open-piston worker control adapter
    porous_sheet_control.*  typed coupled porous-sheet control adapter
    diagnostics.*           conservation, residuals, events, profiling
    checkpoint.*            versioned restart state
    scenarios.*             tunnel, glide, inflation, collapse definitions
    fluid/
        grid.*              block hierarchy and field storage
        scene_surface_geometry.* accepted triangle/cell broad phase
        scene_surface_intersection.* exact triangle/cell narrow phase
        scene_surface_clipping.* barycentric per-cell surface polygons
        scene_surface_ownership.* unique internal cell/face area owners
        scene_surface_crossings.* paired oriented internal-face segments
        scene_surface_face_topology.* sparse multi-sheet face-local index
        scene_surface_face_graph.* provenance-keyed face-local connectivity
        scene_surface_face_chains.* winding-directed chains and loops
        scene_surface_face_loops.* oriented simple-loop geometry and sides
        scene_surface_face_partition.* nested-loop per-region face areas
        geometry.*          exact surface and region reconstruction
        advection.*         bounded donor/limited-MC MAC transport
        projected_advection.* projected nonlinear SSPRK2 transport
        diffusion.*         bounded periodic MAC viscosity verification
        evolution.*         transactional flow steps and subcycling
        projection.*
        interface_jump.*
        planar_face_topology.* generic axis-aligned periodic plane epochs
        planar_pressure_jump.* closed moving multi-layer pressure chains
        planar_region_sweep.* bounded two-epoch regional geometry/GCL
        planar_region_flux.* minimum impermeability/relative-flux screen
        planar_region_opening_flow.* bounded opening-graph feasibility oracle
        planar_region_opening_power.* midpoint pressure-power audit
        porous_interface.* calibrated flux-driven porous jump and ledger
        porous_flow.*       midpoint pressure-driven plug and ledgers
        porous_topology.*   compatibility API over planar face topology
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
`simwing_scene_structure`, `simwing_scene_fluid_surface`,
`simwing_scene_fluid_surface_transfer`,
`simwing_scene_fluid_geometry`,
`simwing_scene_fluid_quadrature`,
`simwing_scene_fluid_grid_epoch`,
`simwing_scene_fluid_opening_cap`,
`simwing_scene_fluid_cell_volume`,
`simwing_transfer`, `simwing_fluid`, `simwing_coupling`,
`simwing_coupled_state`,
`simwing_fluid_structure_bridge`, `simwing_piston_case`,
`simwing_projected_flag_case`,
`simwing_ram_air_cell_case`,
`simwing_porous_sheet_case`, `simwing_open_piston_case`, `simwing_fluid_frame`,
`simwing_periodic_flow_case`, `simwing_pressure_jump_case`,
`simwing_porous_flow_case`,
`simwing_worker_control_protocol`,
`simwing_worker_control_stream`,
`simwing_worker_control_session`, `simwing_periodic_flow_control`,
`simwing_open_piston_control`, `simwing_porous_sheet_control`,
`simwing_viewer_geometry`, `simwing-fsi`,
`simwing_viewer_protocol`,
`simwing-viewer`, and focused test executables. Keep Qt out of the numerical
targets; only the viewer/UI targets link it. Backend interfaces must not be so
abstract that they hide grid layout or force extra full-field copies.

## Delivery phases and gates

Current foundation status: the Qt-free scene-v2 core, XPBD structural adapter,
deterministic scene-to-structure assembly, stable-ID-sorted two-sided fluid
surface extraction with accepted Structure-state capture, immutable structural
diagnostic frames, replayable trace protocol, and standalone Qt/OpenGL trace
viewer are implemented with focused tests. The surface adapter preserves
authored winding, side regions, porous material, and opening order, but does
not yet classify grid crossings or cut cells; it rejects opening-only vertices
until their structural motion is defined. Its conservative-transfer binding
validates the scene and Structure identities again before uniform or
face-resolved barycentric traction can reach XPBD loads. The first grid-facing
stage conservatively bins padded current-triangle AABBs into canonical
cell-major candidates bound to the complete accepted-state fingerprint. It
handles internal grid-plane ambiguity by retaining both adjacent cells and
rejects out-of-domain geometry. A normalized separating-axis narrow phase
then removes false positives, retains contact conservatively, and revalidates
the complete expected pair set. Exact pairs are clipped into flattened
barycentric point/segment/area patches with analytic area and centroids.
Coincident shared-plane area remains duplicated and visibly flagged until a
separate ownership stage pairs it into one canonical internal MAC face with
authored side regions and winding sign. Unique cell and face owners then become
stable barycentric quadrature points carrying authored physical identity into
the existing conservative load transfer; shared-plane area is integrated once
while accepted Structure position and velocity are sampled at the identical
stable points for the reciprocal CFD boundary; no traction is invented here. Explicit ordered one-sided CFD
pressure samples now map to one pressure-difference normal traction through
that same conservative path, with exact equal-pressure cancellation and no
polar or pressure stamp. Ordinary owned polygons also
yield boundary segments: independently clipped adjacent-cell copies must agree
inside a fixed machine-roundoff envelope and pair into one canonical lower-cell
transverse crossing with an in-face negative-to-positive direction; unpaired edges
remain contact and coplanar triangle area remains under face ownership. It
rejects unresolved periodic-domain boundary area and grid-edge-aligned crossing
ambiguity. A bounded sparse index groups crossing and coplanar references under
stable grid-bound MAC-face IDs while retaining each sheet separately and making
no false union-coverage claim. A provenance-keyed graph canonically stitches
shared authored edges, retains opening and grid-edge endpoints, and reports
node degree without claiming an open graph is a closed region. Segments then
become deterministic winding-directed open chains or closed loops independently
per authored region pair. Valid multi-region junctions may terminate several
pair-specific chains at one physical higher-degree node; a branch or winding
conflict within one pair rejects. Simple loops gain signed area, centroid, and winding-derived
enclosed/exterior region identity with self-intersection and degenerate-area
rejection. A bounded containment stage rejects touching loops, requires
authored parent/child region continuity, and closes exact per-region areas on
eligible faces. Simple directed open-chain arrangements whose leaves all reach
the rectangular face boundary likewise close exact per-region areas, global
first moments, and centroids while retaining every source chain. A bounded
half-edge traversal enumerates left faces from
authored winding and rejects unstitched crossings or conflicting labels.
Opening-ended chains, coplanar sheets, and boundary-touching loops remain
unresolved.
The implemented stages are now composed by a versioned `SceneFluidGridEpoch`: one
local build carries a single accepted Structure surface through candidates,
exact intersections, clipping, ownership, crossings, face topology, graph,
chains, loops, partitions, and unique quadrature. The result has a complete
chained fingerprint, individual stage limits, and an aggregate owned-payload
byte bound. Validation rejects nested corruption or any foreign surface step.
The moving regression translates an accepted fabric triangle across a MAC
plane and verifies remapped physical IDs, area, force, and moment. A first
closed-region-cycle volume owner consumes that epoch. Every oriented interface
triangle forms a signed tetrahedron against the grid origin; bounded convex
clipping distributes its chain into exact cells, including wholly interior
cells and surfaces whose face-local contours cross tile boundaries. A separate
whole-surface divergence calculation verifies the global region totals. Nested
analytic tetrahedra, a rigid accepted remap, a valid three-region junction,
and a large cavity with 24 full interior cells are regressions. The complete
real-wing region ledger also closes on a centered coarse grid that crosses the
canopy. Material-only pressure-face accounting ignores same-region
internal-sheet chains without discarding their upstream material identity.
Pair-specific open chains at a stitched interior multi-region junction are
assembled together with the rectangular face boundary and resolve to exact
oriented region areas; opening-ended or dangling arrangements remain explicit
unresolved face partitions. Scene-v2.2 can attach one oriented
boundary-vertex cap disk to an opening. A topology-only owner requires one
closed oriented region cycle around every finite-area material-plus-cap edge,
including valid three-region sheet/cap junctions, and derives final winding
from adjacent fabric with the same region pair. Closed directed material loops
may remain cap-free only while both reference and accepted geometry are
collapsed. Without an explicit disk, planar convex loops preserve their
exact fan and planar concave loops use bounded deterministic ear clipping from
reference geometry. An authored disk may be nonplanar and retains each facet's
normal and identity under accepted motion. Folded, intersecting, degenerate,
or unauthored nonplanar caps reject without creating Structure or traction.
The capped open-tetrahedron regression closes analytic volume and exposes its
mouth area. A second immutable opening epoch gives every cap triangle a stable
centroid sample, interpolates the accepted piecewise-linear vertex velocity,
and retains per-point/per-opening normal surface-sweep rates. Its square-mouth
regression closes rigid material-plus-cap sweep exactly. It does not sample the
Eulerian velocity or claim mass flux. The material clipper now exposes its
exact barycentric triangle/box primitive to a bounded opening-patch owner.
Both material and cap clipping derive every lower and upper cell plane directly
as `gridLower + faceIndex * spacing`; they never reach a shared upper plane by
adding one spacing to the preceding lower coordinate. This preserves bit-exact
adjacent-plane identity even when floating-point addition is non-associative.
Positive cap area is retained once per cell, or deduplicated from paired
coincident cell clips into one canonical non-periodic grid face. A square mouth
closes area and sweep both while grid-aligned and after accepted off-face
motion; periodic-boundary image ambiguity rejects. Cell-owned cap polygons
also pair exact positive-length boundary segments from both adjacent cells
into stable, winding-directed transverse face crossings. Face-owned aperture
area is never collapsed into a line, and grid-edge ambiguity remains explicit.
A bounded planar half-edge arrangement now combines these crossings with
directed material chains on every touched Cartesian face. Disconnected signed
cycles are accounted explicitly; a supported arrangement publishes exact
same-region areas, first moments, and centroids, while an unsupported one remains a first-class unresolved
face. All nine of the coarse real wing's cap-crossed faces now close. Exact
coordinates shared by both endpoints of a clipped edge are preserved through
later-axis clipping; this keeps earlier Cartesian-plane identity without a
geometry tolerance and restores the four rib segments previously lost to
one-ulp interpolation drift. The
immutable touched-face record classifies face-owned opening overlap, coplanar
material, invalid source geometry, material- versus opening-owned dangling
endpoints, unstitched intersections, winding/region ambiguity, and area-closure
failure. It also retains the unique failing source stable ID when one exists.
The accepted pressure epoch retains this capped
material-plus-opening partition,
while face-owned aperture area remains separate. A read-only flux epoch now
binds the entire MAC field. Face owners use the exact normal degree of freedom;
cell owners use bounded degree-three quadrature of periodic staggered
interpolation. It retains signed fluid flow, cap sweep, and relative flow per
patch/opening. The same epoch maps those values into equal-and-opposite outward
balances for the negative- and positive-side regions, including cell-to-cell
crossports, and verifies exact global cancellation. Partial-face tiles integrate
analytically, linear off-face flow is recovered, reversed velocity reverses
sign, and co-moving air/mouth motion has zero relative flux.
Unauthored nonplanar, folded, or self-intersecting openings, opening-only cap
vertices, branching junctions, inconsistent winding, and general moving-boundary
fluid equations remain open. A bounded two-epoch continuity owner now binds
consecutive accepted volume and flux products by their exact surface-state and
producer fingerprints. It trapezoidally integrates endpoint outward relative
flow and reports `delta volume + integrated flow` per region. An analytically
driven expanding cell closes with matching intake transport; removing that
transport produces equal-and-opposite local failures despite a zero global
residual. Authored openings now separately form deterministic undirected
pressure-connectivity components, while retaining their directed transport
sign. Components are ordered by their smallest stable region ID, which also
owns the future pressure gauge. Intake and crossport regressions prove the
expected connectivity; a moving pair of sealed cells proves component-wise
incompatibility cannot be hidden by equal-and-opposite global volume change.
Opening source rates also cancel independently inside each component. A
first immutable pressure-control-volume owner now creates one unknown for
every positive sparse cell/region volume. Cell-major indices and hashed
fixed-grid cell/region IDs preserve identity across accepted motion, while
exact clipped first moments provide cell-local pressure centroids and close
independently in every Cartesian cell. Region and component membership provide
exact volume weights and one deterministic gauge volume per component. The
large-tetrahedron regression retains mixed cut cells and 24 full interior
cells, recovers analytic region first moments on original and translated
grids, and closes every cell, region, component, and the full domain. A
complementary immutable face-link
owner now consumes the exact face partitions and connects only matching
same-region pressure volumes. Every Cartesian link retains the exact centroid
of its region subface; aperture complements derive theirs by subtracting exact
opening first moments from the full face. A resolved nested face retains its exact
exterior, annular-cell, and inner-cell areas as separate links. An untouched
face receives one full-area link directly when both adjacent sparse cells have
one unambiguous common region. When those supports overlap, the accepted
material-plus-virtual-cap closed surface may instead prove one region at the
face centroid by oriented solid angle; the classified region must still own a
control in both cells. A transversely cap-crossed face instead consumes
the exact capped partition, superseding the material-only result; unsupported
capped arrangements retain a distinct unresolved status. The analytic open
tetrahedron reaches pressure links with `0.105 m²` of Cell and `0.895 m²` of
Outside area, and the coarse real wing reaches ten resolved partition faces
rather than one. On the centered 4-by-4-by-4 audit, ten formerly unresolved
material faces now close as stitched multi-region arrangements: 15
material-only plus 43 cap-touched faces produce 58 resolved pressure
partitions, with no unresolved pressure-active Cartesian face. The remaining
294 rejected embedded-opening patches are the distinct projected-centroid
stencil limitation described below. Opening-ended, dangling, coplanar, or
otherwise incomplete material arrangements remain explicitly unresolved. A
face-owned authored opening instead contributes an
oriented cross-region link for each exact cap patch and, when unambiguous, one
same-region link over the complementary face area. The analytic triangular
intake closes its unit Cartesian face as `0.18 m²` of Cell-to-Outside aperture
plus `0.82 m²` of Outside-to-Outside area. The published geometry weight is
area over center distance. A cell-owned patch instead publishes an embedded
same-cell cross-region link along the full authored normal. Its center distance
is the positive normal projection between the two exact cell-region centroids;
when that projection is not positive, the patch count and exact area remain
explicitly unresolved and no conductance link is fabricated. The pressure
operator rejects such incomplete topology. For the fully resolved
subset, a bounded immutable pressure-operator owner now expands every link into
the two directed rows of a symmetric integrated graph Laplacian. Its action is
the conservative sum of `area/distance * (p_row - p_neighbour)`, so constants
are exact null modes and quadratic energy is the positive sum over links. The
operator retains each existing gauge control-volume owner but does not apply a
gauge or solve. It additionally proves that each authored pressure component
is exactly one link-connected graph. The face-aligned intake now passes as one
Cell-plus-Outside component and its unit region jump has exactly `0.18 m` of
graph energy. A tilted off-face intake with otherwise complete Cartesian faces
now joins its authored component through the embedded patch instead of becoming
an accidental no-flow boundary. A
transactional component-wise conjugate-gradient solve now accepts an integrated
RHS only when every component sum is inside an explicit absolute tolerance,
removes just that admitted roundoff, and commits only after recomputing the true
residual. The committed result shifts each existing gauge control volume to
exact zero. Three disconnected manufactured regions and the connected
face-aligned intake recover their prescribed pressure fields; incompatible and
truncated solves preserve the warm start exactly. A subsequent bounded
fixed-epoch adapter now reads one predicted spatial MAC flow per same-region
link and the already accepted fluid-minus-cap-sweep flow per oriented
face-aligned or embedded opening patch. Its integrated RHS is `-rho/dt` times
each control
volume's net outward predicted flow. After the existing transactional solve,
every link receives `dt/rho * area/distance * (p_minus - p_plus)` and the result
is published only when the independently accumulated corrected control-volume
rates meet an explicit absolute/relative bound. The open analytic face keeps
its `0.18 m²` aperture and `0.82 m²` complement as separate projected flows.
Rejected solves retain the predicted diagnostic ledger but expose neither
pressure nor corrected flow. A consecutive-epoch pressure-volume-rate owner
now matches stable cell/region pressure IDs and publishes exact geometry
`dV/dt` for every current unknown plus component/global ledgers. A newly
positive row is marked and receives an exact zero-volume previous endpoint. A
disappeared row can retire its complete previous volume to one unique retained
same-region neighbour; missing or ambiguous retirement rejects. One bounded,
versioned topology-transition product pairs retained rows and publishes both
appearance donors and retirement recipients. Geometry-volume rates,
transported region state, and pressure warm starts all consume its fingerprint
instead of deriving three potentially different mappings. The moving
projection overload adds that rate to
the predicted net outward link flow before RHS assembly and accepts only when
`dV/dt + corrected net outward flow` closes locally. Starting from zero air
velocity, the expanding open tetrahedron develops nonzero pressure and draws
flow inward through its authored intake; omitting the rate remains the exact
frozen zero-flow baseline. General conservative swept-volume remap and unique
corrected MAC/region-momentum continuation for partitioned faces and embedded
openings remain absent. A
bounded return-path adapter now samples accepted sparse pressure on the
authoritative material quadrature. Quadrature-v2 retains exact negative- and
positive-side cell owners for both cell-owned and face-owned patches. Each
side resolves to one current cell/region pressure unknown, and the sampler
requires a shared pressure component before forming a gauge-invariant jump.
That ordered one-sided field delegates to the existing conservative
pressure-traction transfer and then the real Structure load accumulator. The
moving open tetra produces a nonzero load with closed force and moment
ledgers; a uniform pressure warm-start shift is bit-identical after gauge
normalization, and sealed independently gauged sides reject explicitly. This
does not add shear, a polar force, or a second aerodynamic load path. A
versioned bounded pressure epoch now composes the grid remap, opening cap and
patch/crossing topology, capped face partitions, sparse volumes, pressure
controls, face links, and graph operator under one accepted Structure-state
fingerprint. Construction uses
local candidates and publishes only the complete fully resolved chain; nested
corruption, a foreign accepted state, and an excessive aggregate payload all
reject. This is the atomic geometry/operator input for the next transactional
coupling owner and still samples no velocity, solves no pressure, and applies
no load. A
strong feedback owner now uses that atomic input in a real
load-based fixed point. Each iteration rewinds Structure to the accepted
macro-step baseline, advances XPBD under the trapezoidal average of the
accepted start pressure and relaxed end-pressure load, rebuilds the moving
pressure epoch, projects `dV/dt + net flow = 0`, and returns the existing
conservative pressure transfer as the next Aitken candidate. Convergence
requires displacement, velocity, and the difference between the actually
applied load guess and newly solved nodal pressure load. Only that converged
physical iterate commits; exhaustion, an unsupported topology rebase,
projection failure, or an exception restores the exact Structure baseline and
leaves accepted
pressure ownership unchanged. The analytic open cell deforms less than its
no-pressure control and repeats its next accepted macro-step bit-for-bit. This
first feedback owner holds the predicted MAC field fixed across nonlinear
iterations. Accepted corrected link flows can now be collapsed onto one
absolute bulk MAC velocity per Cartesian face, restoring oriented cap sweep at
authored openings and reporting the maximum mixed-region subface velocity
spread. Embedded openings have no unique Cartesian face: collapse validates
and reports them without smearing their flow, while the accepted region state
retains their oriented momentum. The pressure-cell worker uses that derived
field as the next accepted
predictor. This collapse is pressure-correction continuation, not general
cut-cell advection, viscosity, or a topology-rebasing CFD step. Its
topology-stable inverse bookkeeping is now explicit as a separate link-flow
continuation product. It restores previous opening-cap sweep, retains the
absolute velocity deviation of every link from its Cartesian face mean,
recentres those deviations under current link areas, reapplies current cap
sweep, and closes the current bulk face total to roundoff. Pressure projection
has a fingerprinted overload that consumes those exact relative link flows.
This path is not used by the strong worker: in a static open-cell regression,
reusing corrected subface flow without convecting it strongly reduces the next
pressure solve, and repeated use drains the physical pressure load. That
negative oracle makes region-resolved momentum transport—not static state
carry—the required next owner. The first input product for that owner is now
explicit: accepted corrected absolute link velocities are reconstructed into
one immutable momentum vector per positive cell/region control volume.
Controls with Cartesian links only retain the exact component-wise area
average. A control incident to an embedded opening instead solves a bounded
three-dimensional normal equation over all incident explicit link normals;
unconstrained directions retain the cell-centred value of the exact MAC
predictor used by the projection. The product is bound to the pressure,
volume, face-link, opening-patch, and predictor fingerprints and reports total
momentum, kinetic energy, fallback coverage, and link-to-collocated
reconstruction residual. The reconstruction itself performs no time advance.
The first region-momentum time advance is also isolated. Corrected relative link flow
uses donor-cell selection to carry the complete collocated momentum vector
between its two cell/region owners, then graph viscosity exchanges
equal-and-opposite impulse across the same fluid connection. A deterministic
substep count bounds each control volume's outgoing-volume Courant number and
explicit viscous row number. Every post-forcing stage checks energy, and the
complete step checks global internal momentum before publishing. Moving
control volumes advance linearly through the projection's accepted geometry
rate, so uniform flow satisfies the discrete GCL exactly. A bound cell-centred
delta between successive bulk MAC predictors can be applied first; its impulse
and work are separate from the conservative internal ledger. The next
topology-stable boundary retains transported
cell/region velocity while recomputing momentum from current physical volume,
projects the endpoint vectors onto each current link's explicit normal, while
preserving the old exact component arithmetic for Cartesian links, and
subtracts the exact current opening-cap sweep. The moving pressure projection
consumes that fingerprinted first-order predictor together with its dV/dt product and
verifies every predicted link flow exactly. A downstream material-wall product
now remaps the immutable transport to each current nonlinear geometry and
exchanges tangential momentum on both sides of every authoritative material
quadrature point. It uses a local half-volume/incident-area distance, bounded
explicit subcycling, and separate wall-work/dissipation ledgers. The fluid
impulse and equal-and-opposite Structure traction close before publication;
normal traction remains pressure-owned. A bounded first one-ring crossing
adapter sits immediately upstream. Its shared transition pairs retained rows,
records current same-region donors for every appearance, and records the unique
previous same-region recipient for every supported disappearance. Retained
controls keep transported velocity; appeared velocity is donor-area weighted;
and a disappeared source transfers its complete volume and momentum to that
recorded recipient. The mapped source ledger closes before current-volume
geometric correction. Pressure warm state consumes the identical transition,
seeding appearances, preserving retained rows, and dropping retired values.
The strong owner composes that result through wall exchange and pressure
projection transactionally. Cross-material donation and missing or ambiguous
one-ring ownership reject; this is not a general swept-volume remap. Pressure
projection then
consumes the wall-adjusted link predictor, and the existing conservative
transfer applies the combined pressure-plus-shear load to XPBD. This is a
local cut-region wall closure, not a resolved immersed-boundary boundary
layer. The
first composite in-memory checkpoint stores Structure and the accepted sparse
pressure projection. Restore uses a temporary Structure, rebuilds the entire
pressure epoch, resamples the validated projection, and reconstructs the
conservative baseline load before committing either owner. Both initial and
accepted checkpoints replay the exact next result in the original or an
equivalent owner; foreign solver settings, corrupt pressure, and a missing
noninitial projection leave the current accepted run untouched. A
selectable `pressure-cell` worker now exercises that owner through the normal
trace/viewer path. Its soft three-panel tetrahedral cell has one face-aligned
triangular intake and three fixed mouth vertices. Two runs are byte-
deterministic, its in-memory checkpoint reproduces the exact next frame, and
the default 600-step/10-second headless run remains topology-stable. The
mechanical apex drive is gone. A uniform correction maintains
a prescribed `-0.85 m/s` periodic mean wind. The existing symmetric SSPRK2
viscosity/projected-nonlinear-advection/viscosity operator advances that bulk
MAC field on a private candidate, and only an accepted scene pressure solve
commits the next predictor and loads the free apex. After bootstrap, accepted
region momentum receives the bulk-MAC delta, advances through moving-volume
GCL transport, exchanges material-wall momentum, and predicts every current
link in each strong iterate. At two seconds the reference transient reports
about `2.12 Pa` and `17.8 mm`; the molecular-viscosity wall reaction is about
`1.4e-6 N`. The coarse state still relaxes by ten seconds to about `0.00309 Pa`
and `0.0260 mm`, demonstrating that this local exchange is not a resolved wake
or boundary layer. Frames expose the mean-flow pump force plus separate
pressure/wall/total-fluid loads, bulk-flow change, final divergence, viscous
loss, region loss/GCL change, region momentum residual, embedded-opening
collapse count, and wall
loss/momentum residual. It is deliberately a visible bootstrap
diagnostic, not aerodynamic truth. The bulk operator's two intermediate
pressure projections have one velocity per Cartesian face and cannot retain
distinct velocities on multiple cut-region links sharing a face; its viscosity
is periodic bulk viscosity. Material-wall viscosity is applied later on the
region-resolved quadrature and does not turn those bulk stages into a general
immersed-boundary operator.
The final scene projection is sparse and region-aware; its transported fluid
state remains fixed while each strong iterate remaps that state to the current
geometry. The area-collapsed pressure-corrected field remains the bound input
to the next bulk-MAC delta, not the cut-region transport state. Its
bounded canonical checkpoint persists the complete accepted sparse pressure
projection beside the trusted nested Structure payload. Deterministic initial
and accepted round trips, exact next-frame replay, CLI autosave/resume, and
transactional corruption/foreign-file rejection are covered. Restore derives
the exact bulk MAC continuation from that projection instead of duplicating it
or the transient bulk pressure on the wire. `SWPCELL10` additionally persists
the accepted region-momentum state, material-wall traction endpoint, and
transport/wall projection provenance, including embedded-opening normal-
equation diagnostics. Decode bounds every momentum and wall
traction record and validates their complete accepted-epoch binding before
publication. A
canonical Qt-free structural
worker now launches the viewer by default and
publishes accepted steps through a bounded
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
periodic-flow, moving-porous-flow, open-piston, strong-piston, and porous-sheet
adapters bind numerical advance, absolute state, and their distinct complete checkpoint
payloads. Each publishes
immutable
accepted frames, delegates checkpoint persistence, and makes stop terminal
without putting output or file policy in the numerical worker. A bounded
self-framing stream adapter now reads the
envelope payload length without a host-native prefix and flushes every response.
`simwing-fsi --control-stdio` binds that stream to binary stdin/stdout for
all five adapters, suppresses viewer launch and textual stdout, and completes
the trace before acknowledging stop. A restored worker reports its checkpoint's
absolute step/time in Ready and writes only newly accepted frames to its trace.
The moving porous regression crosses the second wrap, checkpoints step 101,
verifies exact step-102 continuation, then restores and publishes only that
same next accepted frame.
Because the moving-flow file codec validates stored state by regenerating a
bounded canonical history, batch runs configured to write that checkpoint
fail before trace creation when restored plus requested steps would exceed the
10,000-step replay ceiling.
The porous-sheet stdio regression advances through its first topology rebase,
persists accepted step 330, reproduces step 331 in the original run, and emits
that same single next frame after a restored Ready response. The original run
then reaches the prescribed-pump collision, reports a numerical failure at the
last accepted absolute state, publishes no rejected frame, and still stops
cleanly. A second process persists that terminal state; its restored Ready
response exposes the same absolute step/time, the repeated advance returns the
same numerical failure, and the completed continuation trace remains empty.
The exact-model capture now exports validated scene-v2.2 skins, authored open
intakes with deterministic nonplanar boundary-vertex cap disks, triangulated
holed ribs, internal sheets, explicit suspension junctions, and the uncollapsed
line graph when supplied explicit physical material and pilot settings. Intake
spanwise lips reuse skin vertices and their side chains follow actual rib-mesh
boundary edges; bounded lip canonicalization removes otherwise motionless
opening-only vertices. The real 3.28 regression verifies oriented nondegenerate
cap facets and captures every opening vertex through the live Structure-to-fluid
surface adapter. The cap owner validates the three-region skin/rib mouth cycles
and consumes the complete real-wing opening set. Pairwise cell-volume measures
close its complete region ledger on a centered coarse grid. One simple
boundary-to-boundary interface resolves exactly, while grid-face junctions
remain pair-specific open chains in the material-only partition. Exact
multi-axis clip-plane identity lets the capped arrangement close all nine
cap-touched faces without snapping or fabricating a junction. Opening
quadrature and grid patches
preserve the full cap area, while authored connectivity and sparse pressure
control volumes assemble. Face-link construction now retains every coarse
embedded opening without admissible projected cell-region centroid separation
as a typed rejection with its patch/opening identity, exact signed
centroid-to-cap-plane distances, count, and area, without publishing a
conductance link. All 24 patches across the two mirrored intake openings are
non-positive in the current coarse fixture: both pressure centroids lie on the
positive side of each local cap plane. This is a non-admissible two-point
stencil, not a near-zero tolerance case, so taking an absolute value is not an
accepted closure. Each rejection also retains a bounded deterministic range of
every same-region Cartesian one-ring neighbor reached from its two side
controls. The provenance binds the exact Cartesian link and root/donor
controls, stores the periodic-image donor offset from the cap centroid, and
classifies signed sidedness without changing operator arithmetic. All 24
coarse patches, and all 294 patches in the refined 4-by-4-by-4 audit, have at
least one correctly sided neighbor on both authored sides. Direct donor
replacement would nevertheless move aperture flux into a neighboring
control-volume row; whole-ring or pairwise volume agglomeration also does not
restore admissible geometry consistently. Pressure-operator assembly therefore
continues to reject the incomplete topology. Resolving it requires a derived
locally conservative symmetric multipoint or hybrid formulation.

The first isolated mixed-hybrid mimetic kernel now exists under
`src/fsi/fluid/mimetic_local_cell.*`; it is not wired into that operator. For
one three-dimensional cell it defines rows
`R_f = area_f (x_f - x_cell)^T` and `N_f = n_f^T`, requires the exact closed-cell
identities `sum(area_f n_f) = 0` and `N^T R = volume I`, and constructs

```text
W = N (N^T R)^-1 N^T + gamma [I - R (R^T R)^-1 R^T],
gamma = trace(N (N^T R)^-1 N^T) / (face_count - 3).
```

`W` is symmetric positive definite and satisfies `W R = N`. The local relation
is `u = -W diag(area) (lambda - p_cell 1)`; eliminating `p_cell` against
`sum(area_f u_f) = source` supplies a conservative trace residual. Manufactured
tests prove translation invariance, exact linear flux on a skew tetrahedron,
roundoff-only source closure, and exact Cartesian reduction: condensing the two
traces on a shared orthogonal face recovers `area / centre_distance`. The
operator stores the exact formula as seven doubles per half-face plus fixed
three-by-three factors and applies it matrix-free; no dense local `n x n`
matrix is retained.

`src/fsi/fluid/mimetic_wall_condensation.*` exactly eliminates a selected
material-wall block after cell-scalar condensation. With `A = diag(area)`,
`g = A W area`, and `d = area^T W area`, the local trace matrix is

```text
H = gamma A^2 + U K U^T,
U = [A N, A R, g],
K = blockdiag((N^T R)^-1, -gamma (R^T R)^-1, -1/d).
```

The wall principal block is positive definite whenever at least one active
trace remains. Its inverse normally uses one equilibrated `7 x 7` Woodbury
core and linear wall-face work. When that auxiliary core is numerically
singular but a wall block of at most eight faces is independently invertible,
the kernel reconstructs the bounded principal block, retains its exact
seven-mode Schur metric, and rebuilds the small solve from the immutable local
operator on demand; no dense wall matrix is stored. The kernel
publishes the exact active Schur action and diagonal, condenses a full local
right-hand side, and reconstructs the eliminated wall traces. With
`Q = U_w^T H_ww^-1 U_w`, repeated action is fused directly as
`D_a + U_a (K - K Q K) U_a^T`; it no longer evaluates two full local balances
and a wall solve. Independent dense tetrahedral oracles verify all four
operations, the active constant null mode, the no-wall identity path, bounds,
and fingerprinted corruption rejection. A four-face refinement sliver with
two walls is the direct-fallback regression. An all-wall cell rejects because
its constant mode requires a global gauge rather than an invertible wall block.

`src/fsi/scene_fluid_mimetic_condensed_trace_system.*` composes those local
Schur products into one immutable global field containing only shared
Cartesian and authored-opening traces. It retains a full-to-reduced mapping,
one deterministic shared gauge per pressure component, exact summed condensed
diagonals, and every nested local fingerprint. Matrix-free application gathers
shared values into each cell and sums the active local Schur rows. A separate
RHS path starts from the full trace source and adds each local wall correction;
post-solve reconstruction recovers each unique wall trace from its owning
cell. Manufactured global fields prove symmetric positive-semidefinite action,
component null modes, diagonal agreement, RHS equivalence, full-trace recovery,
and closure of every original shared and wall equation.

`src/fsi/scene_fluid_mimetic_trace_system.*` now assembles the first global
audit operator over those local kernels. Stable Cartesian and
authored-opening identities create exactly one two-incidence shared trace;
each impermeable material half-face creates a unique one-incidence wall trace
whose equation is zero normal flux. The cell scalar is eliminated by exact
integrated conservation during every matrix-free application. The published
positive operator is symmetric positive semidefinite, retains one exact
component-constant null mode, stores a deterministic trace gauge and Jacobi
diagonal per component/row, and derives a compatible trace right-hand side
from integrated cell sources. All topology, local kernels, limits, and nested
fingerprints remain immutable and independently validated. A pressure
component whose controls are disconnected through shared traces rejects rather
than retaining more null modes than its one authored gauge. No global matrix is
stored.

`src/fsi/scene_fluid_mimetic_trace_solve.*` supplies the first bounded reduced
solve without changing the production pressure path. One solver core accepts
either the complete trace system or the globally material-wall-condensed shared
system. It fixes each retained trace gauge exactly, admits only a declared
component-sum roundoff defect, closes that defect deterministically, and applies
Jacobi-preconditioned conjugate gradients using the stored diagonal. Candidate
traces stay private until an explicitly recomputed residual converges.
Incompatible sources, non-finite arithmetic, or iteration exhaustion leave the
caller's warm start bit-for-bit unchanged. Manufactured reduced fields recover
their gauge-normalized solution; exact wall reconstruction then closes every
row of the original full operator. Compatible cell sources likewise close
local conservation plus all shared and material-wall trace equations. The
public action validates its immutable product on every independent call; the
solver validates once and then uses a private assuming-validated action during
PCG, so integrity scans are not multiplied by the iteration count.

`src/fsi/scene_fluid_mimetic_pressure_solve.*` composes the first atomic
integrated-source pressure transaction. It builds the full trace RHS, condenses
it to shared traces, runs the reduced solve, reconstructs every material-wall
trace, and reevaluates cell scalars, half-face fluxes, local conservation, and
all original trace rows. Reduced fields are only an internal warm-start copy;
failure or a reconstructed RMS residual outside the larger declared solver and
component-compatibility tolerance returns diagnostics with no published state.
The coarse real wing accepts a balanced `+0.02/-0.02` cell-source pair after
307 iterations: reconstructed full residual is `6.17e-9` RMS and `2.59e-6`
maximum, with `5.05e-10` maximum local conservation residual. The transaction
results remain bound to the full, condensed, and optional source fingerprints.

`src/fsi/scene_fluid_mimetic_pressure_state.*` now supplies the accepted-state
ownership boundary. It accepts only a source-bound successful atomic solve,
rebuilds the full trace RHS, independently reconstructs every eliminated wall
trace, and reevaluates all cell scalars and fluxes before requiring exact
agreement with the published result. Shared gauges must already be exactly
zero. Only then are stable control pressures and reduced shared-trace pressures
copied into a bounded fingerprinted state carrying control, full-system,
condensed-system, source, structure, and accepted-epoch provenance. Corruption,
a raw-source solve, a foreign topology, or an unnormalized gauge cannot become
persistent state. The 191,579-row real balanced-source solve crosses this
capture boundary as well.

`src/fsi/scene_fluid_mimetic_pressure_state_persistence.*` makes that endpoint
restartable without weakening its ownership boundary. The `SWMP` envelope is
bounded, checksummed with FNV-64, protocol/versioned, little-endian, and stores
the complete stable control and shared-trace rows plus every state/source/
topology fingerprint and derived summary. Decode is given the trusted rebuilt
control/full/condensed topology, checks counts and identity before publication,
then invokes the full accepted-state validator so stable rows, zero gauges,
epoch, summaries, and fingerprint must all agree. Repeated serialization and
decode/re-encode are byte-identical. Bad magic/version/reserved bits, payload
corruption, truncation, trailing data, byte/record limits, and foreign topology
all preserve the caller's prior state. The coarse real-wing 42,927-shared-trace
endpoint round-trips exactly through the same codec.

`src/fsi/scene_fluid_mimetic_pressure_warm_start.*` now owns the next
consecutive-epoch lifecycle stage. It accepts one immutable pressure state,
the previous/current full and condensed mimetic topologies, and the same
fingerprinted pressure-topology transition already used for moving volume and
transported momentum. Retained control pressures follow stable identity;
appeared controls use the transition's area-weighted retained same-region
donors, and retired controls vanish. A retained shared trace first keeps its
old value, while a genuinely new trace starts at the arithmetic mean of its
rebased endpoint-control pressures. The complete current component is then
shifted by its deterministic gauge value, preserving pressure differences and
making every published gauge exactly zero. Count/working/publication limits,
complete previous/current provenance, explicit trace appearance/retirement
counts and gauge shifts, and fingerprinted validation keep this initialization
separate from accepted state. A moving analytic epoch with one appeared
control exercises the new-trace path and feeds the warm result directly to the
atomic solver. The opt-in pressure-cell shadow now consumes it after bootstrap;
default production selection remains unchanged.

`src/fsi/scene_fluid_mimetic_pressure_sampling.*` crosses the accepted-state
boundary into the existing conservative load path. Each material quadrature
side resolves to its exact cell/region pressure-control row, whose stable
identity and index must agree across the pressure-volume, mimetic-control, and
accepted-state products. Two sheet sides may form a physical pressure jump
only when they share one pressure component; independently gauged absolute
values reject. The bounded fingerprinted output retains both one-sided values,
their exact difference, all source topology/state provenance, and the complete
quadrature binding. It then calls the unchanged pressure-traction quadrature,
so there is no new force path. The moving analytic fixture closes conservative
force and moment transfer, and the coarse real-wing accepted state now samples
every material quadrature point through the same boundary. The opt-in shadow
retains these samples without applying them; the production worker still
selects the graph pressure operator.

`src/fsi/scene_fluid_mimetic_pressure_epoch.*` composes the isolated acceptance
path into one transaction. It starts from an already bounded, fingerprinted
physical source product so predictor-flow and consecutive `dV/dt` ownership
remain visible upstream. Bootstrap uses a bounded zero reduced field; the
moving overload derives the reduced field from the prior accepted state and
the shared topology transition. The transaction then runs the atomic
condensed solve, independently captures accepted pressure state, and samples
all material sides. Only that complete state-plus-sample pair is published.
An incompatible or deliberately exhausted solve returns its diagnostics with
both accepted payloads empty. Analytic regressions require the bootstrap result
to equal the independently assembled state/samples exactly, require the moving
result to reproduce the explicit warm-remap path across a control appearance,
and require nested corruption and bootstrap allocation limits to reject. This
is the rollback-safe boundary consumed by the opt-in shadow experiment below;
the default worker remains on the graph operator.

`src/fsi/scene_fluid_mimetic_pressure_audit.*` owns the complete shadow
endpoint needed to cross that boundary without leaking transient candidates.
One fingerprinted immutable product retains the rebuilt control shells, full
and wall-condensed trace systems, fixed-MAC or transported-wall trace-flow
prediction, physical `rho/dt` source, complete nested solve diagnostics,
accepted pressure state, and material pressure samples. Its aggregate dynamic
storage is independently bounded, and corruption of a nested solve diagnostic
is detected in addition to the constituent products' own fingerprints. The
fixed pre-operator entry point crosses the coarse real wing without first
building the rejected graph operator: it reproduces the independently built
138 controls, 42,927 shared traces, predictor/source rows, and every material
sample exactly. This is significant because some real embedded openings are
geometrically valid mimetic traces even though their old two-point graph
coefficient is inadmissible.

`SceneFluidPressureCoupling` can now construct this endpoint after—and only
after—the normal graph strong iteration converges. The first audited step uses
the same immutable MAC/opening predictor with a bounded zero trace warm field;
later steps use the accepted transported material-wall predictor, the shared
moving-volume topology transition, and the prior audit endpoint's consecutive
warm remap. The candidate remains private until its pressure state and every
material sample accept. A limit or solve failure rewinds the exact Structure
baseline and leaves the graph owner unchanged. Success stores the shadow but
does not apply its pressure samples or alter graph pressure, loads, stepping,
or frames. Four analytic worker steps prove default and audited frames plus
the graph-owned portion of their composite checkpoints are byte-identical
while the shadow changes from fixed bootstrap to transported-wall consecutive
mode.

`src/fsi/scene_fluid_pressure_shadow_comparison.*` makes the next transition
criterion explicit without changing selection. One bounded fingerprinted
product pairs every reference/graph and shadow/mimetic pressure jump by exact
material sample plus cell/region identity, evaluates both through the existing
conservative quadrature transfer, and retains every sample and nodal-force
delta. Diagnostics report reference/shadow transfer closure, pressure
L2/RMS/maximum and relative delta, nodal L2/maximum delta, net force and moment
vectors/norms/relative deltas, and power difference. Corruption or count/byte
failure prevents the endpoint transaction from committing. An independent
coarse real-wing comparison covers all 74,326 material samples, closes both
force and moment transfers below `1e-8`, and has an exact zero-delta self
oracle. The live pressure-cell disagreement is not small: relative pressure
and force delta are about `0.605` at step 4 and `0.601` after 600 steps. This
is evidence against enabling mimetic loads yet, not a tolerance to waive.

`src/fsi/scene_fluid_pressure_owner_transition.*` now makes the selection
decision itself immutable without changing the selected production field. A
fingerprinted policy requires source comparison by default and independently
bounds source-row deltas, material pressure magnitude and fitted shape, nodal
force magnitude and fitted shape, net force, moment, power, and both
conservative-transfer closures. The decision retains the comparison and policy
fingerprints, selected owner, rejection count, and one stable bit per failed
criterion. Any bit selects the reference graph owner; only a zero-bit decision
names the mimetic candidate. The 74,326-sample real-wing exact self-comparison
passes the default numerical thresholds when absence of source rows is
explicitly allowed. The live pressure-cell fails pressure-difference,
pressure-scale, nodal-force-scale, and net-force checks while its source rows
pass, so the typed result retains graph loads. This closes the software
decision boundary but does not turn graph agreement into a physical continuum
oracle. `SceneFluidPressureCoupling` now constructs and owns this decision in
the same accepted-only transaction as the shadow endpoint and comparison,
publishes it through step diagnostics and the pressure-cell accessor, and
clears it with the other transient diagnostics on checkpoint restore. The
headless audit line reports `owner=graph, owner-rejections=0xeb00` at step 4;
those bits identify pressure magnitude/scale, nodal-force scale, net force,
moment, and power. After 600 steps it reports `0x6b00`: the transient power
delta has fallen below policy while the other five blockers remain. A
deliberately permissive policy produces a zero-rejection mimetic candidate
decision while the resulting frame remains byte-identical to the graph-only
worker. No worker consumes the decision to apply the candidate field.

The comparison also retains every upstream control source row. On the live
cell, graph and mimetic geometry rates are exact and their predicted-flow,
continuity, and integrated-source vectors agree to summation roundoff: relative
integrated-source differences are `5.29e-16` at step 4 and `1.66e-16` after
600 steps. A gauge-safe least-squares decomposition then shows that the shadow
pressure jumps and conservative nodal loads are almost pure rescalings of the
graph result. The shadow gains are `2.53035` and `2.50693`, respectively, and
the relative residual after removing that one gain is around `1e-16` to
`2e-16`. Thus the present discrepancy is downstream of forcing and upstream
of sampling/transfer: it is the response of the two spatial pressure
operators on the cut-cell geometry. The Cartesian mimetic kernel already has
an exact area-over-centre-distance oracle; a richer cut-cell manufactured-mode
comparison is required before deciding whether the graph two-point response
or the mixed-hybrid response should own production loads.

`src/fsi/scene_fluid_pressure_operator_response_audit.*` provides that next
offline layer without entering the worker loop. A bounded fingerprinted audit
uses one accepted integrated source when supplied, then constructs six
deterministic component-compatible graph-pressure modes from control centroids,
stable IDs, and authored regions. It applies the graph operator to manufacture
each source,
solves both inverse problems independently, removes one arithmetic pressure
gauge per component, and retains every source, graph response, shadow response,
post-fit residual, and gauge-invariant source work. Solver iterations and final
residuals are part of the
fingerprint; count and byte limits reject before publication.

The first 65-control cut-cell spectrum rules out a global scale error. The
accepted source produces a `2.539` best-fit shadow gain with `2.44%` residual
in the complete gauge-aligned control field. The `x`, `y`, `z`, mixed, and
stable-ID modes produce gains `0.9991`, `1.0049`, `1.0072`, `0.9997`, and
`1.0047`; their relative shape residuals are `8.66%`, `1.23%`, `4.44%`,
`3.78%`, and `16.46%`. Thus the surface pressure-jump result looks almost like
a pure scale only because the current symmetric pump excites a special
operator direction. The sixth manufactured pressure is constant within each
authored region, so its graph source is identically zero on all same-region
links and supported only by the intake. It yields `2.56237` shadow gain,
`1.982%` relative residual, and `0.999804` cosine similarity, reproducing the
accepted outlier and localizing it to authored-opening/interface conductance.
Because this topology has one pressure component and two authored regions, the
same records define an energy-equivalent two-terminal measurement without an
arbitrary pressure average: squared integrated inter-region transfer divided by
source work. The graph intake conductance is exactly `0.18 m`; the mixed-hybrid
response is `0.0700820848 m`, giving a graph-to-shadow conductance ratio of
`2.56841674`. That independent ratio closes the fitted pressure gain while
showing that the transition blocker is the opening-adjacent cut-cell closure,
not a missing global coefficient. The existing orthogonal Cartesian local-cell
oracle still recovers area over centre distance exactly. No global correction
is admissible, and neither non-orthogonal opening response is yet declared
physically authoritative.

`src/fsi/scene_pressure_cell_operator_refinement_audit.*` then removes the
visible tetrahedron's face-aligned special position. Its separately checksummed
rest geometry has a deliberately skew, embedded intake whose vertices avoid
the sampled Cartesian planes. Each bounded isotropic resolution rebuilds the
complete surface ownership, pressure graph, mimetic control shells, full and
wall-condensed trace systems, and the six manufactured response modes. The
product retains nested integrity, exact counts, physical intake area, cell
spacing, dimensional conductances, and conductance multiplied by nominal cell
width and divided by intake area.

The first accepted spectrum is:

| Grid | Controls | Full/reduced traces | Graph conductance | Shadow conductance | Graph/shadow | Normalized graph/shadow |
|---|---:|---:|---:|---:|---:|---:|
| `2^3` | 10 | 38 / 26 | `0.306328 m` | `0.0109311 m` | `28.0237` | `2.85745 / 0.101965` |
| `4^3` | 66 | 206 / 194 | `2.58639 m` | `0.0768124 m` | `33.6715` | `12.0630 / 0.358255` |
| `8^3` | 527 | 1625 / 1565 | `1.90640 m` | `0.215905 m` | `8.82980` | `4.44575 / 0.503493` |

The shadow sequence is smoother under this refinement, whereas the embedded
two-point graph coefficient is strongly grid-phase sensitive. Three coarse
samples do not establish a continuum limit or make the mixed-hybrid response
physically authoritative. They do rule out using the graph result as the
reference merely because it is the current production path. Live arithmetic
therefore remains unchanged while a stronger manufactured continuum oracle is
developed.

The companion fixed-resolution phase audit holds the `4^3` spacing and skew
tetrahedron fixed, then translates the Cartesian lower corner by every
combination of zero and negative half a cell. Its signed phase range is
`[-0.5, 0.5)` so the complete diagnostic geometry stays inside the `4 m`
domain even on the coarsest supported grid. Each accepted phase owns a full
one-sample refinement/response product; an incomplete graph instead retains a
typed `SceneFluidPressureIncompleteFaceOwnership` record, never an exception
message classification.

| Fixed-`4^3` phase result | Count / range | Mean | Population coefficient of variation |
|---|---:|---:|---:|
| Complete graph + shadow products | `6 / 8` | - | - |
| Incomplete face ownership | `2 / 8` | one and two unresolved embedded intake patches | - |
| Normalized graph conductance | `2.55750` - `13.9854` | `6.25249` | `0.77175` |
| Normalized shadow conductance | `0.136954` - `0.408781` | `0.269088` | `0.34030` |

Both rejected phases have all 192 Cartesian faces assigned (`188 + 4` and
`180 + 12` full/partition faces), but one or two cell-owned embedded intake
patches remain unresolved, so pressure-operator assembly correctly refuses to
treat them as sealed walls. The graph coefficient is more than twice as
variable as the shadow coefficient in this deliberately small ensemble. This
strengthens the grid-placement diagnosis and quantifies topology yield; it is
not a statistical convergence study, does not repair unresolved ownership,
and does not promote the shadow operator into live coupling.

`src/fsi/scene_pressure_cell_operator_phase_refinement_audit.*` repeats the
same ordered eight phases at every strictly increasing isotropic resolution
and owns each complete nested phase product under resolution, aggregate-sample,
grid-cell, response, and byte limits. Its topology fraction is derived from
the typed nested statuses and revalidated rather than accepted as telemetry.

| Grid | Complete paired phases | Conditional graph mean / CV | Conditional shadow mean / CV |
|---|---:|---:|---:|
| `2^3` | `4 / 8` | `3.90249 / 0.28334` | `0.101681 / 0.00236` |
| `4^3` | `6 / 8` | `6.25249 / 0.77175` | `0.269088 / 0.34030` |
| `8^3` | `2 / 8` | `4.51872 / 0.01615` | `0.457532 / 0.10045` |

Every rejection still has complete Cartesian face assignment and only typed
unresolved embedded-opening patches. The non-monotone `50% -> 75% -> 25%`
topology yield is itself a blocker. In particular, the small fine-grid graph
CV is based on only two surviving placements and must not be interpreted as
improved robustness. The conditional shadow mean is monotone, but the paired
audit currently builds the graph before the response comparison, so graph
failure censors the shadow sample. An independent shadow-only phase spectrum
is required before assessing convergence; unresolved graph phases remain
explicit and are never converted into sealed walls.

`src/fsi/scene_fluid_mimetic_region_conductance_audit.*` now provides the
graph-independent terminal primitive needed for that spectrum. On a trusted
one-component/two-region mixed-hybrid topology it pairs every permeable
cross-region trace, rejects material-wall traces, and recognizes both
face-aligned Cartesian and embedded authored-opening representations of the
same authored aperture. It allocates a fixed balanced Neumann transfer at
uniform source per aperture area, solves the condensed system from a zero warm
start, removes the component pressure gauge, and defines conductance as
squared achieved transfer divided by source work. Every opening, source,
gauge-aligned pressure, topology fingerprint, solve diagnostic, and aggregate
is retained in one bounded immutable product. The face-aligned pressure-cell
result is `0.0700820848335194 m`, only `3.6e-14 m` from the earlier shadow
response driven by the graph-manufactured source; the embedded two-point
fixture gives `0.0608388978079475 m`. These are deterministic compatibility
oracles for the diagnostic, not Dirichlet data, a continuum truth claim, or a
new live pressure/load owner.

`src/fsi/scene_pressure_cell_mimetic_conductance_phase_refinement_audit.*`
applies that primitive without ever constructing the graph pressure operator.
The bounded immutable product retains every requested placement, not just a
successful response, and distinguishes accepted nested terminal solves from a
typed local-cell linear-consistency rejection carrying the exact control,
region, residual, and tolerance.

The same boundary owns the eight canonical half-cell placements and an
explicit offline numeric profile. `simwing-mimetic-conductance-audit` exposes
that profile without adding fine-grid samples to normal CTest: callers must
name a resolution in `2..64` and either one phase index or `--all-phases`.
Its stable text report retains fingerprints, topology sizes, opening count,
solve residual, normalized response, and any typed local-cell rejection. It
constructs neither the graph pressure owner nor a structural load path.
All-phase runs build, validate, report, and discard one heavy nested product at
a time, retaining only eight normalized values for deterministic compensated
mean/CV aggregation. A selection containing only local-cell rejections remains
a valid immutable audit level with zero conditional statistics, so a
single-phase diagnostic cannot lose its typed failure merely because no other
phase was accepted.

| Grid | Shadow terminal solves | Normalized range | Conditional mean / CV |
|---|---:|---:|---:|
| `2^3` | `8 / 8` | `0.099348` - `0.101965` | `0.100660 / 0.01043` |
| `4^3` | `8 / 8` | `0.129253` - `0.408781` | `0.240930 / 0.39050` |
| `8^3` | `8 / 8` | `0.293734` - `0.888671` | `0.521248 / 0.38828` |
| `16^3` | `8 / 8` | `0.624544` - `1.141234` | `0.902570 / 0.20104` |
| `32^3` | `8 / 8` | `0.956232` - `1.207995` | `1.091977 / 0.07436` |
| `64^3` | `8 / 8` | `1.215742` - `1.321663` | `1.262883 / 0.03130` |

The opt-in streamed runner completed the full `64^3` ensemble in about 29
minutes on the primary Windows development machine:

| Canonical phase | Opening traces | Normalized response | `32 -> 64` increment | Latest contraction | Latest blocker |
|---:|---:|---:|---:|---:|---|
| `0` | `96` | `1.231094375584` | `0.073353040231` | `0.137572` | none |
| `1` | `92` | `1.287760771417` | `0.161243330268` | `1.266587` | noncontracting |
| `2` | `93` | `1.304594221258` | `0.096599035578` | `0.189079` | none |
| `3` | `98` | `1.230172177719` | `0.136386149416` | `0.440183` | none |
| `4` | `94` | `1.215742233756` | `0.232851769643` | `12.096400` | direction + noncontraction |
| `5` | `96` | `1.289829440326` | `0.223719214104` | `2.977992` | direction + noncontraction |
| `6` | `95` | `1.321662620053` | `0.177122935798` | `0.592195` | none |
| `7` | `98` | `1.222206883469` | `0.265975175579` | `1.558100` | direction + noncontraction |

For `16^3 -> 32^3 -> 64^3`, the mean-increment contraction is `0.902325`,
apparent order only `0.148280`, and the corresponding extrapolation is
`2.841728` with a `0.555593` relative fine gap. Those fail the same aggregate
contraction, apparent-order, and extrapolation-gap thresholds even though CV
contracts again by `0.420907`. Three phases reverse direction and four are
noncontracting; phases 0, 2, 3, and 6 pass both phase screens. Applying the
existing screen arithmetic therefore adds aggregate blockers rather than
removing the phase blockers. The streamed report deliberately is not an
invented replacement for the immutable three-source assessment product, so
the published `8^3/16^3/32^3` strict decision remains the production gate and
stays `InsufficientEvidence`.

Thus the earlier `4/8` and `6/8` coarse graph yields really were censoring the
shadow placement spectrum. The five former fine-grid rejections diagnosed a
geometric precision defect rather than an ill-conditioned factorization:
absolute-coordinate shoelace moments lost precision on tiny Cartesian
subfaces, while grid-edge graph nodes depended on sequential clipping order.
The face partitioners now integrate each polygon in a face-local chart and
translate only the published centroid/first moment; the face graph
canonicalizes a grid-edge node from the authored triangle plane and its two
exact Cartesian planes. All `8/8` fine phases now satisfy the unchanged
`1e-10` algebraic-consistency tolerance, without closure fitting. Extending
the same uncensored ensemble to `16^3` exposed two still smaller open-chain
complements. Their closure area remains the exact full-minus-positive
difference, but their centroids now use the independently integrated reverse
boundary polygon whenever subtraction loses more than a fixed coordinate-ULP
envelope. This reduced the representative sliver's divergence-moment defect
from `6.80e-19` to `2.14e-21 m^3`. Its valid four-face operator then exposed a
numerically singular low-rank wall auxiliary core; the bounded direct wall
fallback above accepts it without changing the normal Woodbury path. All
`8/8` phases now solve at `16^3`. At `32^3`, one real complementary region is
only `2.83e-15 m^3`: it passed the existing decomposition-closure envelope but
was erased when that same tolerance also controlled sparse publication.
Cell-volume closure and publication therefore have separate fingerprinted
settings. The offline audit uses a `1e-16 m^3` absolute / `1e-13` relative
publication envelope, cell-local first-moment accumulation, and bounded
centroid repair for cancellation-scale complements; defaults preserve worker
arithmetic. Two valid four-face local operators then report algebraic errors
`5.52e-10` and `1.31e-10`, so the audit explicitly fingerprints a `1e-9`
algebraic-consistency tolerance while production retains the default `1e-10`.
All `8/8` phases solve at `32^3`; phase CV contracts from `0.38828` to
`0.20104` and then to `0.07436`, and the mean drift contracts on the second
interval. Two intervals are stronger continuum evidence but are still not a
convergence result. Because one authored aperture can become
several embedded traces, this audit's uniform area-weighted Neumann source
also differs from the older graph-manufactured source on multi-opening
placements. The live 600-step trace and audited checkpoint remain
byte-identical, so no production pressure or load arithmetic changes.

`src/fsi/scene_pressure_cell_mimetic_conductance_convergence_assessment.*`
turns the last three complete ensembles into a separate bounded immutable
screen instead of treating lower aggregate CV as convergence. The aggregate
mean increments contract by `0.4967106`, giving apparent order `1.00952`, a
Richardson-style extrapolation of `1.278907`, and a `0.146164` relative
fine-to-extrapolated gap. The fine CV and its `0.369850` contraction also pass
the explicit default policy. The phase trajectories do not: phases 2, 4, 5,
and 7 reverse increment direction, while phases 0, 2, 3, and 6 exceed the
`0.75` contraction limit. Only phase 1 passes both requirements. The product
therefore reports `InsufficientEvidence`, typed rejection mask `0x300`, and
the complete adverse counts/trajectories. Its policy, all three source-audit
fingerprints, diagnostics, outcome, storage, and integrity are bound. Turning
off the two phase requirements can produce a read-only trend candidate while
retaining the adverse counts, but this assessment has no worker or load path.

The phase/refinement audit also owns an optional common translation of its
canonical scene and periodic grid, making coordinate-origin sensitivity observable
without changing relative geometry. Translating the `8^3`, phase
`[0,-0.5,0]` sample by `[256,-512,1024] m` preserves all 524 controls,
full/reduced trace counts, and four opening traces. Its intake-area delta is
`5.93e-14 m^2` and its normalized-conductance delta is `5.59e-12`. This
required geometric pairing to compare redundant clipped barycentrics by
zero/nonzero provenance, canonical face-node checks to include a fixed ULP
envelope, and capped arrangements to enter a local chart only when the
coordinate roundoff envelope exceeds their declared minimum tolerance.
Non-finite translations reject before geometry assembly. The normal-origin
worker arithmetic and its two production hashes remain byte-identical.

The experiment is exposed as `simwing-fsi --case pressure-cell
--mimetic-pressure-audit` and reports control/trace counts plus iteration and
predictor mode separately from the established worker summary. The flag is
off by default. `SWPCELL10` composes its settings fingerprint and compact
accepted `SWMP` state with the graph checkpoint. Decode first restores trusted
Structure geometry, rebuilds the mimetic control/full/condensed topology, and
only then lets the independently bounded `SWMP` codec publish pressure rows.
The rebuilt in-memory warm owner deliberately has no invented predictor/source
diagnostics; the next step consumes it through the normal topology transition
and reproduces the uninterrupted consecutive wall-predicted endpoint exactly.
Initial audited checkpoints carry mode identity without inventing pressure,
cross-mode restore is rejected, and nested corruption or limits leave both the
destination checkpoint and live owners unchanged.

`src/fsi/scene_fluid_mimetic_pressure_source.*` owns the physical-unit source
conversion shared with the production projection convention. For each mimetic
control it accepts predicted net-outward volume flow and optional moving-volume
rate, then stores

```text
continuity = dV/dt + predicted_net_outward_flow
integrated_source = -(density / time_step) * continuity
```

in `Pa*m`. The immutable product is count/byte bounded, fingerprinted to the
complete control-cell epoch, and retains compensated component sums so a
physically incompatible predictor remains diagnosable before the trace solver
rejects it. Analytic tests split one continuity residual across geometry and
flow terms and recover the same atomic solution as a direct integrated source;
the real balanced-source audit enters through this product as well.

`src/fsi/scene_fluid_mimetic_trace_flow.*` now supplies both predictor stages.
Every shared trace receives one positive
MinusOrNegative-to-PlusOrPositive relative volume flow. Exact Cartesian
partitions sample their owning MAC face. Face-owned openings and cell-owned
embedded openings instead consume the accepted patch ledger's
fluid-minus-cap-sweep flux, so an opening never needs the graph operator's
rejected two-point coefficient. Material-wall traces remain exactly absent.
The product is bounded and fingerprinted to the shells, trace system, face
links, opening ledger, MAC field, and accepted epoch. Its oriented flows are
accumulated with compensation into per-control net-outward flow, and the
pressure-source adapter maps the accepted consecutive-epoch `dV/dt` field by
stable control identity while requiring its duration to equal the pressure
timestep. The coarse real wing samples all 42,927 shared traces: 101 Cartesian
traces and 42,826 cell-owned opening traces, including the apertures rejected
by the old two-point graph. The moving overload replaces the fixed MAC fluid
sample with the arithmetic mean of accepted material-wall-adjusted endpoint
region velocities projected onto each exact trace normal, then subtracts the
same accepted oriented cap sweep. It matches every existing graph-link
predictor exactly, retains wall-exchange and density provenance, and continues
to sample an opening when its two-point graph coefficient is deliberately made
inadmissible. The pressure-source adapter requires that density to match its
`rho/dt` conversion. The accepted pressure state and its consecutive-epoch
warm remap now own that lifecycle. The pressure-cell can exercise it as a
read-only opt-in shadow; production load selection remains a future,
separately gated stage.

The first immutable scene adapter now assembles each sparse cell/region's
unwrapped half-face shell from exact Cartesian region subfaces, cell- and
face-owned material quadrature, and embedded or face-aligned opening-cap
patches. Material facets become paired zero-normal-flow boundary half-faces;
opening facets retain shared trace/flux identity even when their old two-point
coefficient is rejected. Source identity, the other control, periodic image,
area, centroid, outward normal, per-cell closure matrices, counts, limits, and
fingerprints remain explicit. Nested, face-aligned-opening, and deliberately
rejected embedded-opening fixtures close every control and build the generic
local SPD kernel.

The coarse real wing now closes all 138 sparse controls. Ten untouched faces
whose adjacent sparse cells contain many authored region IDs are independently
proven Outside by the accepted material-plus-cap closed surface. The result
uses neither dominant volume nor shell closure as a label, and every control
is topology-complete. Separately, all 42,826 cell-owned opening patches become
85,652 paired half-faces, proving that the 24 non-admissible embedded two-point
links are geometrically present. A manual refined `4 x 4 x 4` audit likewise
classifies its six formerly ambiguous untouched faces and closes all 358
controls. The adapter safely omits 240 material quadrature sides whose
corresponding positive cell/region volume does not exist; those are
impermeable-wall sides, not missing shared traces. No authored-opening side is
missing. The refined shell set retains 95,984 paired opening half-faces and at
most 2,947 total half-faces in one control, exposing the future dense-local-
matrix cost without paying it: every coarse real-wing shell builds the compact
linear-storage local factorization. The global coarse audit contains 191,579
trace unknowns and 13,132,336 bytes of compact local factor data, and its full
component-constant matrix-free action is roundoff-null. A bounded gauge-fixed
Jacobi-PCG step now runs across the complete real system, reduces its residual,
and rolls back exactly when deliberately truncated. Every one of the 138
coarse controls now builds the exact local wall Schur data, and the global
adapter eliminates 148,652 wall traces into a separate 42,927-row shared
system. It retains 3,986,602 bytes of linear condensation storage, positive
assembled reduced diagonals, a roundoff-null component-constant action, and
full-system wall reconstruction. The same transactional Jacobi-PCG core now
runs on that reduced real system as well. It rolls back exactly after a
deliberately truncated iteration, and a separate manufactured solve reaches
`1e-5` relative RMS within 300 iterations before reconstructing the complete
191,579-row operator below `2e-4` maximum residual. Further local or multilevel
preconditioning and production integration remain open. The fixed MAC or
transported wall-adjusted predictor plus accepted
GCL volume-rate product now assemble a fully fingerprinted physical source.
The atomic source-bound transaction additionally converges a balanced real
source pair and publishes only after reconstructed full-space closure. Its
accepted state now remaps through one moving pressure-topology transaction to
a bounded exact-gauge shared-trace warm start.
The current graph operator and all worker arithmetic are unchanged.
Scene assembly adds per-sheet bending and preserves the junction graph.
It now orients one pilot's line forest toward its harness
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
subcell pocket creates no pressure or flow. The projected nonlinear SSPRK2
operator applies one immutable jump field at both pressure stages; the composed
first-order, Strang, and bounded subcycled paths retain it through every private
step and retry. Empty-jump overloads delegate bit-exactly to the original APIs,
and foreign topology is rejected before mutation. This is a static prescribed
discontinuity under the existing no-added-energy periodic ledger; moving or
powered jump work still requires the coupled interface energy path. A separate
physically grid-bound face-aligned interface now
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
failure. That layer also now contains the first strong-iteration primitive: a
definition-bound vector Aitken delta-squared state with explicit initial,
minimum, and maximum relaxation factors. It uses a deterministic fixed factor
on the first iteration, retains the prior factor for numerically degenerate
residual changes, clips finite dynamic updates, and checkpoints the exact prior
residual/factor/iteration state for rollback. Failed updates and foreign
checkpoint restores are transactional. No worker repeats the fluid/structure
solve around this primitive yet. The adjacent convergence decision is likewise
worker-independent: displacement, velocity, and traction must each satisfy
both explicit absolute and floor-stabilized relative budgets after the minimum
iteration count. An unconverged maximum iteration is reported as exhaustion,
not acceptance. A macro-step iteration owner composes the convergence decision,
Aitken history, and current relaxed interface into one checkpointable algorithm
state. Advance is transactional, convergence and exhaustion are terminal, and
restore reproduces the exact next relaxation/decision. Fluid and Structure
checkpoints deliberately remain solver-owned and must be composed beside this
state. The first composite owner now does that for one real Structure and one
accepted moving-interface fluid epoch: it preconstructs and validates complete
replacement owners, then atomically commits Structure, fluid, and algorithm
state through no-throw moves. Foreign/corrupt members leave all live owners
unchanged, while a valid restore reproduces the next mutation exactly. The
adjacent macro-step retry policy is likewise definition-bound and
checkpointable. An exhausted iteration first enters an explicit pending state,
so the caller restores that composite baseline before acknowledging and
activating the reduced step. Reduction is bounded by a hard minimum step and
retry count, and convergence, retry, and terminal failure remain distinct.
A one-macro-step owner now enforces that ordering: it accepts only a fresh
iteration state, privately retains the composite baseline, and combines the
three-owner restore with retry activation. Replaying an attempt therefore
starts from identical Structure, fluid, and relaxation state. Repeating a
fixed-point iteration within that attempt uses a separate transactional
solver-only checkpoint: Structure and fluid return to the same macro-step
baseline while the advanced relaxed interface and Aitken history are retained.
A generic Qt-free loop now drives a typed physical-solver callback through
those states. It validates each returned Structure/fluid epoch, advances the
relaxation/convergence owner, rewinds only the solvers between iterations, and
performs full rollback before a reduced retry. Converged output remains live;
terminal retry failure, invalid callback output, or an exception restores the
accepted baseline. Its bounded result records each attempt's terminal retry
decision, terminal iteration result, time step, and solver-run count, retaining
the exhausted large-step evidence after a reduced attempt succeeds. The
strong-piston canonical is the first real fluid/Structure callback; a general
wing callback remains future work.
The topology-bound
macro-step surface computes the decision's
inputs from saved-baseline, previous, and current stable-ID kinematics plus two
immutable nodal traction transfers. Maximum displacement and velocity updates
are referenced to macro-step changes, making the reduction invariant to a
world-origin shift or added bulk velocity; traction uses physical nodal-force
magnitudes. Foreign order, topology, Structure, or non-finite data are rejected
before a norm is returned. The first stable-ID fluid-to-structure bridge is present for the exact
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
mode requires the original Eulerian face geometry exactly and accepts either a
separating region pair or a same-region finite panel with a resolved path around
it. Moving-interface samples retain pressure-only transfer for existing callers
and add an explicit complete-reaction selector. The latter validates pressure,
direct constrained-face impulse, their sum, and the selected aggregate before
mapping; diagnostics keep pressure and selected load ledgers separate. Its
rigid-normal-translation mode instead keeps transverse tiling and material
patches immutable while explicitly separating the current Eulerian grid plane
from the unwrapped physical plane. It accepts only a matching rigid structural
translation and uniform normal velocity, with independent position and velocity
correspondence ledgers, and is verified for X, Y, and Z normals and periodic
grid-plane wrap. Accepted porous traction can cross the same stable-ID planar
mapping without reinterpreting the fluid load as a structural load: the bridge
selects the equal-and-opposite sheet reaction, excludes separately prescribed
pressure sources, and independently closes source-to-structure impulse and
work while carrying porous dissipation as its own energy ledger. Fixed and
rigid-normal translated correspondence are supported; general deforming or
nonplanar porous-sheet mapping is not. The `--case piston` worker now evaluates compatible
start/end face-resolved fluid samples, trapezoidally integrates them, accepts
the impulse through XPBD, and publishes the resulting moving two-triangle
surface and CFD ledger fields to the standalone viewer. A distinct headless
light-piston canonical now exercises the actual strong-coupling driver. Its
`6 kg` plate is lighter than the `28.8 kg` projected fluid, and every iteration
repeats moving-interface projection, face-resolved pressure mapping, temporal
integration, and XPBD acceptance from one physical baseline. Aitken relaxes the
scalar terminal speed, topology-bound motion/traction residuals decide
convergence, and only the converged Structure/fluid state persists. The
`--case strong-piston` worker publishes only those accepted epochs, including
the true iteration/retry counts, interface closure, residuals, traction, and
transferred impulse/work. Its endpoint pressure-force slope is the analytic
`10.8 kg` discrete added mass. With trapezoidal start/end force integration,
the fixed-point speed therefore uses `6 kg + 10.8 kg/2` in its denominator;
the worker publishes and tests both the recovered mass and speed residual. Its
in-memory checkpoint composes only the accepted Structure and fixed-topology
moving-interface fluid epoch, validates their velocity closure before commit,
and replays the exact next coupled result and immutable frame. Iteration and
rejected-attempt state remain deliberately outside that accepted restart
boundary. The deterministic bounded `SWSPCKP1` envelope nests the existing
Structure and moving-interface fluid codecs, revalidates their coupled
canonical owner, and enforces outer checksum plus nested byte limits before
transactional decode. The headless batch worker routes that accepted boundary
through atomic `--checkpoint-in`/`--checkpoint-out` replacement and absolute
`--checkpoint-every` autosaves. A typed control adapter likewise exposes only
accepted macro-step safe points. Its separate-process regression checkpoints
step two, matches the exact third frame in the original trace, and restores a
one-frame continuation with identical bytes. The
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
projection, an isolated fingerprinted SPD mixed-hybrid mimetic local-cell
kernel with exact linear consistency/conservation and Cartesian two-point
equivalence, and fixed-topology face-aligned moving
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
  that layered field replayable without involving Qt. The same immutable field
  now survives both projected SSPRK2 pressure stages and every private Strang
  substep without changing the balanced slab velocity; empty fields preserve the
  legacy arithmetic exactly. The disconnected moving-interface projection now
  accepts that field on unconstrained faces as well: sharp sources close within
  each connected fluid component, constrained velocities remain exact, and an
  overlapping moving/jump face is rejected transactionally. The translating
  sealed-slab canonical reconstructs a second analytic pressure slab inside one
  component without spurious flow. The nonlinear Darcy-Forchheimer iteration
  can now use that combined disconnected solve directly. A translating
  two-region slab closes endpoint and midpoint porous slip on unconstrained
  interior faces while retaining exact wall velocity, final moving-reaction
  diagnostics, separate porous jump/work/dissipation ledgers, deterministic
  replay, empty-topology delegation, and transactional failure. This remains
  fixed-grid coexistence, not a general moving porous-sheet or cut-cell solve.
  The accepted nonlinear wrapper now crosses the planar control-volume and
  physical cut-surface boundaries only as a whole. A porous resolved opening
  matches summed tile flux to swept chamber volume and independently sampled
  GCL transport while retaining the final moving reaction. Unaccepted wrappers
  and mismatched nested projection diagnostics are rejected before either
  downstream ledger is exposed. Accepted moving/porous state now also has an
  immutable in-memory checkpoint with a combined moving, porous, and prescribed
  jump topology fingerprint. It canonicalizes porous definitions and retains
  the pre-projection velocity needed to reconstruct endpoint or midpoint
  constitutive samples and prescribed-source power exactly. Capture and restore
  revalidate the calibrated law, combined crossing ownership, nested moving
  projection, jump/dissipation ledgers, and committed fields. A distinct
  bounded, checksummed, deterministic little-endian codec now nests the complete
  moving-only envelope and adds that provenance, canonical porous definitions,
  prescribed jumps, nonlinear diagnostics, and outer acceptance. Decode is
  transactional, reuses both accepted-state validators, and enforces independent
  byte, field, interface, porous-crossing, pressure-jump, region, and diagnostic
  limits before publishing output. A conservative porous-traction adapter now
  reconstructs every accepted calibrated tile at its subcell position, groups
  faces deterministically by stable surface ID, and exposes equal-and-opposite
  fluid/sheet force and impulse. Fluid pressure work uses constitutive-time
  fluid velocity, sheet work uses authored sheet velocity, and their deficit
  closes against porous dissipation per face, per surface, and globally.
  Separately prescribed pressure sources are deliberately not attributed to
  porous-sheet traction. The stable-ID planar fluid/structure bridge consumes
  that accepted boundary directly, transfers only the sheet reaction through
  the existing overlap quadrature, and requires both force-to-impulse and
  power-to-work closure over the recorded projection time step.
  The regression budgets are:
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
A new thin multi-layer planar topology oracle supplies the next moving-interface
prerequisite without changing the production solver. It canonicalizes stable
layer entities into a closed periodic region chain, expands each complete plane
over every transverse tile on X, Y, or Z, and rigidly translates the whole chain
through retained, adjacent, and periodic-wrapped face epochs. Exact-boundary and
skipped-segment motion, broken or nonperiodic region chains, mixed axes,
duplicate surfaces, and invalid coordinates reject transactionally. When the
canonical `+70/-70 Pa` pocket enters one face segment, both ordered layers and
their distinct fractions remain inspectable but the current dense stencil sums
them to zero. A separate static regional profile partitions the full unwrapped
period into exact layer-bounded volumes, requires the pressure-jump cycle and
every repeated region potential to close exactly, and applies one declared
volume-mean gauge. The canonical pocket then retains its physical `2.4 m^3`
interior volume and 70 Pa regional pressure difference even while both layers
share a face; separated layers produce the same regional profile while the
ordinary sharp projection also recovers its cell pressure. Non-closing cycles,
inconsistent repeated regions, and non-finite gauges reject transactionally.
The profile owns no regional velocity degrees of freedom and does not enter the
production projection. It therefore exposes rather than closes the remaining
regional-subcell momentum/flux requirement; it is not a leakage closure,
deforming folded-surface tracker, or moving folded-interface solve.
A separate bounded `planar_region_sweep.*` product now binds two such profiles
by stable surface, side-region, pressure-jump, order, and one-segment topology
continuity. It derives each layer displacement and velocity, compares every
interval's exact volume change with the independently signed lower/upper
boundary sweep, and aggregates that ledger per region and globally. Rigid
translation has zero regional change; a breathing pocket assigns equal and
opposite interior/exterior changes and closes on X/Y/Z. Both signed periodic
rebases retain identity. Crossed layers, exact-boundary or skipped motion,
foreign identities, and count/region/byte-limit violations reject before a
candidate is published. This is the geometric-conservation half of moving
regional flow only: no Eulerian subcell velocity or relative flux exists yet,
so the no-leakage gate remains open.
`planar_region_flux.*` now quantifies that open gate without selecting a solver.
For each interval it fits the least-squares uniform axial fluid velocity to the
two layer velocities, retains both signed one-sided outward-relative flows, and
checks `delta volume + integrated outward relative flow = 0` independently of
impermeability. Rigid motion has zero slip and zero relative flow, including a
periodic rebase. The canonical breathing pocket closes continuity only with
`0.2 m/s` minimum interface slip and `1.6 m^3/s` total absolute one-sided flow;
the exterior ledger is equal and opposite. All-axis area scaling, input-ledger
revalidation, tolerances, and count/region/byte bounds are explicit. This is a
compatibility screen, not a regional projection or permission to let fluid
cross fabric.
The compatibility artifact carries a deterministic nonzero semantic
fingerprint. Validation reconstructs all interval quantities from the retained
primitive volumes, boundary velocities, duration, area, identities, and
tolerance policy, then rebuilds stable-ID chain continuity, sorted regional
summaries, global aggregates, flags, and owned-storage accounting. Mutated or
foreign diagnostic payloads therefore reject transactionally; this hardens an
offline measurement and does not promote it into accepted fluid state.
`planar_region_opening_flow.*` then asks the complementary topology question:
whether explicit oriented openings can carry the relative flow required by the
same moving regional volumes without attributing it to fabric. Each connected
opening component must balance independently. Compatible components solve a
gauge-fixed area-weighted graph Laplacian, minimizing `sum(flow^2 / area)`;
parallel openings therefore share one normal velocity. The closed breathing
pocket reports two incompatible sealed components even though its global
volume change cancels. A `0.5 m^2` exterior-to-pocket opening instead carries
`1.6 m^3/s` at `3.2 m/s`, and a three-region serial graph closes its zero-net
intermediate region. Orientation, all-axis swept-volume scaling, deterministic
opening-ID ordering, dense-factorization work/storage bounds, semantic
fingerprinting, and source-bound reconstruction are explicit. These analytic
opening links are not yet scene-v2 aperture geometry, a subcell velocity basis,
or permission to wire the result into production pressure arithmetic.
`planar_region_opening_power.*` binds a feasible allocation back to the two
endpoint regional pressure profiles. Each oriented opening retains midpoint
pressure drop and `flow * pressure drop`; positive power follows a passive
high-to-low pressure release, while negative power reports the exact local
external-power deficit. The aggregate independently closes opening pressure
power plus regional midpoint `p * dV/dt`. Canonical inflation is kinematically
feasible through one opening but moves `1.6 m^3/s` up a 70 Pa rise, requiring
`112 W`; deflation releases `112 W`. A serial graph can use `40 W` from one
downhill link against a `24 W` uphill link, so local deficits and the zero net
external deficit remain separate. The audit is immutable, source-bound, and
bounded; it supplies neither dynamic pressure nor a constitutive opening law.
A calibrated Darcy-Forchheimer adapter now samples resolved X/Y/Z MAC normal
velocity relative to each authored sheet tile and emits the corresponding
signed sharp jump. It retains tile area, volume flow, and nonnegative pressure
dissipation, and its exact monotone inverse is stable across the pure-linear,
pure-quadratic, and combined fits. A periodic prescribed-flux plane with an
explicit balancing pressure source preserves a `195 Pa` porous loss without
changing its uniform velocity. A bounded transactional Picard solve now closes
nonuniform porous flux and jump against the full periodic projection at either
the endpoint or temporal midpoint. It reprojects the original predicted field
on every iterate, permits separately prescribed pressure-source jumps, requires
both normal-velocity and constitutive jump residuals, retains tile samples and
dissipation at the selected constitutive time, and rolls back both fields if its
bound is exhausted. The uniform endpoint canonical reaches the analytic
`0.4 m/s` under a `20 Pa` drive, while the nonlinear, orientation-reversed,
moving-sheet-relative, and heterogeneous-tile cases remain deterministic and
divergence-free. Midpoint mode closes pressure work minus porous dissipation to
the kinetic-energy change and reproduces the exact nonlinear scalar plug oracle.
The accepted diagnostic separately sums every porous and prescribed oriented
jump into fluid force/impulse and power/work, evaluated at the selected
constitutive time, while porous dissipation remains its own nonnegative ledger.
The first complete periodic macro-step now consumes those ledgers explicitly:
aggregate momentum subtracts jump impulse and its energy ceiling includes jump
work. Its uniform midpoint driven acceleration and unforced porous decay close
analytically, an exhausted nonlinear solve rolls back advection and diffusion as
well as projection, and empty topology preserves the original arithmetic
bit-for-bit. A complete fixed-grid second-order path now applies an
implicit-midpoint porous half-step on each side of the existing full
viscosity/transport/viscosity Strang step. The symmetric wrapper exactly matches
that three-operator composition, sums both jump impulse/work ledgers, keeps
porous dissipation separate, rolls back earlier operators after a later failure,
and observes second-order temporal refinement for the driven uniform canonical.
Empty topology delegates to the original bulk path without changing its fields
or nested diagnostics. A stage-resolved planar overload accepts complete sheet
definitions at the two porous half-step midpoints. It validates immutable
stable identity, region orientation, and resistance before assembly, then
requires the second physical sample to retain or advance exactly one topology
segment from the first. It also requires physical displacement to match the
trapezoidal integral of the two sampled normal velocities over their `dt/2`
separation, using explicit finite nonnegative absolute/relative tolerances. The
positive-wrap regression crosses wrapped face `3 -> 0` and signed image
`0 -> 1` inside one symmetric macro-step; diagnostics own both unwrapped epochs
and the kinematic residual, while metadata and solver failure remain
transactional. General moving cut-cell topology remains open.
`MovingPorousFlowCase` is the first owning-frame harness for that path. Its
prescribed `0.4 m/s` plane starts at unwrapped `x=3.48 m`, crosses
`face=3,image=0` to `face=0,image=1` between the first pair of porous midpoint
samples, and advances through a second wrap while retaining the accumulated
fluid field. The final constitutive sheet jump and fixed pump are published as
separate quads with both unwrapped stage epochs and the closed kinematic,
dissipation, work, momentum, and energy diagnostics. The Qt-free
`simwing-fsi --case moving-porous-flow` path publishes the accepted frames and
completed trace headlessly. Its immutable in-memory checkpoint binds public
case/grid/step/kinematic/topology metadata to a private owning payload and
replays the initial, ordinary, and second-wrap epochs bit-identically; rejected
restore is transactional. The bounded/checksummed `SWMF` persistence format
stores all fluid fields and sharp crossings explicitly and regenerates the
accepted diagnostic through bounded canonical replay before publication. CLI
checkpoint flags provide additional-step resume, absolute accepted-step
autosave cadence, final-write deduplication, and atomic same-path replacement.
Coupled structure remains open.
A pressure-driven uniform-plug companion uses the exact nonlinear implicit
midpoint to close pressure impulse and the
driving-work/porous-dissipation/kinetic-energy identity on every step. Its
endpoint Strang/SSPRK2 evolution retains the accepted plug and all crossings at
both internal pressure stages while exposing the porous and zero-net-jump
periodic gauge-closure planes in owning frames; the physical driving rise stays
separate in the impulse/work ledger. The case converges to the analytic
`1.74165739 m/s`, `250 Pa` steady state. This closes a one-degree-of-freedom
second-order porous-flow oracle; general moving porous/cut-cell topology remains
open. A separate `--case porous-sheet` canonical now couples the nonuniform
midpoint projection directly to the planar stable-ID sheet-reaction bridge and
XPBD. Its analytic linear-resistance midpoint relation independently closes
fluid and sheet momentum against prescribed pump impulse, and closes their
kinetic-energy changes plus porous dissipation against pump work on every
accepted step. The immutable trace shows the translated two-sided sheet and
keeps pressure jump, sheet impulse, pump work, and porous loss separate. At the
first accepted midpoint sample in each of six consecutive ordinary dual-cell
segments it explicitly rebinds the crossing to the next MAC face while
preserving unwrapped physical sheet and fluid state. The fifth event wraps the
Eulerian face from `7` to `0` and advances the signed periodic image to `1`.
That decision uses a pure topology
selector whose epoch owns axis, wrapped face coordinate, and signed periodic
image. It supports X/Y/Z, one adjacent segment in either direction, and both
periodic wraps while rejecting exact MAC-plane ambiguity and skipped segments.
The adjacent complete-plane assembler consumes that epoch plus unwrapped
physical position, rigid normal velocity, two-sided stable IDs, and calibrated
resistance. It emits exactly one canonically ordered porous crossing per
transverse MAC tile for X, Y, or Z only after validating the complete
definition. Both the static porous-flow oracle and this moving-sheet oracle use
that common assembly boundary; neither keeps a private X-only tiling loop.
Each accepted pre-pump epoch survives persistent checkpoint round trips and
continues bit-identically. A later collision with the next periodic image of
the prescribed pump surface is rejected by the case after selection; this
remains a planar oracle, not a
general moving cut-cell remap or a general strong-coupling solve.
A direction-specific companion instantiates the same boundary with negative
pump pressure, reverses both fluid and XPBD momentum, crosses six negative
rebases including `0 -> 7` into signed image `-1`, and persists that wrapped
epoch with distinct provenance and case fingerprint. It then rejects the next
negative pump image transactionally. This reverse mode is a focused symmetry
oracle; the worker CLI continues to expose the positive canonical.
Its immutable in-memory checkpoint owns the nested Structure state, MAC
velocity, pressure, and last accepted coupled diagnostics. Restore validates
the case identity, exact step/time epoch, rigid sheet state, uniform fluid
state, field energy, and cumulative pump momentum before one transactional
commit; initial and accepted checkpoints reproduce the exact next frame in an
equivalent rebuilt worker.
The distinct `SWPS` persistent envelope reuses the validated Structure codec
and stores the topology version, axis, wrapped face, signed periodic image, all
three MAC velocity components, and pressure explicitly. The in-memory owner and
diagnostic frame carry that same complete topology epoch rather than rebuilding
it from a face index. Motion direction is immutable construction metadata with
its own case fingerprint; there are no runtime-varying controls. Decode
regenerates its coupled diagnostic in the owner's direction across topology
epochs by bounded replay and
requires the decoded Structure and every field sample to match that replay
bit-for-bit. Magic,
protocol, reserved bits, payload size, checksum, nested Structure state, total
bytes, scalar samples, and replay steps are bounded and rejected before the
destination checkpoint changes. The standard worker checkpoint flags now write
and restore this envelope before trace creation. Same-path resume is atomic,
steps are additional, autosaves use absolute accepted-step multiples, and the
final accepted state is not written twice.
The case checkpoint must restore the initial state and a later accepted state,
then reproduce the next frame bit-for-bit in both the original and an equivalent
rebuilt worker. Version, case fingerprint, grid, sample count, step, time, and
payload presence are independently rejected before mutation. Its persistent
encoding must be byte-deterministic across repeats and decode/re-encode, retain
all nested diagnostics exactly, and continue bit-for-bit after either initial
or accepted-state decode. The default untrusted-input limits are `256 MiB`,
5,000,000 scalar samples, and 10,000 replay steps. Bad magic/version/reserved
bits,
truncation, trailing data, checksum corruption, and every configured limit must
leave the destination checkpoint unchanged. The porous-sheet CLI integration
writes through the first topology rebase at step 330, restores it, advances
three additional steps, atomically replaces the same restart file, and decodes
it again without launching Qt. The initial interval of 165 writes at steps 165
and 330; the resumed interval of two writes at absolute step 332 and final step
333. A second process chain writes the ordinary wrapped epoch at step 900
(`face=0`, signed image `1`), resumes three additional steps, and verifies the
same wrapped epoch at step 903 with only new trace frames. Interval mode without
an output path is rejected before a worker or trace is created.
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
The periodic, open-piston, and strong-piston stdio integrations each send
advance-two, checkpoint, advance-one, and stop to a separate worker process.
They require
an exact ready/advanced/checkpointed/
advanced/stopped response stream with no trailing stdout, a completed
three-frame trace, and a step-two checkpoint whose next frame is byte-identical
to trace step three. A second process restores that checkpoint, reports Ready at
step two, advances once, and produces a completed one-frame trace whose step
three is byte-identical to independent replay. The in-process open-piston
adapter additionally advances through the first topology rebase in one command,
delegates that composite checkpoint, and reproduces the next frame exactly.
Control mode is rejected for the structural and weak sealed-piston cases.
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
deterministic translation without leaving the initial viewer scale. Its
headless strong-coupling counterpart instead uses `6 kg`; its deterministic
multi-iteration solution differs from the weak one-pass response and closes
accepted interface/structural velocity before retaining persistent Structure
and fluid fields. The open
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
- flow through a calibrated porous sheet (the prescribed-flux constitutive,
  transactional fixed-grid endpoint/midpoint iteration, and pressure-driven
  uniform-plug subsets are complete; full-flow second-order moving-grid coupling
  remains);
- moving piston/membrane with analytic volume and work balance;
- flexible flag and a thin shell with one interface per cell (the one-way,
  fixed-reference normal-gust flag subset is complete; displaced geometry and
  tangential-flow interaction remain open; the five-panel open-cell subset now
  exercises shared-node multi-surface transfer, but its CFD geometry is also
  fixed);
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
The coupled porous-sheet oracle now also exercises the accepted
fluid-to-structure boundary end to end: a periodic pump drives midpoint flow,
the porous adapter excludes that pump from material traction, the bridge maps
only the sheet reaction, and XPBD receives the same impulse and work. Its six
positive and six negative MAC-face topology rebases, both signed wrap
directions, per-epoch checkpoint replay, and later explicit pump-collision
rejection remain intentionally smaller gates than general moving porous
topology. Checkpoint topology belongs to the accepted
constitutive midpoint rather than the possibly farther endpoint; this permits
the terminal accepted state to persist and reproduce the same collision
transactionally after restore.
General interpolated cut-cell pressure metrics, curved or transversely deforming
correspondence, nonplanar or opening topology events, sealed deforming chambers,
a complete coupled energy ledger for those general cases, and richer
case-specific control messages remain open, so the
full Phase 2 piston gate is not yet closed.

The first flexible CFD-to-XPBD canonical is now available as `--case flag`.
A finite 4-by-4 same-region MAC-face panel sits inside a periodic 12-by-8-by-8
domain. Sinusoidal mean-flow increments create a persistent accelerating
cross-flow; every step projects that field around the stationary panel, then
maps the complete face-resolved constraint reaction into a matching 5-by-5
orthotropic membrane grid. Two fixed node rows encode both position and slope
at the edge clamp, eliminating the otherwise admissible rigid hinge rotation.
Accepted frames expose pressure traction, direct constraint traction, complete
reaction traction, conservative nodal loads, divergence, and normal
displacement. Deterministic regressions require signed load reversal, visible
non-rigid free-edge motion, exact anchors, converged projection, and force and
moment transfer closure. This is deliberately one-way: the displaced membrane
does not alter the CFD surface, the stationary reference has zero interface
power, and the case does not claim moving cut cells, classical tangential flag
aerodynamics, two-way energy closure, or wing flow.

The next fixed-reference canonical is available as `--case ram-cell`. It uses
five 4-by-4 MAC-face surfaces for the back, left, right, bottom, and top of a
one-metre open rectangular cavity. Each stable-ID surface independently clips
to a matching panel and transfers its complete pressure-plus-direct constraint
reaction into one 89-node shared-edge XPBD shell with 160 membranes and 200
sheet-local dihedrals. Two perimeter rows clamp mouth position and slope while
the other nodes deform. The regression requires deterministic accepted frames,
converged projection, exact fixed nodes, local mouth exchange with near-zero
net incompressible flux, visible outward bulging, and per-panel plus aggregate
force/moment closure. This verifies multi-panel load assembly and structural
continuity; it does not model changing cavity volume in the fluid, feed fabric
motion back to the fixed grid, or claim physical ram-air inflation.

The standalone worker also includes a deliberately non-CFD curved structural
canonical: a one-metre soft fabric hemisphere held by three equatorial anchors
and a compliant rim. Intrinsic per-triangle charts, signed rest-shape bending
hinges, and a time-varying four-lobe analytic follower-pressure mode exercise
hundreds of membrane elements and immutable viewer topology without exciting
the free rigid rotation admitted by two positional pins. This `--case
hemisphere` path is a
visual and structural regression only; it must not be cited as curved-interface
fluid coupling or aerodynamic validation.

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

- connect a real fluid/Structure solve to the generic rollback-safe repeated
  macro-step driver, then add IQN-ILS;
- drive adaptive macro-step control from the existing retry policy and coupled
  residual budgets;
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
