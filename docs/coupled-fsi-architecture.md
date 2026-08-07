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

This document is the target architecture.
`docs/legacy/leparagliding/CONTINUE.md` remains the record of the imported
Playground implementation and its measured regression guards.

## Imported baseline

The initial SimWing tree was copied from the working tree at
`C:\CODE\LeParagliding` on 2026-08-07:

- Git commit: `497b5632154fb09ef14c4b47943d3fa31a184309`.
- Imported tracked and non-ignored files: 1,182.
- Excluded: `.git`, generated build trees, and ignored generated outputs.
- Preserved alongside the import: `ChatGPT.txt`.
- The source working tree was clean for tracked files at final synchronization;
  its non-ignored design and fixture files were included.
- SimWing has no Git metadata yet. Repository initialization and any history
  relationship to LEparagliding are an explicit project decision, not an
  accidental side effect of the copy.

Before implementation starts, freeze this baseline in version control and run
the existing test suite. The source repository changed during the import, so
future transfers must be deliberate commits or patches rather than live-folder
synchronization.

One inherited inconsistency must be resolved before changing solver arithmetic:
current code and some documents use 96 targeted free-flight cable sweep pairs,
while older handoff text describes three. Record the executable baseline and
update the stale document; do not silently choose either value.

### Initial verification snapshot

The imported tree was configured on Windows with CMake 4.1.1, Visual Studio
2022/MSVC 19.44, Qt 6.11.1, and Open CASCADE 8.0. The ten focused Release
targets for Playground/soft-wing physics all compiled successfully.

All nine focused CTest cases passed in 0.24 seconds:
`playground_metrics`, `playground_pressure_solve`, `playground_cell_air`,
`softwing_material`, `softwing_cell_volume`, `playground_cells`,
`playground_material`, `playground_contact`, and
`playground_contact_integration`.

A first run appeared to hang because the narrowly built executables could not
find `Qt6Core.dll`; Windows displayed a modal loader-error dialog outside
CTest's output. Prepending `C:\Qt\6.11.1\msvc2022_64\bin` to `PATH` corrected
the test environment. This was not a simulation failure. Only focused targets
were built, so the remaining configured tests were not executed in this
snapshot.

## Executive decision

Build two complementary simulation modes behind one scene format and one
structural model:

1. **Fast flight mode** retains reduced-order aerodynamics and a lumped cell-air
   network. It supports interactive design work, deterministic regressions, and
   early validation of the structural and suspension model. The unshipped
   SoftWingLab aerodynamic implementation is a useful source, but it must be
   ported and revalidated, not copied wholesale.
2. **Full FSI mode** uses an adaptive Cartesian finite-volume flow solver with a
   projection step and a sharp, two-sided immersed-interface treatment of the
   zero-thickness canopy and ribs. It owns internal and external air and uses
   strong partitioned coupling to XPBD. This is the path to inflation, collapse,
   and reinflation.

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
| `src/softwing` XPBD constraints, contact, pneumatics, and suspension | Evolve behind a structural adapter | It provides the best current large-deformation foundation and deterministic tests. |
| `playground_metrics.*` | Extract and generalize | Shape, load, collapse, and settling metrics are essential validation instruments. |
| Playground prescribed pressure, bounded Cp projection, polar force, and cell pressure stamping | Keep only as legacy/fast-mode backends | CFD must own pressure in full mode; applying both would double-count fluid loads. |
| Playground GUI and monolithic worker | Replace incrementally with Flight Lab controller/view | The CFD worker must be isolated from the Qt GUI and support checkpointed long runs. |
| `lep-sim.json` version 1 | Retain a compatibility reader; supersede for new runs | It lacks authoritative triangles, two-sided region labels, material coordinates, seam topology, and full rib/intake geometry. |
| SoftWingLab reduced-order VLM/wake implementation | Selectively port into fast mode | It contains useful conservative transfer and coupling diagnostics, but its own records show unresolved circulation and planform gates. |
| XFLR5 | Keep as a polar and aerodynamic comparison tool | It is not a dynamic FSI solver. |

The existing deterministic tunnel pressure path remains an oracle until the new
scene and structural adapter reproduce it. It must not constrain the final CFD
physics after equivalent regressions exist.

## Target process architecture

```text
LEparagliding design + airfoils + Studio-only material data
                         |
               child calculation engine
                         |
             exact OCCT model + scene-v2
                         |
         +---------------+----------------+
         |                                |
  Studio / Flight Lab              simwing-fsi worker
  controls, viewer, plots         checkpointable, no Qt UI
         |                                |
         |            +-------------------+-------------------+
         |            |                                       |
         |     structure adapter                       fluid backend
         |     XPBD/contact/lines               adaptive FV + sharp IIM
         |            |                                       |
         |            +------- strong coupling driver --------+
         |                                |
         +--------- versioned frames, metrics, logs -----------+
```

The worker is a separate executable, following the successful calculation
engine boundary. A crash, out-of-memory condition, or long CFD step must not
take down the editor. The same executable must run locally, from a command line,
and later through WSL or a remote scheduler without changing simulation logic.

Use a small versioned control protocol. Commands include load scene, initialize,
advance, pause, checkpoint, resume, cancel, and request fields. Responses include
progress, convergence, metrics, diagnostics, preview frames, and fatal errors.
Do not stream full CFD fields to the GUI every step.

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
and diagonals from sampled curves in the simulator. Retain the coarse v1 writer
only for compatibility while v2 is established.

Scene validation is a standalone library and command-line tool. It rejects:

- missing or inconsistent units and transforms;
- zero-area, duplicate, non-finite, or inconsistently oriented faces;
- invalid side-region assignments or a disconnected/unexpected region graph;
- accidental faces over intakes and missing rib-hole boundaries;
- missing material or mass assignments;
- dangling seams, suspension edges, or attachments;
- duplicate stable IDs and incompatible checkpoint/scene checksums.

## Structural model

The first structural backend wraps the existing `softwing_core`; it does not
rewrite XPBD during CFD development. Refactor the current Playground body
construction into a Qt-free adapter with these responsibilities:

- build nodes, membrane/bending/seam constraints, ribs, diagonals, lines,
  payload, and contact from scene-v2;
- expose a contiguous coupling surface with stable IDs;
- accept time-integrated nodal forces or impulses from conservative traction
  transfer;
- predict and correct interface position and velocity;
- save and restore all state for implicit-coupling iterations;
- expose surface normals, region labels, contact state, cell volume, and
  topology without leaking Playground widget concepts;
- compute structural energy, constraint residuals, penetration, line loads,
  and material strain in every accepted step.

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
- Cell pressure is measured from the fluid solution. The existing 0D gas model
  may provide a fast backend or an explicitly documented unresolved-flow
  closure, but cannot stamp pressure in full CFD mode.
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
    diagnostics.*           conservation, residuals, events, profiling
    checkpoint.*            versioned restart state
    scenarios.*             tunnel, glide, inflation, collapse definitions
    fluid/
        grid.*              block hierarchy and field storage
        geometry.*          discrete surface and region reconstruction
        advection.*
        diffusion.*
        projection.*
        interface_jump.*
        turbulence.*
        boundaries.*

src/gui/
    flight_lab_page.*       controls and results, no solver arithmetic
    flight_lab_controller.* worker lifecycle and protocol
    flight_lab_view.*       sampled structure/flow visualization

tools/
    simwing_fsi_main.cpp    headless worker executable
    simwing_scene_main.cpp  inspect/validate/convert scenes
    simwing_fsi_bench.cpp   canonical deterministic/throughput cases
```

Likely CMake targets are `simwing_scene`, `simwing_structure`, `simwing_fluid`,
`simwing_coupling`, `simwing-fsi`, and focused test executables. Keep Qt out of
the numerical targets. Backend interfaces must not be so abstract that they
hide grid layout or force extra full-field copies.

## Delivery phases and gates

### Phase 0 — Freeze and measure the inherited baseline

Work:

- choose repository history strategy and commit the import;
- build Release and run all inherited tests;
- archive the current deterministic tunnel CSV and named glide/tuck scenarios;
- reconcile the cable-sweep documentation/code discrepancy;
- record representative runtime, memory, shape metrics, line loads, and known
  failures;
- define scenario files instead of relying on UI defaults.

Gate: every inherited regression has an explicit pass/fail result; the baseline
can be reproduced from a clean checkout; no new solver code has landed yet.

### Phase 1 — Scene-v2 and structural extraction

Work:

- extend the exact model capture/export with authoritative simulation
  triangles, material coordinates, regions, rib holes, seams, and attachments;
- implement binary payload, reader, validator, and v1 compatibility adapter;
- extract Playground body construction and stepping into the structural adapter;
- move reusable metrics behind solver-independent views;
- add rollback/checkpoint support and structural conservation diagnostics.

Gate: v2 reproduces the v1 tunnel baseline within declared tolerances, all
scene validation fixtures pass, and a headless structural run has no Qt/GUI
dependency.

### Phase 2 — Calibratable fast Flight Lab

Work:

- port only the reviewed SoftWingLab reduced-order aerodynamic pieces;
- preserve its conservative force/moment/power transfer tests;
- use the 0D cell-air network only in this mode;
- expose backend identity and validation level in every output;
- replace the Playground UI with Flight Lab only after feature and regression
  parity; keep the legacy launch option during migration.

Gate: deterministic tunnel and free-flight cases meet conservation budgets and
measured polar/trim tolerances. Unresolved circulation or planform gates remain
visible failures, not waived checks.

This phase is strongly recommended because it matures the scene, structure,
protocol, metrics, and coupling instrumentation before expensive CFD is added.
It is not a substitute for full inflation CFD.

### Phase 3 — CFD verification kernel

Work:

- evaluate AMReX as grid/AMR/linear-solver infrastructure with a custom sharp
  discrete-surface layer;
- in parallel conceptually, use IBAMR or OpenFOAM/preCICE as an external
  reference for selected canonical cases, not as a required GUI dependency;
- implement projection, moving-interface jump conditions, conservative
  transfer, refinement, and checkpoints;
- verify CPU first, then add GPU kernels behind identical numerical tests.

Required canonical cases:

- Taylor-Green vortex and manufactured divergence/pressure solutions;
- channel/flat-plate flow and a static pressure jump across a membrane;
- flow through a calibrated porous sheet;
- moving piston/membrane with analytic volume and work balance;
- flexible flag and a thin shell with one interface per cell;
- two closely spaced and folded sheets with multiple crossings per cell;
- resolved opening that closes below grid scale and reopens;
- force, moment, and power transfer under rigid translation and rotation.

Gate: observed convergence is consistent with the intended order away from
interface singularities; mass, force, moment, and power errors satisfy written
budgets; folded-interface and unresolved-gap cases do not leak or change fluid
regions spuriously.

An LBM prototype proceeds beyond this phase only if it demonstrates adaptive
refinement, sharp two-sided pressure, folded-interface handling, acceptable
memory, a compatible license, and better measured cost/accuracy than the
projection solver.

### Phase 4 — One-way CFD on an inflated wing

Work:

- run frozen designed and XPBD-equilibrated shapes;
- perform mesh/refinement/domain/time-step studies;
- compare pressure maps, force/moment polars, wake, and intake mass flow with
  wind-tunnel data, XFLR5 where applicable, and an independent CFD reference;
- establish wall/turbulence model uncertainty.

Gate: integrated loads and pressure distributions are grid-converged within a
declared band and agree with validation data closely enough to justify two-way
coupling. A matching global lift coefficient alone is insufficient.

### Phase 5 — Two-way inflated-wing FSI

Work:

- implement rollback-capable Aitken coupling, then IQN-ILS;
- add adaptive macro-step control and coupled residual budgets;
- validate static inflation and steady deformed trim before maneuvers;
- compare structural strain, pressure, cell shape, line loads, and pilot motion.

Gate: results are insensitive to further coupling iterations, structural
substeps, and fluid time-step reduction within written tolerances; no hidden
force path is needed to keep the canopy inflated.

### Phase 6 — Resolved ram-air inflation

Work:

- enable internal/external region connectivity through intakes and crossports;
- start with restrained symmetric inflation, then moving-payload launch cases;
- calibrate permeability and inlet loss only against independent measurements;
- validate fill sequence, cell pressures, projected area, opening shock, line
  loads, and inflation time from synchronized experiment data.

Gate: symmetric inflation repeats across grid/time-step refinement and matches
the validation envelope without artificial pressure injection or geometry
stabilization.

### Phase 7 — Collapse and reinflation

Work:

- add fluid/contact sealing and reconnection events;
- validate controlled asymmetric disturbances in a tunnel configuration;
- progress through frontal collapse, asymmetric collapse, cravat-like contact,
  and recovery inputs;
- use ensembles for chaotic cases and compare distributions, not one trace.

Gate: mass does not cross sealed contact, topology events converge under
refinement, and measured collapse/reopening thresholds and time histories fall
within declared uncertainty bands.

### Phase 8 — Productization and performance

Work:

- GPU-port proven hotspots, add distributed block execution if required;
- provide progressive preview fields and robust cancellation/checkpoint resume;
- create presets for interactive fast mode and offline validation mode;
- add result provenance: commit, scene checksum, backend, grid, time steps,
  tolerances, material set, validation level, and hardware;
- retire the legacy Playground only after saved scenarios and workflows have a
  migration path.

Gate: production cases are restartable, reproducible within the declared
parallelism contract, bounded in memory, diagnosable, and covered by automated
and experimental regression suites.

## First implementation slice

Do not begin with a whole-wing CFD solver. The first mergeable vertical slice
should be:

1. `scene-v2` schema plus a deterministic exporter for one representative
   design;
2. scene validator and stable-ID round-trip tests;
3. Qt-free structural adapter that loads v2 and reproduces the legacy pinned
   tunnel shape/load metrics;
4. generic surface coupling view plus a synthetic, analytically known pressure
   jump;
5. conservative traction transfer tests for total force, moment, and work;
6. a minimal worker protocol that advances and checkpoints this structural
   case;
7. Flight Lab GUI displaying its geometry, convergence, and conservation
   diagnostics.

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
