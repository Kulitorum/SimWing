# LEparagliding Studio repository guide

This is the repository-level working guide for Codex and other coding agents.
It records where behavior lives, how the pieces communicate, which compatibility
contracts are load-bearing, and how to verify changes. Keep it current when a
change moves ownership or invalidates a command or invariant.

## Project identity

- SimWing is a separate GPL-3.0 project for coupled XPBD-CFD paraglider
  simulation, bootstrapped as an independent root from LEparagliding Studio.
  `UPSTREAM.md` records the source commit without importing upstream history. The
  imported C++20/Qt 6 application is a desktop port of Pere Casellas'
  LEparagliding 3.28 engineering program and retains its application/target
  names during the staged migration.
- `project(... VERSION ...)` in `CMakeLists.txt` is the Studio/release version
  (currently 0.4.1). "3.28" is the calculation core and input-format version;
  it is intentionally a different number.
- Qt Widgets provides the application shell and editors. Open CASCADE (OCCT)
  8.0 owns exact NURBS construction, STEP/XCAF I/O, tessellation, selection,
  and the primary 3D viewport.
- The numerical/drawing oracle remains the translated Fortran core. New code
  wraps it; it does not casually replace or "clean up" it.
- XFLR5 v6.62 and the SoftWingLab-derived XPBD core are vendored. Both have
  local changes and are no longer drop-in upstream copies.

### SimWing remake boundary

- Reuse the isolated XPBD structure primitives in `src/softwing`: membranes,
  constraints, cables, suspension, contact, and their canonical tests.
- Do not treat the inherited Playground simulation as a baseline or migration
  target. Its aerodynamic loads, bounded pressure solve, cell-air network,
  flight dynamics, metrics, scenarios, and results are discarded for SimWing.
- The new FSI path must not link `playground_sim`, `playground_pressure_solve`,
  `playground_cell_air`, or `playground_analysis`, and it does not read
  `lep-sim.json` v1.
- Build scene-v2 directly from authoritative model geometry, assemble a new
  Qt-free structural adapter from XPBD primitives, and let CFD be the sole
  owner of internal/external pressure and aerodynamic traction.
- Archived Playground records explain inherited code only. They are not
  SimWing physics acceptance criteria.
- A standalone `simwing-viewer` is an early diagnostic requirement. Interactive
  runs launch it by default; `--no-viewer` is the explicit CI/headless path.
  Rendering consumes immutable, droppable snapshots and must never block or
  alter solver arithmetic. Closing the viewer cannot stop the worker unless the
  user sends an explicit stop command.

## Start every task here

1. Run `git status --short --branch` and preserve all user changes, including
   untracked design/airfoil files. Generated build directories may coexist with
   a running Studio instance.
2. Read the subsystem-specific design record before changing that subsystem.
   Inherited LEparagliding records are under `docs/legacy/leparagliding/`;
   consult its `CONTINUE.md` for Playground work and `BACKLOG.md` for
   deliberately deferred editor/Print work.
3. Inspect the relevant public header and its test before editing the large
   implementation file. Several `.cpp` files are thousands of lines long but
   have deliberately small boundaries.
4. Do not kill a `LEparagliding.exe` process you did not start. It may contain
   unsaved user work and it locks the link target. See the build-lock note under
   Verification.

## Build, test, and run

Windows is the primary local environment:

```powershell
cmake --preset windows
cmake --build --preset release --parallel
ctest --preset release
```

The preset uses Visual Studio 2022 x64 and writes executables to
`build/bin/Release`. CMake auto-detects the newest compatible
`C:/Qt/6.*/msvc2022_64` and checks the usual local OCCT roots. Override OCCT
with `-DLEP_OCCT_ROOT=C:\path\to\occt` when needed. Required packages are Qt
6.5+ components Core/Widgets/OpenGL/OpenGLWidgets/XML and OCCT 8.0 toolkits.

Useful narrower commands:

```powershell
cmake --build build --config Release --target leparagliding-engine
ctest --test-dir build -C Release --output-on-failure -R "engine|preset"
ctest --test-dir build -C Release --output-on-failure -R "playground"
```

When only narrow test targets have been built in a fresh Windows build tree,
prepend the selected Qt `bin` directory to `PATH` before running CTest. Without
it, Qt-linked test executables open a modal `Qt6Core.dll was not found` loader
dialog that looks like a stalled test. A full deployed GUI build may already
have copied the required runtime next to the executables.

Run the products with:

```powershell
.\build\bin\Release\LEparagliding.exe
.\build\bin\Release\leparagliding-engine.exe <design-file> <output-directory>
.\build\bin\Release\LEparagliding.exe --headless <design-file> <output-directory>
.\build\bin\Release\simwing-fsi.exe [--case structural|hemisphere|flag|ram-cell|pressure-cell|piston|strong-piston|open-piston|periodic-flow|porous-flow|moving-porous-flow|porous-sheet|pressure-jump] [--steps N] [--trace <file>] [--checkpoint-in <file>] [--checkpoint-out <file>] [--checkpoint-every N] [--mimetic-pressure-audit] [--control-stdio] [--no-viewer]
.\build\bin\Release\simwing-viewer.exe [--follow] <trace-file>
```

Engine flags are `--preview`, `--no-construction-curves`, and
`--resource-dir <directory>`. GUI/developer entry points include
`--smoke-test`, `--studio-self-test <design>`, `--validate-presets <dir>`,
`--xflr5`, and `--playground <lep-sim.json>`.

Linux/macOS use a Ninja build and an OCCT 8.0 install prefix; the exact recipe
is in `README.md`. `docs/legacy/leparagliding/CLAUDE.md` additionally records
the inherited machine's WSL setup, case-sensitive verification, release
procedure, and signed Windows installer.

## System map

```text
design.txt + relative airfoil files
        |
        +--> Studio: DesignDocument -> per-section editors
        |       |                         |
        |       +---- assembled current editor text -----+
        |                                                  |
        +--------------------------------------------------v
                                    temporary input + resource directory
                                                      |
                                            QProcess child engine
                                                      |
                       PreparedInput -> path/C ABI -> translated 3.28 core
                                                      |
             +----------------------+-----------------+------------------+
             |                      |                                    |
       legacy files          NURBS capture callbacks             2D capture hooks
   DXF/reports/lines        -> OCCT XCAF or AP242 STEP          -> flat-part JSON
                                    |                                    |
                              lep-sim.json                         Print/Cut tab
                                    |
                              Playground XPBD

Preview: temporary output + binary `lep-3d.xbf`; the directory is discarded
after the GUI has loaded XCAF/simulation/flat-part data.
Export: user output directory + `lep-3d.step`; files remain on disk.
```

The GUI never calls the translated core in-process. The child-process boundary
keeps the UI responsive and contains legacy parser aborts/access violations.

## Build targets and dependencies

- `f2c_runtime`: C90 runtime from `third_party/libf2c`.
- `lep_occt`: interface target collecting the required OCCT headers/toolkits.
- `leparagliding-engine`: engine boundary, flat capture, translated core, and
  exact model builder.
- `flatparts`: structured part model, deterministic nester, PDF/DXF writers.
- `softwing_core`: dependency-free XPBD/cloth/contact/pneumatics/suspension
  core and topology-bound persistent SoftBody and suspension checkpoint codecs
  under `src/softwing`.
- `lep_nurbs_model`: exact OCCT model builder and direct authoritative-capture
  to scene-v2 exporter, shared by the engine and focused exporter tests.
- `leparagliding_engine_runtime`: translated calculation core and its path,
  migration, flat-capture, and exact-model dependencies, shared by the thin
  `leparagliding-engine` CLI and the representative real-capture exporter test.
- `simwing_scene`: Qt-free scene-v2 data model, deterministic validation, and
  bounded binary serialization under `src/fsi`.
- `simwing_structure`: Qt-free SimWing-facing adapter around the retained XPBD
  primitives, with a bounded persistent composite checkpoint codec. It links
  `softwing_core` and no Playground library.
- `simwing_transfer`: Qt-free topology-bound coupling surface, exact uniform
  triangle and barycentric patch traction integration, independent
  force/moment/power ledgers, and validated additive application to
  `simwing_structure`.
- `simwing_coupling`: Qt-free trapezoidal macro-step integration of immutable
  transfer samples into nodal impulse/angular-impulse/work ledgers, plus
  checkpoint-transactional XPBD acceptance through equivalent average loads.
  It also owns definition-bound, bounded vector Aitken relaxation with exact
  in-memory iteration checkpoints and the absolute-plus-relative
  displacement/velocity/traction convergence decision. Its stable-ID surface
  reduces those norms from saved-baseline and consecutive kinematics/loads
  without coordinate-origin or bulk-velocity dependence. A checkpointable
  macro-step iteration owner composes current vector, Aitken history,
  convergence, and terminal exhaustion.
- `simwing_coupled_state`: Qt-free composite rollback owner for one real
  `Structure`, one accepted moving-interface fluid epoch, and the strong
  iteration state. It validates complete replacements before committing all
  three through no-throw moves. It also owns the bounded, checkpointable
  macro-step reduction policy. A one-macro-step state owner captures a fresh
  three-owner baseline and makes restore-and-begin the only integrated retry
  transition. Its separate solver-only rewind restores Structure and fluid
  between fixed-point iterations without discarding the advanced Aitken state.
  Its generic loop driver validates each callback-produced physical epoch,
  advances convergence, performs both rollback modes, retains only converged
  output, and restores the baseline on callback or terminal retry failure. Its
  bounded result retains one terminal decision and solver-run count per
  attempted time step. No general wing worker supplies that callback yet; the
  strong-piston canonical is the first real solver-chain callback and
  regression owner.
- `simwing_fluid_structure_bridge`: Qt-free stable-ID bridges from one
  face-aligned fluid-pressure surface to one structural coupling surface. The
  first is the exact uniform subset; the planar face-resolved path clips fixed
  MAC tiles against reference triangles and transfers nonuniform tile traction
  through conservative barycentric patches. Its rigid-normal mode retains
  those material patches while the physical plane moves and the Eulerian grid
  plane rebases, with independent position/velocity/force/moment/power checks.
  The porous overload transfers only the accepted equal-and-opposite sheet
  reaction, excludes separately prescribed pressure sources, and closes the
  source and mapped impulse/work ledgers while retaining porous dissipation as
  a separate energy term.
- `simwing_projected_flag_case`: Qt-free one-way flexible-fabric canonical. A
  finite nonseparating MAC-face panel diverts a persistent accelerating gust;
  the complete pressure-plus-direct constraint reaction crosses the
  face-resolved conservative bridge into an edge-clamped XPBD membrane. Its
  geometry remains a fixed CFD reference and it is not moving cut-cell or
  two-way FSI.
- `simwing_ram_air_cell_case`: Qt-free fixed-reference open-cell canonical.
  Five face-aligned surfaces independently map their complete CFD constraint
  reactions into one shared-node five-panel XPBD shell. Its mouth is open and
  slope-clamped, but displaced panels do not alter the fluid geometry, so it is
  a multi-panel transfer/deformation gate rather than physical ram-air
  inflation or two-way FSI.
- `simwing_piston_case`: Qt-free visible verification harness crossing fluid
  projection, the face-resolved bridge, temporal coupling, XPBD acceptance, and
  immutable viewer frames. It also owns the headless light-piston strong-
  coupling canonical that supplies a real moving-interface/transfer/XPBD
  callback to `simwing_coupled_state`. `--case strong-piston` publishes its
  accepted immutable diagnostic frames. Its bounded
  persistence codec is routed through the batch checkpoint flags with atomic
  same-file replacement, and its typed standard-I/O control adapter exposes
  only complete accepted macro-steps.
- `simwing_porous_sheet_case`: Qt-free midpoint oracle coupling a
  pressure-driven porous sheet to XPBD. It closes fluid/sheet/pump momentum,
  pump work, porous dissipation, and kinetic-energy ledgers before publishing
  the translated sheet. Its immutable composite checkpoint restores Structure,
  fluid fields, and the last accepted ledger with exact next-frame replay. The
  bounded checksummed persistent envelope stores the existing Structure codec
  plus the topology version, axis, wrapped face, signed periodic image, and MAC
  fields, then regenerates the ledger by bounded deterministic replay before
  publishing a decode. It explicitly rebinds the
  sheet across six consecutive ordinary MAC faces, including a `7 -> 0` wrap
  to signed periodic image `1`, and rejects the later collision with the next
  periodic image of the prescribed pump surface. Every accepted pre-pump
  topology epoch persists and continues bit-identically. A focused reverse
  instance has distinct provenance/fingerprint, flips pump and coupled momentum,
  crosses `0 -> 7` into signed image `-1`, resumes that epoch, and rejects the
  next negative pump image. The CLI remains the positive canonical. The reusable
  `fluid/planar_face_topology.*` selector owns strict open-segment placement,
  bidirectional adjacent-face rebases, X/Y/Z axes, and unwrapped periodic-image
  tracking; `fluid/porous_topology.*` preserves the porous compatibility API.
  `fluid/planar_porous_sheet.*` validates an authored epoch and
  expands it into the complete deterministic X/Y/Z transverse crossing plane
  shared with `simwing_porous_flow_case`; the one-way oracle rejects its pump
  collision after selection.
- `simwing_moving_porous_flow_case`: Qt-free owning-frame canonical for the
  stage-resolved full-flow operator. A prescribed `0.4 m/s` plane crosses
  `face=3,image=0` to `face=0,image=1` in its first macro-step, continues through
  a second wrap without resetting fluid state, and publishes the final sheet
  and pump crossings with both stage epochs and conservation ledgers. The
  headless `moving-porous-flow` CLI case is exposed, and an immutable in-memory
  checkpoint restores the complete fluid, diagnostic, kinematic, and topology
  epoch transactionally. Its bounded/checksummed `SWMF` persistence stores the
  fluid fields and sharp crossings explicitly, regenerates diagnostics by
  bounded deterministic replay, and accepts only a bit-identical canonical
  epoch. The headless checkpoint flags provide additional-step resume,
  absolute-step autosave cadence, final-write deduplication, and atomic
  same-path replacement. Structure coupling remains open.
- `simwing_open_piston_case`: Qt-free driven open-piston harness crossing a
  connected-fluid projection, partial-cell control-volume geometry, opening
  transport, exact planar topology rebasing, CFD reaction, XPBD acceptance,
  composite structure/fluid checkpoint replay, deterministic bounded
  persistent restart encoding exposed through the headless checkpoint flags,
  and immutable viewer frames.
- `simwing_scene_structure`: deterministic scene-v2 membrane, per-sheet
  bending, junction, and cable assembly into `simwing_structure`.
- `simwing_scene_fluid_surface`: deterministic compact ownership of the
  authoritative scene-v2 fluid regions, porous fabric, oriented surface
  triangles, and openings, plus immutable capture of accepted Structure
  vertex motion. It does not classify grid crossings or cut cells.
- `simwing_scene_fluid_surface_transfer`: binds that authoritative surface and
  its accepted motion to the existing conservative triangle/quadrature
  traction transfer without weakening scene, Structure, or topology identity.
- `simwing_scene_fluid_geometry`: builds a bounded deterministic cell-major
  broad phase from the fingerprinted accepted surface state and one uniform
  grid, applies an exact normalized-SAT triangle/cell narrow phase, and clips
  barycentric surface polygons into each exact cell. It then resolves internal
  shared-plane duplicates into one oriented face owner and pairs transverse
  clipped-polygon segments into unique oriented internal-face crossings. A
  bounded sparse face-local index retains multiple crossing and coplanar-sheet
  references without treating their summed measures as unions. A provenance-
  keyed graph then stitches adjacent triangle segments through stable authored
  vertices/edges and explicitly marks opening and grid-edge endpoints. Its
  segments become winding-directed open chains or closed loops independently
  per authored region pair. Valid multi-region junctions may terminate several
  pair-specific chains at one physical higher-degree node; branching within
  one region pair still rejects. Simple closed loops then gain
  signed area, centroid, and enclosed/exterior region identity without merging
  nested loops or sealing open chains. Eligible non-touching loop nests and
  simple directed open-chain arrangements whose leaves all reach the
  rectangular face boundary are finally resolved into exact per-region
  MAC-face areas. The latter retain every source chain and enumerate bounded
  left faces from authored winding. Opening-ended and other incomplete faces
  stay explicit unresolved. None of these products is a cut-cell volume or
  complete grid-region classification.
- `simwing_scene_fluid_quadrature`: converts those unique cell/face area owners
  into stable barycentric quadrature points with authored material, sheet, and
  two-sided region identity plus exact side-cell provenance, then delegates
  ordered traction to the existing conservative Structure transfer. Its
  pressure adapter converts explicit one-sided CFD samples into
  `(p_negative-p_positive)*normal` exactly once; it does not calculate
  pressure, shear, or any polar force.
- `simwing_scene_fluid_grid_epoch`: composes one accepted Structure surface
  through every grid-geometry stage and unique quadrature under one versioned
  fingerprint. Per-stage bounds remain active and the completed vector payload
  has an aggregate byte ceiling. It does not infer cut-cell volumes or solve
  moving-boundary fluid equations.
- `simwing_scene_fluid_opening_cap`: builds topology-only virtual caps for
  simple authored openings. Material and cap incidences must form one closed
  oriented region cycle around every finite-area edge, including valid
  three-region sheet/cap junctions; adjacent fabric with the same region pair
  fixes negative-to-positive winding. A closed material boundary may remain
  cap-free only while both reference and accepted geometry are collapsed.
  Planar convex loops retain their exact fan and planar concave loops use
  bounded deterministic reference-geometry ear clipping. Scene-v2.2 may
  instead provide one oriented boundary-vertex disk, allowing nonplanar facets
  and their individual normals to remain stable under accepted motion.
  Folded, self-intersecting, degenerate, or unauthored nonplanar caps reject.
  Caps follow accepted Structure vertices but never
  become fabric or traction. Its companion one-point triangle
  quadrature exactly integrates accepted piecewise-linear cap velocity into a
  bounded surface-sweep ledger. Exact grid patches then own positive cap area
  once per cell or non-periodic face and retain their clipped polygons.
  Cell-owned polygons also pair their exact shared-face boundary segments into
  bounded, directed virtual-cap crossings. Both adjacent cell clips must agree
  before one crossing is published; face-owned aperture area is counted but
  never converted into a line, and grid-edge ambiguity remains explicit. The
  accepted pressure epoch owns this crossing product. A bounded planar
   half-edge arrangement then combines those crossings with directed material
   chains, including disconnected signed cycles, and publishes exact
   same-region area, global first moment, and centroid plus a
   resolved-or-unresolved record for every touched
  Cartesian face. Its polygon moments use a local chart before the public
  global centroid/first moment is reconstructed. Face-owned aperture area
  remains with its existing owner. When a distant coordinate origin makes a
  fixed ULP envelope exceed the minimum arrangement tolerance, topology also
  moves into that chart; ordinary-origin arithmetic remains unchanged.
  Material-only face accounting omits same-region internal-sheet chains from
  pressure-region boundaries while retaining their audited count upstream. A
  bounded region-separating chain arrangement accepts stitched interior
  multi-region junctions. On the 4-by-4-by-4 real-wing audit this removes all
  ten pressure-active material failures: 15 material-only plus 43 cap-touched
  faces publish resolved pressure partitions, leaving only the separately
  typed embedded-opening stencil rejections.
  Exact coordinates shared by a clipped edge's endpoints remain canonical
  through later-axis clipping, preventing one-ulp drift from erasing an
  earlier face owner. On the coarse real wing this closes all nine cap-crossed
  faces and reaches ten resolved pressure-link partition faces. Rejection
  statuses distinguish face-owned
  aperture overlap, coplanar material, invalid source geometry, opening-owned
  dangling endpoints, unstitched intersections, winding/region ambiguity, and
  area-closure failure. Neither product makes the cap fabric. A
  read-only evaluator samples exact face-normal MAC values or bounded
  staggered interpolation, subtracts cap motion, and publishes oriented volume
  flow. Each opening contributes equal-and-opposite outward fluid, surface,
  and relative-flow balances to its negative- and positive-side regions, so
  intakes and crossports share one canonical conservative source convention.
  It does not modify projection or fluid connectivity.
- `simwing_scene_fluid_cell_volume`: reconstructs bounded deterministic sparse
  per-cell region volumes for the first capped-manifold subset. It requires one
  Outside region and closed oriented material-plus-cap region cycles. Pairwise
  signed measures support valid multi-region junctions. Grid-face crossings
  split into pair-specific open chains at those junctions. Stitched
  region-separating arrangements now close; opening-ended, dangling, or
  otherwise incomplete arrangements remain explicit unresolved face
  partitions.
  Each oriented triangle becomes a signed tetrahedron against the grid origin;
  bounded convex clipping distributes it across exact cells, including full
  interior cells and open face-local tile chains. The same clipped tetrahedra
  retain exact first moments and cell-local centroids. Each tetrahedron, cell
  volume, and cell first moment closes independently before the result agrees
  with whole-surface region volumes. General swept-volume remap, unauthored
  nonplanar or folded/self-intersecting mouths, opening-only cap vertices,
  general junction-aware swept-volume remap, and periodic-boundary ambiguity
  remain unsupported.
- `simwing_scene_fluid_region_continuity`: binds two consecutive accepted
  volume and opening-flux epochs. It trapezoidally integrates each region's
  outward relative opening flow and reports `delta volume + integrated flow`
  without modifying either solver. Inputs retain exact surface-state and
  producer fingerprints; skipped epochs and mutated products reject.
- `simwing_scene_fluid_region_connectivity`: treats authored intakes and
  crossports as pressure-connectivity edges, canonicalizes components and one
  stable gauge region per component, then reduces the two-epoch continuity
  ledger component-wise. Sealed components with nonzero volume change reject
  pressure compatibility even when equal-and-opposite changes cancel globally.
- `simwing_scene_fluid_pressure_control_volume`: turns every positive sparse
  cell-region volume into one immutable volume-weighted pressure unknown. It
  preserves the exact cell-region centroid and cell/region/component ownership,
  derives a stable ID from the fixed-grid cell/region key, and selects a
  deterministic gauge control volume in each authored-opening component.
- `scene_fluid_mimetic_control_cell.{h,cpp}` is an audit-only immutable adapter
  from those controls plus Cartesian links, material quadrature, and opening
  patches to complete periodic-image-unwrapped half-face shells. It retains
  paired trace/wall identity, other-control provenance, exact area, centroid,
  outward normal, area-vector and `N^T R` closure, limits, and fingerprints.
  Analytic nested and open fixtures close and build the generic SPD kernel.
  On the coarse real wing, all 138 controls close and are topology-complete;
  the accepted closed material-plus-cap surface proves that ten otherwise
  ambiguous untouched faces are Outside without using dominant volume or shell
  closure as a label. Its 42,826 cell-owned opening patches become 85,652
  paired half-faces, so rejected two-point apertures are not missing geometry.
  The manual `4 x 4 x 4` audit likewise classifies six untouched faces and
  closes all 358 controls. It safely omits 240 zero-volume material sides,
  misses no opening side, retains 95,984 paired opening half-faces, and reaches
  2,947 half-faces in its largest control. It does not replace the production
  graph operator.
- `scene_fluid_mimetic_trace_system.{h,cpp}` is the audit-only immutable
  global mixed-hybrid topology and matrix-free action over those completed
  shells. Shared Cartesian and authored-opening identities form exactly one
  two-incidence trace, while every material half-face owns a unique
  one-incidence zero-flux wall trace. It retains one deterministic gauge per
  component, positive Jacobi diagonals, compact local kernels, bounded storage,
  and nested fingerprints; disconnected controls within one component reject.
  Cell-scalar elimination is locally conservative;
  the resulting trace action is symmetric positive semidefinite with one
  component-constant null mode and a compatible source-derived right-hand
  side. The coarse real wing has 191,579 trace unknowns and 13,132,336 bytes of
  compact local factors.
- `scene_fluid_mimetic_trace_flow.{h,cpp}` samples one immutable oriented
  relative predictor flow for every shared mimetic trace. The fixed bootstrap
  reads exact Cartesian MAC partitions and accepted per-patch relative opening
  flux. The moving overload instead projects accepted material-wall-adjusted
  endpoint region velocities onto every exact trace normal and subtracts the
  same accepted cap sweep. It exactly matches existing graph predictors where
  those links exist and still covers a deliberately rejected embedded opening.
  Material-wall traces stay absent and impermeable. The bounded product retains
  complete topology/opening/MAC or wall-exchange provenance and source density.
  On the coarse real wing the fixed sampler covers all 42,927 shared traces:
  101 Cartesian traces plus 42,826 cell-owned opening traces, including every
  aperture omitted by the rejected two-point graph.
- `scene_fluid_mimetic_trace_solve.{h,cpp}` is the audit-only transactional
  gauge-fixed Jacobi-PCG solve shared by the full and material-wall-condensed
  matrix-free systems. It admits only bounded component-sum roundoff, uses the
  exact stored diagonal, publishes only after a fresh residual check, and
  preserves the warm start exactly on incompatibility, non-finite arithmetic,
  or exhaustion. Manufactured reduced fields recover deterministically;
  reconstructed wall traces close every original full row, and source-driven
  local/trace conservation closes. One complete real-wing iteration over both
  the 191,579-row full system and its 42,927-row condensed system remains
  finite, reduces the residual, and rolls back when truncated. The condensed
  real system additionally reaches `1e-5` relative RMS within 300 iterations
  and reconstructs the full operator below `2e-4` maximum row residual. Its
  trusted PCG loop validates the immutable product once, then uses the internal
  assuming-validated action; public one-shot actions still validate every
  call. Jacobi is not practical on the uncondensed full system. Stronger
  local/multilevel preconditioning and production integration remain open.
- `scene_fluid_mimetic_condensed_trace_system.{h,cpp}` composes the local wall
  Schur products into an immutable shared Cartesian/opening trace topology. It
  owns exact full/reduced mapping, deterministic shared component gauges,
  summed positive diagonals, bounded nested storage, matrix-free action, full
  RHS condensation, and unique wall-trace reconstruction. Manufactured global
  fields prove symmetric positive-semidefinite action, component null modes,
  diagonal agreement, reduced-RHS equivalence, full-field recovery, and closure
  of every original row. The coarse product has 42,927 reduced traces and
  3,986,602 bytes of local Schur data.
- `fluid/mimetic_wall_condensation_detail.h` and
  `scene_fluid_mimetic_condensed_trace_system_detail.h` expose only the
  assuming-validated repeated actions used inside the already validated PCG
  transaction. The local action is fused as the exact diagonal-plus-seven-mode
  active Schur form `D_a + U_a (K - K Q K) U_a^T`; it does not repeat two full
  balances or a wall solve. Do not call either detail entry point from an
  unvalidated public boundary.
- `scene_fluid_mimetic_pressure_solve.{h,cpp}` is the audit-only atomic
  integrated-source transaction. It assembles the full and reduced RHS, solves
  the condensed traces, reconstructs every wall trace, reevaluates cell
  scalars/half-face fluxes/all original trace rows, and publishes fields only
  when both reduced and reconstructed RMS residuals fit the declared residual
  and component-compatibility envelope. Failure returns diagnostics with empty
  state. A balanced coarse real-wing source pair converges in 307 iterations,
  with `6.17e-9` reconstructed RMS, `2.59e-6` maximum original-row residual,
  and `5.05e-10` maximum cell-conservation residual. Results retain full,
  condensed, and optional source fingerprints.
- `scene_fluid_mimetic_pressure_state.{h,cpp}` captures the immutable accepted
  endpoint of a source-bound atomic solve. Before publication it independently
  rebuilds the full RHS, reconstructs every material-wall trace from the
  reduced field, reevaluates all cell scalars/fluxes, requires exact agreement
  with the solve result, and requires zero shared gauges. It persists only
  stable control pressures and shared trace pressures with complete
  control/full/condensed/source provenance, count/byte bounds, and fingerprinted
  corruption rejection. The coarse real source solve now crosses this boundary.
- `scene_fluid_mimetic_pressure_state_persistence.{h,cpp}` owns the deterministic
  `SWMP` little-endian restart envelope for that endpoint. It encodes the
  complete stable control/shared-trace rows, state/source/topology provenance,
  derived summaries, FNV-64 payload checksum, protocol/state versions, and
  reserved fields under independent byte/control/trace limits. Decode receives
  trusted rebuilt mimetic control/full/condensed topology, rejects foreign
  identity before publication, and then runs the complete state validator.
  Repeated and decode/re-encode bytes are identical; magic/version/reserved,
  checksum, truncation, trailing data, record limits, and a foreign topology
  leave the caller's destination unchanged. The coarse real-wing state also
  round-trips exactly.
- `scene_fluid_mimetic_pressure_warm_start.{h,cpp}` maps that accepted state
  across exactly one fingerprinted pressure-topology transition without
  mutating either epoch. Retained controls follow stable identity, appeared
  controls use the transition's same-region area-weighted donors, and retired
  controls disappear. Retained shared traces preserve their old pressure
  difference; new traces start from the mean of their rebased endpoint-control
  pressures. One current-component offset then fixes every deterministic gauge
  to exact zero. The product owns bounded working/publication storage, complete
  previous/current topology provenance, trace appearance/retirement counts,
  gauge shifts, and corruption-detecting fingerprints. A moving analytic case
  exercises a real control/trace appearance and feeds the result directly to
  the atomic solve. The opt-in pressure-cell shadow now consumes it after the
  first audited endpoint; the default graph worker remains unchanged.
- `scene_fluid_mimetic_pressure_sampling.{h,cpp}` is the accepted-state load
  adapter. Every material quadrature side resolves through the exact
  cell/region pressure-control row shared by the mimetic topology and state.
  The two authored sides must belong to one pressure component before their
  gauge-safe difference is published. The bounded fingerprinted sample set
  then reuses the existing conservative pressure quadrature and Structure-load
  transfer without another force path. Analytic and coarse real-wing tests
  require exact one-sided values plus force/moment closure. The opt-in worker
  shadow retains these samples but does not apply them; production load
  selection is still the graph pressure operator.
- `scene_fluid_pressure_shadow_comparison.{h,cpp}` is the immutable evidence
  gate at that unchanged load boundary. It pairs every graph/reference and
  mimetic material pressure jump by exact control/region identity, evaluates
  both through the existing conservative transfer, retains every sample and
  nodal-force delta, and reports L2/RMS/maximum/relative pressure plus net,
  nodal, moment, and power differences under independent count/byte bounds.
  It never applies loads. The pressure-cell reports about 60% relative
  pressure/force difference both at step 4 and after 600 steps, so production
  selection must remain on the graph path. The same product now binds every
  graph projection source row to its mimetic control/source row: the source
  vectors agree within `5.3e-16` relative at step 4 and `1.7e-16` after 600
  steps. In this symmetric case the remaining pressure and nodal-load fields
  are almost pure scalar multiples (shadow gains `2.53035` and `2.50693`, with
  post-fit relative shape residuals near `2e-16`). The unresolved transition
  question is therefore the cut-cell graph-versus-mixed-hybrid operator
  response, not source units, gauge-safe sampling, or conservative transfer.
- `scene_fluid_pressure_owner_transition.{h,cpp}` turns that immutable
  evidence into a separate typed decision without applying either field. Its
  fingerprinted policy requires source agreement by default and independently
  bounds pressure magnitude/shape, nodal-force magnitude/shape, net force,
  moment, power, and both conservative-transfer closures. Any failed criterion
  selects the reference graph owner and records a stable rejection bit; only a
  zero-rejection decision names the mimetic candidate. The exact real-wing
  self-comparison exercises the positive branch with source evidence explicitly
  waived, while the live pressure-cell records pressure magnitude, pressure
  scale, nodal-force scale, and net-force rejections and therefore retains graph
  loads. The opt-in coupling transaction now owns that decision only after the
  comparison succeeds, publishes it through diagnostics and the pressure-cell
  accessor, and reports `owner=graph` with rejection mask `0xeb00` at step 4
  and `0x6b00` after 600 steps, where the transient power-delta rejection has
  cleared. A deliberately permissive zero-rejection policy selects the mimetic
  candidate while leaving the production frame byte-identical, proving that
  selection is still diagnostic. This is an agreement gate, not a claim that
  the graph operator is a continuum oracle, and no worker consumes the
  decision to apply loads yet.
- `scene_fluid_pressure_operator_response_audit.{h,cpp}` is the offline
  bounded inverse-response discriminator for that remaining question. It runs
  the accepted source plus six deterministic component-compatible coordinate,
  high-frequency, and authored-region modes through independent graph and
  mimetic solves, aligns
  control pressures by component gauge, retains every source/response/residual
  row, and fingerprints solver residuals, fitted gain/shape diagnostics, and
  gauge-invariant source work. A one-component/two-region mode additionally
  reports its energy-equivalent two-terminal conductance.
  On the 65-control pressure cell the accepted direction has about `2.54` gain
  with `2.4%` control-field shape residual, while the five bulk modes have
  `0.999`-`1.007` gains and `1.2%`-`16.5%` residuals. A piecewise-constant
  authored-region mode, whose graph source exists only on intake links,
  reproduces the outlier at `2.562` gain and `1.98%` residual. The live
  pressure-cell graph conductance is exactly `0.18 m`, versus `0.0700821 m`
  for the mixed-hybrid response; their `2.5684` energy ratio independently
  closes the observed pressure gain. The mismatch is therefore opening/cut-
  cell closure, not a global pressure-unit correction or bulk mode. This audit
  is never run inside live coupling and never applies loads.
- `scene_fluid_mimetic_region_conductance_audit.{h,cpp}` is the graph-
  independent one-component/two-region terminal primitive for the next
  discriminator. It identifies every permeable cross-region trace in the
  trusted mixed-hybrid topology, including Cartesian traces for face-aligned
  authored openings and `AuthoredOpeningTrace` rows for embedded openings,
  while excluding material walls. A fixed balanced transfer is distributed
  at uniform areal source density, solved directly through the condensed
  mixed-hybrid system, and reduced to gauge-invariant source work and
  conductance. The product owns every paired opening/source/pressure row,
  topology and solve provenance, count/byte bounds, and a complete integrity
  fingerprint. On the face-aligned pressure cell it gives
  `0.0700820848335194 m`, agreeing with the older graph-manufactured shadow
  response within `3.6e-14 m`; the rejected-two-point embedded fixture gives
  `0.0608388978079475 m`. This is an offline Neumann terminal experiment, not
  a continuum oracle or load path.
- `scene_pressure_cell_geometry.{h,cpp}` owns both the unchanged visible
  worker tetrahedron/default `4^3` grid and a separate fixed, deliberately
  skew intake used only by refinement diagnostics. Keeping construction here
  prevents the worker and offline oracle from drifting while preserving their
  distinct scene checksums.
- `scene_pressure_cell_operator_refinement_audit.{h,cpp}` rebuilds the skew
  rest geometry on bounded isotropic grids and nests the complete six-mode
  response audit for every sample. Its `2^3`, `4^3`, and `8^3` conductance
  ratios are `28.024`, `33.672`, and `8.830`; cell-width/area-normalized shadow
  conductance rises smoothly `0.102 -> 0.358 -> 0.503`, while the graph values
  jump `2.857 -> 12.063 -> 4.446`. This proves severe grid sensitivity in the
  embedded two-point graph response, but does not yet make the shadow operator
  a production or physical oracle.
- `scene_pressure_cell_operator_phase_audit.h` and the same implementation
  run a bounded fixed-`4^3` placement ensemble over the eight zero/negative-
  half-cell corner phases. Six placements build complete operators; two retain
  typed pressure-operator failures with exactly one and two unresolved
  embedded-opening patches. Across accepted placements the normalized graph
  conductance has population coefficient of variation `0.77175`, versus
  `0.34030` for the shadow response. This is grid-phase sensitivity evidence,
  not a convergence result or permission to change the live pressure owner.
- `scene_pressure_cell_operator_phase_refinement_audit.{h,cpp}` nests that
  identical eight-phase ensemble at `2^3`, `4^3`, and `8^3`. Complete graph
  and shadow products exist for only `4/8`, `6/8`, and `2/8` placements. Their
  conditional normalized graph means are `3.9025`, `6.2525`, and `4.5187`;
  conditional shadow means are `0.10168`, `0.26909`, and `0.45753`. The latter
  trend is not an unbiased convergence result because graph rejection censors
  the phase set before the paired response audit.
- `scene_pressure_cell_mimetic_conductance_phase_refinement_audit.{h,cpp}`
  runs the graph-independent area-weighted terminal source across the complete
  phase/refinement matrix without constructing the graph operator. It owns all
  40 requested sample records, nested terminal responses, conditional
  statistics, bounds, and exact typed local-cell linear-consistency failures.
  All `8/8` phases solve at every level, with normalized mean/CV
  `0.100660/0.01043`, `0.240930/0.39050`, `0.521248/0.38828`, and
  `0.902570/0.20104`, and `1.091977/0.07436` at
  `2^3`/`4^3`/`8^3`/`16^3`/`32^3`; the `32^3` range is
  `0.956232`-`1.207995`. The former five `8^3` rejections were
  removed by local-coordinate polygon moments and canonical authored-
  triangle/grid-edge intersections, without closure fitting or relaxing the
  default `1e-10` algebraic tolerance. Tiny `16^3` complementary boundary
  regions additionally retain the full-minus-positive closure area while
  selecting an independently integrated reverse-polygon centroid only beyond
  a fixed coordinate-ULP envelope. The resulting four-face/two-wall sliver
  uses the bounded direct wall-Schur fallback when its normal `7 x 7`
  auxiliary core is numerically singular. At `32^3`, independently
  fingerprinted sparse-publication tolerances retain a real `2.83e-15 m^3`
  complementary region, cell-local first moments bound its centroid, and an
  explicit offline `1e-9` local algebraic tolerance accepts the two otherwise
  valid operators whose errors are `5.52e-10` and `1.31e-10`; production keeps
  `1e-10`. Phase CV contracts on both `8^3 -> 16^3` and `16^3 -> 32^3`, but
  this is not a convergence claim. This removes graph censoring from the
  shadow-only matrix, and the area-weighted
  multi-opening source is deliberately not the older graph-manufactured
  source. The production 600-step trace and audited checkpoint remain
  byte-identical. Its common scene/grid translation setting additionally
  preserves the selected `8^3`, phase `[0,-0.5,0]` topology under a
  `[256,-512,1024] m` shift: intake area and normalized conductance change by
  only `5.93e-14 m^2` and `5.59e-12`, respectively. Non-finite translations
  reject before assembly.
  The complete streamed opt-in `64^3` ensemble accepts `8/8`, with normalized
  mean/CV `1.262883/0.03130` and range `1.215742`-`1.321663`. Over
  `16^3 -> 32^3 -> 64^3`, mean-increment contraction is only `0.902325`,
  apparent order is `0.148280`, and extrapolation gap is `0.555593`; all fail
  the existing screen despite another CV contraction. Three phases reverse
  direction and four are noncontracting. The compact streamed report does not
  replace the immutable three-source assessment, and the production decision
  remains `InsufficientEvidence`.
- `scene_pressure_cell_mimetic_conductance_convergence_assessment.{h,cpp}`
  consumes the finest levels of three compatible immutable phase audits and
  screens aggregate Richardson evidence separately from every same-phase
  trajectory. For `8^3`/`16^3`/`32^3`, the mean contraction is `0.49671`,
  apparent order is `1.00952`, extrapolated normalized conductance is
  `1.27891`, and the fine gap is `0.14616`. Those aggregate checks pass, but
  four phases reverse increment direction and four exceed the `0.75`
  contraction bound; only one phase passes both. The strict typed outcome is
  `InsufficientEvidence` with mask `0x300`. Policy, three source fingerprints,
  trajectories, diagnostics, bounds, and outcome are immutable. A permissive
  policy can name only a read-only continuum-trend candidate; neither outcome
  enters worker pressure selection or load application.
- `scene_fluid_mimetic_pressure_epoch.{h,cpp}` is the first atomic acceptance
  boundary over the isolated mimetic pressure path. Bootstrap allocates a
  bounded zero reduced field; the consecutive overload consumes the exact
  fingerprinted warm remap. Both solve the source-bound condensed system,
  independently capture accepted state, and sample every material pressure
  side before publication. Exhaustion or incompatibility returns solve
  diagnostics with no state or sample payload. The result retains complete
  current topology/source/quadrature provenance plus optional warm/transition
  fingerprints and nested bounds. It accepts an already assembled physical
  source product; trace-flow and `dV/dt` ownership remain explicit upstream.
- `scene_fluid_mimetic_pressure_audit.{h,cpp}` owns the complete immutable
  shadow endpoint: rebuilt control/full/condensed topology, exact fixed-MAC or
  transported-wall trace predictor, physical source, accepted state, solve
  diagnostics, and material samples. Its endpoint fingerprint covers the
  complete nested solve diagnostics and every nested product fingerprint;
  aggregate dynamic storage is independently bounded. A fixed pre-operator
  entry point crosses the coarse real wing even though the graph operator
  rejects some embedded openings, reproducing all 138 controls, 42,927 shared
  traces, source rows, and material samples. `SceneFluidPressureCoupling`
  invokes the endpoint only after graph convergence when explicitly enabled:
  bootstrap samples the unchanged MAC predictor, later epochs consume the
  accepted transported wall predictor and consecutive warm remap. Audit
  failure rolls Structure and the graph owner back; success never applies the
  mimetic samples. `simwing-fsi --case pressure-cell
  --mimetic-pressure-audit` exposes this read-only experiment. The default
  frame and graph-checkpoint portion remain identical. `SWPCELL10` composes
  the compact accepted `SWMP` rows and audit-settings fingerprint with that
  graph restart, rebuilds trusted topology before decode, and resumes the exact
  consecutive wall-predicted endpoint without persisting transient topology.
  Each live endpoint also publishes the bounded graph-versus-mimetic pressure/
  load comparison; its failure rolls back with the endpoint.
- `scene_fluid_mimetic_pressure_source.{h,cpp}` is the immutable physical-unit
  bridge from per-control predicted net-outward flow plus optional `dV/dt` to
  the integrated source `-(rho/dt)*(dV/dt + net_outward)` in `Pa*m`. It is
  bounded and fingerprinted to the complete mimetic control topology, retains
  compensated per-component compatibility sums, and can feed the atomic solve
  without a raw field. It now also accumulates the oriented shared-trace flow
  product into exact per-control net-outward flow and maps an accepted
  consecutive-epoch pressure-volume-rate product by stable control identity.
  Source provenance retains both optional upstream fingerprints and rejects a
  volume-rate duration inconsistent with the pressure timestep. A transported
  wall-adjusted flow additionally requires exact source-density agreement.
- `simwing_scene_fluid_pressure_face_link`: binds those pressure unknowns
   across Cartesian faces. Exact resolved partitions create one same-region
   link per positive region area and retain its exact subface centroid. Capped material-plus-opening partitions
  supersede the material-only result on transversely crossed faces, while an
  unresolved capped arrangement retains its own status and publishes no link.
  Provably untouched faces use an exact common region when available; when
  sparse supports overlap, the accepted material-plus-cap closed surface may
  classify the face-centroid region by oriented solid angle, and that region
  must still own controls in both cells. A face-aligned authored
  opening creates oriented cross-region links over its exact patch area and an
  unambiguous same-region link over the complementary face area. A cell-owned
  opening patch creates an embedded same-cell cross-region link along its
  authored normal, with conductance distance from the exact projected pressure
   centroids. A patch without admissible projected centroid separation retains
   a typed rejection with patch/opening identity, signed centroid-to-cap-plane
   distances, exact count and area, and publishes no fabricated link. The
   rejection also owns bounded deterministic ranges of same-region Cartesian
   one-ring neighbors for both side controls, exact link/control provenance,
   periodic-image offsets, and signed sidedness. This is reconstruction
   evidence only and does not reroute aperture flux.
  Material/open-chain/coplanar overlap and any surface classification that is
  non-integral, non-unique, or absent from either cell remain explicit and
  unlinked.
- `simwing_scene_fluid_pressure_operator`: assembles those links into a
  symmetric integrated graph Laplacian with one retained gauge owner per
  pressure component. It requires complete face ownership and exactly one
  link-connected graph per authored component. Resolved face-aligned and
  embedded off-face openings join their regions; unresolved Cartesian faces or
  embedded patches reject instead of becoming implicit sealed walls. That
  rejection is a typed `SceneFluidPressureIncompleteFaceOwnershipError` with
  exact resolved/unresolved face and embedded-patch counts, so offline studies
  do not classify failures by diagnostic text. Its
  deterministic component-wise
  CG solve
  admits only bounded compatibility roundoff, commits only after an explicitly
  recomputed residual converges, and normalizes every retained gauge to exact
  zero. It deliberately owns no physical RHS or link-flow correction; the next
  adapter owns those stages.
- `scene_fluid_pressure_projection.{h,cpp}` supplies the first physical
  link-flow adapter on that topology. It maps predicted MAC normals to one
  oriented flow per same-region link, substitutes the accepted relative
  opening flux on face-aligned or embedded aperture links, assembles `-rho/dt`
  times net
  outward control-volume flow, solves transactionally, and corrects every link
  independently. Pressure and corrected flow publish only after an explicit
  continuity check. Its consecutive-epoch overload adds the exact geometry
  volume rate and instead closes `dV/dt + net flow = 0` locally. It does not
  reconstruct one ambiguous corrected MAC value across partitioned faces and
  embedded openings, and it does not itself remap pressure/momentum state.
- `scene_fluid_pressure_volume_rate.{h,cpp}` matches two consecutive accepted
  pressure-volume epochs by stable cell/region identity and publishes exact
  per-control-volume, component, and global geometry-volume rates. A newly
  positive row is explicitly marked and receives an exact zero-volume previous
  endpoint. The topology-aware overload retires a disappeared row's complete
  previous volume to its unique retained same-region neighbour; ambiguous or
  donorless retirement rejects.
- `scene_fluid_pressure_topology_transition.{h,cpp}` is the single versioned,
  bounded owner of consecutive pressure-topology identity. It pairs retained
  rows, publishes current-topology retained donors for appearances, and
  publishes the unique previous-topology retained recipient for each supported
  disappearance. Volume rates, region-state rebase, and pressure warm starts
  consume this same fingerprinted mapping; they do not rediscover ownership.
- `scene_fluid_pressure_sampling.{h,cpp}` resolves both authored sides of each
  quadrature-v2 patch to exact cell/region pressure unknowns from an accepted
  projection. Both controls must share one pressure component before their
  gauge-invariant difference can enter the existing pressure-traction and
  conservative Structure transfer. Rejected projections and independently
  gauged sealed sides do not publish loads.
- `scene_fluid_pressure_epoch.{h,cpp}` atomically composes one accepted grid
  remap, opening cap/patch/crossing topology, capped Cartesian-face
  partitions, sparse cell volumes, pressure controls, face links, and ungauged
  operator. Its fingerprint binds the complete chain to one Structure surface
  state and its aggregate payload is bounded. It owns geometry/operator
  preparation only, not velocity, a pressure solution, or Structure mutation.
- `scene_fluid_pressure_coupling.{h,cpp}` owns the first strong
  feedback step for that pressure path. It iterates canonical end-pressure
  nodal loads with Aitken relaxation, rewinds Structure before every solve,
  applies start/end loads trapezoidally through `simwing_coupling`, rebuilds
  the moving pressure epoch, and commits only when displacement, velocity, and
  applied-versus-solved load residuals converge. Exhaustion and all failures
  restore the exact structural baseline and preserve accepted pressure state.
  Its in-memory checkpoint composes Structure with the accepted sparse
  projection; restore rebuilds the epoch and conservative baseline load before
  committing, and exact-next-step replay is required from initial and accepted
  checkpoints in both the original and an equivalent owner.
  The bounded first crossing subset admits current pressure controls that
  appear through current one-ring donors and accepted controls that disappear
  into one unique previous one-ring recipient. One shared topology-transition
  product drives geometry-volume retirement, pressure warm-state mapping, and,
  when region transport is active, the explicit region rebase below before
  wall exchange. Ambiguous or unsupported donor topology still rolls back the
  entire macro step.
  Its predicted MAC field is fixed during each nonlinear iteration. Accepted
  corrected link flows can be conservatively area-collapsed back to one
  absolute bulk MAC velocity per Cartesian face, including oriented intake-cap
  sweep, with the lost mixed-subface velocity spread reported explicitly.
  Embedded opening links have no unique Cartesian face: collapse validates and
  reports them without smearing their flow, while accepted region momentum
  retains their explicit-normal state. This is a continuation primitive; it
  does not yet own cut-cell advection, viscosity, or topology rebase/remap.
- `scene_fluid_pressure_link_flow.{h,cpp}` owns the explicit topology-stable
  inverse of that collapse. It restores previous opening-cap sweep, carries
  absolute per-link velocity deviations, recentres them under current link
  areas, reapplies current sweep, and exactly preserves the current bulk face
  total. Pressure projection can consume the fingerprinted result explicitly.
  Do not enable static carry in the worker: without region-resolved convection
  it strongly drains repeated pressure and fabric load. The regression that
  demonstrates this is a load-bearing negative oracle.
- `scene_fluid_region_momentum.{h,cpp}` reconstructs accepted corrected
  absolute links into one immutable momentum vector per positive cell/region
  pressure volume. Cartesian-only controls preserve the exact component-wise
  area average. Controls incident to embedded openings use a bounded 3x3
  normal-equation reconstruction over every incident explicit link normal;
  its unconstrained nullspace retains the cell-centred value from the exact
  MAC predictor fingerprinted by the projection. Its momentum/energy/fallback
  and reconstruction-residual ledgers are the input boundary for conservative
  region transport, not a time advance or wall model.
- `scene_fluid_region_transport.{h,cpp}` owns the first region-momentum step.
  Corrected relative link flow carries donor-cell vector momentum; graph
  viscosity applies equal-and-opposite impulses. Moving pressure volumes
  advance linearly through their accepted geometry rates, giving the uniform
  flow an exact discrete GCL update. An optional cell-centred delta between two
  bound bulk MAC predictors transmits the worker's pump/advection/viscosity
  split without erasing cut-region velocity differences. Bounded deterministic
  subcycling, stage energy, bulk impulse/work, and global momentum ledgers gate
  publication. Material-wall viscosity is a separate downstream product; this
  transport owner still does not perform topology rebase or the next pressure
  projection.
- `scene_fluid_region_rebase.{h,cpp}` owns the bounded first current-topology
  rebase and consumes the shared pressure-topology transition rather than
  deriving ownership independently. Retained controls keep transported
  velocity; each appeared control receives the transition's area-weighted
  retained same-region velocity. Each disappeared control transfers its
  complete transported volume and momentum to the transition's unique retained
  same-region recipient. The mapped source ledger closes before current-volume
  geometric correction. Pressure warm starts use that same transition to
  retain stable rows, seed appeared rows, and drop retired values.
  Cross-material donation, missing/ambiguous one-ring ownership, corruption,
  and configured limits reject before publication. This is a bounded
  conservative state rebase, not a general swept-volume remap.
- `scene_fluid_region_wall.{h,cpp}` remaps one immutable accepted transport to
  each current strong-coupling geometry, or consumes the explicit current-
  topology region rebase, and applies two-sided tangential
  viscous exchange at the authoritative material quadrature. Its local wall
  distance is the greater of a configured floor and half control-volume
  volume per incident wall area. Deterministic explicit subcycling bounds each
  aggregate viscous row. Fluid impulse, equal-and-opposite Structure traction,
  wall work, nonnegative dissipation, and exact accepted endpoint provenance
  gate publication; pressure remains the sole owner of normal traction. This
  is a coarse wall closure, not a resolved boundary layer.
- `scene_fluid_region_link_flow.{h,cpp}` remaps an accepted fixed-epoch
  transport candidate onto the next topology-stable moving pressure epoch. It
  preserves stable cell/region velocity while current physical volume changes,
  projects endpoint velocity onto each current link's explicit normal while
  preserving exact component arithmetic on Cartesian links, and subtracts
  exact current opening-cap sweep. Pressure projection can consume this
  fingerprinted first-order GCL predictor together with the consecutive
  volume-rate product. A second overload consumes the current material-wall-
  adjusted region state and preserves that exchange as explicit projection
  provenance. It consumes an already rebased wall state but does not perform a
  topology rebase itself.
- `scene_pressure_cell_case.{h,cpp}` exposes the shared visible geometry as the selectable
  `simwing-fsi --case pressure-cell` diagnostic. The soft open tetrahedral cell
  has three fixed intake vertices and no mechanical actuator. It maintains a
  prescribed periodic mean wind through a uniform pump, advances a private
  collapsed-MAC candidate with the symmetric SSPRK2 viscosity/projected-
  nonlinear-advection/viscosity bulk operator, and commits it only after
  applying the final scene pressure load to Structure. Immutable frames publish
  displacement, triangle pressure jump, separate pressure/wall/total-fluid and
  pump vectors, bulk-flow change/divergence/viscous loss, wall
  loss/conservation, iteration count, bulk MAC speed, and subface-collapse
  spread, plus the embedded-opening count retained outside the Cartesian bulk
  field. The final projection is cut-region-aware, but the intermediate bulk
  projections do not retain distinct per-region subface velocities. After
  bootstrap it
  owns accepted region momentum, applies the bulk-MAC predictor delta,
  transports that state through moving-volume GCL, exchanges material-wall
  momentum, and uses the adjusted current-link predictor in every strong
  iterate. Its
  in-memory checkpoint must reconstruct the accepted derived predictor and
  replay the exact next frame. It remains an airflow-bootstrap canonical, not
  a wing or general immersed-boundary CFD claim.
- `scene_pressure_cell_checkpoint_persistence.{h,cpp}` owns the bounded
  `SWPCELL10` canonical restart envelope. It persists the trusted nested
  Structure payload and complete accepted sparse pressure projection, rebuilds
  canonical identity/settings instead of trusting them from the wire, and
  validates exact next-frame replay before decoded state becomes visible. The
  complete accepted region-momentum state, accepted material-wall traction,
  and projection provenance are serialized; audited checkpoints additionally
  nest only the compact accepted `SWMP` state and rebuild its control/full/
  condensed topology from the restored Structure before decode. The derived
  accepted MAC field, mimetic topology, and transient bulk pressure are not.
- `simwing_fluid`: Qt-free periodic staggered-grid field operators and
  transactional pressure projection verification kernel. Its isolated
  `fluid/mimetic_local_cell.*` mixed-hybrid building block constructs a
  fingerprinted SPD inverse flux inner product from complete oriented
  half-face geometry, enforces area-vector/divergence-theorem closure and
  linear consistency, eliminates the cell scalar conservatively, and reduces
  exactly to `area / centre-distance` on condensed Cartesian neighbors. Its
  exact diagonal-plus-two-rank-three representation stores seven doubles per
  half-face and applies matrix-free; every closed coarse real-wing scene shell
  builds it without dense local storage. A linear-consistency rejection now
  carries its geometry scale, face-area range, closure/moment residuals,
  condition estimates, stabilization scale, measured algebraic residual, and
  unchanged tolerance as a typed error; the global trace builder adds the
  exact control/grid-cell/region identity before propagating it to offline
  diagnostics. Its
  `fluid/mimetic_wall_condensation.*` companion exactly eliminates selected
  zero-flux wall rows through the trace operator's diagonal-plus-seven-low-rank
  form, one equilibrated `7 x 7` Woodbury inverse, and linear face storage. A
  numerically singular auxiliary core conditionally falls back only when the
  independently assembled principal wall block has at most eight faces and is
  invertible; it retains the exact seven-mode Schur metric and reconstructs
  later small solves without persistent dense storage. It
  publishes an exact active Schur action/diagonal, RHS condensation, and wall
  reconstruction without a dense wall matrix. All 138 coarse controls accept
  it, reducing 148,652 wall traces to 42,927 shared traces with 3,986,602 bytes
  of Schur data and positive assembled reduced diagonals. It does not alter the production
  pressure graph. The target also includes a
  validated ordered multi-crossing sharp pressure-jump field and
  a calibrated flux-driven Darcy-Forchheimer crossing adapter with explicit
  per-tile volume-flow and dissipation ledgers. A bounded transactional Picard
  solve closes that nonlinear endpoint or midpoint law against the periodic
  sharp projection with optional separately prescribed jumps. It also owns
  fixed-topology, face-aligned moving-interface constraints. Its first open planar
  control-volume operator closes partial-cell geometry, surface sweep, and
  resolved-opening transport, then transactionally rebases by one MAC face at
  an exact crossing. A separate bounded planar cut-surface operator places the
  projection's face-resolved complete constraint reaction on the congruent
  physical plane, can resample that macro-step-average reaction at endpoint
  kinematics, and closes area, force, moment, power, and periodic-image ledgers.
  Its in-memory checkpoints capture accepted moving-only pressure, velocity,
  interface topology, and diagnostics, and separately capture accepted
  moving/porous state including the predicted-field provenance required to
  reconstruct midpoint samples, canonical resistance definitions, prescribed
  jumps, and complete nested diagnostics. Both are immutable and bound to exact
  grid/topology fingerprints. Distinct deterministic bounded/checksummed
  little-endian codecs persist the moving-only state and the complete
  moving/porous composition; the latter nests and revalidates the former rather
  than duplicating its interface/field wire contract. Accepted porous samples
  also map to deterministic stable-ID face/surface traction diagnostics with
  equal-and-opposite fluid/sheet force and impulse, distinct fluid/sheet power
  and work, and explicit porous-dissipation energy closure; separately
  prescribed jumps remain outside that sheet-load ownership.
  Its viscosity operators are a bounded explicit periodic MAC Euler oracle and
  a two-stage SSPRK2 path, both with the sharp
  `nu*dt*sum(1/h^2) <= 0.5` per-stage stability contract, exact zero/uniform
  modes, and accepted momentum/energy ledgers. The latter has observed
  second-order temporal convergence for the discrete Fourier eigenproblem.
  Its transport baseline advects all three MAC components with conservative
  donor-cell updates. The prescribed-uniform oracle is
  maximum-principle bounded for total absolute CFL at or below one, preserves
  periodic momentum and solenoidal modes, and delegates bit-exactly from a
  variable-flow operator. That operator reconstructs shared staggered
  control-volume fluxes from a divergence-free MAC field, supports nonlinear
  self-advection, and preserves periodic component momentum, bounds, and
  non-increasing energy under its local outgoing-CFL limit. Both donor forms
  remain first order. A selectable monotonized-central (MC) reconstruction
  supplies the conservative limited higher-order path. Its Euler stages may add
  the expected
  `O(dt^2)` energy only inside an enclosing SSPRK2 transaction; the committed
  update still enforces old-time bounds, momentum, and non-increasing energy.
  Smooth uniform-flow L1 refinement approaches second order.
  A projected SSPRK2 nonlinear operator composes two such self-advection stages
  with intermediate/final pressure projections. It preserves transactional
  pressure/velocity commit, stage eligibility, momentum, and non-increasing
  energy, has observed second-order temporal refinement on a fixed grid, and
  accepts either donor-cell or limited MC reconstruction. A Galilean-translated
  Taylor-Green vortex with `dt` proportional to `h^2` observes first-order donor
  and near-second-order limited-MC L1 spatial refinement on 16/32/64 grids. It
  does not establish cut-cell or moving-interface accuracy.
  The composed periodic evolution step selects uniform or nonlinear transport
  and Euler or SSPRK2 viscosity, then runs the zero-mean projection on
  candidates and commits velocity and pressure together only after every stage
  and the aggregate conservation ledger pass.
  A second-order temporal path instead uses symmetric Strang splitting:
  half-step SSPRK2 viscosity, full projected nonlinear SSPRK2 transport, and
  the matching viscous half step. Its exact sub-integrator ledger and
  fixed-grid refinement are tested; an independent continuously viscous,
  Galilean-translated Taylor-Green solution verifies near-second-order L1
  refinement of the complete limited-MC path with `dt` proportional to `h^2`.
  Donor-cell remains the default and limited MC reconstruction is selectable.
  A bounded subcycling wrapper advances one requested outer interval through
  equal Strang steps. It sizes the first schedule from the viscous limit and
  restarts the whole private interval only after an explicit CFL, limited-MC
  maximum-principle, or viscous stability rejection. Projection and ledger
  failures are fatal; at most 4096 substeps are allowed, and the caller fields
  commit only after an independent outer momentum/energy ledger passes.
  Projected SSPRK2, the composed first-order step, Strang splitting, and its
  retrying subcycler all accept one immutable sharp-jump field for the complete
  interval. Every internal projection sees the same ordered crossings; an
  empty field is bit-exact to the original path. The balanced slab remains
  pressure-sharp without spurious velocity through every stage and substep.
  Powered/moving jump work is not hidden inside the no-added-energy periodic
  ledger. It has no Playground dependency and is not yet a dynamic arbitrary-
  interface or whole-wing flow solver.
- `simwing_fluid_frame`: Qt-free owning adapter from one accepted periodic MAC
  state to immutable cell-centre diagnostic points with pressure, speed,
  velocity, exact finite-volume divergence, and diagnostic centred-curl
  vorticity fields. It never aliases the solver arrays.
- `simwing_pressure_jump_case`: Qt-free visible static split-slab worker. Its
  frame adapter retains one oriented quad for every ordered subcell crossing,
  including multiple layers on one grid face, alongside owning cell pressure,
  velocity, speed, and divergence samples. It is a diagnostic oracle, not a
  moving folded-interface or cut-cell evolution case.
- `simwing_porous_flow_case`: Qt-free pressure-driven porous plug worker. An
  implicit-midpoint Darcy-Forchheimer momentum solve closes pressure impulse,
  driving work, porous dissipation, and kinetic energy before the accepted
  endpoint pressure circuit crosses a complete Strang/SSPRK2 fluid step and
  immutable layered frames. It remains a uniform one-degree-of-freedom
  verification oracle.
- `simwing_periodic_flow_case`: Qt-free Taylor-Green verification worker that
  advances the bounded Strang/SSPRK2 subcycler and publishes accepted fluid
  frames. Its immutable in-memory checkpoint binds the exact grid, case
  definition, fields, last accepted diagnostics, step, and time. Its bounded,
  checksummed little-endian file codec round-trips that complete payload for
  `--checkpoint-in`/`--checkpoint-out`; resumed `--steps` are additional
  intervals. Optional `--checkpoint-every N` writes at absolute accepted-step
  multiples and the final accepted state, always through the same atomic
  replacement path. It is not a canopy, cut-cell, or aerodynamic-truth case.
- `simwing_worker_control_protocol`: Qt-free transport-neutral binary commands
  and responses for safe-point worker orchestration. Version 1 bounds complete
  messages, advance counts, and error text; request IDs correlate advance,
  checkpoint, and stop commands with ready, advanced, checkpointed, stopped,
  or error responses. Process transport and command execution remain separate
  integration layers.
- `simwing_worker_control_stream`: Qt-free self-framing stream adapter for the
  control messages. It reads the envelope's bounded little-endian payload size,
  adds no native prefix, rejects partial/oversized input transactionally, and
  flushes every complete output message for pipe use.
- `simwing_worker_control_session`: Qt-free case-neutral safe-point executor.
  Accepted advances publish immutable frames through an injected sink,
  checkpoint persistence is an injected action, and stop is terminal. It owns
  no pipe, socket, file, viewer, or scheduling policy.
- `simwing_periodic_flow_control`, `simwing_moving_porous_flow_control`,
  `simwing_open_piston_control`, `simwing_strong_piston_control`, and
  `simwing_porous_sheet_control`: typed
  adapters binding their numerical owners and complete checkpoint payloads to
  the shared control session.
- `simwing_viewer_geometry`: Qt-free deterministic diagnostic-vector glyph
  builder. It owns bounded arrow-segment output and links only the viewer
  protocol; the Qt/OpenGL viewer only uploads the resulting geometry.
- `simwing_viewer_protocol`: Qt-free immutable diagnostic-frame and trace
  protocol shared by future workers and the standalone viewer.
- `playground_contact`: Qt-free bounded Playground cloth-contact features,
  topology exclusions, projection and coverage diagnostics.
- `playground_sim`: widget-independent mesh parser/body builder/step and shape
  metrics; linked by both GUI and headless tests/tools.
- `xfoil` and `xflr5core`: vendored aerodynamic application libraries.
- `LEparagliding`: Qt GUI, linked to OCCT, XFLR5, Playground, and flat parts;
  it depends on the engine executable because previews launch it.
- `softwing-bench` and `nesting-bench`: developer tools, not shipped.
- The remaining small executables in `CMakeLists.txt` are test drivers.

When adding a production source file, add it to the appropriate explicit CMake
target. XFLR5 is the exception: its source list is generated from upstream
qmake manifests into `cmake/xflr5_sources.cmake` by
`cmake/gen_xflr5_sources.py`; do not replace that list with a glob.

## Directory and ownership map

### Calculation engine

- `src/engine/main.cpp`: CLI parsing, path validation, output cleanup, engine
  orchestration, OCCT model write, and companion `lep-sim.json` write.
- `src/engine/input_migration.{h,cpp}`: non-destructive calculation copy. It
  removes the Studio history trailer and blank records, and appends disabled
  sections 33-37 to an older design that has section 32. It rejects gaps.
- `src/engine/engine_paths.{h,cpp}`: UTF-8/native path bridge, safe output path
  handling, Fortran-shaped string helpers, `fopen` wrappers, and the planform
  row parser used by the translated core.
- `src/legacy/leparagliding_core.cpp`: mechanically translated 3.28 numerical
  and drawing core (about 49k lines). It is compiled as permissive C++ by MSVC
  and as C on GCC/Clang. Treat Fortran indexing, formatted I/O, execution order,
  and C89 declaration rules as compatibility constraints.
- Root `leparagliding3.28.f` is the active original-source oracle;
  `leparagliding3.17.f` is retained historical input. `Data/Originals` contains
  upstream/source artifacts, not application code.

The core produces the traditional calculations and drawings, but modern
outputs enter through narrow C callbacks. Keep those callbacks C-clean because
the core is a C translation unit on non-MSVC platforms.

### Exact model and viewport

- `src/model/nurbs_model.{h,cpp}` owns the capture state and all OCCT model
  construction. Legacy callbacks provide transformed profiles, analytical
  ballooning data, legacy tessellation for comparison, ribs/holes, diagonals,
  mini-ribs, and labelled line segments.
- `NurbsModel` constructs exact rational span curves and lofted
  Extrados/Vents/Intrados faces, mirrors the half wing, sews topology while
  preserving designed intake openings, adds ribs/holes/diagonals/lines to a
  named XCAF assembly, validates source/legacy agreement, and writes AP242 STEP
  or binary XCAF.
- `writeSimMesh` serializes a deliberately coarse full-wing JSON mesh in
  millimetres, Z-up: welded skin quads with surface tags, rib/hole loops,
  diagonal straps, and labelled suspension lines.
- `src/gui/paraglider_view.{h,cpp}` imports STEP or XCAF, triangulates with
  OCCT, hosts the native OCCT OpenGL viewer, maps XCAF labels to stable colour
  roles, and owns camera/navigation, selection, visibility, clipping, X-ray,
  and measurement. There is no second application-side 3D mesh builder.
- `ParagliderView::ColorRole` values are persisted by ordinal. Append roles;
  never reorder existing values or their `colorSettingsKeys` in
  `mainwindow.cpp`.

Model units are millimetres. The Playground converts the simulation JSON to
metres at its boundary.

### Studio shell, document, and editors

- `src/gui/main.cpp`: application bootstrap and command-line/self-test modes.
  `--headless` is detected before creating the full GUI application path and
  delegates to the sibling engine executable.
- `src/gui/mainwindow.{h,cpp}`: application composition and workflow. It owns
  file/preset opening, editor pages, QProcess runs, preview/export/XFLR5 modes,
  model and output refresh, part-to-source navigation, preferences, settings,
  and the major tabs.
- `src/gui/design_document.{h,cpp}`: parsing and persistence of numbered
  design sections plus embedded version history and Studio-only spline data.
- `src/gui/section_help.*` reads generated manual chapters; `section_specs.*`
  describes the dominant tables and the few generic curve-enabled sections.
- `src/gui/section1_curves.*`, `spline_fit.*`, and `section_grid.*` are pure
  C++ parser/math helpers with focused tests. Keep reusable parsing out of
  widgets when practical.
- `curve_editor.*` is the generic interactive plot. `section1_curve_panel.*`
  handles section 1 B-spline truth; `airfoil_panel.*` handles section 2 airfoil
  files/splines; `holes_panel.*` handles section 4 hole geometry;
  `grid_curve_panel.*` currently exposes spec-declared curves (section 3 and
  section 30). `section_grid_panel.*` remains built/tested but is intentionally
  not instantiated: the generic value-grid UI was retired.
- `geometry_preprocessor.*` is a pure C++ port of the separate v1.6 Canigo
  preprocessor. `geometry_preprocessor_dialog.*` is its UI and patches the
  resulting rib rows into section 1.

Each section has its own `QPlainTextEdit` undo stack. Once that live stack is
exhausted, MainWindow continues through the saved per-section history. Enter
builds the current in-memory design; Shift+Enter inserts a record. Preview and
export do not implicitly save.

`DesignDocument` recognizes headers with `* <number>. <title>`, preserves the
preamble and final-newline choice, rejects duplicate sections and interior
blank records, and assembles text in section order. Save uses `QSaveFile` and
appends a comment trailer delimited by:

```text
* >>> LEPARAGLIDING STUDIO HISTORY V1 >>>
...
* <<< LEPARAGLIDING STUDIO HISTORY V1 <<<
```

The trailer is base64 compact JSON. It contains a SHA-256-linked chain of
qCompressed full-wing snapshots (UTC time, parent, changed sections, summary)
and optional Studio spline JSON. Restoring an old version does not delete newer
history; saving the restored payload appends a new latest revision. The engine
always strips this trailer from its temporary calculation input.

Shipped presets are read-only masters. Saving one asks for a new design path
outside `resources/presets` and copies adjacent resource/airfoil files with it.

### Preview, export, and XFLR5 flow

`MainWindow::startCalculation` always writes the current editors to a temporary
`leparagliding.txt` and passes the real design directory as `--resource-dir`, so
relative airfoil references resolve correctly.

- Preview writes to a temporary output directory and adds `--preview`; the
  engine writes `lep-3d.xbf`, which the viewport loads faster than STEP.
- Export writes to the selected output directory and keeps the standard STEP.
  The construction-curve preference affects export only.
- XFLR5 transfer creates a temporary design with section 36 forced on, runs the
  engine, and calls the locally added `MainFrame::lepImportLepWing` on the
  generated `xflr5` folder. Transfers are hash-deduplicated and queued behind
  an in-flight preview.
- After a successful run, `lep-sim.json` is copied into Playground state and
  `lep-2d-parts.json` is loaded into Print state before the temporary directory
  disappears.

The six user-facing export files are:

- `leparagliding.dxf`: 2D manufacturing plans.
- `lep-3d.step`: exact named OCCT model.
- `lep-3d.dxf`: legacy reference wireframe.
- `lep-out.txt`: calculation report.
- `lines.txt`: suspension data.
- `run-log.txt`: calculation progress/diagnostics.

Optional legacy settings can additionally create `stl/`, `xflr5/`, and surface
DXFs. `lep-3d.xbf`, `lep-sim.json`, and `lep-2d-parts.json` are internal Studio
handoffs, although they are useful developer fixtures.

### Flat parts and Print/Cut

- Hooks declared in `src/flatparts/flat_capture.h` are inserted around the
  translated core's 2D drawing loops. `flat_capture.cpp` attributes otherwise
  anonymous DXF primitives to stable categories/pieces, chains line segments,
  classifies cut/seam/mark roles, converts to local millimetres/y-up, and writes
  `lep-2d-parts.json`.
- `flat_parts.*` parses that JSON into `FlatPartSet`/`FlatPiece`.
- `nesting.*` performs deterministic irregular strip packing. It rasterizes
  concave outlines to bitmasks, tries rotations and orders, minimizes tiled
  page count, reports only improving layouts, and cooperatively checks cancel.
- `nest_worker.*` runs nesting on a `QThread` with generation IDs so stale
  queued results cannot overwrite a newer run.
- `flat_parts_view.*` draws either the selected parts or the exact packed
  layout. Preview and exporters share `frameFor`; do not independently
  reimplement placement transforms.
- `sheet_export.*` tiles paper with overlap/furniture into a true-size PDF or
  writes R12 DXF (continuous canvas or one machine-bed file per sheet).
- `print_page.*` owns selection, scale/target area, paper/bed/rotation options,
  the live pack, and export controls. Export uses the captured packed options,
  not whatever widgets have changed since the pack.

Known Print correctness gaps are documented in
`docs/legacy/leparagliding/BACKLOG.md` and
`docs/legacy/leparagliding/flat-part-orientation.md`: `grainAngleDeg=90` with
`grainSource=assumed-from-plan-layout` is a placeholder, not engineering data;
and `outerBoundary()` misses small extents of some cut geometry. Quarter-turns
avoid bias but do not prove the intended grain direction. Never enable free
rotation for direct fabric cutting. Seam allowances remain fixed millimetres
when the wing itself is scaled.

### Playground and soft-wing simulation

- `src/gui/playground_sim.{h,cpp}` is the widget-independent adapter: parse and
  refine `lep-sim.json`, build the XPBD body/constraints/ribs/lines/pilot/cells,
  stamp pressure, impose the wing polar, step, contact, brakes, grabbing, and
  free-flight recentring.
- `playground_pressure_solve.*` is the Qt-free deterministic bounded weighted
  equality projection for the final exterior-Cp field. Force, pitch and the
  two half-differential rows enter in priority stages; it reports physical
  saturation separately from numerical failure.
- `playground_metrics.*` captures the rest-shape baseline and measures section
  deviation, span/area/volume, slack, asymmetry, dents/twist, line-row loads,
  cell/kink/collapse diagnostics, settling, heatmap fields, and CSV.
- `playground_page.*` is the lazy-created OpenGL tab, controls, live solver,
  HUD/session log, settle run, and analysis entry point.
- `playground_analysis.*` runs alpha sweeps on a worker using a fresh body per
  point and the same settle/measure path as the bench.
- `tools/softwing_bench.cpp` is the authoritative headless driver for the same
  `playground_sim` code. `tools/softwing_gpu.*` is an experimental benchmark
  backend, not the GUI simulation path.
- `src/softwing` supplies generic XPBD distance/cable constraints,
  deterministic coloured parallel sweeps, membrane/contact/pneumatics/canopy/
  suspension machinery. The Playground defaults to its calibrated legacy
  distance-truss skin and also offers an explicitly experimental orthotropic
  membrane skin with compression softening and four-node dihedral hinges.
  Free flight additionally requests three deterministic cable-only
  reverse/forward sweep pairs to propagate payload load through the deep
  suspension graph; zero is the exact historical/tunnel path.

The Playground is an instrumented relative/structural design sandbox, not a
certified absolute aerodynamic solver. Wind-tunnel mode is the default claim;
free flight is explicitly experimental. Simulation time is fixed at 60 Hz and
can advance more slowly than wall time. Worker count is explicit because
reproducibility is per worker count.

Before changing Playground forces, pressure, flight trim, or controls, read all
of `docs/legacy/leparagliding/CONTINUE.md` and
`docs/legacy/leparagliding/playground-shape-analysis.md`. The essential
invariants are:

1. System-level forces enter through the actual suspension/load path; prior
   point/body/couple shortcuts visibly damaged the canopy.
2. The wing polar fully cancels the pressure field's resultant. A new pressure
   request belongs in `wingForce`, or cancellation erases it. The one explicit
   post-solve exception is bounded-free-flight viscous skin traction: only the
   positive `q*S*Cd_polar - achievedPressureDrag` deficit, never missing lift
   or the excess-frontal fabric heuristic.
3. Do not feed per-node breathing velocity back into pressure. Rotational
   feedback uses the canopy rigid-body fit.
4. Verify angle-of-attack sign dynamically: sinking must raise alpha.
5. Free-flight damping is relative to bulk velocity or it becomes fake drag.
6. Drag must decelerate, never propel.
7. Solver left is negative mesh X (the viewer's right). The UI crossover lives
   in `setBrakePull`; bench flags and `SimBody::ribHalf` use solver sense.
8. Constructor-time widget values do not guarantee controls reached the view;
   add new defaults to `ensureView`'s explicit push list.
9. Production retrim solves the FINAL exterior Cp inside `[-3, 1]`, then
   reconstructs `p_inside-q*Cp`; never reintroduce an unbounded increment plus
   a post-solve clamp. q≈0 faces have no authority and stay at their prior.
10. The objective prior is the same-frame legacy 4x4 load-path proposal, but
    every intermediate pressure is clamped through the per-face physical
    interval `[p_inside-q, p_inside+3q]`. Never reuse the legacy absolute clamp
    here: its interval can invert during the pre-inflated soft start.
11. Active faces stay within +/-0.05 Cp of that preferred topology. Production
    authority is cached and re-verified on live geometry: halve an infeasible
    hint, then try one +0.02 probe. Force starts at 0; normally feasible pitch
    and differential start at 1. Preserve the cached normalised dual warm
    starts. Exact bisection belongs to pure tests/offline audits, not every
    full-canopy frame.
12. Force rows freeze before pitch, which freezes before the common L-R
    lift/drag authority. Preserve the prior half difference and the explicit
    factor two for +D/-D half loading.
13. Pinned bounded lift/drag telemetry reports the achieved pressure force.
    Free-flight telemetry adds separately diagnosed polar-only skin traction
    and pilot drag. Traction must be q-gated, area-weighted, dissipative, zero
    in pinned/legacy modes, and must never reinject `lastFabricDragNewtons`.
    The old symmetric tunnel row is an explicit regression oracle. Use
    `--legacy-pressure`; never make the legacy path the production default.
14. TrimmedGlide q comes from achieved vertical wing support versus total
    dynamic-node weight, followed by eight bounded-only co-moving structural
    relaxation frames and a state-clean recalibration. DropFromRest and legacy
    launch do not receive this relaxation. Keep the 30/90/250 kg support/cap
    regression and report calibration residuals/iterations.
15. Three free-flight cable sweep pairs are targeted load-path conditioning,
    interleaved before general structural sweeps. Tunnel and zero-pair paths
    preserve historical arithmetic exactly; do not replace this with hidden
    line ballast or global extra cloth iterations.

Weathercock incidence and the bounded final-Cp retrim are now implemented.
Remaining free-flight stability and material-model limitations are recorded in
`docs/legacy/leparagliding/CONTINUE.md`; do not infer that bounded authority
makes this a certified aerodynamic solver.

### SimWing FSI foundations

- `src/fsi/scene.{h,cpp}` owns scene-v2.2's SI/right-handed/Z-up contract,
  stable-ID entities, authored fabric-sheet identity, two-sided fluid regions,
  materials, openings with optional oriented boundary-vertex cap disks, paired
  seam chains, pilot, suspension junctions, attachments, suspension lines,
  validation, and bounded deterministic binary round trips. Binary v3 writes
  cap topology and still reads cap-free binary v2. Diagonal/mini-rib sheets may
  retain one connected cell on both sides; skin and ribs may not.
- `src/model/nurbs_model.{h,cpp}` can export scene-v2 directly from analytical
  captures without reading `lep-sim.json`. It preserves authored open intakes,
  exact triangulated rib/crossport faces, internal sheets, and the segmented
  suspension graph. Open intakes now carry deterministic nonplanar
  boundary-vertex cap disks. Their spanwise lips reuse the coarse skin grid;
  their side chains reuse the actual OCCT rib-triangulation boundary, with a
  maximum 0.25 mm sub-mesh canonicalization when a captured lip corner falls
  between rib mesh samples; larger disagreement rejects export. The
  real-design regression verifies complete cap construction through the
  three-region skin/rib junctions, cap orientation, exact fabric boundary-edge
  reuse, and live Structure motion for every opening vertex. Physical
  fabric/line/pilot settings are mandatory and no engine CLI writes this scene
  yet because the design format has no authoritative source for them. Its
  present local intrinsic charts are not manufacturing flat-pattern UVs, and
  it does not yet author paired seam chains. The cap owner now consumes the
  complete real-wing opening set. Its full multi-region volume ledger closes
  on a centered coarse grid that crosses the canopy. One simple
  boundary-to-boundary interface is now an exact oriented partition; grid-face
  junctions remain preserved as pair-specific chains in the material-only
  partition, while exact multi-axis clip-plane identity lets the capped
  arrangement close all nine touched faces.
  Opening quadrature and grid patches retain the full cap area, and authored
  connectivity plus sparse pressure control volumes assemble. Pressure-link
  construction now retains the 24 coarse embedded patches across two mirrored
  intake openings as typed non-positive-distance rejections with source
  identity and signed centroid-to-cap-plane distances. Both pressure centroids
  lie on the positive side of each local cap plane, so never hide this
   non-admissible two-point stencil with an absolute value. It publishes no
   fabricated conductance link and is rejected by pressure-operator assembly.
   All 24 coarse rejections, and all 294 in the refined 4-by-4-by-4 audit,
   retain correctly sided same-region Cartesian one-ring support on both
   authored sides with periodic images unwrapped. Donor substitution would
   violate local row conservation. The isolated mixed-hybrid mimetic
   local-cell kernel now proves the SPD/conservative formulation and exact
   Cartesian equivalence on manufactured closed cells. The immutable scene
   shell audit now combines Cartesian region subfaces, material-wall
   quadrature, and opening-cap patches and rejects area-vector or divergence-
   theorem nonclosure. Oriented closed-surface classification resolves ten
   coarse and six refined untouched faces as Outside, closing all 138 coarse
   and all 358 refined controls without dominant-volume or closure fitting.
   The global audit now pairs all shared traces, retains one zero-flux wall
   trace per material half-face and one gauge per pressure component, and
   applies the symmetric positive-semidefinite condensed operator matrix-free.
   The coarse real-wing system has 191,579 trace unknowns and 13,132,336 bytes
   of compact local factors. Its component-constant action is roundoff-null.
   The bounded gauge-fixed Jacobi-PCG solve now recovers manufactured fields
   on both full and material-wall-condensed systems, reconstructs a solution
   that closes every original row, closes compatible source-driven fluxes, and
   executes finite residual-reducing iterations across both real systems with
   exact rollback when truncated. The condensed real solve also reaches
   `1e-5` relative RMS within 300 iterations, reconstructing all 191,579 rows
   below `2e-4` maximum residual. Exact local wall-trace condensation accepts
   all coarse controls; its action is fused into the exact seven-mode Schur
   form, and its global 42,927-row shared system is assembled with full wall
   reconstruction. The atomic integrated-source transaction also converges a
   balanced real source pair and publishes only after reconstructed full-space
   closure. The physical-unit source product now applies the production
   `-(rho/dt)` sign to sampled fixed or material-wall-adjusted trace flow plus
   accepted GCL volume rates. Its accepted source-bound result now captures an
   independently reconstructed immutable pressure state, and that state maps
   across one accepted topology transition into an exact-gauge reduced-trace
   warm start. Stronger preconditioning and production integration remain open.
   The production operator and subsequent worker path are unchanged.
- `src/fsi/structure.{h,cpp}` is the new Qt-free XPBD boundary. It owns nodal
  loads/state, trusted constraint/membrane/bending assembly, explicit optional
  fabric self-contact, rigid-pilot suspension, diagnostics, and composite
  rollback checkpoints without including any Playground header.
  `structure_checkpoint_persistence.{h,cpp}` composes the SoftBody and optional
  suspension codecs with the public step/time/load state. Both encode and
  decode validate through an equivalent rebuilt `Structure`; decode also
  rejects disagreement between public node state and the opaque body payload
  before publishing transactional output.
- `src/fsi/scene_structure.{h,cpp}` deterministically assembles supported
  scene-v2 fabric, suspension junctions, cable segments, and one rigid
  pilot/harness tree into that XPBD boundary. It derives fabric mass from SI
  material-chart area, creates signed dihedrals only across consistently
  oriented edges within one authored sheet, and enables contact only from
  explicit worker-supplied settings. It still rejects seams until their stitch
  and tributary load-sharing model is implemented.
- `src/fsi/scene_fluid_surface.{h,cpp}` extracts stable-ID-sorted compact fluid
  regions, used porous materials, oriented two-sided triangles, ordered
  openings, and optional authored cap disks directly from validated scene-v2,
  then captures one immutable
  fingerprinted accepted Structure epoch in that surface order. Opening-only
  construction vertices are rejected until they have trusted structural motion; grid
  crossing, cut-cell classification, and pressure ownership remain downstream.
- `src/fsi/scene_fluid_surface_transfer.{h,cpp}` converts that exact oriented
  triangle set into the topology-bound `ConservativeSurfaceTransfer`, validates
  every accepted state against both scene-surface and Structure fingerprints,
  and delegates uniform or barycentric-quadrature traction and load application.
  It does not choose traction or discard two-sided fluid ownership.
- `src/fsi/fluid/scene_surface_geometry.{h,cpp}` conservatively bins every
  current authoritative triangle AABB into canonical cell-major candidates,
  including both cells when zero-thickness geometry lies on an internal grid
  plane. The result is byte/count bounded and bound to the complete accepted
  surface-state fingerprint. It rejects out-of-domain surfaces; periodic-image
  selection, face crossings, volume fractions, and region reconstruction
  remain explicit future stages.
- `src/fsi/fluid/scene_surface_intersection.{h,cpp}` filters those candidates
  with normalized separating-axis triangle/AABB tests. Exact contact is retained
  on both adjacent cells, optional separation tolerance cannot exceed the broad-
  phase padding, results are count/byte bounded and fully epoch-bound, and
  validation recomputes every expected pair. It still assigns no clipped
  polygon, volume fraction, face crossing, or fluid region.
- `src/fsi/fluid/scene_surface_clipping.{h,cpp}` clips each exact pair against
  all six cell planes while carrying original-triangle barycentric coordinates.
  It publishes flattened point/segment/area patches, exact area and centroids,
  and coincident-boundary masks. Coordinates already shared by both endpoints
  remain bit-exact through later-axis interpolation so prior plane identity is
  not lost. Material and opening-cap clippers derive both sides of every cell
  directly from `gridLower + faceIndex * spacing`; never form an upper plane as
  `cellLower + spacing`, because adjacent cells can then disagree through
  non-associative floating-point addition. A triangle lying on a shared grid
  plane stays duplicated and explicitly flagged on both cells; no traction
  integration may consume both before a unique face-ownership stage. Volume fractions, face
  crossings, and region labels remain unassigned.
- `src/fsi/fluid/scene_surface_ownership.{h,cpp}` retains ordinary positive
  area as cell-owned patches, pairs exact lower/upper duplicates into one
  canonical MAC-face owner, carries authored negative/positive region IDs and
  the triangle-normal axis sign, and excludes point/segment contact from area.
  Unpaired positive area on the periodic domain boundary is rejected until an
  explicit periodic-image owner exists.
- `src/fsi/fluid/scene_surface_crossings.{h,cpp}` extracts positive-length
  boundary segments from ordinary owned cell patches, pairs independently
  clipped adjacent-cell positions only within a fixed machine-roundoff
  envelope, requires matching barycentric zero/nonzero provenance rather than
  redundant nonzero-weight equality, publishes one stable lower-cell
  crossing, and carries the in-face
  negative-to-positive direction plus authored region/material/sheet identity.
  Unpaired triangle edges remain contact diagnostics, coplanar area remains
  face-owned, and grid-edge-aligned crossings await an explicit edge owner.
- `src/fsi/fluid/scene_surface_face_topology.{h,cpp}` groups those transverse
  crossings and coplanar owners by canonical internal MAC face with stable
  grid-bound face IDs and bounded flattened references. It preserves multiple
  sheets independently; summed length and area are diagnostics, not union
  coverage or inferred region labels.
- `src/fsi/fluid/scene_surface_face_graph.{h,cpp}` stitches face-local segment
  endpoints by stable scene-vertex/edge provenance, recomputes shared-edge
  intersections canonically, derives grid-edge endpoints from the authored
  triangle plus the two exact Cartesian planes, admits only a fixed
  coordinate-ULP envelope when checking canonical nodes against clipped
  inputs, retains authored opening boundaries, and publishes bounded node
  degrees/connectivity. It does
  not require closed curves or infer a fluid region from an open graph.
- `src/fsi/fluid/scene_surface_face_chains.{h,cpp}` orients graph segments from
  triangle winding, extracts deterministic open chains and closed loops, keeps
  opening/grid-boundary endpoints explicit, and requires one authored region
  pair per chain. Branches, conflicting winding, and region changes are
  rejected before any face partition is inferred.
- `src/fsi/fluid/scene_surface_face_loops.{h,cpp}` measures simple closed
  chains in a right-handed face chart, rejects self-intersection and degenerate
  area, and uses winding to identify the enclosed versus exterior authored
  region. Open chains stay unresolved and separate loop areas are not treated
  as a union or final face partition.
 - `src/fsi/fluid/scene_surface_face_partition.{h,cpp}` rejects touching or
   intersecting loops, builds smallest-parent containment, requires authored
   region continuity through nesting, and closes exact per-region area and
   global first moment on faces
  containing only interior closed loops. It also closes simple directed
  open-chain arrangements whose leaves reach the rectangular face boundary,
  retains every source chain, rejects unstitched crossings/conflicting region
   winding, and enumerates exact positive-area faces and centroids through a
   bounded half-edge traversal. Polygon areas and moments are evaluated about
   local references before global centroids/first moments are published.
   Opening-ended chains, coplanar sheets,
  boundary-touching loops, and empty active faces remain explicitly unresolved.
- `src/fsi/scene_fluid_quadrature.{h,cpp}` turns each unique positive-area
  owner into one stable-ID barycentric quadrature point, preserves authored
  region/material/sheet/role metadata and the exact grid cell owning each
  side, binds the accepted surface epoch and transfer topology, and feeds
  ordered finite traction through the existing conservative load adapter. The
  same barycentric points sample accepted Structure position and velocity for
  the reciprocal CFD boundary. Shared-plane area is integrated exactly once;
  this boundary does not invent pressure or traction.
- `src/fsi/scene_fluid_pressure_traction.{h,cpp}` validates ordered one-sided
  pressure samples against that accepted quadrature epoch and converts their
  difference along the current oriented triangle normal before delegating to
  conservative transfer. Equal pressures cancel exactly; no shear, polar, or
  additional aerodynamic load is introduced.
- `src/fsi/scene_fluid_grid_epoch.{h,cpp}` transactionally composes every
  accepted surface/grid geometry stage and conservative quadrature into one
  immutable, fingerprinted epoch. Nested stage identities cannot cross
  Structure steps, each stage retains its own limits, and the completed vector
  payload has an aggregate byte ceiling. Explicit unresolved face partitions
  remain unresolved; this owner does not infer cell volumes or pressure.
- `src/fsi/scene_fluid_opening_cap.{h,cpp}` matches every boundary edge of the
  separating surface to exactly one authored opening and derives cap winding
  from the adjacent triangle's region order. Planar loops preserve the exact
  convex fan or receive bounded deterministic concave ear clipping. An authored
  boundary-vertex disk may instead retain nonplanar facets and their individual
  normals. Accepted motion retains triangle identity and rejects folds,
  intersections, and degeneracy. This immutable bounded fluid topology enters
  volume, opening flux, pressure connectivity, and momentum, but never
  Structure, material quadrature, or traction.
- `src/fsi/scene_fluid_opening_quadrature.{h,cpp}` gives every virtual cap
  triangle a topology-stable centroid sample of accepted position and linear
  vertex velocity. Per-point and per-opening `area*(velocity dot normal)`
  ledgers exactly close rigid material-plus-cap surface sweep for both planar
  and explicitly faceted nonplanar caps.
- `src/fsi/scene_fluid_opening_patch.{h,cpp}` clips those triangles with the
  same exact barycentric triangle/box primitive used by material geometry.
  Positive area is owned once: off-face pieces belong to cells, while paired
  coincident copies become one canonical non-periodic grid-face owner. Exact
  clipped polygons and area/sweep closure remain bounded; periodic-boundary
  image ambiguity rejects.
- `src/fsi/scene_fluid_opening_flux.{h,cpp}` evaluates one immutable MAC field
  over those owners. Face-owned patches read their exact normal degree of
  freedom; cell-owned polygons use deterministic degree-three triangle
  quadrature of periodic staggered interpolation. Per-patch/opening ledgers
  retain fluid flow, cap sweep, and signed negative-to-positive relative flow.
  This diagnostic never mutates projection, pressure, connectivity, or loads.
- `src/fsi/scene_fluid_cell_volume.{h,cpp}` consumes one validated grid epoch
  and clips the signed tetrahedral chain induced by every oriented interface
  triangle and the grid origin into exact Cartesian cells. It publishes
  positive sparse region volumes, includes wholly interior cells without
  requiring closed face-local contours, enforces every cell's full volume, and
  cross-checks their global sum against a separate closed-surface calculation.
  Its nested opening-cap epoch can close the supported authored-mouth subset
  without treating a cap as fabric. It rejects unsupported opening geometry
  and non-manifold/inconsistently wound surfaces rather than guessing closure.
  Decomposition closure and sparse positive-region publication have separate
  tolerances. Defaults preserve production arithmetic; the offline refinement
  audit opts into smaller published slivers, cell-local first moments, and a
  bounded centroid repair for cancellation-scale complements.
- `softwing::SuspensionSystem` now exposes a validated transactional production
  checkpoint for its complete suspension and rigid-payload state, bound to the
  registered body/contact topology. `softwing::SoftBodyCheckpoint` captures the
  complementary node/constraint/membrane/bending/contact state, including
  contact warm starts; `simwing_structure` composes both transactionally.
  `softwing/checkpoint_persistence.h` serializes the opaque body checkpoint's
  mutable state through a deterministic, checksummed, byte-bounded
  little-endian envelope. Decode overlays only onto an equivalent live
  topology template, so wire data cannot redefine masses, connectivity,
  materials, or contact registration.
  `softwing/suspension_checkpoint_persistence.h` independently preserves the
  complete suspension/rigid-payload checkpoint, including stable identities,
  controls, multipliers, ledgers, diagnostics, and bounded UTF-8 text. Decode
  binds counts and identities to an equivalent checkpoint template, verifies
  the complete state fingerprint, and leaves final semantic validation to the
  live owner's transactional `restore()`. `simwing_structure` composes both
  payloads into its enclosing persistent format.
- `src/viewer/viewer_protocol.{h,cpp}` owns immutable sampled diagnostic frames
  and replayable traces. It uses nonzero 64-bit stable entity/region IDs,
  transactional decoding, configurable limits, and a 256 MiB default encoded
  frame ceiling.
- `src/viewer/structure_frame.{h,cpp}` maps committed structure state to those
  immutable frames, including rigid harness vertices and suspension segments,
  without inventing unavailable per-element quantities.
- `src/viewer/fluid_frame.{h,cpp}` copies one accepted periodic MAC state into
  owning cell-centre diagnostic points. Pressure is exact cell data; velocity
  is the arithmetic average of the two bounding MAC faces per component;
  divergence reuses the projection operator; vorticity is the centred periodic
  curl of that published velocity.
- `src/viewer/viewer_window.{h,cpp}` and `tools/simwing_viewer_main.cpp` form
  the separate Qt/OpenGL trace viewer. Replay and growing-trace follow remain
  bounded and asynchronous. Unreferenced diagnostic vertices render as
  scalar-coloured points without changing the trace protocol. Vertex vector
  fields have an independent selector and render through bounded arrow glyphs,
  so scalar colouring and vector direction can be inspected together.
- `src/viewer/vector_glyphs.{h,cpp}` builds those arrows outside Qt/OpenGL.
  It normalizes length by the field maximum, derives automatic scale from the
  dimensional point spacing, and caps large fields with a deterministic stride.
- `src/viewer/viewer_camera.{h,cpp}` owns the tight scene-relative perspective
  depth range. The OpenGL widget requests and verifies a 24-bit depth buffer,
  restores depth state after every HUD paint, and offsets filled triangles so
  coplanar diagnostic lines do not z-fight.
- `src/fsi/canonical_case.{h,cpp}` and `tools/simwing_fsi_main.cpp` are the
  first end-to-end worker slice. The case is an analytic structural harness,
  not aerodynamic truth; it writes only accepted steps and launches the
  sibling viewer by default. `--no-viewer` must remain Qt-free and unthrottled.
- `tools/simwing_mimetic_conductance_audit_main.cpp` is the opt-in Qt-free
  high-resolution evidence runner. It consumes the canonical offline profile
  and grid phases owned by the mimetic conductance phase-audit boundary. An
  all-phase run validates and discards one nested product at a time before
  compact deterministic aggregation; it must never select worker pressure
  ownership or apply structural loads.
- `src/fsi/hemisphere_case.{h,cpp}` is the larger structural/viewer canonical:
  a soft triangulated fabric hemisphere held at three equatorial points with a
  compliant rim, intrinsic membrane charts, signed rest-shape hinges, and a
  time-varying four-lobe analytic follower-pressure mode. The third
  non-collinear support removes the rigid rotation left by two positional pins.
  It deliberately tests curved loaded structure, not CFD truth.
- `src/fsi/projected_flag_case.{h,cpp}` is the first flexible CFD-to-XPBD
  worker canonical, selected with `--case flag`. It preserves one persistent
  periodic MAC field, applies sinusoidal cross-flow increments, projects them
  around a finite stationary same-region panel, and conservatively transfers
  the complete constraint reaction to a 5-by-5 membrane grid. Two fixed node
  rows encode clamp position and slope so the result bends instead of rotating
  as a rigid hinge. Deformed geometry is intentionally not returned to CFD;
  reference-surface power is zero and no two-way energy closure is claimed.
- `src/fsi/ram_air_cell_case.{h,cpp}` extends the fixed-reference diagnostic to
  one open rectangular fabric cell. Five stable-ID 4-by-4 wall surfaces map
  complete reactions through independent bridges into an 89-node shared-edge
  shell with per-sheet membranes and bending. Two mouth-perimeter rows encode
  clamp position and slope. Net mouth flux remains zero for the incompressible
  dead-ended cavity while the trace reports local RMS exchange. It deliberately
  does not return displaced wall geometry to CFD or claim moving-volume
  inflation.
  `--checkpoint-in`, `--checkpoint-out`, and `--checkpoint-every` are supported
  by periodic-flow, moving-porous-flow, open-piston, and porous-sheet workers;
  all restore before
  trace creation,
  advance additional steps, autosave only accepted absolute-step multiples,
  and atomically replace the final file without duplicating the last write.
  Periodic-flow, moving-porous-flow, open-piston, and porous-sheet
  `--control-stdio` are also
  Qt-free: stdin
  accepts only bounded binary commands, stdout contains only flushed binary
  responses, diagnostics stay on stderr, and the trace is completed before a
  successful stopped response.
  With `--checkpoint-in`, Ready reports the restored absolute step/time and the
  new trace contains only subsequently accepted frames.
  The real-export regression additionally runs the 3.28 fixture through direct
  scene export, structural assembly, a coupled accepted step, composite replay,
  and a completed viewer trace using explicitly synthetic physical settings.
- `src/fsi/piston_case.{h,cpp}` is the second worker case selected with
  `--case piston`. Its synthetic heavy piston makes the accepted motion visible
  while preserving the analytic tributary-mass translation. It is an
  end-to-end fixed-topology verification harness, not general moving-grid FSI.
  The same files contain the distinct `6 kg` `--case strong-piston` added-mass
  canonical. Each
  strong iteration restores the accepted physical epoch, projects the real
  `28.8 kg` fluid at the relaxed interface speed, transfers its pressure through
  the planar bridge and temporal adapter, advances XPBD, and returns the actual
  structural speed plus topology-bound residuals. Only convergence commits the
  new persistent Structure/fluid state. The endpoint pressure-force slope must
  recover the canonical `10.8 kg` discrete added mass; trapezoidal force
  integration then gives the fixed point with `10.8/2 kg` added to the `6 kg`
  structural denominator. Its in-memory checkpoint contains only accepted
  Structure and moving-interface fluid state, validates their shared velocity
  closure before no-throw commit, and replays the exact next result/frame.
  `strong_piston_checkpoint_persistence.*` wraps those accepted nested codecs
  in the deterministic checksummed `SWSPCKP1` envelope; it never stores an
  in-progress iteration or rejected attempt. The batch CLI accepts the normal
  checkpoint flags for this case. Its typed standard-I/O control binding
  publishes and checkpoints only completed accepted macro-steps.
- `src/fsi/porous_sheet_case.{h,cpp}` is selected with `--case porous-sheet`.
  Its analytic linear-resistance midpoint relation drives the same accepted
  nonuniform porous projection, sheet-reaction bridge, temporal transfer, and
  XPBD boundary used by later FSI. Every step closes the prescribed pump
  impulse/work against fluid and sheet momentum/kinetic energy plus porous
  dissipation. The plane may translate rigidly across six consecutive ordinary
  dual-cell boundaries, including one periodic wrap; each accepted midpoint
  sample commits the next complete MAC-face epoch. The following next-image
  pump-surface collision is a hard topology rejection. The constructor also
  supports the direction-reversed focused oracle with distinct case identity;
  persistence rebuilds and replays in the owner's direction. Its in-memory checkpoint
  keeps the nested Structure
  checkpoint, velocity, pressure, and accepted coupled diagnostics private,
  revalidates the canonical epoch before commit, and reproduces the exact next
  frame in the same or an equivalent rebuilt worker.
  `porous_sheet_checkpoint_persistence.{h,cpp}` wraps that state in the distinct
  `SWPS` envelope. It reuses Structure's persistent codec, stores the four MAC
  scalar fields and the complete versioned axis/face/periodic-image topology
  epoch explicitly, and uses bounded replay to
  regenerate diagnostics and require bit-identical structure/fluid state
  before transactional decode.
  The standard checkpoint CLI flags restore before trace creation, atomically
  replace output, use absolute autosave cadence, avoid a duplicate final write,
  and count `--steps` as additional intervals.
- `src/fsi/open_piston_case.{h,cpp}` is the third worker case selected with
  `--case open-piston`. An explicit actuator accelerates and then drives a
  `6000 kg` planar plate at `0.05 m/s`; CFD supplies the independently exposed
  resisting pressure. The connected fluid is routed around a resolved opening,
  and only accepted geometry/transport/structure states reach the trace. At an
  exact cell crossing the worker verifies old/new volume continuity, remaps the
  constraint within a written velocity budget, and commits the next epoch with
  the complete frame. Pressure reaches XPBD through the moving face-resolved
  bridge, not the older uniform-only subset, and only after the fluid-side
  planar cut-surface reaction is accepted. Grid and unwrapped physical planes,
  periodic-image closure, plus correspondence residuals remain explicit in
  every frame. The projection reaction is treated as a macro-step average and
  sampled at both endpoint velocities; structure, fluid, actuator, and total
  momentum/kinetic-energy ledgers must close before commit. Its composite
  checkpoint restores the full Structure state, accepted fluid fields and
  topology epoch, partial-cell offset, and committed ledgers; continuation must
  replay bit-for-bit in the same or an equivalent rebuilt worker.
  `open_piston_checkpoint_persistence.{h,cpp}` stores that complete composite
  in a deterministic, checksummed, byte-bounded envelope, validates both
  nested solver payloads through a rebuilt worker, then reconstructs the
  canonical epoch, stable identities, geometry, acceptance state, and
  force/moment/power/conservation relationships for every diagnostic ledger.
  It preserves the caller's prior output on any failed decode. The worker CLI
  restores it before opening the trace and uses the shared atomic output path,
  absolute autosave cadence, final-save deduplication, and additional-step
  resume semantics.
- `src/fsi/periodic_flow_case.{h,cpp}` and
  `src/fsi/periodic_flow_checkpoint.cpp` are selected with
  `--case periodic-flow`. It advances a Galilean-shifted Taylor-Green field on
  the periodic grid, writes only accepted subcycled states, and exposes cell
  pressure/velocity snapshots through the standalone viewer. It is a smooth
  CFD verification harness, not whole-wing aerodynamics. Advance remains
  transactional through frame validation, and its private-payload checkpoint
  restores the initial or any committed state into the same or an equivalent
  rebuilt worker with bit-identical continuation. The persistent codec is
  versioned, byte-bounded, checksummed, deterministic, transactional on decode,
  and validates complete nested diagnostics before restore. The CLI writes it
  by atomic replacement. Periodic autosave occurs only after accepted frame
  publication, uses absolute accepted-step multiples so cadence survives
  restart, skips a duplicate final write, and reports its write count;
  checkpoint flags are also supported by open-piston, strong-piston,
  moving-porous-flow, and porous-sheet but rejected for the structural and weak
  sealed-piston cases; interval mode requires an output path.
- `src/fsi/worker_control_protocol.{h,cpp}` owns the first
  transport-independent safe-point message contract. Command and response
  envelopes have distinct magic, one protocol version, a bounded payload
  length, and a checksum. Decoding is transactional; every response carries
  accepted step and time, while `Advanced` reports produced frames and `Error`
  carries one bounded failure code/message. It does not select stdin, pipes,
  sockets, or a scheduler transport.
- `src/fsi/worker_control_stream.{h,cpp}` owns incremental command/response
  framing over standard C++ streams. Physical EOF is clean only between
  envelopes; the worker CLI still requires an explicit stop. The adapter reads
  the bounded payload length before allocation, delegates checksum/type checks
  to the protocol codec, preserves caller output on every failure, and flushes
  each successful write.
- `src/fsi/worker_control_session.{h,cpp}` executes validated protocol commands
  synchronously on a numerical owner thread. Advance commits and sends each
  immutable accepted frame individually; a later sink failure leaves that
  accepted safe point visible in the absolute error response. Checkpoint calls
  an injected action without changing solver state. Stop rejects later
  advance/checkpoint commands and repeated stop remains idempotent.
  `periodic_flow_control.{h,cpp}` and `open_piston_control.{h,cpp}` bind the
  corresponding case, absolute step/time queries, and typed checkpoint sink.
  `strong_piston_control.{h,cpp}` applies the same boundary to whole accepted
  strong-coupling macro-steps; active iterations and rejected retries remain
  invisible.
- `src/fsi/fluid/grid.{h,cpp}` owns the uniform periodic Cartesian grid,
  cell-centred scalar fields, unique periodic MAC face velocities, and the
  paired finite-volume divergence/gradient operators.
- `src/fsi/fluid/projection.{h,cpp}` owns the deterministic zero-mean periodic
  pressure solve. It commits pressure and velocity together only on
  convergence; failure leaves both inputs bit-identical for later step retry.
- `src/fsi/fluid/diffusion.{h,cpp}` owns the first laminar velocity-evolution
  operator: forward Euler viscosity on each translated periodic MAC component
  with the centred seven-point Laplacian. It rejects an excessive diffusion
  number without mutation, preserves component momentum, cannot increase
  kinetic energy inside its declared stability interval, and is the exact
  stage oracle for a two-stage SSPRK2 integrator. SSPRK2 composes two full Euler
  candidates and the old field through a convex average, preserves the same
  stability/conservation contract, and is second order in time for the
  discrete viscous eigenproblem. Coupled advection/viscosity splitting is not
  yet the intended second-order production time integrator.
- `src/fsi/fluid/advection.{h,cpp}` owns the bounded donor-cell and limited MC
  transport operators for periodic MAC components. The prescribed-uniform oracle is a
  convex combination when `sum(abs(U_i)*dt/h_i) <= 1` and commutes with the
  discrete divergence. The variable-flow path averages a divergence-free MAC
  advector onto each translated component control volume, uses one shared
  upwind flux per periodic face, and supports safe self-advection aliasing.
  Accepted donor steps preserve component momentum, bounds, and non-increasing
  kinetic energy under the local outgoing-CFL limit. MC reconstructs
  monotonized slopes on the same shared faces. Its individual Euler stages may
  bypass only the energy check inside `advectVelocityByMacFlowSspRk2`; the
  convex aggregate rechecks original component bounds, momentum, and energy.
  A discontinuous pulse stays bounded, and smooth periodic uniform transport
  shows near-second-order L1 refinement. Nonlinear spatial behavior is verified
  through the projected operator below, not this fixed-advector enclosure.
- `src/fsi/fluid/projected_advection.{h,cpp}` owns pressure-projected nonlinear
  SSPRK2 transport. Stage one is self-advected and projected before becoming
  stage two's advector; the old field and twice-advanced prediction are then
  convexly averaged and projected again. Velocity and pressure commit together
  only after both selected transport stages, both projections, and the aggregate
  momentum/energy ledger pass. The fixed-grid vortical regression observes
  second-order temporal refinement. Donor or limited MC reconstruction is
  selected explicitly; intermediate MC Euler energy is accepted only here and
  the projected aggregate must still be non-increasing. The analytic
  Galilean-translated Taylor-Green regression suppresses temporal contamination
  with `dt` proportional to `h^2`; donor L1 ratios stay in `[1.7, 2.2]` and MC
  ratios in `[3.0, 5.0]` across 16/32/64 grids.
- `src/fsi/fluid/evolution.{h,cpp}` composes selectable uniform or nonlinear
  self-advection, Euler or SSPRK2 viscosity, and periodic projection into the
  first complete fluid macro-step. Its second path composes half-step SSPRK2
  viscosity, projected nonlinear SSPRK2 transport, and a symmetric viscous half
  step into the first complete second-order temporal periodic flow integrator;
  it exposes donor-cell or limited MC reconstruction explicitly.
  Its bounded subcycling wrapper uses equal steps, pre-sizes from the explicit
  half-viscosity limit, and restarts the full private interval with a finer
  schedule only for reported advection/diffusion stability or limited-MC bound
  rejection. The final schedule is retained exactly in diagnostics. Projection
  and conservation failures never trigger retry.
  Stage diagnostics remain intact, the aggregate momentum/energy/divergence
  ledger is checked independently, and failure at advection, diffusion,
  projection, either Strang sub-integrator, or final conservation commits
  neither velocity nor pressure. Its analytic viscous translating-Taylor-Green
  L1 regression uses 12/24/48 grids with `dt <= 0.12 h^2` and requires both
  successive error ratios in `[3.0, 5.0]`; this is a smooth periodic full-flow
  oracle, not a cut-cell or moving-interface accuracy claim.
  The first macro-step also has a porous-projection overload. Flow settings
  remain sole owner of density, time step, and linear tolerances; the porous
  input supplies only nonlinear controls. Its aggregate momentum residual
  subtracts diagnosed jump impulse, and its energy ceiling includes diagnosed
  jump work. Midpoint driven acceleration and unforced porous decay close those
  ledgers; nonlinear failure rolls back every prior stage and empty topology
  delegates bit-for-bit to the original path. A separate second-order
  fixed-grid wrapper symmetrically composes midpoint porous half-step, complete
  bulk Strang/SSPRK2 step, and midpoint porous half-step. It sums both interface
  impulse/work ledgers, retains porous dissipation, rolls back all operators on
  failure, and delegates empty topology to the original bulk path without
  changing its fields or nested diagnostics. Its stage-resolved planar overload
  accepts distinct complete sheet definitions at the two porous midpoints,
  preserves identity/regions/resistance, admits only a retained or adjacent
  topology epoch, and requires their displacement to equal trapezoidal
  integration of sampled normal velocity over `dt/2` within explicit
  tolerances. It retains both unwrapped epochs and that kinematic residual in
  diagnostics, and crosses a `3 -> 0`, image `0 -> 1` wrap transactionally
  inside one macro-step. This is moving planar porous-source sampling, not a
  cut-cell remap.
- `src/fsi/fluid/checkpoint.{h,cpp}` owns the versioned in-memory checkpoint for
  one accepted moving-interface fluid epoch. Its immutable payload includes
  pressure, velocity, exact interface kinematics, and projection diagnostics;
  public grid metadata and a stable topology fingerprint are rebound before a
  restore candidate is returned. `checkpoint_persistence.cpp` round-trips that
  complete state through a deterministic, checksummed, byte-bounded
  little-endian envelope. Decode rebuilds and validates the grid/interface and
  calls the accepted-state checkpoint validator before committing its output.
- `src/fsi/fluid/interface_jump.{h,cpp}` owns canonical grid-face crossings,
  signed two-region pressure jumps, the jump-corrected gradient, and its paired
  Poisson source. Crossings on the same face-normal cell segment retain stable
  surface identity, have distinct open-interval positions, and must form an
  ordered continuous region chain. The dense stencil uses their deterministic
  signed sum; this static pressure-jump subset is not moving folded topology.
- `src/fsi/fluid/planar_face_topology.{h,cpp}` is the generic pure
  axis-aligned moving-plane selector. `porous_topology.*` aliases and delegates
  to it so the established porous contract and arithmetic remain unchanged.
- `src/fsi/fluid/planar_pressure_jump.{h,cpp}` canonicalizes a closed periodic
  chain of stable planar pressure layers, expands every layer over all
  transverse X/Y/Z tiles, and transactionally translates each layer through at
  most one adjacent or wrapped face epoch. Same-face layers remain distinct in
  the sparse field, but the present dense stencil still aggregates their signed
  jumps. Its static regional profile exactly partitions the unwrapped period,
  requires a closed pressure potential, and volume-gauges each distinct region;
  it owns no regional velocity degree of freedom and does not enter production
  projection arithmetic.
- `src/fsi/fluid/planar_region_sweep.{h,cpp}` binds two static profiles into a
  bounded two-epoch geometric-conservation ledger. Stable layer identity,
  region sides, jumps, order, and one-segment topology motion are immutable;
  per-interval boundary sweep closes against volume change before per-region
  and global aggregation. It owns no Eulerian regional flux and is not a
  leakage solver.
- `src/fsi/fluid/planar_region_flux.{h,cpp}` revalidates that sweep and fits
  one least-squares uniform axial fluid velocity per interval. It reports both
  one-sided outward-relative flows, impermeability and continuity separately,
  and the minimum unavoidable slip for incompatible boundary motion. Its
  nonzero semantic fingerprint and deep validator bind primitive inputs,
  derived intervals, stable-ID chain order, sorted region/global ledgers,
  tolerance flags, and storage accounting. It is an offline screen and never
  advances or applies fluid state.
- `src/fsi/fluid/planar_region_opening_flow.{h,cpp}` consumes the same sweep
  plus explicit oriented analytic opening links. It requires independent
  connected-component source balance and solves a bounded gauge-fixed
  area-weighted graph Laplacian for the minimum-norm compatible flow. Sealed
  breathing remains infeasible; the output is fingerprinted, source-bound,
  and never becomes a pressure, velocity, or scene-aperture state.
- `src/fsi/fluid/planar_region_opening_power.{h,cpp}` binds only feasible
  opening allocations to the previous/current static regional pressures. It
  separates local uphill pressure-power deficits from net graph demand and
  closes midpoint opening power against regional `p*dV/dt`. It does not invent
  an aerodynamic energy source or opening constitutive law.
- `src/fsi/fluid/planar_region_fragment.{h,cpp}` splits the current planar
  pressure profile at Cartesian faces and expands every axial segment over the
  transverse tiles. It retains layer/grid boundary provenance, periodic image,
  stable fragment identity, physical volume/centroid, regional pressure, and
  exact cell/region/domain closure. It owns no connectivity or velocity state.
- `src/fsi/fluid/planar_region_fragment_topology.{h,cpp}` pairs all six faces
  of those fragments. Same-region Cartesian pairs receive a geometric weight;
  pressure-layer pairs retain two-sided wall geometry and signed static jump
  with zero conductance. The graph assigns deterministic connected components
  but owns no velocity, momentum, opening link, or pressure solve.
- `src/fsi/fluid/planar_region_fragment_pressure_operator.{h,cpp}` assembles
  an ungauged symmetric integrated graph Laplacian from only the same-region
  fragment links. Layer walls create no entries. It retains one deterministic
  gauge per component but owns no RHS, solve, velocity correction, or
  production pressure state.
- `src/fsi/fluid/planar_region_fragment_pressure_solve.{h,cpp}` solves only a
  component-compatible correction on that operator. It removes admitted RHS
  roundoff, commits roundoff-zero volume-mean corrections, keeps the static
  regional potential separate, and rolls back the warm start on failure. It
  owns no divergence source or production pressure update.
- `src/fsi/fluid/planar_region_fragment_pressure_projection.{h,cpp}` supplies
  per-link static and topology-stable moving projection over that solve.
  Same-region links own oriented area flow and receive pressure correction;
  pressure-layer walls require exact zero relative flow. The moving overload
  composes the bound fragment `dV/dt`. Velocity and pressure commit together
  only after recomputed continuity closes. Rebases, face momentum mass,
  opening conductance, and production integration remain outside this slice.
- `src/fsi/fluid/planar_region_fragment_volume_rate.{h,cpp}` reconstructs
  topology-stable previous volume and constant geometry `dV/dt` for every
  current regional fragment from its layer-boundary displacements. It closes
  cell, region, component, and global ledgers, rejects Cartesian rebases, and
  is consumed only by the opt-in moving projection, never production.
- `src/fsi/fluid/planar_region_fragment_velocity_metric.{h,cpp}` assigns the
  diagonal dual-volume geometry for regional normal velocities. Same-region
  grid links own one shared degree of freedom; each pressure-layer wall owns
  two independent one-sided traces. Fragment, component, axis, and domain
  volumes close without carrying velocity, density, energy, or production
  state.
- `src/fsi/fluid/planar_region_fragment_velocity_state.{h,cpp}` applies one
  positive density and one finite normal velocity to every regional metric
  degree. It retains per-DOF momentum/energy and closes fragment, component,
  axis, and global mass/momentum/kinetic-energy ledgers without averaging the
  independent pressure-wall traces. It owns no wall prescription, pressure-
  projection certificate, transport, pressure work, or production state.
- `src/fsi/fluid/planar_region_fragment_projection_energy.{h,cpp}` certifies a
  static or topology-stable moving correction against before/after regional
  velocity states. It closes gradient velocity change, pressure impulse versus
  diagonal momentum, midpoint pressure work versus kinetic-energy change,
  component gauges, and corrected continuity. Moving audits bind volume rates,
  material-wall trace velocity, and the affine geometry-work identity. It
  excludes authored jump work, topology rebase, transport, and production.
- `src/fsi/fluid/porous_interface.{h,cpp}` applies a calibrated normal
  Darcy-Forchheimer resistance to resolved MAC velocity relative to an authored
  sheet. It emits canonical signed sharp jumps and retains per-tile area,
  volume flux, and nonnegative dissipation. Its bounded Picard solver
  repeatedly samples that law and projects the original predicted MAC field
  until both normal-velocity and pressure-jump residuals close. Constitutive
  sampling is selectable at the endpoint or the midpoint between the original
  prediction and candidate. Optional
  prescribed jumps remain separately owned in every trial; empty topology
  delegates to the original projection. It supports heterogeneous fixed-grid
  tiles and prescribed sheet-normal motion. Midpoint mode matches the scalar
  nonlinear plug oracle and closes its pressure-work/dissipation/kinetic-energy
  identity. Accepted diagnostics sum all porous and prescribed oriented jumps
  into explicit fluid force/impulse and power/work ledgers while keeping porous
  dissipation separate. It is composed into both the first periodic macro-step
  and the symmetric porous-half/bulk-Strang/porous-half second-order fixed-grid
  integrator. Its nonlinear iteration can also select the disconnected
  moving-interface/jump projector as the inner solve, retaining the final
  moving reaction diagnostic while porous/prescribed jump ledgers remain
  separate. Impermeable and porous faces cannot overlap. This supports
  fixed-grid moving boundaries around porous flow, not moving porous cut-cell
  topology.
- `src/fsi/fluid/porous_topology.{h,cpp}` preserves the moving-sheet API over
  the generic planar selector. A topology epoch carries axis, wrapped face
  coordinate, and signed periodic image. Selection retains the current strict
  `(0,1)` crossing or advances exactly one face in either direction, including
  X/Y/Z wraps; exact MAC-plane placement, skipped segments, invalid
  axes/versions, and non-finite positions are rejected before caller state
  changes.
- `src/fsi/fluid/planar_porous_sheet.{h,cpp}` expands one validated authored
  topology epoch into exactly one deterministic crossing for every transverse
  MAC tile on X, Y, or Z. It owns the common stable-ID, two-sided-region,
  material, rigid-normal-kinematics, and strict physical-placement checks used
  by both static porous flow and the translating porous-sheet oracle.
- `src/fsi/fluid/porous_flow.{h,cpp}` owns the pressure-driven uniform-plug
  midpoint oracle. It solves the nonlinear Darcy-Forchheimer response exactly
  for each time step and independently closes pressure impulse and
  work/dissipation/kinetic-energy ledgers before committing its scalar state.
- `src/fsi/porous_flow_case.{h,cpp}` maps that accepted plug velocity to a
  periodic MAC field, a porous loss plane and zero-net-jump gauge-closure plane,
  verifies both internal pressure stages of a complete sharp-jump Strang step,
  and publishes an owning trace frame. The physical driving rise stays explicit
  in the plug ledger. It is not general nonuniform or moving porous coupling.
- `src/viewer/pressure_jump_frame.{h,cpp}` maps accepted sharp-jump state to an
  immutable owning frame. Cell-centre samples retain pressure and velocity
  diagnostics; every authored crossing retains its fraction, signed jump,
  oriented region pair, normal, and separate quad, including periodic wrapping
  and multiple layers on one face.
- `src/fsi/pressure_jump_case.{h,cpp}` repeatedly projects a static periodic
  split slab from a fresh zero state and publishes only the accepted frame. Its
  48 crossings resolve two ordered transitions at each slab boundary and
  preserve the analytic `-125/+125 Pa` field without spurious flow.
- `src/fsi/fluid/moving_interface.{h,cpp}` owns physically grid-bound MAC-face
  velocity constraints and their stable connected fluid regions. Its
  disconnected projection preserves one prior pressure mean per region,
  rejects nonzero regional volume rate before mutation, and reports canonical
  MAC-tile bounds, traction, force, velocity, and power as well as aggregate
  pressure impulse/work per surface. This pressure sampling is exact for the
  piecewise-constant slab test, not general surface reconstruction. Distinct
  side IDs require a separating topology; equal IDs describe a nonseparating
  surface connected around a resolved grid path.
  Adjacent-cell pressure traction remains separately diagnosed. The complete
  coupled load additionally includes the direct reaction required to replace a
  predicted constrained-face velocity: on the uniform MAC grid,
  `F_reaction = A*(p_minus-p_plus)*n - rho*V_face*(U-u*)/dt`. The direct term is
  exactly zero for an already compatible prescribed velocity.
  The same disconnected solve may retain an immutable sharp pressure-jump
  field on unconstrained faces. Its paired source is gauged per connected fluid
  region, the jump-corrected gradient never modifies constrained MAC values,
  and a face claimed by both boundary types is rejected before mutation. Empty
  jump fields preserve the moving-only arithmetic bit-for-bit. This fixed-grid
  coexistence is not a moving porous or folded-interface coupling model.
- `src/fsi/fluid/moving_control_volume.{h,cpp}` binds one complete
  nonseparating planar surface and one complete open MAC plane. Within the
  first partial cell it independently closes geometric volume change, surface
  sweep, projected opening transport, velocity, and pressure-power ledgers. An
  overload accepts only a fully accepted moving-porous wrapper whose nested
  base projections remain identical, preventing a converged inner solve from
  escaping an unaccepted nonlinear iteration. The porous-opening canonical
  matches its resolved tile flux to both the swept volume and GCL transport. An
  explicit candidate rebase converts the completed partial cell to a full
  reference layer on the next positive-axis MAC plane while preserving stable
  IDs and chamber volume; it rejects skipped planes, changed regions, broken
  ledgers, and collision with the opening. General cut-cell pressure metrics,
  nonplanar topology events, folded moving surfaces, and multiple crossings in
  moving/cut-cell topology remain future work.
- `src/fsi/fluid/planar_cut_surface.{h,cpp}` owns the bounded physical-plane
  pressure-reaction geometry for that control volume. It retains each canonical
  MAC tile as the Eulerian source, translates its application rectangle to the
  unwrapped physical plane only when grid-plane plus partial-cell offset is a
  matching periodic image, and independently closes aggregate area, force,
  moment, and power. It transfers the projection constraint reaction; it does
  not interpolate a new cell pressure or claim general cut-cell reconstruction.
  Its moving-porous overload enforces the same complete outer acceptance and
  nested-projection identity before exposing the physical reaction.
  Its temporal resampler preserves that accepted reaction force while changing
  only congruent physical geometry, rigid normal velocity, moment, and power so
  a macro-step-average projection load can be integrated at endpoint kinematics.
- `src/fsi/transfer.{h,cpp}` owns canonical stable-ID coupling topology and
  immutable transfer results. It integrates either current triangle
  area/centroid or explicit area/barycentric quadrature patches, distributes
  traction resultants barycentrically, checks independent surface/nodal force,
  moment, and power ledgers, and binds application to the exact Structure
  definition/surface fingerprints.
- `src/fsi/coupling.{h,cpp}` owns the versioned time-integrated interface
  exchange. Macro-step-local samples must begin at zero, increase strictly,
  retain one moment reference and one topology binding, and are integrated by
  the trapezoidal rule into independent surface/nodal impulse, angular impulse,
  and work ledgers. Acceptance requires the exact Structure step duration,
  applies equivalent average nodal loads across its internal substeps, and
  restores the checkpoint from before load application on any failure. It also
  owns the definition-bound vector Aitken delta-squared state: first-iteration
  fixed relaxation, bounded dynamic factors, deterministic degenerate fallback,
  and transactional exact iteration checkpoint/restore. Its stateless
  convergence gate requires displacement, velocity, and traction to each pass
  absolute and floor-stabilized relative budgets, after a minimum iteration
  count, and distinguishes convergence from maximum-iteration exhaustion. Its
  topology-bound residual reducer takes baseline/previous/current kinematics
  and consecutive immutable traction transfers, uses maximum nodal updates,
  references motion to macro-step changes, and rejects foreign bindings before
  returning origin- and Galilean-invariant motion norms. Its macro-step
  iteration owner composes the current relaxed vector, Aitken state, and
  convergence decision into an exactly replayable checkpoint. Failed advances
  and restores are transactional; converged/exhausted states are terminal.
  Fluid and Structure rollback payloads remain in their solver owners.
- `src/fsi/coupled_state.{h,cpp}` composes those three rollback domains without
  flattening their ownership. It owns a real `Structure`, accepted
  moving-interface fluid state, and iteration owner, prevalidates complete
  replacement instances, then commits all three with no-throw moves. Invalid
  composite identity, Structure, fluid topology, or iteration state mutates
  none of them. The same boundary owns a definition-bound macro-step retry
  policy: iteration exhaustion creates a pending reduced step, the caller must
  restore the composite baseline before activating it, and minimum-step or
  retry-budget exhaustion is terminal and exactly checkpointable. Its
  one-macro-step coordinator accepts only a fresh iteration, privately retains
  that baseline, and composes restoration with retry activation so discarded
  Structure, fluid, or Aitken state cannot leak into the reduced attempt. A
  distinct transactional solver checkpoint restores only Structure and fluid
  between iterations, preserving the new relaxed interface and Aitken history.
  `runStrongCouplingMacroStep` drives this state machine around a typed callback
  that may mutate only the physical solver owners. It validates every returned
  epoch before Aitken advance, automatically retries exhausted steps, retains
  converged output, and restores the baseline before propagating callback or
  validation failure and before returning terminal retry failure.
- `src/fsi/fluid_structure_bridge.{h,cpp}` owns the first deliberately narrow
  uniform grid-to-structure bridge. One stable fluid surface can drive one
  structural coupling surface only when adjacent-cell pressure traction is
  uniform and independent area, force, and power ledgers close.
- `src/fsi/face_resolved_bridge.{h,cpp}` owns two explicit planar
  correspondence modes. Both clip nonoverlapping, fully covering axis-aligned
  reference MAC tiles against consistently oriented structural triangles once,
  then map each current tile traction through stable overlap-area/centroid
  barycentric quadrature. `FixedMaterial` preserves the original exact face
  binding. `RigidNormalTranslation` accepts equal-sided nonseparating surfaces,
  normal grid-plane rebasing, and an explicitly unwrapped physical plane only
  while all transverse tile geometry and structural node coordinates remain
  fixed and all fluid/structural normal velocities agree. Per-face and
  aggregate area, force, moment, and power ledgers must close. Curved,
  transverse-deforming, or general Eulerian remap and strong-coupling
  convergence remain future work.
- New numerical targets must not link the inherited `playground_*` libraries.
  Cross the scene/structure/viewer boundaries through explicit adapters and
  versioned data only.

### Vendored and generated content

- `third_party/libf2c`: vendored runtime. Avoid edits unless fixing a proven
  compatibility/runtime issue; preserve `Notice`.
- `third_party/xflr5`: vendored v6.62 application and XFoil. Read
  `third_party/xflr5/PROVENANCE.md` before changing it and record any new local
  embedding patch there. The upstream app bootstrap is intentionally excluded.
- `src/softwing`: originated from SoftWingLab but has intentionally diverged.
  Read `src/softwing/README.md`; future syncs are merges, never blind copies.
- `resources/manual`: generated offline manual and qrc. Do not hand-edit its
  HTML/images; run `python tools/extract_manual.py`. Delete the gitignored
  `tools/manual.en.html` first only when an upstream refresh is intended.
- `resources/presets`: shipped, self-contained designs plus
  `presets.json`. Update its README when adapting upstream files. All presets
  must pass Studio validation, not just the more tolerant engine parser.
- `cmake/xflr5_sources.cmake`, `installer/version.iss`, build trees, deployed
  DLLs, installer output, and most calculation output are generated. The first
  has a generator script; `installer/version.iss` is regenerated at configure
  time and is gitignored.

## Verification matrix

There are 164 configured tests on Windows. The Fortran-reference test is
Windows-only; local `gui_smoke` and `studio_model_smoke` exercise display/model
paths that release CI deliberately excludes from its offscreen test command.
The light added-mass path is covered both inside `simwing_piston_case` and by
`simwing_fsi_strong_piston_headless`; preserve multiple real coupling solves,
weak-response separation, accepted interface/Structure velocity closure,
deterministic persistent Structure/fluid state, immutable-frame provenance and
residual fields, the `10.8 kg` discrete pressure added-mass slope and matching
trapezoidal fixed-point speed, accepted checkpoint next-step/frame replay,
transactional public-version/interface/Structure/fluid-topology rejection, and
deterministic persistent decode/re-encode and continuation, bounded outer and
nested sizes, magic/protocol/reserved/checksum/truncation/trailing/identity
rejection without output mutation, and accepted-only CLI trace publication.
For second-order fixed-grid porous evolution, `simwing_fluid_evolution` must
also preserve exact porous-half/bulk-Strang/porous-half composition, summed
interface ledgers, later-stage rollback, bit-exact empty-topology bulk
delegation, and the driven-flow temporal-refinement ratios in `[3.8, 4.2]`.
For porous flow inside fixed-topology moving regions,
`simwing_fluid_porous_interface` must preserve endpoint/midpoint constitutive
closure, exact constrained velocity, separate moving-reaction and porous-jump
ledgers, bit-exact empty delegation, overlap rejection, and nonlinear rollback.
`simwing_fluid_control_volume` must additionally keep porous opening-tile flux,
moving-volume change, projected opening transport, and physical cut reaction
consistent, while rejecting incomplete or internally inconsistent wrappers.
Mimetic pressure-state or consecutive-epoch warm-remap changes require both
`simwing_scene_fluid_mimetic_trace_system` and
`simwing_scene_fluid_pressure_volume_rate`; accepted-state surface sampling
also requires `simwing_model_scene_real_export`. Preserve source-bound accepted
capture, stable control/trace retention, same-region appearance donors,
new-trace endpoint initialization, trace retirement, exact current gauges,
bounded fingerprinted storage, direct atomic-solve consumption, exact
cell/region side sampling, shared-component gauge safety, and conservative
force/moment transfer. Atomic epoch changes must additionally preserve
bootstrap/consecutive identity, nested limits, and diagnostics-only rollback.
Persistent-state changes must preserve `SWMP` deterministic round trips,
transactional corruption/limit rejection, and trusted-topology rebinding.

| Change area | Minimum relevant checks |
|---|---|
| engine paths/migration/legacy I/O | `engine_help`, `engine_rejects_missing_input`, `engine_tolerates_blank_lines`, `engine_matches_fortran_reference` |
| translated numerical/drawing core | full suite; the Windows Fortran reference is mandatory, plus `preset_gnua7_calculates` |
| design parsing/history/editor orchestration | `studio_model_smoke`, `presets_validate_in_studio`, `gui_headless_help` |
| Section 1/splines/airfoils/grid helpers | `section1_curves`, `spline_fit`, `airfoil_file`, `section_grid` as applicable |
| geometry preprocessor | all four `preprocessor_matches_fortran_*` tests |
| NURBS/XCAF/viewport/engine callback changes | `engine_matches_fortran_reference`, `studio_model_smoke`, and `gui_smoke`; inspect engine model statistics/warnings |
| flat capture/nesting/PDF/DXF | `flatparts_export`; run `nesting-bench` when placement quality/performance can move |
| Playground body/pressure/cells/material | `playground_pressure_solve`, `playground_cells`, `playground_metrics`, `playground_material`, `softwing_material`, and the deterministic bench guards below |
| Playground contact | `playground_contact`, `playground_contact_integration`, plus a relevant bench/GUI scenario |
| SimWing mimetic local cell/wall condensation | `simwing_fluid_mimetic_local_cell`, `simwing_model_scene_real_export`; preserve complete unwrapped half-face geometry, exact area-vector and `N^T R = volume I` closure, fingerprinted SPD flux inner products in the exact seven-doubles-per-half-face diagonal-plus-two-rank-three factorization, matrix-free application, linear-field exactness, roundoff-only integrated source conservation, translation invariance, skew-tetrahedron manufacture, Cartesian `area / centre-distance` equivalence, bounded allocation, and transactional invalid/corrupt input rejection; preserve exact diagonal-plus-seven-low-rank material-wall Schur action/diagonal/RHS/reconstruction against independent dense oracles, fused `D_a + U_a (K - K Q K) U_a^T` action, active constant nullspace, no-wall identity, all-wall rejection, normal equilibrated `7 x 7` Woodbury inversion, the conditional at-most-eight-wall direct principal-block fallback with no persistent dense storage, the four-face/two-wall sliver oracle, linear bounded storage, fingerprinted corruption rejection, and all 138 coarse real controls reducing 148,652 walls to 42,927 shared traces with 3,986,602 bytes and positive assembled condensed diagonals |
| SimWing mimetic scene shells | `simwing_scene_fluid_mimetic_control_cell`, `simwing_model_scene_real_export`; preserve control-major stable half-face identity, periodic-image unwrapping, exact paired Cartesian/material/opening source ownership, per-control area-vector and `N^T R` closure, analytic nested/face-opening/rejected-embedded completeness, compact local-kernel acceptance for every ready shell including all coarse real-wing controls, conservative unresolved-face incident marking, authoritative oriented-surface classification of multiply supported untouched faces, all 138 coarse and all 358 manually refined real-wing controls closed without dominant-volume or closure fitting, bounded storage, and fingerprinted corruption rejection without changing production pressure arithmetic |
| SimWing mimetic trace system/condensation/solve | `simwing_scene_fluid_mimetic_trace_system`, `simwing_model_scene_real_export`; preserve exact two-incidence Cartesian/opening trace pairing, unique one-incidence zero-flux material-wall traces, deterministic component gauges, rejection of disconnected within-component trace topology, compact local-kernel binding, positive Jacobi diagonals, matrix-free symmetric positive-semidefinite action, one component-constant null mode, compatible source-derived right-hand sides, bounded storage and nested fingerprint/corruption rejection; preserve exact immutable full-to-shared topology, shared component gauges, summed condensed diagonals, symmetric positive-semidefinite/null matrix-free action, full-RHS condensation, unique wall reconstruction, manufactured full-row closure, bounded nested storage/corruption rejection, and the coarse 148,652-wall-to-42,927-shared/3,986,602-byte audit; preserve exact fixed-MAC and transported wall-adjusted oriented trace-flow sampling, accepted cap-sweep subtraction, coverage of rejected embedded graph apertures, source-density/GCL-duration provenance, shared full/condensed gauge-fixed Jacobi-PCG manufacture, declared component-roundoff correction only, exact gauge publication, freshly recomputed residual convergence, reconstructed original-row closure, local/source/trace closure, validate-once trusted inner application, transactional rollback on incompatibility/non-finite arithmetic/exhaustion, finite residual-reducing truncated audits of the real 191,579-trace full and 42,927-trace condensed systems, condensed convergence to `1e-5` relative RMS within 300 iterations with reconstructed maximum full-row residual below `2e-4`, bounded fingerprinted `-(rho/dt)*(dV/dt + net_outward)` physical source conversion with compensated component sums, atomic source-bound publication only after reconstructed full-space closure including the 307-iteration balanced real-source audit, and independently reconstructed bounded accepted-state capture with exact zero shared gauges and control/full/condensed/source provenance, without changing production pressure arithmetic |
| SimWing capped cell volumes/opening kinematics | `simwing_scene_fluid_opening_cap`, `simwing_scene_fluid_opening_quadrature`, `simwing_scene_fluid_opening_patch`, `simwing_scene_fluid_opening_face_crossing`, `simwing_scene_fluid_capped_face_partition`, `simwing_scene_fluid_opening_flux`, `simwing_scene_fluid_cell_volume`, `simwing_scene_fluid_region_continuity`, `simwing_scene_fluid_region_connectivity`, `simwing_scene_fluid_pressure_control_volume`, `simwing_scene_fluid_pressure_face_link`, `simwing_scene_fluid_pressure_operator`, `simwing_scene_fluid_pressure_epoch`, `simwing_scene_fluid_pressure_solve`, `simwing_scene_fluid_pressure_projection`, `simwing_scene_fluid_pressure_volume_rate`, `simwing_scene_fluid_pressure_sampling` plus the authoritative scene/grid geometry checks when their inputs change; preserve accepted-state/epoch binding, atomic full pressure-geometry/operator composition with nested validation and aggregate storage limits, scene-v2.2 optional oriented boundary-vertex cap disks with bounded topology/wire storage and cap-free binary-v2 read compatibility, boundary-edge-complete topology-only cap winding, exact convex-fan compatibility, bounded deterministic reference-geometry ear clipping for planar concave loops, authored nonplanar facet identity and individual normals across accepted motion, rejection of folded/self-intersecting/degenerate or unauthored nonplanar caps and opening-only cap vertices, topology-stable one-point triangle kinematics, exact piecewise-linear normal surface sweep and rigid material-plus-cap GCL closure, exact polygon partition with unique cell/non-periodic-face ownership across accepted motion, exact paired transverse cap crossings from both adjacent clipped-cell owners, retained face-area ownership and explicit grid-edge ambiguity, bounded material-plus-cap planar arrangements with disconnected signed-cycle accounting, one resolved-or-unresolved record per touched face, exact analytic capped pressure-link area handoff, complete MAC-field binding, exact partial-face flow, deterministic off-face staggered quadrature, negative-to-positive sign, moving-cap subtraction and co-moving zero-relative-flow invariance, equal-and-opposite per-region intake/crossport sources with exact global cancellation, consecutive-epoch identity and endpoint-trapezoidal `delta volume + outward relative flow` closure, explicit local mismatch reporting even when the global ledger cancels, deterministic intake/crossport connected components with one stable gauge each, exact component source cancellation, rejection of changing sealed pressure components, cell-major positive-volume pressure unknowns with stable fixed-grid cell/region IDs, exact cell/region/domain volume and per-cell first-moment closure, cell-local analytic centroids preserved into pressure controls, component membership and deterministic gauge-volume ownership, exact full-face and nested mixed-region pressure-link area closure including oriented closed-surface classification of multiply supported untouched faces, same-region/component endpoint binding, oriented face-aligned cross-region aperture links with exact same-region complementary area, oriented cell-owned embedded aperture links using projected centroid separation, typed unresolved embedded-patch provenance/reason/signed-plane-distance/count/area plus bounded periodic-aware same-region one-ring support/sidedness without a fabricated link, no absolute-value repair or donor rerouting of a non-admissible stencil, no dominant-cell/interface smearing, symmetric positive-semidefinite integrated graph energy, exact constant nullspace, component row-sum conservation, one connected link graph and retained gauge per authored component, exact face-aligned and embedded aperture graph energy, deterministic manufactured multi-component and open-intake CG recovery, explicit component compatibility with roundoff-only correction, exact gauge normalization, explicitly recomputed final residual, fixed-epoch `-rho/dt` physical RHS assembly, exact MAC-versus-relative-aperture predicted flow ownership including off-face staggered samples, independent per-link pressure correction, corrected control-volume continuity, rejected-attempt non-publication, stable consecutive pressure-volume identity, exact per-control/component geometry `dV/dt`, one shared fingerprinted topology-transition mapping, explicit appearance marking with zero previous volume, conservative unique-neighbour disappearance retirement, local `dV/dt + net flow` closure, zero-velocity pressure-driven intake inflation, exact side-cell pressure sampling, shared-component gauge safety, gauge-shift-invariant traction, conservative pressure load delivery, independent-gauge rejection, transactional incompatibility/non-convergence rollback, operator rejection of unresolved Cartesian or embedded connectivity, periodic-boundary ambiguity rejection, no cap-to-Structure/traction path, deterministic signed-tetrahedron clipping and sparse ordering, exact analytic nested/open/large-cavity volumes including full interior cells, local/global closure, bounded work/publication, and explicit rejection of non-manifold/inconsistent winding |
| SimWing scene pressure feedback | `simwing_scene_fluid_pressure_coupling`, `simwing_coupling`; preserve exact macro-step baseline rewind before every nonlinear solve, stable or bounded one-ring appearance/unique-retirement consecutive pressure epochs, trapezoidal accepted-start/relaxed-end total fluid-load application, the existing conservative pressure/shear-to-Structure path, displacement/velocity/applied-versus-solved-load convergence, deterministic Aitken continuation, accepted pressure warm-state continuation, exact next-step replay, and full Structure/pressure/wall-owner rollback on exhaustion, unsupported topology, projection failure, or exception; accepted corrected-flow MAC collapse must restore oriented opening-cap sweep, close every Cartesian face's total volume flow, quantify mixed-subface velocity loss, validate/report embedded links without smearing them onto a face, and remain a continuation primitive rather than a claim of cut-cell advection, viscosity, or general topology rebase; the explicit inverse link-flow continuation must bind consecutive topology and velocity/opening identities, carry absolute subface deviations, current-area-recentre them, reapply current cap sweep, close the bulk face flow, reject corruption/limits, and remain disabled in the worker while the static-carry negative oracle proves that missing region-resolved transport drains pressure load; reconstructed region momentum must bind the accepted pressure/volume/link/opening/predictor chain, preserve exact Cartesian component averages, solve bounded explicit-normal equations for embedded-incident controls with exact MAC fallback in their nullspace, and preserve deterministic momentum/energy/fallback ledgers; the transported coupling overload must bind one accepted prior projection/face epoch, keep that fluid transport immutable until convergence, build one bounded fingerprinted transition shared by volume rates, momentum rebase, and pressure warm state, use exact stable-ID velocity retention plus same-region one-ring area-weighted donors for appeared controls, conservatively retire each disappeared source only through one unique previous same-region neighbour, reject cross-material/donorless/ambiguous ownership, remap independently for every nonlinear geometry iterate, perform bounded two-sided tangential material-wall exchange with equal-and-opposite Structure traction and closed momentum/work/dissipation ledgers, feed its adjusted link predictor into pressure projection, and retain a minimal accepted wall endpoint for exact restart |
| SimWing mimetic pressure shadow | `simwing_scene_pressure_cell_case`, `simwing_model_scene_real_export`, `simwing_fsi_pressure_cell_mimetic_audit_headless`, and `simwing_fsi_pressure_cell_mimetic_audit_checkpoint_resume`; require a complete fingerprinted endpoint over rebuilt control/full/condensed topology, predictor, physical source, nested solve diagnostics, accepted pressure state, and material samples; fixed pre-operator real-wing capture must cover all 138 controls and 42,927 shared traces without requiring the rejected graph operator; live audit must run only after graph convergence, bootstrap from fixed MAC flow, then use transported wall prediction plus consecutive warm remap, retain but never apply sampled loads, build one bounded fingerprinted comparison over every source row, pressure jump, and conservative nodal load with exact identity, source/pressure/force/moment/power deltas and zero-delta real-geometry oracle, preserve the observed roundoff-equivalent source vectors plus nearly pure `2.53`/`2.51` shadow pressure/load gain as evidence that the remaining roughly 60% disagreement is an operator-response transition blocker, use the separate bounded seven-mode offline inverse-response audit to prove that the accepted direction is about `2.54` gain, five bulk modes stay near unit gain, and the intake-only authored-region mode reproduces `2.562` fitted gain plus a `0.18 m` graph versus `0.0700821 m` shadow energy-equivalent conductance (`2.5684` ratio) rather than applying a global scale repair; preserve the graph-independent bounded mixed-hybrid terminal audit with uniform area-weighted balanced transfer, explicit face-aligned Cartesian versus embedded authored-opening trace identity, `0.0700820848335194 m` face-aligned compatibility closure and `0.0608388978079475 m` embedded fixture oracle; preserve the separately checksummed skew-intake `2^3`/`4^3`/`8^3` refinement spectrum, nested audit integrity and limits, graph/shadow ratios `28.024`/`33.672`/`8.830`, normalized shadow trend `0.102`/`0.358`/`0.503`, the fixed-`4^3` eight-phase yield of six accepted plus two typed incomplete-face rejections with one/two unresolved embedded patches, accepted graph/shadow normalized-conductance coefficients of variation `0.77175`/`0.34030`, the multi-resolution eight-phase topology yields `4/8`/`6/8`/`2/8` and conditional graph means `3.9025`/`6.2525`/`4.5187` versus shadow means `0.10168`/`0.26909`/`0.45753`, and explicit graph-censoring/non-authoritative interpretation; preserve the graph-free 40-sample area-weighted terminal matrix with `8/8` shadow yield at every level, normalized mean/CV `0.100660/0.01043`, `0.240930/0.39050`, `0.521248/0.38828`, `0.902570/0.20104`, and `1.091977/0.07436`, `32^3` range `0.956232`-`1.207995`, local-coordinate positive/complementary face moments, canonical triangle/grid-edge nodes, bounded at-most-eight-wall direct Schur fallback for a singular auxiliary core, separate closure/publication tolerance identity, cell-local cut-volume first moments, bounded cancellation-scale centroid repair, no closure fitting, and the explicit offline `1e-9` algebraic tolerance while production keeps `1e-10`; preserve the three-level immutable convergence screen with mean contraction `0.49671`, apparent order `1.00952`, extrapolation `1.27891`, relative fine gap `0.14616`, four phase direction reversals, four noncontracting phases, only one phase passing both, and strict `InsufficientEvidence` mask `0x300`, while a permissive zero-mask result remains a read-only trend candidate; preserve the `[256,-512,1024] m` common-translation oracle on the selected fine phase with exact topology and intake-area/normalized-response deltas below `1e-12 m^2`/`5e-11`; leave default frames and the graph portion of composite checkpoints byte-identical, rollback graph/Structure on endpoint/comparison-limit failure, persist only independently bounded compact `SWMP` rows plus settings identity, rebuild trusted topology before decode, reject cross-mode/corrupt/oversized restart transactionally, and reproduce the exact uninterrupted next endpoint and comparison after CLI resume |
| SimWing scene region momentum | `simwing_scene_fluid_pressure_projection`; preserve immutable region-momentum reconstruction identity and exact MAC fallback ownership, donor-cell vector-momentum transfer over corrected relative links, equal-and-opposite graph-viscous impulse, deterministic local outgoing-Courant/viscous-number subcycling, non-increasing post-forcing stage energy, global internal-momentum conservation, exact moving-volume GCL for uniform flow, explicit bound bulk-MAC delta impulse/work, bounded transactional non-publication, exact zero fixed points, topology-stable first-order current-epoch remap at retained region velocity, bounded shared-transition rebase with stable-ID retention, same-region one-ring area-weighted appearance donors, unique-neighbour disappearance retirement, exact source-volume/momentum closure before geometric correction, pressure warm-start mapping parity, cross-material/donorless/ambiguous rejection, current-volume momentum accounting, endpoint-averaged current-link prediction on each explicit normal with exact Cartesian fast-path arithmetic, exact authored-opening cap-sweep subtraction, explicit transport/rebase/wall projection provenance, two-sided tangential wall exchange at exact cell/region-owned material quadrature, normal-pressure ownership, bounded aggregate viscous subcycling, nonnegative wall dissipation, equal-and-opposite fluid/Structure impulse, exact zero-viscosity predictor preservation, accepted wall-traction capture/integrity, and rejection of general swept-volume-remap or resolved-boundary-layer claims |
| SimWing visible scene pressure cell | `simwing_scene_pressure_cell_case`, `simwing_fsi_pressure_cell_headless`, and the `simwing_fsi_pressure_cell_checkpoint_*` chain; preserve byte-deterministic accepted frames, three fixed intake vertices with no mechanical actuator, the prescribed periodic mean-wind pump, private-candidate symmetric SSPRK2 viscosity/projected-nonlinear-advection/viscosity bulk prediction with finite acceptance, bounded final divergence, nonnegative viscous loss, and no predictor mutation before scene acceptance, accepted cut-region momentum reconstruction and moving-volume transport after bootstrap, explicit-normal embedded-opening continuation without Cartesian smearing, bound bulk-MAC increment impulse/work, current-geometry material-wall exchange and adjusted current-link prediction in every strong iterate, topology-stable moving pressure epochs, final cut-region-aware pressure feedback, pressure-corrected bulk MAC continuation between accepted steps, reported maximum MAC speed/subface collapse spread plus region loss/GCL/momentum residual and wall loss/momentum residual, separate conservative pressure/wall/total-fluid load fields, exact derived-MAC reconstruction after restore without persisting transient bulk pressure, bounded deterministic initial/accepted `SWPCELL10` checkpoint round trips with complete accepted region momentum, accepted wall traction, transported-region/wall predictor provenance, and optional compact mimetic state, exact split-run next-frame and audited-endpoint replay, transactional foreign/cross-mode/corrupt/truncated/trailing/oversized momentum/wall/SWMP rejection, sustained bounded pressure-driven motion, per-triangle area-weighted pressure jump, global pump vector field, bulk-flow telemetry, strong-iteration telemetry, and explicit airflow-bootstrap diagnostic scope without a general immersed-boundary advection, resolved boundary-layer, or aerodynamic-truth claim |
| SimWing authoritative scene/grid geometry | `simwing_scene_fluid_surface`, `simwing_scene_fluid_surface_transfer`, `simwing_scene_surface_geometry`, `simwing_scene_surface_intersection`, `simwing_scene_surface_clipping`, `simwing_scene_surface_ownership`, `simwing_scene_surface_crossings`, `simwing_scene_surface_face_topology`, `simwing_scene_surface_face_graph`, `simwing_scene_surface_face_chains`, `simwing_scene_surface_face_loops`, `simwing_scene_surface_face_partition`, `simwing_scene_fluid_quadrature`, `simwing_scene_fluid_pressure_traction`, `simwing_scene_fluid_grid_epoch`; preserve complete accepted-state fingerprints, canonical bounded broad/narrow phases, exact barycentric clipping and area partition, explicit point/segment contact, duplicate shared-plane masks, unique internal face ownership with authored sides/winding, roundoff-bounded adjacent transverse-segment pairing with one canonical lower-cell result, stable identity and in-face authored-side direction, contact/coplanar separation, periodic-boundary and grid-edge ambiguity rejection before explicit ownership exists, sparse stable grid-bound active-face identity, separate multiple-sheet crossing/coplanar references without false union claims, provenance-keyed canonical shared-edge stitching, explicit opening/grid-edge endpoints and bounded degree connectivity without false closure or region claims, deterministic winding-directed open/closed chains with consistent authored region pairs and rejection of branches/conflicts, right-handed simple-loop signed area/centroid/enclosed-region classification with self-intersection and degenerate rejection, bounded non-touching smallest-parent nesting with authored-region continuity and exact per-region face-area closure while unresolved open/coplanar/boundary cases remain explicit, stable unique-area quadrature IDs, physical metadata, and exact negative/positive side-cell provenance, one immutable aggregate remap with nested epoch identity and per-stage plus aggregate bounds, ordered one-sided pressure-difference normal traction with exact equal-pressure cancellation, conservative load binding, and no shared-plane force double counting or secondary aerodynamic path |
| SimWing scene/structure/viewer foundations | `simwing_scene`, `simwing_model_scene_export`, `simwing_model_scene_real_export`, `simwing_structure`, `simwing_scene_structure`, `simwing_scene_fluid_surface`, `simwing_scene_fluid_surface_transfer`, `simwing_fluid_structure_bridge`, `simwing_projected_flag_case`, `simwing_ram_air_cell_case`, `simwing_viewer_protocol`, `simwing_vector_glyphs`, `simwing_structure_frame`, plus `softwing_material`/`softwing_cell_volume` when core primitives change, `softwing_suspension_checkpoint` for payload/suspension state, and `softwing_checkpoint_persistence` for opaque body-state wire changes; preserve canonical stable-ID surface ownership, authored winding/region sides/material/opening order, accepted Structure-epoch binding, exact topology-bound uniform/quadrature transfer and load application, deterministic multi-panel load addition on shared nodes, exact fixed mouth clamps, transactional bounded rejection, owning deterministic vector geometry, relative magnitude/direction, dimensional auto-scaling, bounded deterministic sampling, topology-template binding, and transactional decode |
| SimWing fluid grid/projection/interface | `simwing_scene_surface_geometry`, `simwing_scene_surface_intersection`, `simwing_scene_surface_clipping`, `simwing_fluid_projection`, `simwing_fluid_interface_jump`, `simwing_fluid_porous_interface`, `simwing_fluid_porous_topology`, `simwing_fluid_moving_interface`, `simwing_fluid_control_volume`, `simwing_fluid_cut_surface`, `simwing_fluid_checkpoint`, `simwing_fluid_diffusion`, `simwing_fluid_advection`, `simwing_fluid_variable_advection`, `simwing_fluid_projected_advection`, `simwing_fluid_evolution`; preserve fingerprinted accepted surface epochs, deterministic cell-major padded-AABB candidates, normalized-SAT pruning of false positives, conservative contact/two-cell grid-plane retention, tolerance bounded by broad-phase padding, complete recomputed narrow-phase validation, exact barycentric point/segment/area clipping, exact shared-coordinate preservation across later-axis clipping, analytic area partition away from coincident boundaries, explicit duplicate boundary-plane masks before ownership selection, explicit out-of-domain rejection, count/byte bounds without cut-cell claims, discrete gradient/divergence adjointness, periodic momentum, non-increasing projection/diffusion/committed-transport energy, transactional failure and composed all-stage rollback, deterministic replay, Taylor-Green invariance, manufactured second-order pressure and viscous spatial-eigenvalue convergence, exact zero/uniform/Nyquist viscous modes and the `0.5` per-stage stability boundary, exact Euler-stage SSPRK2 composition and observed second-order viscous temporal convergence, bounded donor-cell uniform transport with exact CFL-one shift and divergence commutation, exact uniform delegation from the variable path, divergence-free staggered control-volume flux closure, safe nonlinear self-advection aliasing, observed first-order uniform and variable-shear refinement, exact limited-MC SSPRK2 stage composition, discontinuous-pulse bounds, confined intermediate Euler energy exception, and near-second-order smooth-wave L1 refinement, exact four-stage projected nonlinear SSPRK2 composition with donor/MC selection, repeated stage eligibility, observed fixed-grid second-order nonlinear temporal refinement, translating Taylor-Green donor/MC spatial refinement with time error suppressed by `dt` proportional to `h^2`, exact symmetric half-viscosity/projected-transport/half-viscosity Strang composition with donor/MC selection, closed sub-integrator energy ledger, full-flow rollback, observed second-order temporal refinement, analytic viscous translating-Taylor-Green limited-MC L1 ratios in `[3.0, 5.0]` across 12/24/48 grids with `dt <= 0.12 h^2`, exact viscous pre-sizing and CFL/limited-bound equal-step subcycling retries, manual final-schedule equivalence, bounded substep-limit rollback, and fatal projection failure without retry, exact standalone-versus-composed uniform/nonlinear and Euler/SSPRK2 stage equivalence, stable region/interface orientation, static sharp-jump balance through both projected SSPRK2 stages and every Strang substep, bit-exact empty-jump delegation, foreign-topology rejection before mutation, monotone sign-correct Darcy-Forchheimer forward/inverse evaluation, canonical X/Y/Z porous sampling, relative sheet velocity, per-tile area/flux/nonnegative dissipation, strict open-segment porous placement, deterministic bidirectional one-face selection on X/Y/Z, signed periodic-image wraps, exact-boundary/skipped-segment rejection, complete deterministic X/Y/Z planar crossing assembly with validated identity/regions/material/kinematics and strict unwrapped placement, stage-resolved symmetric planar porous motion with retained/adjacent topology continuity, immutable two-stage identity/material, owned unwrapped epochs, in-step `3 -> 0` positive wrap, and transactional metadata/numerical failure, compatible flux-driven porous pressure loss without spurious velocity, bounded deterministic endpoint/midpoint Picard closure against prescribed jumps, analytic linear and nonlinear uniform flow, scalar midpoint-oracle and work/dissipation/energy agreement, explicit oriented jump force/impulse and power/work ledgers, orientation symmetry, heterogeneous tile response, independent velocity/jump residuals, and nonlinear rollback, exact X/Y/Z face constraints, separate adjacent-pressure and complete direct-enforcement reaction ledgers, canonical face-tile geometry/traction ledgers, per-region compatibility/gauges, analytic translating-slab pressure impulse/work, open-piston partial-cell/surface-sweep/opening-flux GCL closure, exact X/Y/Z one-plane rebase volume continuity, accepted physical cut-plane area/force/moment/power, macro-step-average endpoint resampling, periodic-image closure, immutable grid/topology-bound accepted-state checkpoint replay, deterministic bounded persistent field/interface/diagnostic round trips, rebased topology persistence, checksum enforcement, and transactional corruption/truncation/trailing-data/sample/face/region/surface-limit rejection |
| SimWing moving porous full-flow case | `simwing_moving_porous_flow_case`, `simwing_fsi_moving_porous_flow_headless`, `simwing_fsi_moving_porous_flow_checkpoint_write`, `simwing_fsi_moving_porous_flow_checkpoint_resume`, `simwing_fsi_moving_porous_flow_checkpoint_verify`, `simwing_fsi_moving_porous_flow_checkpoint_replay_limit`, `simwing_fsi_moving_porous_flow_rejects_foreign_checkpoint`; preserve deterministic owning frames, the first-step `face=3,image=0` to `face=0,image=1` staged wrap, continuous accumulated fluid state, five topology rebases through the second wrap at step 101, complete final sheet/pump crossing planes, both unwrapped stage epochs, closed kinematic residual, positive porous dissipation, bounded momentum residual, immutable earlier frames, initial/ordinary/second-wrap in-memory checkpoint replay, private-payload and case/grid/step/kinematic/topology metadata binding, deterministic bounded/checksummed `SWMF` round trips with explicit fluid fields and sharp crossings, bounded canonical diagnostic regeneration, corruption/truncation/trailing-data/version/reserved/definition/grid/step/topology/field/crossing/byte/sample/replay-limit rejection without output mutation, transactional rejected restore, same-path additional-step CLI resume after the second wrap, absolute autosave cadence, final-write deduplication, exact checkpoint-write telemetry, fail-fast rejection before trace creation when batch checkpoint output would exceed the replay ceiling, foreign-format rejection before trace creation, completed traces, and Qt-free headless execution |
| SimWing periodic fluid worker/snapshots | `simwing_periodic_flow_case`, `simwing_fsi_periodic_flow_headless`, `simwing_fsi_periodic_flow_checkpoint_write`, `simwing_fsi_periodic_flow_checkpoint_resume`, `simwing_fsi_periodic_flow_checkpoint_verify`, `simwing_fsi_checkpoint_rejects_foreign_case`, `simwing_fsi_checkpoint_interval_requires_output`, `simwing_viewer_protocol`; preserve transactional advance through frame validation, accepted-only publishing, owning cell-centre pressure/speed/velocity fields, exact MAC divergence, diagnostic centred-curl vorticity, the discrete Taylor-Green vorticity oracle, stable IDs, immutable grid/definition-bound checkpoint metadata and payload, initial/same/rebuilt-worker bit-identical replay, deterministic bounded little-endian persistent round trips, payload checksum, corruption/truncation/trailing-data/limit rejection without output mutation, additional-step CLI resume, same-path atomic checkpoint replacement followed by successful decode, absolute accepted-step autosave cadence, duplicate-final-write suppression, exact checkpoint-write telemetry, required autosave output, foreign-case flag rejection, rejected-restore non-mutation, completed traces, and Qt-free headless execution |
| SimWing pressure-jump snapshots/worker | `simwing_pressure_jump_frame`, `simwing_pressure_jump_case`, `simwing_fsi_pressure_jump_headless`, `simwing_viewer_protocol`; preserve owning deterministic cell samples, a separate oriented quad and triangle fields for every ordered crossing, periodic X/Y/Z placement, region-sided normals, diagnostic layered pressure reconstruction, the 48-crossing `-125/+125 Pa` split slab, zero spurious flow, bit-identical repeated static projection, completed traces, and Qt-free headless execution |
| SimWing porous-flow worker | `simwing_fluid_porous_interface`, `simwing_porous_flow_case`, `simwing_fsi_porous_flow_headless`, `simwing_pressure_jump_frame`; preserve the implicit-midpoint nonlinear plug solve, orientation symmetry, unforced dissipation, transactional invalid-input rejection, independently closed pressure impulse and driving-work/porous-loss/kinetic-energy ledgers, analytic steady Darcy-Forchheimer speed and pressure loss, complete endpoint Strang/SSPRK2 evolution with both pressure stages retaining all crossings and no spurious velocity, owning porous-layer/global diagnostics, deterministic replay, completed traces, and Qt-free headless execution |
| SimWing worker control | `simwing_worker_control_protocol`, `simwing_worker_control_stream`, `simwing_periodic_flow_control`, `simwing_moving_porous_flow_control`, `simwing_open_piston_control`, `simwing_strong_piston_control`, `simwing_porous_sheet_control`, `simwing_fsi_control_stdio`, `simwing_fsi_moving_porous_flow_control_stdio`, `simwing_fsi_open_piston_control_stdio`, `simwing_fsi_strong_piston_control_stdio`, `simwing_fsi_porous_sheet_control_stdio`, `simwing_fsi_porous_sheet_collision_control_stdio`, `simwing_fsi_control_rejects_unsupported_case`; preserve distinct command/response magic, versioned bounded little-endian envelopes, checksums, nonzero request correlation outside the ready response, positive bounded advances, exact absolute step/time responses, produced-frame counts only on advance, bounded coded error text only on error, byte-deterministic round trips, cross-type rejection, transactional corruption/truncation/trailing-data/limit failures, self-framing without a host prefix, pre-allocation stream bounds, clean EOF only between envelopes, per-message flush, case-neutral execution with typed checkpoint adapters, accepted-frame publication at individual safe points, exact periodic, moving-porous-flow, open-piston, strong-piston, and porous-sheet checkpoint delegation, moving porous control advance through both wraps and exact next-frame replay, visible absolute state after output or numerical failure, no solver mutation on checkpoint/protocol failure, terminal stop, idempotent repeated stop, binary-only protocol stdout, explicit-stop trace completion, no viewer launch, exact end-to-end response sequences for all five CLI workers, moving-porous step-101 checkpoint and exact original/resumed step-102 traces, typed open-piston and porous-sheet replay through their first topology rebases, periodic/open-piston/strong-piston three-frame traces with step-two checkpoint replay into step three, porous-sheet rebased step-330 checkpoint replay into step 331 followed by pump-collision failure with no rejected frame, process-level terminal checkpoint persistence, restored terminal Ready and repeated failure with an empty completed trace, restored ordinary Ready absolute state, and resumed ordinary traces containing only the exact next accepted frame |
| SimWing conservative transfer | `simwing_transfer`; preserve stable topology/Structure binding, analytic uniform and barycentric-quadrature area/force/moment, rigid translation/rotation power, independent ledger closure, additive nodal load application, and rejection before mutation for foreign results/structures |
| SimWing macro-step coupling | `simwing_coupling`; preserve strictly ordered local sample time, topology/duration binding, analytic moving-piston impulse and pressure-volume work, independent temporal ledger closure, momentum delivery through XPBD, deterministic replay, pre-load checkpoint rollback on failure, definition-bound vector Aitken relaxation with analytic affine convergence, explicit factor bounds, deterministic degenerate fallback, exact iteration replay, and transactional update/restore rejection, topology-bound maximum-nodal residual reduction from saved-baseline/previous/current kinematics and consecutive tractions with world-origin/bulk-velocity invariance and foreign-binding rejection, deterministic convergence decisions that require displacement/velocity/traction to each pass absolute and floor-stabilized relative tolerances, enforce minimum iterations, and distinguish convergence from iteration-budget exhaustion, and a checkpointable macro-step iteration owner with exact replay, transactional failure, and terminal convergence/exhaustion |
| SimWing composite coupling rollback | `simwing_coupled_state`; preserve interface-identity binding, accepted moving-interface fluid capture, complete Structure/fluid/iteration replacement validation before commit, no-throw three-owner commit, exact baseline and next-mutation replay, current unaccepted-fluid checkpoint rejection, transactional foreign composite/Structure/fluid-topology/iteration rejection, definition-bound bounded macro-step retries with hard minimum-step/retry-count failure, exact replay, and transactional sequencing/checkpoint rejection, fresh-baseline macro-step ownership whose integrated retry transition restores all three owners before activating the reduced step, a distinct transactional solver-only rewind that returns Structure/fluid to baseline while preserving advanced iteration state, and the bounded generic loop's large-step exhaustion/reduced-step convergence, per-run physical baseline, accepted final state, terminal-failure rollback, unaccepted-epoch rejection, and exception rollback |
| SimWing fluid/structure bridge and piston workers | `simwing_fluid_structure_bridge`, `simwing_piston_case`, `simwing_porous_sheet_case`, `simwing_open_piston_case`, `simwing_fsi_piston_headless`, `simwing_fsi_porous_sheet_headless`, `simwing_fsi_porous_sheet_checkpoint_write`, `simwing_fsi_porous_sheet_checkpoint_resume`, `simwing_fsi_porous_sheet_checkpoint_verify`, `simwing_fsi_porous_sheet_wrapped_checkpoint_write`, `simwing_fsi_porous_sheet_wrapped_checkpoint_resume`, `simwing_fsi_porous_sheet_wrapped_checkpoint_verify`, `simwing_fsi_porous_sheet_rejects_foreign_checkpoint`, `simwing_fsi_open_piston_headless`, `simwing_fsi_open_piston_rebase_headless`, `simwing_fsi_open_piston_checkpoint_write`, `simwing_fsi_open_piston_checkpoint_resume`, `simwing_fsi_open_piston_checkpoint_verify`, `simwing_fsi_open_piston_rejects_foreign_checkpoint`; preserve the strict uniform subset, planar face-resolved nonuniform transfer, stable surface/geometry binding, complete nonoverlapping coverage, per-face and aggregate area/force/moment/power closure, porous sheet-reaction ownership with prescribed-source exclusion and closed source/mapped impulse/work/dissipation ledgers, analytic porous-sheet midpoint momentum and pump-work/porous-loss/kinetic-energy closure across all six explicit ordinary MAC-face rebases including the `7 -> 0` signed-image wrap, direction-reversed pump/momentum closure through six negative rebases and the `0 -> 7` signed-image `-1` wrap, distinct directional provenance/fingerprint and transactional cross-direction checkpoint rejection, exact persistent continuation from every pre-pump topology epoch, process-level same-path additional-step restart at the ordinary wrapped `face=0,image=1` epoch, later next-image pump-collision rollback, immutable initial/ordinary/rebased/terminal porous-sheet checkpoint restore with topology bound to the accepted constitutive midpoint, deterministic bounded/checksummed `SWPS` round trips, complete topology version/axis/wrapped-face/signed-image ownership, nested Structure validation, explicit field and bounded-replay identity, exact next-frame or repeated terminal-collision replay, transactional public-metadata/magic/version/reserved/checksum/truncation/trailing/topology-version/topology-axis/topology-face/topology-image/byte/sample/replay/nested-state rejection, same-path additional-step porous-sheet CLI resume from the rebased epoch, absolute autosave cadence, final-write deduplication, and foreign-format rejection before trace creation, rigid-normal X/Y/Z grid/physical-plane correspondence and velocity binding, analytic impulse delivery, explicit actuator-versus-complete-CFD reaction, bit-identical replay through periodic topology crossings and composite checkpoint restore, deterministic bounded/checksummed composite persistence, ordinary/rebased decode-reencode and next-frame equivalence, transactional magic/version/checksum/truncation/trailing/topology/limit rejection including recomputed-checksum diagnostic identity/geometry/acceptance corruption, atomic same-path additional-step resume from a rebase epoch, absolute autosave cadence, final-write deduplication, cross-format rejection before trace creation, open-piston structure/fluid/actuator/system momentum residual below `1e-8 N*s` and energy residual below `2e-9 J`, accepted-only frames, and Qt-free headless execution |
| packaging/resources/CMake | configure from clean metadata, build Release, and run the full suite |

For `simwing_fluid_evolution`, the stage-resolved planar porous check also
requires `dt/2` trapezoidal position/normal-velocity consistency, explicit
finite nonnegative tolerances, an owned residual within that bound, and
transactional rejection of otherwise topology-valid motion that would
teleport the sheet.

For `simwing_fluid_interface_jump`, static sharp-jump balance includes ordered
same-face crossings, distinct open-interval positions, continuous region-chain
validation, deterministic signed aggregation, bit-identical split/compact slab
projection, and zero spurious pressure/flow for a balanced folded subcell
pocket. This does not close the moving folded-interface gate.

Also run `simwing_fluid_planar_pressure_jump` for changes to generic planar
epochs or layered pressure motion. It requires complete transverse X/Y/Z
planes, stable closed region chains, deterministic authored-order
canonicalization, adjacent and periodic one-face translation, and transactional
exact-boundary/skipped-segment rejection. The canonical same-face thin pocket
must retain both layer fractions while its current dense jump remains exactly
zero. Its static profile must retain the physical interval volumes, recover the
70 Pa regional difference under zero and nonzero volume-mean gauges, remain
rigid-translation invariant, and reject non-closing or inconsistent pressure
potentials. Only the separated-face state enters the current projection as an
intermediate cell pressure. The profile has no regional velocity unknown; this
exposes the unresolved subcell momentum/flux limitation and is not acceptance
of moving folded-interface physics.

The same test covers `planar_region_sweep.*`: require rigid and breathing
geometry/sweep closure on X/Y/Z, positive and negative periodic rebases, stable
endpoint identity, deterministic authored-order canonicalization, finite
positive duration, and count/region/byte bounds. Exact-boundary, skipped,
crossed, or foreign layer motion must reject transactionally. This is only the
regional GCL boundary ledger; no Eulerian relative-flux claim is permitted.

Also require the regional flux compatibility screen to keep rigid motion
exactly impermeable and continuous, expose breathing motion as non-impermeable
while closing `delta volume + integrated outward relative flow`, scale physical
flow with X/Y/Z area, retain periodic-rebase sealing, and reject corrupted
source profiles/ledgers plus invalid tolerance or storage policies. It remains
diagnostic; never interpret its least-squares interval velocity as a production
regional projection. Require deterministic nonzero fingerprints, and mutate
primitive/derived interval fields, surface identity, region order, aggregates,
finiteness, tolerance policy, and validator limits to prove deep rejection.

For `planar_region_opening_flow.*`, require sealed rigid motion to remain
feasible with zero flow and sealed breathing to fail per connected component
even when the global volume change cancels. One oriented opening must close the
canonical pocket with the correct sign; reversing it reverses flow. Parallel
openings split by area at equal normal velocity, a serial three-region graph
balances its intermediate region, and X/Y/Z follow physical swept volume.
Reject foreign/same-region/duplicate/zero-area openings, corrupted immutable
outputs, invalid settings, and interval/region/opening/owned-byte/dense-byte/
factorization-work limit violations. Do not treat this analytic graph oracle as
scene-v2 opening geometry or production fluid state.

For `planar_region_opening_power.*`, require rigid zero power, `112 W` external
pressure-power demand for canonical 70 Pa inflation, passive `112 W` release
under reversed volume motion, and invariance to authored opening orientation.
Parallel-link power must follow the area-weighted flow split. The serial graph
must retain its `+40/-24 W` local links, `+16 W` net opening power, `-16 W`
regional volume power, and zero net external demand while still flagging the
uphill link. Reject infeasible allocation sources, mutated fingerprints/
opening/region/aggregate ledgers, invalid tolerances, and opening/region/byte or
nested source-limit violations. Never treat the audit as dynamic-pressure or
momentum closure.

For `planar_region_fragment.*`, require the canonical 16-cell grid to publish
24 controls, with four axial cells split into three exterior/pocket/exterior
fragments and the pocket retained as four `0.6 m^3` controls bounded by layer
IDs 10/20. Require unique deterministic fragment IDs, 13.6/2.4 `m^3` regional
volume and 70 Pa pressure difference, exact per-cell volume and first-moment
closure, X/Y/Z physical scaling, rigid translation, and signed periodic
rebases. Reject mutated fragment/region/cell/fingerprint ledgers and interval/
region/cell/fragment/byte limit violations. Do not infer connectivity, a
regional velocity basis, mimetic pressure acceptance, or production wiring.

For `planar_region_fragment_topology.*`, require the canonical 24 controls to
close as 72 paired links: 64 same-region grid links, eight nonconductive layer
walls, and 28 periodic grid links. Require one 20-fragment/13.6 `m^3` exterior
component and one 4-fragment/2.4 `m^3` pocket component. Each fragment must
close exactly six incidences and its analytic boundary area; the 56 `m^2`
unique area must equal half the 112 `m^2` incident area. Each layer must retain
four 1 `m^2` wall tiles, 0.4 m center distance, authored orientation, and the
signed +/-70 Pa jump while publishing zero conductance. Require X/Y/Z closure,
stable graph/component identity during within-segment rigid motion, signed
periodic rebase closure, and deep rejection of mutated link/fragment/component/
fingerprint products plus link/component/byte/nested-source limits. Do not
infer a regional face velocity, opening conductance, momentum operator,
mimetic pressure acceptance, or production wiring.

For `planar_region_fragment_pressure_operator.*`, require 24 rows and 128
directed entries from the 64 same-region links, with all eight layer walls
excluded. Require two deterministic component gauges, exact regional-pressure
null action, zero integrated action per component, symmetric arbitrary-vector
action, and nonnegative energy equal to the link sum. The canonical unique and
diagonal geometry weights must be `160/3 m` and `320/3 m`; exterior/pocket
component weights must be `728/15 m` and `4.8 m`. Require X/Y/Z closure and
stable row/entry/component identity during within-segment motion while metric
fingerprints update. Reject mutated operator/row/entry/component/member/source
products, wrong-sized or non-finite pressure fields, and row/entry/component/
byte/nested-topology limit violations. Do not infer a pressure RHS or solve,
regional velocity correction, opening link, mimetic acceptance, or production
ownership.

For `planar_region_fragment_pressure_solve.*`, require deterministic recovery
of a manufactured two-component, zero-volume-mean correction within `3e-11
Pa`, with recomputed maximum residual below `2e-12 Pa*m` and correction
volume-mean residual below `3e-16 Pa`. A zero RHS must remove only correction
gauges and preserve every static +/-70 Pa wall jump. Require X/Y/Z zero-RHS
closure, independent component compatibility diagnostics, and bit-exact warm-
start rollback for incompatible or iteration-truncated solves. Reject mutated
operator/source products, wrong-sized or non-finite fields, non-finite/empty
tolerance policies, and a zero iteration bound. Do not infer a regional
velocity divergence, physical pressure RHS, projection acceptance, or
production ownership.

For `planar_region_fragment_pressure_projection.*`, require deterministic
manufactured-divergence cancellation below `3e-14 m3/s` through X/Y/Z, exact
preservation of uniform wall-tangential flow, and exact zero velocity on every
pressure-layer wall. Require the physical `-rho/dt` integrated RHS, the
matching `dt/rho` same-region link correction, independent corrected-
continuity recomputation, source fingerprints, and bounded working storage.
Pressure non-convergence must preserve both link velocity and pressure warm
start bit-for-bit. Reject moving geometry, nonzero wall flow, mutated sources,
wrong-sized/non-finite fields, invalid density/time/tolerances, and working or
nested-operator limit violations. Do not infer topology-rebase ownership,
face momentum mass, kinetic-energy acceptance,
opening conductance, or production ownership.
For its moving overload, require exact duration/rate binding and the physical
`-rho/dt * (dV/dt + net outward flow)` RHS. A topology-stable rigid `0.1 m`
translation must close moving continuity below `1e-11 m3/s` on X/Y/Z while
preserving zero layer-wall relative flow; reprojecting its balanced field must
remain stable. A sealed breathing pocket must expose its `1.6 m3/s` component
deficit and corresponding `3.84 Pa*m` incompatibility, rolling back both
fields. Reject corrupt or foreign volume rates, duration mismatch, nested rate
limits, and truncated solves. Do not infer opening flux, topology-transition
ownership, momentum/energy acceptance, or production integration.

For `planar_region_fragment_volume_rate.*`, require canonical breathing to
publish `+1.6/-1.6 m3/s` pocket/exterior component rates while every fixed
Cartesian cell and the global domain close. Require bit-exact
`change = current - previous`, strictly positive reconstructed previous
volumes, exact boundary displacement/velocity provenance, static zero rates,
rigid-motion zero component change with nonzero local exchange, and X/Y/Z
cross-section scaling. Require deterministic source fingerprints and fragment,
cell, region, component, global, and storage ledgers. Reject any layer crossing
a Cartesian segment, mutated products/sources, and fragment/cell/region/
component/byte/nested-topology limit violations. Do not infer topology-rebase
appearance or retirement, pressure RHS composition beyond the topology-stable
moving overload, face momentum mass, opening flow, or production ownership.

For `planar_region_fragment_velocity_metric.*`, require one shared normal-
velocity degree of freedom per same-region Cartesian link and two independent
one-sided trace degrees per pressure-layer wall. The canonical X case must
publish 64 shared plus 16 trace degrees, with `44.8 m3` shared, `3.2 m3`
trace, and `48 m3` total dual volume. Require every fragment's six incidences,
every exterior/pocket component, and the domain to close physical volume
independently on X/Y/Z; rigid motion must preserve stable degree identity while
metric fingerprints update, and breathing must update one-sided wall volume.
Reject mutated products/sources, inconsistent half distances, and DOF,
fragment, component, byte, or nested-topology limit violations. Do not infer
velocity values, density, kinetic-energy acceptance, wall prescription,
advection, pressure acceptance, or production ownership.

For `planar_region_fragment_velocity_state.*`, require exact binding to every
metric degree, positive finite density, finite normal velocity/momentum/energy,
and half-volume accumulation back to every fragment. At `1.25 kg/m3`, require
the canonical uniform `[2,-0.5,0.25] m/s` state to close `20 kg` per axis,
`[40,-10,5] kg*m/s`, and `43.125 J`; distinct `+1/-2 m/s` traces on one wall
must remain exterior/pocket-owned and close `-2.5 kg*m/s` plus `3.25 J`
globally. Require X/Y/Z profile-axis coverage and topology-stable breathing to
preserve uniform global momentum/energy while component mass follows current
volume. Reject mutated products/sources, wrong-sized or non-finite samples,
non-positive density, aggregate overflow, and sample/fragment/component/owned-
byte/working-byte/nested-metric limit violations. Do not infer a wall velocity
prescription, pressure-projection certificate, transport acceptance, pressure
work, or production ownership.

For `planar_region_fragment_projection_energy.*`, require a source-bound
before/after state pair with projection density, one correction per shared
degree, `dt/rho * delta-p / distance` velocity change, `dt * area * delta-p`
pressure impulse, diagonal momentum closure, midpoint pressure-work/kinetic-
energy closure per degree/component/global ledger, and roundoff-zero volume-
weighted correction gauges. Static audits require exact-zero one-sided wall
traces, non-increasing kinetic energy, X/Y/Z continuity below `3e-14 m3/s`,
pure-gradient after-energy below `1e-26 J`, work residual below `3e-18 J`, and
bit-exact tangential null flow. Moving audits require exact volume-rate source
and duration binding, both traces at material-wall velocity, local
`dV/dt + flow` closure, and `delta-K = geometry-pressure work - correction
kinetic energy`. Rigid `0.1 m/s` translation from zero grid flow must receive
twice the correction energy as final geometry work, retain half as new kinetic
energy, and close affine/final-work residuals below `4e-13 J` on X/Y/Z; sealed
breathing must reject. Reject mismatched wall traces, gauge shifts, inconsistent
states, corrupt sources, wrong-sized/non-finite pressure, invalid settings, and
all correction/pressure/component/owned/working/nested state or volume-rate
limits. Do not infer authored static-jump work, topology rebase, momentum
transport, or production ownership.

The Windows reference fixture is `tests/fixtures/3.28/leparagliding.txt` with
adjacent `gnuC2.txt`; expected reports are under `tests/reference/3.28`.
`leparagliding.dxf` and `lep-3d.dxf` must match the native gfortran run by SHA-256;
reports require identical structure/fields with only the documented 0.00015
display tolerance. Do not regenerate the oracle from the translated engine.

Geometry-preprocessor fixtures are under `tests/fixtures/pre-1.6` and compare
against native Fortran output within 0.03 cm/0.06 degrees.

The current deterministic Playground guards (details and interpretation in
`docs/legacy/leparagliding/CONTINUE.md`) are:

```powershell
.\build\bin\Release\softwing-bench.exe build\aero\gnuC2\lep-sim.json --shape --csv --legacy-pressure
.\build\bin\Release\softwing-bench.exe build\aero\gnuC2\lep-sim.json --pressure-acceptance
.\build\bin\Release\softwing-bench.exe build\aero\gnuC2\lep-sim.json --glide 300
.\build\bin\Release\softwing-bench.exe build\aero\Swoop\lep-sim.json --tuck 250 --substeps 60 --iterations 4
```

The `--dive` scenario is chaotic and non-monotonic; never call one run a
regression. Sweep disturbances and compare distributions. The tunnel CSV guard
is the cheap bit-comparable signal.

### Verifying GUI changes without disturbing the user

A running Studio locks `LEparagliding.exe`. Do not terminate it. Windows allows
renaming a running executable aside before relinking, but preserve and report
any existing renamed copies. Prefer a narrower non-GUI target when sufficient.

For screenshots, do not steal focus or inject global input. Background Qt
windows can be captured with `PrintWindow(..., PW_RENDERFULLCONTENT)`, but a
native/OpenGL child can still be blank or flaky. A blank capture is not proof
that a widget failed to paint. Use an offscreen Qt harness and `widget->grab()`
for ordinary widgets; verify native GL behavior through its status/metrics and
an explicitly launched instance. The detailed machine recipe lives in
`.claude/skills/verify/SKILL.md` even though it is not a global Codex skill.

## Code and change conventions

- Match the existing four-space style and local brace/Qt conventions; there is
  no repository-wide formatter configuration to run blindly.
- Use `QStringLiteral`/`QLatin1String` for fixed Qt text. Keep reusable math and
  parsers in standard C++ when that permits small non-widget tests.
- Preserve exact byte layout where a helper promises format-preserving token
  replacement or Fortran parity. Do not normalize whole design sections as a
  side effect of one edit.
- Prefer atomic file writes (`QSaveFile`) for user documents/exports where the
  existing layer supports them.
- Keep long calculations, nesting, and sweeps off the GUI thread. Worker
  callbacks execute on the worker thread unless explicitly queued; never touch
  widgets from them.
- Preserve deterministic results. The nester deliberately fixes tie-breaking;
  XPBD parallel sweeps deliberately define reproducibility by worker count.
- Comments in the model, flat-part, and Playground code often record rejected
  physical/geometric approaches with measurements. Treat them as design
  constraints, not cleanup targets. Update the matching design document when
  changing the conclusion.
- The app embeds XFLR5 in its own `QApplication`. Do not reintroduce XFLR5's
  bootstrap or application-wide style mutations.
- New outputs loaded from a temporary preview must be consumed before
  `calculationDirectory_` is reset.

## Documentation and release map

- `README.md`: user-visible features, output contract, portable build recipe,
  architecture summary, and release overview.
- `UPSTREAM.md`: independent-project boundary, bootstrap commit, licensing,
  and upstream remote.
- `docs/legacy/leparagliding/`: inherited planning, machine/release, Playground,
  performance, Print/Cut, and claim-boundary records. These remain load-bearing
  for imported subsystems but are not the active SimWing roadmap.
- `docs/coupled-fsi-architecture.md`: target XPBD-CFD architecture, scene and
  process boundaries, strong-coupling algorithm, validation requirements, and
  phased delivery gates for inflation, collapse, and reinflation.
- `third_party/xflr5/PROVENANCE.md` and `src/softwing/README.md`: vendor origins
  and local divergence contracts.

Release version comes only from `CMakeLists.txt`. GitHub CI runs on `v*` tags
or manual dispatch, builds Windows/Linux/macOS, and packages a Windows zip,
Linux AppImage, and macOS DMG. Every release needs a hand-written description
with a first-visitor summary, version changes, and platform install notes; see
`docs/legacy/leparagliding/CLAUDE.md`. The signed Windows installer is built locally with
`installer/build_installer.ps1`; commit first because configure bakes the HEAD
hash into `installer/version.iss` and the installer filename.
