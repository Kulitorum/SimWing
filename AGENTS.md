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
.\build\bin\Release\simwing-fsi.exe [--case structural|piston|open-piston|periodic-flow|porous-flow|porous-sheet|pressure-jump] [--steps N] [--trace <file>] [--checkpoint-in <file>] [--checkpoint-out <file>] [--checkpoint-every N] [--control-stdio] [--no-viewer]
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
- `simwing_piston_case`: Qt-free visible verification harness crossing fluid
  projection, the face-resolved bridge, temporal coupling, XPBD acceptance, and
  immutable viewer frames.
- `simwing_porous_sheet_case`: Qt-free midpoint oracle coupling a
  pressure-driven porous sheet to XPBD. It closes fluid/sheet/pump momentum,
  pump work, porous dissipation, and kinetic-energy ledgers before publishing
  the translated sheet. Its immutable composite checkpoint restores Structure,
  fluid fields, and the last accepted ledger with exact next-frame replay. The
  bounded checksummed persistent envelope stores the existing Structure codec
  plus explicit topology and MAC fields, then regenerates the ledger by bounded
  deterministic replay before publishing a decode. It explicitly rebinds the
  sheet to the next MAC face at its first dual-cell crossing and rejects the
  later collision with the prescribed pump surface.
- `simwing_open_piston_case`: Qt-free driven open-piston harness crossing a
  connected-fluid projection, partial-cell control-volume geometry, opening
  transport, exact planar topology rebasing, CFD reaction, XPBD acceptance,
  composite structure/fluid checkpoint replay, deterministic bounded
  persistent restart encoding exposed through the headless checkpoint flags,
  and immutable viewer frames.
- `simwing_scene_structure`: deterministic scene-v2 membrane, per-sheet
  bending, junction, and cable assembly into `simwing_structure`.
- `simwing_fluid`: Qt-free periodic staggered-grid field operators and
  transactional pressure projection verification kernel, including a
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
- `simwing_periodic_flow_control`, `simwing_open_piston_control`, and
  `simwing_porous_sheet_control`: typed adapters binding their numerical owners
  and complete checkpoint payloads to the shared control session.
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

- `src/fsi/scene.{h,cpp}` owns scene-v2.1's SI/right-handed/Z-up contract,
  stable-ID entities, authored fabric-sheet identity, two-sided fluid regions,
  materials, openings, paired seam chains, pilot, suspension junctions,
  attachments, suspension lines, validation, and bounded deterministic binary
  round trips. Diagonal/mini-rib sheets may retain one connected cell on both
  sides; skin and ribs may not.
- `src/model/nurbs_model.{h,cpp}` can export scene-v2 directly from analytical
  captures without reading `lep-sim.json`. It preserves authored open intakes,
  exact triangulated rib/crossport faces, internal sheets, and the segmented
  suspension graph. Physical fabric/line/pilot settings are mandatory and no
  engine CLI writes this scene yet because the design format has no authoritative
  source for them. Its present local intrinsic charts are not manufacturing
  flat-pattern UVs, and it does not yet author paired seam chains.
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
- `src/fsi/canonical_case.{h,cpp}` and `tools/simwing_fsi_main.cpp` are the
  first end-to-end worker slice. The case is an analytic structural harness,
  not aerodynamic truth; it writes only accepted steps and launches the
  sibling viewer by default. `--no-viewer` must remain Qt-free and unthrottled.
  `--checkpoint-in`, `--checkpoint-out`, and `--checkpoint-every` are supported
  by periodic-flow, open-piston, and porous-sheet workers; all restore before
  trace creation,
  advance additional steps, autosave only accepted absolute-step multiples,
  and atomically replace the final file without duplicating the last write.
  Periodic-flow, open-piston, and porous-sheet `--control-stdio` are also
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
- `src/fsi/porous_sheet_case.{h,cpp}` is selected with `--case porous-sheet`.
  Its analytic linear-resistance midpoint relation drives the same accepted
  nonuniform porous projection, sheet-reaction bridge, temporal transfer, and
  XPBD boundary used by later FSI. Every step closes the prescribed pump
  impulse/work against fluid and sheet momentum/kinetic energy plus porous
  dissipation. The plane may translate rigidly across one dual-cell boundary;
  the accepted midpoint sample commits the next MAC-face binding and exposes
  its coordinate and rebase count. A later pump-surface collision is a hard
  topology rejection. Its in-memory checkpoint keeps the nested Structure
  checkpoint, velocity, pressure, and accepted coupled diagnostics private,
  revalidates the canonical epoch before commit, and reproduces the exact next
  frame in the same or an equivalent rebuilt worker.
  `porous_sheet_checkpoint_persistence.{h,cpp}` wraps that state in the distinct
  `SWPS` envelope. It reuses Structure's persistent codec, stores the four MAC
  scalar fields and topology epoch explicitly, and uses bounded replay to
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
  checkpoint flags are also supported by open-piston but rejected for the
  structural and sealed-piston cases, and interval mode requires an output
  path.
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
  changing its fields or nested diagnostics.
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
  restores the checkpoint from before load application on any failure.
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

There are 91 configured tests on Windows. The Fortran-reference test is
Windows-only; local `gui_smoke` and `studio_model_smoke` exercise display/model
paths that release CI deliberately excludes from its offscreen test command.
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
| SimWing scene/structure/viewer foundations | `simwing_scene`, `simwing_model_scene_export`, `simwing_model_scene_real_export`, `simwing_structure`, `simwing_scene_structure`, `simwing_viewer_protocol`, `simwing_vector_glyphs`, `simwing_structure_frame`, plus `softwing_material`/`softwing_cell_volume` when core primitives change, `softwing_suspension_checkpoint` for payload/suspension state, and `softwing_checkpoint_persistence` for opaque body-state wire changes; preserve owning deterministic vector geometry, relative magnitude/direction, dimensional auto-scaling, bounded deterministic sampling, topology-template binding, and transactional decode |
| SimWing fluid grid/projection/interface | `simwing_fluid_projection`, `simwing_fluid_interface_jump`, `simwing_fluid_porous_interface`, `simwing_fluid_moving_interface`, `simwing_fluid_control_volume`, `simwing_fluid_cut_surface`, `simwing_fluid_checkpoint`, `simwing_fluid_diffusion`, `simwing_fluid_advection`, `simwing_fluid_variable_advection`, `simwing_fluid_projected_advection`, `simwing_fluid_evolution`; preserve discrete gradient/divergence adjointness, periodic momentum, non-increasing projection/diffusion/committed-transport energy, transactional failure and composed all-stage rollback, deterministic replay, Taylor-Green invariance, manufactured second-order pressure and viscous spatial-eigenvalue convergence, exact zero/uniform/Nyquist viscous modes and the `0.5` per-stage stability boundary, exact Euler-stage SSPRK2 composition and observed second-order viscous temporal convergence, bounded donor-cell uniform transport with exact CFL-one shift and divergence commutation, exact uniform delegation from the variable path, divergence-free staggered control-volume flux closure, safe nonlinear self-advection aliasing, observed first-order uniform and variable-shear refinement, exact limited-MC SSPRK2 stage composition, discontinuous-pulse bounds, confined intermediate Euler energy exception, and near-second-order smooth-wave L1 refinement, exact four-stage projected nonlinear SSPRK2 composition with donor/MC selection, repeated stage eligibility, observed fixed-grid second-order nonlinear temporal refinement, translating Taylor-Green donor/MC spatial refinement with time error suppressed by `dt` proportional to `h^2`, exact symmetric half-viscosity/projected-transport/half-viscosity Strang composition with donor/MC selection, closed sub-integrator energy ledger, full-flow rollback, observed second-order temporal refinement, analytic viscous translating-Taylor-Green limited-MC L1 ratios in `[3.0, 5.0]` across 12/24/48 grids with `dt <= 0.12 h^2`, exact viscous pre-sizing and CFL/limited-bound equal-step subcycling retries, manual final-schedule equivalence, bounded substep-limit rollback, and fatal projection failure without retry, exact standalone-versus-composed uniform/nonlinear and Euler/SSPRK2 stage equivalence, stable region/interface orientation, static sharp-jump balance through both projected SSPRK2 stages and every Strang substep, bit-exact empty-jump delegation, foreign-topology rejection before mutation, monotone sign-correct Darcy-Forchheimer forward/inverse evaluation, canonical X/Y/Z porous sampling, relative sheet velocity, per-tile area/flux/nonnegative dissipation, compatible flux-driven porous pressure loss without spurious velocity, bounded deterministic endpoint/midpoint Picard closure against prescribed jumps, analytic linear and nonlinear uniform flow, scalar midpoint-oracle and work/dissipation/energy agreement, explicit oriented jump force/impulse and power/work ledgers, orientation symmetry, heterogeneous tile response, independent velocity/jump residuals, and nonlinear rollback, exact X/Y/Z face constraints, separate adjacent-pressure and complete direct-enforcement reaction ledgers, canonical face-tile geometry/traction ledgers, per-region compatibility/gauges, analytic translating-slab pressure impulse/work, open-piston partial-cell/surface-sweep/opening-flux GCL closure, exact X/Y/Z one-plane rebase volume continuity, accepted physical cut-plane area/force/moment/power, macro-step-average endpoint resampling, periodic-image closure, immutable grid/topology-bound accepted-state checkpoint replay, deterministic bounded persistent field/interface/diagnostic round trips, rebased topology persistence, checksum enforcement, and transactional corruption/truncation/trailing-data/sample/face/region/surface-limit rejection |
| SimWing periodic fluid worker/snapshots | `simwing_periodic_flow_case`, `simwing_fsi_periodic_flow_headless`, `simwing_fsi_periodic_flow_checkpoint_write`, `simwing_fsi_periodic_flow_checkpoint_resume`, `simwing_fsi_periodic_flow_checkpoint_verify`, `simwing_fsi_checkpoint_rejects_foreign_case`, `simwing_fsi_checkpoint_interval_requires_output`, `simwing_viewer_protocol`; preserve transactional advance through frame validation, accepted-only publishing, owning cell-centre pressure/speed/velocity fields, exact MAC divergence, diagnostic centred-curl vorticity, the discrete Taylor-Green vorticity oracle, stable IDs, immutable grid/definition-bound checkpoint metadata and payload, initial/same/rebuilt-worker bit-identical replay, deterministic bounded little-endian persistent round trips, payload checksum, corruption/truncation/trailing-data/limit rejection without output mutation, additional-step CLI resume, same-path atomic checkpoint replacement followed by successful decode, absolute accepted-step autosave cadence, duplicate-final-write suppression, exact checkpoint-write telemetry, required autosave output, foreign-case flag rejection, rejected-restore non-mutation, completed traces, and Qt-free headless execution |
| SimWing pressure-jump snapshots/worker | `simwing_pressure_jump_frame`, `simwing_pressure_jump_case`, `simwing_fsi_pressure_jump_headless`, `simwing_viewer_protocol`; preserve owning deterministic cell samples, a separate oriented quad and triangle fields for every ordered crossing, periodic X/Y/Z placement, region-sided normals, diagnostic layered pressure reconstruction, the 48-crossing `-125/+125 Pa` split slab, zero spurious flow, bit-identical repeated static projection, completed traces, and Qt-free headless execution |
| SimWing porous-flow worker | `simwing_fluid_porous_interface`, `simwing_porous_flow_case`, `simwing_fsi_porous_flow_headless`, `simwing_pressure_jump_frame`; preserve the implicit-midpoint nonlinear plug solve, orientation symmetry, unforced dissipation, transactional invalid-input rejection, independently closed pressure impulse and driving-work/porous-loss/kinetic-energy ledgers, analytic steady Darcy-Forchheimer speed and pressure loss, complete endpoint Strang/SSPRK2 evolution with both pressure stages retaining all crossings and no spurious velocity, owning porous-layer/global diagnostics, deterministic replay, completed traces, and Qt-free headless execution |
| SimWing worker control | `simwing_worker_control_protocol`, `simwing_worker_control_stream`, `simwing_periodic_flow_control`, `simwing_open_piston_control`, `simwing_porous_sheet_control`, `simwing_fsi_control_stdio`, `simwing_fsi_open_piston_control_stdio`, `simwing_fsi_porous_sheet_control_stdio`, `simwing_fsi_control_rejects_unsupported_case`; preserve distinct command/response magic, versioned bounded little-endian envelopes, checksums, nonzero request correlation outside the ready response, positive bounded advances, exact absolute step/time responses, produced-frame counts only on advance, bounded coded error text only on error, byte-deterministic round trips, cross-type rejection, transactional corruption/truncation/trailing-data/limit failures, self-framing without a host prefix, pre-allocation stream bounds, clean EOF only between envelopes, per-message flush, case-neutral execution with typed checkpoint adapters, accepted-frame publication at individual safe points, exact periodic, open-piston, and porous-sheet checkpoint delegation, visible absolute state after output failure, no solver mutation on checkpoint/protocol failure, terminal stop, idempotent repeated stop, binary-only protocol stdout, explicit-stop trace completion, no viewer launch, exact end-to-end response sequences for all three workers, open-piston and porous-sheet replay through their first topology rebases, three accepted trace frames, step-two checkpoint replay into trace step three, restored Ready absolute state, and a resumed trace containing only the exact next accepted frame |
| SimWing conservative transfer | `simwing_transfer`; preserve stable topology/Structure binding, analytic uniform and barycentric-quadrature area/force/moment, rigid translation/rotation power, independent ledger closure, additive nodal load application, and rejection before mutation for foreign results/structures |
| SimWing macro-step coupling | `simwing_coupling`; preserve strictly ordered local sample time, topology/duration binding, analytic moving-piston impulse and pressure-volume work, independent temporal ledger closure, momentum delivery through XPBD, deterministic replay, and pre-load checkpoint rollback on failure |
| SimWing fluid/structure bridge and piston workers | `simwing_fluid_structure_bridge`, `simwing_piston_case`, `simwing_porous_sheet_case`, `simwing_open_piston_case`, `simwing_fsi_piston_headless`, `simwing_fsi_porous_sheet_headless`, `simwing_fsi_porous_sheet_checkpoint_write`, `simwing_fsi_porous_sheet_checkpoint_resume`, `simwing_fsi_porous_sheet_checkpoint_verify`, `simwing_fsi_porous_sheet_rejects_foreign_checkpoint`, `simwing_fsi_open_piston_headless`, `simwing_fsi_open_piston_rebase_headless`, `simwing_fsi_open_piston_checkpoint_write`, `simwing_fsi_open_piston_checkpoint_resume`, `simwing_fsi_open_piston_checkpoint_verify`, `simwing_fsi_open_piston_rejects_foreign_checkpoint`; preserve the strict uniform subset, planar face-resolved nonuniform transfer, stable surface/geometry binding, complete nonoverlapping coverage, per-face and aggregate area/force/moment/power closure, porous sheet-reaction ownership with prescribed-source exclusion and closed source/mapped impulse/work/dissipation ledgers, analytic porous-sheet midpoint momentum and pump-work/porous-loss/kinetic-energy closure across its first explicit MAC-face rebase, later pump-collision rollback, immutable initial/ordinary/rebased porous-sheet checkpoint restore, deterministic bounded/checksummed `SWPS` round trips, explicit topology epoch, nested Structure validation, explicit field and bounded-replay identity, exact next-frame replay, transactional public-metadata/magic/version/reserved/checksum/truncation/trailing/byte/sample/replay/nested-state rejection, same-path additional-step porous-sheet CLI resume from the rebased epoch, absolute autosave cadence, final-write deduplication, and foreign-format rejection before trace creation, rigid-normal X/Y/Z grid/physical-plane correspondence and velocity binding, analytic impulse delivery, explicit actuator-versus-complete-CFD reaction, bit-identical replay through periodic topology crossings and composite checkpoint restore, deterministic bounded/checksummed composite persistence, ordinary/rebased decode-reencode and next-frame equivalence, transactional magic/version/checksum/truncation/trailing/topology/limit rejection including recomputed-checksum diagnostic identity/geometry/acceptance corruption, atomic same-path additional-step resume from a rebase epoch, absolute autosave cadence, final-write deduplication, cross-format rejection before trace creation, open-piston structure/fluid/actuator/system momentum residual below `1e-8 N*s` and energy residual below `2e-9 J`, accepted-only frames, and Qt-free headless execution |
| packaging/resources/CMake | configure from clean metadata, build Release, and run the full suite |

For `simwing_fluid_interface_jump`, static sharp-jump balance includes ordered
same-face crossings, distinct open-interval positions, continuous region-chain
validation, deterministic signed aggregation, bit-identical split/compact slab
projection, and zero spurious pressure/flow for a balanced folded subcell
pocket. This does not close the moving folded-interface gate.

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
