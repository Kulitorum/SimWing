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
fixed-grid temporal refinement.
Both the projected transport and Strang flow paths select donor-cell or limited
monotonized-central reconstruction explicitly; donor-cell remains the default.
The original single-stage transport and Euler/SSPRK2 viscosity modes compose
with the zero-mean
pressure solve in one transactional periodic fluid step. Every stage runs on
candidate fields; velocity and pressure commit together only when transport
bounds, viscous stability, projection convergence, and the final
momentum/energy ledger all pass. Focused rollback cases fail each stage
independently without changing either caller field. This is the first complete
nonlinear fluid evolution path, but donor-cell convection and first-order
operator splitting still keep it short of the intended second-order production
Navier-Stokes scheme.
A validated single-crossing sharp-interface
field now preserves a prescribed two-sided static pressure jump without
smearing or spurious flow, including across the periodic domain boundary.
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
delivers the same total impulse to structural momentum. Alongside the strict
uniform fluid-to-structure subset, a planar face-resolved bridge clips canonical
MAC tiles against structural triangles and conservatively maps nonuniform face
traction while closing per-face and aggregate area, force, moment, and power
ledgers. Its rigid-normal mode keeps those material patches while the physical
plane moves and the Eulerian grid plane rebases, provided transverse geometry
remains fixed and fluid/structural normal velocities agree. The
`simwing-fsi --case piston` harness crosses that
face-resolved bridge, temporal coupling, XPBD acceptance, and replayable viewer
frames with visible deterministic motion. The third
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
both ordinary and periodic topology rebases. Persistent checkpoint files and
worker control messages remain future protocol work. A real 3.28 fixture
also crosses export, structural assembly, one coupled pilot/suspension step,
checkpoint replay, and completed trace using synthetic physical settings.
Export still requires explicit physical material/pilot settings;
manufacturing-pattern UVs, exact authored attachment vertices, structural seam
assembly, curved or transversely deforming grid-to-surface correspondence,
general cut-cell pressure metrics, nonplanar topology events, multiple crossings
per face, and AMR CFD remain open. These worker cases validate the pipeline;
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
```

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
