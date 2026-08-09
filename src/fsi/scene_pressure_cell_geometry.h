#pragma once

#include "fluid/grid.h"
#include "scene_structure.h"

namespace simwing::fsi {

inline constexpr char scenePressureCellCaseChecksum[] =
    "sha256:simwing-scene-pressure-feedback-cell-v9";
inline constexpr char scenePressureCellCaseSolverId[] =
    "simwing-fsi-scene-pressure-feedback-worker-v9";
inline constexpr char scenePressureCellRefinementChecksum[] =
    "sha256:simwing-scene-pressure-refinement-cell-v1";

// Shared immutable construction for the unchanged visible pressure-cell
// tetrahedron and a separately checksummed skew tetrahedron used only by the
// offline refinement audit. They share the 4 m periodic domain and mechanical
// fixing policy; the audit varies only Cartesian counts for its own geometry.
[[nodiscard]] Scene makeScenePressureCellGeometry();

[[nodiscard]] Scene makeScenePressureCellRefinementGeometry();

[[nodiscard]] SceneStructureAssembly makeScenePressureCellAssembly(
    const Scene& scene);

[[nodiscard]] fluid::PeriodicCartesianGrid makeScenePressureCellGrid(
    fluid::GridCellCounts cellCounts = {4, 4, 4},
    fluid::Vector3 lowerMeters = {});

} // namespace simwing::fsi
