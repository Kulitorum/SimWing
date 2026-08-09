#include "scene_fluid_capped_face_partition.h"
#include "scene_structure.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>

namespace {

using namespace simwing::fsi;

int failures = 0;

void check(const bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++failures;
    }
}

void checkNear(const double actual,
               const double expected,
               const double tolerance,
               const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::fprintf(stderr,
                     "FAIL: %s (actual %.17g expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

template<typename Callback>
void expectInvalid(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected, message);
}

template<typename Callback>
void expectLimited(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::length_error&) {
        rejected = true;
    }
    check(rejected, message);
}

Vec3 rotate(const Vec3& point, const fluid::GridFaceAxis axis) {
    if (axis == fluid::GridFaceAxis::X) return point;
    if (axis == fluid::GridFaceAxis::Y) {
        return {point.z, point.x, point.y};
    }
    return {point.y, point.z, point.x};
}

std::array<Vec2, 3> chart(const std::array<Vec3, 4>& positions,
                         const std::array<std::size_t, 3>& vertices) {
    const auto& first = positions[vertices[0]];
    const auto& second = positions[vertices[1]];
    const auto& third = positions[vertices[2]];
    const Vec3 edge{second.x - first.x,
                    second.y - first.y,
                    second.z - first.z};
    const Vec3 diagonal{third.x - first.x,
                        third.y - first.y,
                        third.z - first.z};
    const double edgeLength = std::hypot(edge.x, edge.y, edge.z);
    const double projected = (edge.x * diagonal.x
                              + edge.y * diagonal.y
                              + edge.z * diagonal.z) / edgeLength;
    const double diagonalSquared = diagonal.x * diagonal.x
        + diagonal.y * diagonal.y + diagonal.z * diagonal.z;
    return {{{0.0, 0.0},
             {edgeLength, 0.0},
             {projected, std::sqrt(std::max(
                 0.0, diagonalSquared - projected * projected))}}};
}

Scene openingTetrahedron(const fluid::GridFaceAxis axis) {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-capped-face-partition";
    scene.metadata.exporterVersion =
        "scene-fluid-capped-face-partition-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    std::array<Vec3, 4> positions{{
        rotate({1.0, 1.7, 1.3}, axis),
        rotate({1.5, 1.2, 1.2}, axis),
        rotate({2.5, 1.8, 1.2}, axis),
        rotate({2.5, 1.2, 1.8}, axis),
    }};
    for (std::size_t index = 0; index < positions.size(); ++index) {
        scene.vertices.push_back({10 + index, positions[index]});
    }
    const std::array<std::array<std::size_t, 3>, 3> sides{{
        {{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}},
    }};
    for (std::size_t index = 0; index < sides.size(); ++index) {
        scene.triangles.push_back({
            500 + index,
            {10 + sides[index][0], 10 + sides[index][1],
             10 + sides[index][2]},
            chart(positions, sides[index]),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

fluid::PeriodicCartesianGrid grid() {
    return {{4, 4, 4}, {}, {4.0, 4.0, 4.0}};
}

struct Fixture {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceTransfer transfer;
    SceneFluidSurfaceState state;
    SceneFluidGridEpoch epoch;
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet quadrature;
    SceneFluidOpeningGridPatchSet patches;
    SceneFluidOpeningFaceCrossingSet crossings;

    explicit Fixture(const fluid::GridFaceAxis axis)
        : scene(openingTetrahedron(axis)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          transfer(surface.definition,
                   structureAssembly.mappings, structure),
          state(captureSceneFluidSurfaceState(
              surface.definition, structureAssembly.mappings, structure)),
          epoch(buildSceneFluidGridEpoch(
              surface.definition, state, grid(), transfer)),
          caps(buildSceneFluidOpeningCaps(surface.definition, state)),
          quadrature(buildSceneFluidOpeningQuadrature(
              surface.definition, state, caps)),
          patches(buildSceneFluidOpeningGridPatches(
              surface.definition, state, caps, quadrature, grid())),
          crossings(buildSceneFluidOpeningFaceCrossings(
              surface.definition, state, caps, quadrature, patches,
              grid())) {}

    SceneFluidCappedFacePartitionSet partitions(
        const SceneFluidCappedFacePartitionSettings& settings = {},
        const SceneFluidCappedFacePartitionLimits& limits = {}) const {
        return buildSceneFluidCappedFacePartitions(
            surface.definition, state, grid(), transfer, epoch, caps,
            quadrature, patches, crossings, settings, limits);
    }
};

void checkRegionMoments(
    const SceneFluidCappedFacePartitionSet& partitions,
    const SceneFluidCappedFacePartition& partition) {
    Vec3 summedMoment;
    for (std::size_t offset = 0;
         offset < partition.regionAreaCount; ++offset) {
        const auto& area = partitions.regionAreas[
            partition.firstRegionArea + offset];
        checkNear(area.firstMomentMeters3.x,
                  area.areaSquareMeters * area.centroidMeters.x,
                  3.0e-15,
                  "capped partition x first moment matches centroid");
        checkNear(area.firstMomentMeters3.y,
                  area.areaSquareMeters * area.centroidMeters.y,
                  3.0e-15,
                  "capped partition y first moment matches centroid");
        checkNear(area.firstMomentMeters3.z,
                  area.areaSquareMeters * area.centroidMeters.z,
                  3.0e-15,
                  "capped partition z first moment matches centroid");
        summedMoment.x += area.firstMomentMeters3.x;
        summedMoment.y += area.firstMomentMeters3.y;
        summedMoment.z += area.firstMomentMeters3.z;
    }
    const auto fluidGrid = grid();
    const auto lower = fluidGrid.lowerMeters();
    const auto spacing = fluidGrid.cellSpacingMeters();
    const Vec3 center{
        partition.axis == fluid::GridFaceAxis::X
            ? lower.x + static_cast<double>(partition.i) * spacing.x
            : lower.x + (static_cast<double>(partition.i) + 0.5) * spacing.x,
        partition.axis == fluid::GridFaceAxis::Y
            ? lower.y + static_cast<double>(partition.j) * spacing.y
            : lower.y + (static_cast<double>(partition.j) + 0.5) * spacing.y,
        partition.axis == fluid::GridFaceAxis::Z
            ? lower.z + static_cast<double>(partition.k) * spacing.z
            : lower.z + (static_cast<double>(partition.k) + 0.5) * spacing.z,
    };
    checkNear(summedMoment.x,
              partition.faceAreaSquareMeters * center.x, 3.0e-15,
              "capped partition x first moments close the face");
    checkNear(summedMoment.y,
              partition.faceAreaSquareMeters * center.y, 3.0e-15,
              "capped partition y first moments close the face");
    checkNear(summedMoment.z,
              partition.faceAreaSquareMeters * center.z, 3.0e-15,
              "capped partition z first moments close the face");
}

void testClosedOpeningSection(const fluid::GridFaceAxis axis) {
    Fixture fixture(axis);
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "capped-face opening tetrahedron assembles");
    check(fixture.epoch.facePartitions.partitions.empty()
              && fixture.epoch.facePartitions.unresolvedActiveFaceCount > 0,
          "material-only face partition retains the opening-ended section");
    const auto first = fixture.partitions();
    const auto repeated = fixture.partitions();
    check(first == repeated
              && first.version == sceneFluidCappedFacePartitionVersion
              && first.fingerprint != 0
              && first.touchedFaceCount == 1
              && first.unresolvedTouchedFaceCount == 0
              && first.faces.size() == 1
              && first.faces.front().status
                  == SceneFluidCappedFaceStatus::Resolved
              && first.faces.front().failureSourceStableId
                  == invalidStableId
              && first.faces.front().partitionIndex == 0
              && first.partitions.size() == 1
              && first.materialChainReferences.size() == 1
              && first.openingCrossingReferences.size() == 1,
          "virtual cap closes one deterministic material face section");
    if (first.partitions.empty()) return;
    const auto& partition = first.partitions.front();
    check(partition.axis == axis
              && partition.activeFaceIndex
                  != invalidSceneFluidActiveFaceIndex
              && partition.regionAreaCount == 2
              && partition.faceAreaSquareMeters == 1.0,
          "capped partition retains exact Cartesian face ownership");
    if (axis == fluid::GridFaceAxis::X) {
        check(partition.i == 2 && partition.j == 1 && partition.k == 1,
              "X capped partition uses the expected face coordinate");
    } else if (axis == fluid::GridFaceAxis::Y) {
        check(partition.i == 1 && partition.j == 2 && partition.k == 1,
              "Y capped partition uses the expected face coordinate");
    } else {
        check(partition.i == 1 && partition.j == 1 && partition.k == 2,
              "Z capped partition uses the expected face coordinate");
    }
    std::map<StableId, double> areas;
    for (std::size_t offset = 0;
         offset < partition.regionAreaCount; ++offset) {
        const auto& area = first.regionAreas[
            partition.firstRegionArea + offset];
        areas.emplace(area.regionId, area.areaSquareMeters);
    }
    checkNear(areas.at(2), 0.105, 3.0e-15,
              "capped partition retains analytic cell cross-section area");
    checkNear(areas.at(1), 0.895, 3.0e-15,
              "capped partition retains analytic exterior complement area");
    checkNear(partition.assignedAreaSquareMeters, 1.0, 3.0e-15,
              "capped partition closes assigned area");
    checkNear(partition.areaResidualSquareMeters, 0.0, 3.0e-15,
              "capped partition has zero face-area residual");
    checkRegionMoments(first, partition);
    validateSceneFluidCappedFacePartitions(
        first, fixture.surface.definition, fixture.state, grid(),
        fixture.transfer, fixture.epoch, fixture.caps,
        fixture.quadrature, fixture.patches, fixture.crossings);
}

void testCorruptionSettingsAndLimits() {
    Fixture fixture(fluid::GridFaceAxis::X);
    const auto accepted = fixture.partitions();
    auto corrupt = accepted;
    corrupt.regionAreas.front().areaSquareMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidCappedFacePartitions(
            corrupt, fixture.surface.definition, fixture.state, grid(),
            fixture.transfer, fixture.epoch, fixture.caps,
            fixture.quadrature, fixture.patches, fixture.crossings); },
        "capped-face validation rejects payload corruption");
    corrupt = accepted;
    corrupt.faces.front().status =
        SceneFluidCappedFaceStatus::UnpairedMaterialEndpoint;
    expectInvalid(
        [&] { validateSceneFluidCappedFacePartitions(
            corrupt, fixture.surface.definition, fixture.state, grid(),
            fixture.transfer, fixture.epoch, fixture.caps,
            fixture.quadrature, fixture.patches, fixture.crossings); },
        "capped-face validation rejects status corruption");
    corrupt = accepted;
    corrupt.faces.front().partitionIndex =
        invalidSceneFluidCappedFacePartitionIndex;
    expectInvalid(
        [&] { validateSceneFluidCappedFacePartitions(
            corrupt, fixture.surface.definition, fixture.state, grid(),
            fixture.transfer, fixture.epoch, fixture.caps,
            fixture.quadrature, fixture.patches, fixture.crossings); },
        "capped-face validation rejects ownership corruption");

    SceneFluidCappedFacePartitionSettings invalidSettings;
    invalidSettings.geometryToleranceMeters = 0.0;
    expectInvalid(
        [&] { static_cast<void>(fixture.partitions(invalidSettings)); },
        "capped-face partition rejects zero geometry tolerance");

    SceneFluidCappedFacePartitionLimits limits;
    limits.maximumTouchedFaces = 0;
    expectLimited(
        [&] { static_cast<void>(fixture.partitions({}, limits)); },
        "capped-face partition bounds touched faces");
    limits = {};
    limits.maximumPartitions = 0;
    expectLimited(
        [&] { static_cast<void>(fixture.partitions({}, limits)); },
        "capped-face partition bounds results");
    limits = {};
    limits.maximumReferences = 0;
    expectLimited(
        [&] { static_cast<void>(fixture.partitions({}, limits)); },
        "capped-face partition bounds source and retained references");
    limits = {};
    limits.maximumSegmentPairTests = 0;
    expectLimited(
        [&] { static_cast<void>(fixture.partitions({}, limits)); },
        "capped-face partition bounds planar work");
    limits = {};
    limits.maximumPartitionBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.partitions({}, limits)); },
        "capped-face partition bounds owned storage");
}

} // namespace

int main() {
    try {
        testClosedOpeningSection(fluid::GridFaceAxis::X);
        testClosedOpeningSection(fluid::GridFaceAxis::Y);
        testClosedOpeningSection(fluid::GridFaceAxis::Z);
        testCorruptionSettingsAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d capped-face-partition check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all capped-face-partition checks passed");
    return 0;
}
