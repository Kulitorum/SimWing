#include "scene_fluid_opening_face_crossing.h"
#include "scene_structure.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
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

enum class MouthGeometry {
    Transverse,
    FaceAligned,
    GridEdgeAligned,
};

Scene openingTetrahedron(const MouthGeometry geometry) {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-opening-face-crossing";
    scene.metadata.exporterVersion =
        "scene-fluid-opening-face-crossing-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    std::array<Vec3, 4> positions{{
        {1.0, 1.7, 1.3},
        {1.5, 1.2, 1.2},
        {2.5, 1.8, 1.2},
        {2.5, 1.2, 1.8},
    }};
    if (geometry == MouthGeometry::FaceAligned) {
        positions[1] = {2.0, 1.2, 1.2};
        positions[2] = {2.0, 1.8, 1.2};
        positions[3] = {2.0, 1.5, 1.8};
    } else if (geometry == MouthGeometry::GridEdgeAligned) {
        positions[1].y = 0.8;
        positions[2].y = 1.2;
        positions[3].y = 1.2;
    }
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

SceneFluidSurfaceAssembly checkedSurface(const Scene& scene) {
    auto result = assembleSceneFluidSurface(scene);
    if (!result.ok()) {
        for (const auto& diagnostic : result.diagnostics) {
            std::fprintf(stderr, "fluid surface: %s\n",
                         diagnostic.message.c_str());
        }
    }
    return result;
}

struct Fixture {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;
    SceneFluidSurfaceState state;
    SceneFluidOpeningCapSet caps;
    SceneFluidOpeningQuadratureSet quadrature;
    SceneFluidOpeningGridPatchSet patches;

    explicit Fixture(const MouthGeometry geometry)
        : scene(openingTetrahedron(geometry)),
          surface(checkedSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition),
          state(captureSceneFluidSurfaceState(
              surface.definition, structureAssembly.mappings, structure)),
          caps(buildSceneFluidOpeningCaps(surface.definition, state)),
          quadrature(buildSceneFluidOpeningQuadrature(
              surface.definition, state, caps)),
          patches(buildSceneFluidOpeningGridPatches(
              surface.definition, state, caps, quadrature, grid())) {}

    SceneFluidOpeningFaceCrossingSet crossings(
        const SceneFluidOpeningFaceCrossingLimits& limits = {}) const {
        return buildSceneFluidOpeningFaceCrossings(
            surface.definition, state, caps, quadrature, patches,
            grid(), limits);
    }
};

void testTransverseOpeningCrossing() {
    Fixture fixture(MouthGeometry::Transverse);
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "opening-face transverse tetrahedron assembles");
    const auto first = fixture.crossings();
    const auto repeated = fixture.crossings();
    check(first == repeated
              && first.version == sceneFluidOpeningFaceCrossingVersion
              && first.fingerprint != 0
              && first.candidateSegmentCount == 2
              && first.unpairedContactSegmentCount == 0
              && first.faceOwnedPatchCount == 0
              && first.crossings.size() == 1,
          "transverse cap crossing is paired exactly and deterministically");
    if (first.crossings.empty()) return;
    const auto& crossing = first.crossings.front();
    check(crossing.axis == fluid::GridFaceAxis::X
              && crossing.i == 2 && crossing.j == 1 && crossing.k == 1
              && crossing.openingId == 700
              && crossing.negativeSideRegionId == 2
              && crossing.positiveSideRegionId == 1
              && crossing.role == OpeningRole::Intake
              && crossing.sourcePointStableId
                  == fixture.quadrature.points.front().stableId,
          "cap crossing retains exact face, source, and authored regions");
    const auto& lower = fixture.patches.patches[
        crossing.lowerCellPatchIndex];
    const auto& upper = fixture.patches.patches[
        crossing.upperCellPatchIndex];
    check(lower.ownerKind == SceneFluidOpeningPatchOwnerKind::Cell
              && upper.ownerKind == SceneFluidOpeningPatchOwnerKind::Cell
              && lower.cell.i == 1 && upper.cell.i == 2,
          "cap crossing retains both adjacent clipped-cell owners");
    checkNear(crossing.first.positionMeters.x, 2.0, 0.0,
              "first crossing point lies exactly on the X face");
    checkNear(crossing.second.positionMeters.x, 2.0, 0.0,
              "second crossing point lies exactly on the X face");
    checkNear(crossing.lengthMeters, std::sqrt(0.18), 2.0e-15,
              "cap crossing retains analytic segment length");
    checkNear(first.crossingLengthMeters, crossing.lengthMeters, 0.0,
              "cap crossing aggregate length closes");
    checkNear(crossing.midpointMeters.y, 1.35, 2.0e-15,
              "cap crossing retains analytic midpoint y");
    checkNear(crossing.midpointMeters.z, 1.35, 2.0e-15,
              "cap crossing retains analytic midpoint z");
    const Vec3 capNormal = fixture.quadrature.points.front()
        .unitNormalNegativeToPositive;
    const double projectedMagnitude = std::hypot(capNormal.y, capNormal.z);
    checkNear(crossing.negativeToPositiveDirectionInFace.x, 0.0, 0.0,
              "cap crossing direction remains in its Cartesian face");
    checkNear(crossing.negativeToPositiveDirectionInFace.y,
              capNormal.y / projectedMagnitude, 2.0e-16,
              "cap crossing direction preserves projected normal y");
    checkNear(crossing.negativeToPositiveDirectionInFace.z,
              capNormal.z / projectedMagnitude, 2.0e-16,
              "cap crossing direction preserves projected normal z");
    validateSceneFluidOpeningFaceCrossings(
        first, fixture.surface.definition, fixture.state, fixture.caps,
        fixture.quadrature, fixture.patches, grid());
}

void testFaceAreaAndGridEdgeOwnership() {
    Fixture aligned(MouthGeometry::FaceAligned);
    const auto faceArea = aligned.crossings();
    check(faceArea.crossings.empty()
              && faceArea.candidateSegmentCount == 0
              && faceArea.faceOwnedPatchCount == 1,
          "face-owned cap area is counted but never converted into a crossing");

    Fixture gridEdge(MouthGeometry::GridEdgeAligned);
    bool rejected = false;
    bool retainedAsContact = false;
    try {
        const auto contact = gridEdge.crossings();
        retainedAsContact = contact.crossings.empty()
            && contact.unpairedContactSegmentCount > 0;
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check(rejected || retainedAsContact,
          "grid-edge cap segment rejects or remains explicit contact");
}

void testCorruptionAndLimits() {
    Fixture fixture(MouthGeometry::Transverse);
    const auto accepted = fixture.crossings();
    auto corrupt = accepted;
    corrupt.crossings.front().lengthMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidOpeningFaceCrossings(
            corrupt, fixture.surface.definition, fixture.state,
            fixture.caps, fixture.quadrature, fixture.patches, grid()); },
        "opening-face-crossing validation rejects payload corruption");

    SceneFluidOpeningFaceCrossingLimits limits;
    limits.maximumCandidateSegments = accepted.candidateSegmentCount - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.crossings(limits)); },
        "opening-face crossing bounds candidate work");
    limits = {};
    limits.maximumCrossings = 0;
    expectLimited(
        [&] { static_cast<void>(fixture.crossings(limits)); },
        "opening-face crossing bounds result count");
    limits = {};
    limits.maximumCrossingBytes = accepted.ownedStorageBytes - 1;
    expectLimited(
        [&] { static_cast<void>(fixture.crossings(limits)); },
        "opening-face crossing bounds owned storage");
}

} // namespace

int main() {
    try {
        testTransverseOpeningCrossing();
        testFaceAreaAndGridEdgeOwnership();
        testCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d opening-face-crossing check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all opening-face-crossing checks passed");
    return 0;
}
