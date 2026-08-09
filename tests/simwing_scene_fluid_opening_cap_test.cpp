#include "scene_fluid_opening_cap.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <utility>

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

template<std::size_t VertexCount>
std::array<Vec2, 3> intrinsicChart(
    const std::array<Vec3, VertexCount>& positions,
    const std::array<std::size_t, 3>& vertices) {
    const Vec3& first = positions[vertices[0]];
    const Vec3& second = positions[vertices[1]];
    const Vec3& third = positions[vertices[2]];
    const Vec3 edge{second.x - first.x,
                    second.y - first.y,
                    second.z - first.z};
    const Vec3 diagonal{third.x - first.x,
                        third.y - first.y,
                        third.z - first.z};
    const double edgeLength = std::hypot(edge.x, edge.y, edge.z);
    const double projected = (diagonal.x * edge.x
                              + diagonal.y * edge.y
                              + diagonal.z * edge.z) / edgeLength;
    const double diagonalSquared = diagonal.x * diagonal.x
        + diagonal.y * diagonal.y + diagonal.z * diagonal.z;
    return {{{0.0, 0.0},
             {edgeLength, 0.0},
             {projected, std::sqrt(std::max(
                 0.0, diagonalSquared - projected * projected))}}};
}

Scene openTetraScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-opening-cap";
    scene.metadata.exporterVersion = "scene-fluid-opening-cap-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    const std::array<Vec3, 4> positions{{
        {1.2, 1.5, 1.45},
        {2.8, 1.2, 1.15},
        {2.8, 1.8, 1.15},
        {2.8, 1.5, 1.75},
    }};
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    const std::array<std::array<std::size_t, 3>, 3> faces{{
        {{0, 2, 1}},
        {{0, 1, 3}},
        {{0, 3, 2}},
    }};
    for (std::size_t face = 0; face < faces.size(); ++face) {
        scene.triangles.push_back({
            500 + face,
            {10 + faces[face][0], 10 + faces[face][1],
             10 + faces[face][2]},
            intrinsicChart(positions, faces[face]),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

Scene nonPlanarOpeningScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-fluid-nonplanar-cap";
    scene.metadata.exporterVersion = "scene-fluid-opening-cap-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    scene.vertices = {
        {10, {1.0, 1.0, 1.0}},
        {11, {2.0, 0.5, 0.5}},
        {12, {2.0, 1.5, 0.5}},
        {13, {2.1, 1.5, 1.5}},
        {14, {2.0, 0.5, 1.5}},
    };
    const std::array<std::array<std::size_t, 3>, 4> faces{{
        {{0, 2, 1}}, {{0, 3, 2}}, {{0, 4, 3}}, {{0, 1, 4}},
    }};
    for (std::size_t face = 0; face < faces.size(); ++face) {
        const auto& ids = faces[face];
        const Vec3& first = scene.vertices[ids[0]].positionMeters;
        const Vec3& second = scene.vertices[ids[1]].positionMeters;
        const Vec3& third = scene.vertices[ids[2]].positionMeters;
        const Vec3 edge{second.x - first.x, second.y - first.y,
                        second.z - first.z};
        const Vec3 diagonal{third.x - first.x, third.y - first.y,
                            third.z - first.z};
        const double edgeLength = std::hypot(edge.x, edge.y, edge.z);
        const double projected = (edge.x * diagonal.x
                                  + edge.y * diagonal.y
                                  + edge.z * diagonal.z) / edgeLength;
        const double diagonalSquared = diagonal.x * diagonal.x
            + diagonal.y * diagonal.y + diagonal.z * diagonal.z;
        scene.triangles.push_back({
            500 + face,
            {10 + ids[0], 10 + ids[1], 10 + ids[2]},
            {{{0.0, 0.0}, {edgeLength, 0.0},
              {projected, std::sqrt(std::max(
                  0.0, diagonalSquared - projected * projected))}}},
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13, 14}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

Scene concaveOpeningScene() {
    Scene scene = nonPlanarOpeningScene();
    scene.metadata.designChecksum = "sha256:scene-fluid-concave-cap";
    scene.vertices[3].positionMeters = {2.0, 1.0, 0.8};
    std::array<Vec3, 5> positions;
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        positions[vertex] = scene.vertices[vertex].positionMeters;
    }
    const std::array<std::array<std::size_t, 3>, 4> faces{{
        {{0, 2, 1}}, {{0, 3, 2}}, {{0, 4, 3}}, {{0, 1, 4}},
    }};
    for (std::size_t face = 0; face < faces.size(); ++face) {
        scene.triangles[face].materialCoordinates =
            intrinsicChart(positions, faces[face]);
    }
    return scene;
}

Scene selfIntersectingOpeningScene() {
    Scene scene;
    scene.metadata.designChecksum =
        "sha256:scene-fluid-self-intersecting-cap";
    scene.metadata.exporterVersion = "scene-fluid-opening-cap-test/1";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900.0, 650.0, 220.0, 0.015,
         0.041, 0.02, 0.0125, 2.5e-12},
    };
    constexpr double pi = 3.14159265358979323846;
    std::array<Vec3, 6> positions;
    positions[0] = {1.0, 1.0, 1.0};
    const std::array<std::size_t, 5> starOrder{{0, 2, 4, 1, 3}};
    for (std::size_t index = 0; index < starOrder.size(); ++index) {
        const double angle = 0.5 * pi
            + 2.0 * pi * static_cast<double>(starOrder[index]) / 5.0;
        positions[index + 1] = {
            2.0,
            1.0 + 0.4 * std::cos(angle),
            1.0 + 0.4 * std::sin(angle),
        };
    }
    for (std::size_t vertex = 0; vertex < positions.size(); ++vertex) {
        scene.vertices.push_back({10 + vertex, positions[vertex]});
    }
    for (std::size_t edge = 0; edge < starOrder.size(); ++edge) {
        const std::array<std::size_t, 3> face{{
            0, 1 + (edge + 1) % starOrder.size(), 1 + edge}};
        scene.triangles.push_back({
            500 + edge,
            {10 + face[0], 10 + face[1], 10 + face[2]},
            intrinsicChart(positions, face),
            2, 1, 100, 900, SurfaceRole::Skin,
        });
    }
    scene.openings = {
        {700, {11, 12, 13, 14, 15}, 2, 1, OpeningRole::Intake},
    };
    return scene;
}

struct Fixture {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly structureAssembly;
    Structure structure;

    explicit Fixture(Scene source)
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          structureAssembly(assembleSceneStructure(scene)),
          structure(structureAssembly.definition) {}

    SceneFluidSurfaceState state() const {
        return captureSceneFluidSurfaceState(
            surface.definition, structureAssembly.mappings, structure);
    }
};

void testPlanarCapAndAcceptedMotion() {
    Fixture fixture(openTetraScene());
    check(fixture.surface.ok() && fixture.structureAssembly.ok(),
          "opening-cap tetrahedron assembles");
    const auto state = fixture.state();
    const auto first = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);
    const auto repeated = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);
    check(first == repeated
              && first.version == sceneFluidOpeningCapVersion
              && first.fingerprint != 0
              && first.separatingBoundaryEdgeCount == 3
              && first.caps.size() == 1
              && first.triangles.size() == 1,
          "planar authored opening produces one deterministic virtual cap");
    const auto& cap = first.caps.front();
    const auto& triangle = first.triangles.front();
    check(cap.openingId == 700
              && cap.negativeSideRegionIndex == 1
              && cap.positiveSideRegionIndex == 0
              && cap.role == OpeningRole::Intake
              && triangle.vertexIndices
                  == std::array<std::size_t, 3>({1, 2, 3})
              && cap.unitNormalNegativeToPositive.x > 0.99,
          "cap winding closes the surface from cell to Outside");
    checkNear(cap.areaSquareMeters, 0.18, 1.0e-14,
              "triangular mouth cap has its analytic area");
    validateSceneFluidOpeningCaps(
        first, fixture.surface.definition, state);

    auto reversedScene = openTetraScene();
    std::ranges::reverse(reversedScene.openings.front().orderedVertexIds);
    Fixture reversed(std::move(reversedScene));
    const auto reversedState = reversed.state();
    const auto reversedCaps = buildSceneFluidOpeningCaps(
        reversed.surface.definition, reversedState);
    check(reversedCaps.triangles.front().vertexIndices
              == triangle.vertexIndices
              && reversedCaps.caps.front().unitNormalNegativeToPositive.x
                  == cap.unitNormalNegativeToPositive.x
              && reversedCaps.caps.front().areaSquareMeters
                  == cap.areaSquareMeters,
          "surface boundary fixes cap winding independently of loop direction");

    for (std::size_t node = 0;
         node < fixture.structure.definition().nodes.size(); ++node) {
        fixture.structure.addExternalForce(
            node,
            {fixture.structure.definition().nodes[node].massKg,
             0.0, 0.0});
    }
    StructureStepSettings settings;
    settings.timeStepSeconds = 0.1;
    settings.substeps = 1;
    settings.constraintIterations = 4;
    settings.gravityMetersPerSecondSquared = {};
    settings.velocityDampingPerSecond = 0.0;
    const auto diagnostics = fixture.structure.step(settings);
    const auto movedState = fixture.state();
    const auto moved = buildSceneFluidOpeningCaps(
        fixture.surface.definition, movedState);
    check(diagnostics.finite
              && moved.fingerprint != first.fingerprint
              && moved.surfaceStateFingerprint != first.surfaceStateFingerprint,
          "accepted Structure motion produces a distinct opening-cap epoch");
    checkNear(moved.caps.front().areaSquareMeters,
              cap.areaSquareMeters, 1.0e-14,
              "rigid opening translation preserves cap area");
    checkNear(moved.caps.front().centroidMeters.x,
              cap.centroidMeters.x + 0.01, 1.0e-13,
              "virtual cap follows accepted opening vertices");
}

void testRejectionCorruptionAndLimits() {
    Fixture fixture(openTetraScene());
    const auto state = fixture.state();
    const auto accepted = buildSceneFluidOpeningCaps(
        fixture.surface.definition, state);

    auto missingOpeningScene = openTetraScene();
    missingOpeningScene.openings.clear();
    Fixture missing(std::move(missingOpeningScene));
    const auto missingState = missing.state();
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            missing.surface.definition, missingState)); },
        "unauthored separating boundary remains rejected");

    Fixture nonPlanar(nonPlanarOpeningScene());
    check(nonPlanar.surface.ok() && nonPlanar.structureAssembly.ok(),
          "nonplanar opening fixture assembles before cap rejection");
    const auto nonPlanarState = nonPlanar.state();
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            nonPlanar.surface.definition, nonPlanarState)); },
        "nonplanar boundary without authored cap geometry is rejected");

    auto authoredNonPlanarScene = nonPlanarOpeningScene();
    authoredNonPlanarScene.openings.front().capTriangleVertexIds = {
        {{11, 12, 13}}, {{11, 13, 14}},
    };
    Fixture authoredNonPlanar(std::move(authoredNonPlanarScene));
    check(authoredNonPlanar.surface.ok()
              && authoredNonPlanar.structureAssembly.ok(),
          "authored nonplanar opening fixture assembles");
    if (authoredNonPlanar.surface.ok()
        && authoredNonPlanar.structureAssembly.ok()) {
        const auto authoredState = authoredNonPlanar.state();
        const auto authoredCaps = buildSceneFluidOpeningCaps(
            authoredNonPlanar.surface.definition, authoredState);
        check(authoredCaps.caps.size() == 1
                  && authoredCaps.triangles.size() == 2
                  && authoredCaps.triangles[0].vertexIndices
                      == std::array<std::size_t, 3>({1, 2, 3})
                  && authoredCaps.triangles[1].vertexIndices
                      == std::array<std::size_t, 3>({1, 3, 4}),
              "authored nonplanar disk retains its exact triangle identity");
        checkNear(
            authoredCaps.caps.front().areaSquareMeters,
            std::sqrt(1.01), 1.0e-14,
            "authored nonplanar cap retains the analytic faceted area");
        check(std::abs(
                  authoredCaps.triangles[0]
                      .unitNormalNegativeToPositive.z)
                      > 0.09
                  && std::abs(
                      authoredCaps.triangles[1]
                          .unitNormalNegativeToPositive.y)
                      > 0.09
                  && authoredCaps.triangles[0]
                         .unitNormalNegativeToPositive
                         .z
                      != authoredCaps.triangles[1]
                             .unitNormalNegativeToPositive
                             .z,
              "authored nonplanar cap publishes distinct triangle normals");
        validateSceneFluidOpeningCaps(
            authoredCaps, authoredNonPlanar.surface.definition,
            authoredState);

        for (std::size_t node = 0;
             node < authoredNonPlanar.structure.definition().nodes.size();
             ++node) {
            authoredNonPlanar.structure.addExternalForce(
                node,
                {authoredNonPlanar.structure.definition().nodes[node].massKg,
                 0.0, 0.0});
        }
        StructureStepSettings authoredMotion;
        authoredMotion.timeStepSeconds = 0.1;
        authoredMotion.substeps = 1;
        authoredMotion.constraintIterations = 4;
        authoredMotion.gravityMetersPerSecondSquared = {};
        authoredMotion.velocityDampingPerSecond = 0.0;
        const auto authoredMotionDiagnostics =
            authoredNonPlanar.structure.step(authoredMotion);
        const auto movedAuthoredCaps = buildSceneFluidOpeningCaps(
            authoredNonPlanar.surface.definition,
            authoredNonPlanar.state());
        check(authoredMotionDiagnostics.finite
                  && movedAuthoredCaps.triangles[0].vertexIndices
                      == authoredCaps.triangles[0].vertexIndices
                  && movedAuthoredCaps.triangles[1].vertexIndices
                      == authoredCaps.triangles[1].vertexIndices,
              "accepted motion retains authored nonplanar cap identity");
        checkNear(
            movedAuthoredCaps.caps.front().areaSquareMeters,
            authoredCaps.caps.front().areaSquareMeters, 1.0e-14,
            "rigid authored nonplanar motion preserves faceted area");
        checkNear(
            movedAuthoredCaps.caps.front().centroidMeters.x,
            authoredCaps.caps.front().centroidMeters.x + 0.01,
            1.0e-13,
            "authored nonplanar cap follows accepted boundary motion");

        auto reversedAuthoredScene = nonPlanarOpeningScene();
        reversedAuthoredScene.openings.front().capTriangleVertexIds = {
            {{11, 13, 12}}, {{11, 14, 13}},
        };
        std::ranges::reverse(
            reversedAuthoredScene.openings.front().orderedVertexIds);
        Fixture reversedAuthored(std::move(reversedAuthoredScene));
        const auto reversedAuthoredCaps = buildSceneFluidOpeningCaps(
            reversedAuthored.surface.definition,
            reversedAuthored.state());
        check(reversedAuthoredCaps.triangles[0].vertexIndices
                      == authoredCaps.triangles[0].vertexIndices
                  && reversedAuthoredCaps.triangles[1].vertexIndices
                      == authoredCaps.triangles[1].vertexIndices
                  && reversedAuthoredCaps.caps.front().areaSquareMeters
                      == authoredCaps.caps.front().areaSquareMeters,
              "fabric winding canonicalizes a reversed authored nonplanar disk");
    }

    Fixture concave(concaveOpeningScene());
    check(concave.surface.ok() && concave.structureAssembly.ok(),
          "concave opening fixture assembles");
    if (concave.surface.ok() && concave.structureAssembly.ok()) {
        const auto concaveState = concave.state();
        const auto concaveCaps = buildSceneFluidOpeningCaps(
            concave.surface.definition, concaveState);
        const auto repeatedConcaveCaps = buildSceneFluidOpeningCaps(
            concave.surface.definition, concaveState);
        check(concaveCaps == repeatedConcaveCaps
                  && concaveCaps.caps.size() == 1
                  && concaveCaps.triangles.size() == 2
                  && concaveCaps.triangles[0].vertexIndices
                      == std::array<std::size_t, 3>({1, 2, 3})
                  && concaveCaps.triangles[1].vertexIndices
                      == std::array<std::size_t, 3>({1, 3, 4}),
              "planar concave opening receives deterministic reference-geometry ear clipping");
        checkNear(
            concaveCaps.caps.front().areaSquareMeters, 0.4, 1.0e-14,
            "concave virtual cap has its analytic polygon area");
        checkNear(
            concaveCaps.caps.front().centroidMeters.y,
            19.0 / 24.0, 1.0e-14,
            "concave virtual cap has its analytic area centroid y");
        checkNear(
            concaveCaps.caps.front().centroidMeters.z,
            97.0 / 120.0, 1.0e-14,
            "concave virtual cap has its analytic area centroid z");
        validateSceneFluidOpeningCaps(
            concaveCaps, concave.surface.definition, concaveState);

        auto reversedConcaveScene = concaveOpeningScene();
        std::ranges::reverse(
            reversedConcaveScene.openings.front().orderedVertexIds);
        Fixture reversedConcave(std::move(reversedConcaveScene));
        const auto reversedConcaveCaps = buildSceneFluidOpeningCaps(
            reversedConcave.surface.definition, reversedConcave.state());
        check(reversedConcaveCaps.triangles[0].vertexIndices
                      == concaveCaps.triangles[0].vertexIndices
                  && reversedConcaveCaps.triangles[1].vertexIndices
                      == concaveCaps.triangles[1].vertexIndices,
              "surface boundary fixes concave ear clipping independently of authored loop direction");

        for (std::size_t node = 0;
             node < concave.structure.definition().nodes.size(); ++node) {
            concave.structure.addExternalForce(
                node,
                {concave.structure.definition().nodes[node].massKg,
                 0.0, 0.0});
        }
        StructureStepSettings motion;
        motion.timeStepSeconds = 0.1;
        motion.substeps = 1;
        motion.constraintIterations = 4;
        motion.gravityMetersPerSecondSquared = {};
        motion.velocityDampingPerSecond = 0.0;
        const auto motionDiagnostics = concave.structure.step(motion);
        const auto movedConcaveCaps = buildSceneFluidOpeningCaps(
            concave.surface.definition, concave.state());
        check(motionDiagnostics.finite
                  && movedConcaveCaps.triangles[0].vertexIndices
                      == concaveCaps.triangles[0].vertexIndices
                  && movedConcaveCaps.triangles[1].vertexIndices
                      == concaveCaps.triangles[1].vertexIndices,
              "accepted motion retains concave reference triangulation identity");
        checkNear(
            movedConcaveCaps.caps.front().areaSquareMeters,
            concaveCaps.caps.front().areaSquareMeters, 1.0e-14,
            "rigid concave opening translation preserves cap area");
        checkNear(
            movedConcaveCaps.caps.front().centroidMeters.x,
            concaveCaps.caps.front().centroidMeters.x + 0.01, 1.0e-13,
            "concave virtual cap follows accepted opening motion");
    }

    Fixture selfIntersecting(selfIntersectingOpeningScene());
    check(selfIntersecting.surface.ok()
              && selfIntersecting.structureAssembly.ok(),
          "self-intersecting opening fixture assembles topologically before cap rejection");
    if (selfIntersecting.surface.ok()
        && selfIntersecting.structureAssembly.ok()) {
        const auto selfIntersectingState = selfIntersecting.state();
        expectInvalid(
            [&] { static_cast<void>(buildSceneFluidOpeningCaps(
                selfIntersecting.surface.definition,
                selfIntersectingState)); },
            "self-intersecting planar opening is rejected before triangulation");
    }

    auto duplicateOpeningScene = openTetraScene();
    duplicateOpeningScene.openings.push_back(
        duplicateOpeningScene.openings.front());
    duplicateOpeningScene.openings.back().id = 701;
    Fixture duplicate(std::move(duplicateOpeningScene));
    check(duplicate.surface.ok() && duplicate.structureAssembly.ok(),
          "duplicate opening-claim fixture assembles before cap rejection");
    if (duplicate.surface.ok() && duplicate.structureAssembly.ok()) {
        const auto duplicateState = duplicate.state();
        expectInvalid(
            [&] { static_cast<void>(buildSceneFluidOpeningCaps(
                duplicate.surface.definition, duplicateState)); },
            "one separating boundary edge cannot be claimed twice");
    }

    auto corrupt = accepted;
    corrupt.triangles.front().areaSquareMeters += 0.01;
    expectInvalid(
        [&] { validateSceneFluidOpeningCaps(
            corrupt, fixture.surface.definition, state); },
        "opening-cap validation rejects payload corruption");

    SceneFluidOpeningCapSettings invalidSettings;
    invalidSettings.minimumTriangleAreaSquareMeters = 0.0;
    expectInvalid(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            fixture.surface.definition, state, invalidSettings)); },
        "opening-cap construction rejects invalid tolerances");

    SceneFluidOpeningCapLimits limits;
    limits.maximumCaps = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            fixture.surface.definition, state, {}, limits)); },
        "opening-cap count is bounded");
    limits = {};
    limits.maximumBoundaryEdges = 2;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            fixture.surface.definition, state, {}, limits)); },
        "opening-cap surface-edge work is bounded");
    limits = {};
    limits.maximumCapTriangles = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            fixture.surface.definition, state, {}, limits)); },
        "opening-cap triangle count is bounded");
    limits = {};
    limits.maximumCapBytes = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            fixture.surface.definition, state, {}, limits)); },
        "opening-cap result bytes are bounded");

    Fixture concaveLimits(concaveOpeningScene());
    const auto concaveLimitState = concaveLimits.state();
    limits = {};
    limits.maximumBoundaryIntersectionTests = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            concaveLimits.surface.definition, concaveLimitState,
            {}, limits)); },
        "opening-cap simple-boundary work is bounded");
    limits = {};
    limits.maximumTriangulationPointTests = 0;
    expectLimited(
        [&] { static_cast<void>(buildSceneFluidOpeningCaps(
            concaveLimits.surface.definition, concaveLimitState,
            {}, limits)); },
        "opening-cap ear-clipping work is bounded");
}

} // namespace

int main() {
    try {
        testPlanarCapAndAcceptedMotion();
        testRejectionCorruptionAndLimits();
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "unexpected exception: %s\n", exception.what());
        return 1;
    }
    if (failures != 0) {
        std::fprintf(stderr, "%d opening-cap check(s) failed\n", failures);
        return 1;
    }
    std::puts("all scene fluid opening-cap checks passed");
    return 0;
}
