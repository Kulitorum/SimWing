#include "fluid/scene_surface_face_partition.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <utility>

namespace {
using namespace simwing::fsi;
using namespace simwing::fsi::fluid;
int failures = 0;
void check(bool value, const char* message) {
    if (!value) { std::fprintf(stderr, "FAIL: %s\n", message); ++failures; }
}
void checkNear(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual-expected)>tolerance) {
        std::fprintf(stderr, "FAIL: %s (actual %.17g expected %.17g)\n",
                     message, actual, expected); ++failures;
    }
}
template<class F> void expectInvalid(F&& f, const char* message) {
    bool rejected=false; try { f(); } catch (const std::invalid_argument&) { rejected=true; }
    check(rejected,message);
}
template<class F> void expectLimited(F&& f, const char* message) {
    bool rejected=false; try { f(); } catch (const std::length_error&) { rejected=true; }
    check(rejected,message);
}

void addTetra(Scene& scene, StableId base, StableId triangleBase,
              StableId regionInside, StableId regionOutside,
              double leftX, double rightX, double yRadius, double zRadius,
              StableId sheet) {
    scene.vertices.insert(scene.vertices.end(), {
        {base+0, {leftX, 1.5, 1.45}},
        {base+1, {rightX, 1.5-yRadius, 1.45-zRadius}},
        {base+2, {rightX, 1.5+yRadius, 1.45-zRadius}},
        {base+3, {rightX, 1.5, 1.45+zRadius}},
    });
    const std::array<Vec2,3> uv{{{0,0},{1,0},{0,1}}};
    scene.triangles.insert(scene.triangles.end(), {
        {triangleBase+0,{base+0,base+2,base+1},uv,regionInside,regionOutside,100,sheet,SurfaceRole::Skin},
        {triangleBase+1,{base+0,base+1,base+3},uv,regionInside,regionOutside,100,sheet,SurfaceRole::Skin},
        {triangleBase+2,{base+0,base+3,base+2},uv,regionInside,regionOutside,100,sheet,SurfaceRole::Skin},
        {triangleBase+3,{base+1,base+2,base+3},uv,regionInside,regionOutside,100,sheet,SurfaceRole::Skin},
    });
}

Scene nestedScene() {
    Scene scene;
    scene.metadata.designChecksum="sha256:scene-face-partition";
    scene.metadata.exporterVersion="scene-face-partition-test/1";
    scene.regions={{1,RegionKind::Outside,"outside"},{2,RegionKind::Cell,"outer"},
                   {3,RegionKind::Cell,"inner"}};
    scene.fabricMaterials={{100,"fabric",900,650,220,0.015,0.041,0.02,0.0125,2.5e-12}};
    addTetra(scene,10,500,2,1,1.2,2.8,0.3,0.3,900);
    addTetra(scene,20,600,3,2,1.6,2.4,0.1,0.1,901);
    return scene;
}

Scene boundarySplitScene(
    const GridFaceAxis axis = GridFaceAxis::X) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-face-boundary-split";
    scene.metadata.exporterVersion = "scene-face-partition-test/2";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900, 650, 220, 0.015, 0.041, 0.02,
         0.0125, 2.5e-12},
    };
    scene.vertices = {
        {10, {1.2, 0.8, 1.15}},
        {11, {2.8, 0.8, 1.15}},
        {12, {2.8, 2.2, 1.65}},
        {13, {1.2, 2.2, 1.65}},
    };
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {1.6, 0.0}, {1.6, 1.4}}},
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {10, 12, 13},
         {{{0.0, 0.0}, {1.6, 1.4}, {0.0, 1.4}}},
         1, 2, 100, 900, SurfaceRole::Skin},
    };
    if (axis == GridFaceAxis::Y) {
        scene.metadata.designChecksum += "-y";
        for (auto& vertex : scene.vertices) {
            const Vec3 point = vertex.positionMeters;
            vertex.positionMeters = {point.z, point.x, point.y};
        }
    } else if (axis == GridFaceAxis::Z) {
        scene.metadata.designChecksum += "-z";
        for (auto& vertex : scene.vertices) {
            const Vec3 point = vertex.positionMeters;
            vertex.positionMeters = {point.y, point.z, point.x};
        }
    }
    return scene;
}

Scene boundaryJunctionScene(
    const GridFaceAxis axis = GridFaceAxis::X) {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-face-boundary-junction";
    scene.metadata.exporterVersion = "scene-face-partition-test/3";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell-a"},
        {3, RegionKind::Cell, "cell-b"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900, 650, 220, 0.015, 0.041, 0.02,
         0.0125, 2.5e-12},
    };
    scene.vertices = {
        {10, {1.5, 1.5, 1.5}},
        {11, {2.5, 1.5, 1.5}},
        {12, {2.5, 2.5, 1.5}},
        {13, {2.5, 1.5, 2.5}},
        {14, {2.5, 0.5, 1.5}},
    };
    const std::array<Vec2, 3> chart{{
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}}};
    scene.triangles = {
        {500, {10, 11, 12}, chart,
         1, 2, 100, 900, SurfaceRole::Skin},
        {501, {10, 11, 13}, chart,
         2, 3, 100, 901, SurfaceRole::Rib},
        {502, {10, 11, 14}, chart,
         3, 1, 100, 902, SurfaceRole::Skin},
    };
    if (axis == GridFaceAxis::Y) {
        scene.metadata.designChecksum += "-y";
        for (auto& vertex : scene.vertices) {
            const Vec3 point = vertex.positionMeters;
            vertex.positionMeters = {point.z, point.x, point.y};
        }
    } else if (axis == GridFaceAxis::Z) {
        scene.metadata.designChecksum += "-z";
        for (auto& vertex : scene.vertices) {
            const Vec3 point = vertex.positionMeters;
            vertex.positionMeters = {point.y, point.z, point.x};
        }
    }
    return scene;
}

Scene boundaryJunctionWithSameRegionSheetScene() {
    Scene scene = boundaryJunctionScene();
    scene.metadata.designChecksum += "-same-region-sheet";
    scene.vertices.insert(scene.vertices.end(), {
        {20, {1.5, 1.2, 1.7}},
        {21, {2.5, 1.2, 1.9}},
        {22, {2.5, 1.4, 1.9}},
    });
    scene.triangles.push_back({
        503, {20, 21, 22},
        {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.2}}},
        3, 3, 100, 903, SurfaceRole::Diagonal});
    return scene;
}

Scene sameRegionSheetScene() {
    Scene scene;
    scene.metadata.designChecksum = "sha256:scene-face-same-region-sheet";
    scene.metadata.exporterVersion = "scene-face-partition-test/4";
    scene.regions = {
        {1, RegionKind::Outside, "outside"},
        {2, RegionKind::Cell, "cell"},
    };
    scene.fabricMaterials = {
        {100, "fabric", 900, 650, 220, 0.015, 0.041, 0.02,
         0.0125, 2.5e-12},
    };
    scene.vertices = {
        {10, {1.5, 1.2, 1.7}},
        {11, {2.5, 1.2, 1.9}},
        {12, {2.5, 1.4, 1.9}},
    };
    scene.triangles = {
        {500, {10, 11, 12},
         {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.2}}},
         2, 2, 100, 900, SurfaceRole::Diagonal},
    };
    return scene;
}
PeriodicCartesianGrid grid(){ return {{4,4,4},{},{4,4,4}}; }

struct Pipeline {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly sa;
    Structure structure;
    SceneFluidSurfaceState state;
    SceneFluidGridCandidateSet candidates;
    SceneFluidGridIntersectionSet intersections;
    SceneFluidGridPatchSet patches;
    SceneFluidPatchOwnership ownership;
    SceneFluidFaceCrossingSet crossings;
    SceneFluidFaceTopology topology;
    SceneFluidFaceGraph graph;
    SceneFluidFaceChainSet chains;
    SceneFluidFaceLoopSet loops;

    explicit Pipeline(Scene source = nestedScene())
        : scene(std::move(source)),
          surface(assembleSceneFluidSurface(scene)),
          sa(assembleSceneStructure(scene)),
          structure(sa.definition),
          state(captureSceneFluidSurfaceState(
              surface.definition, sa.mappings, structure)),
          candidates(buildSceneFluidGridCandidates(
              surface.definition, state, grid())),
          intersections(intersectSceneFluidSurfaceWithGrid(
              surface.definition, state, grid(), candidates)),
          patches(clipSceneFluidSurfaceToCells(
              surface.definition, state, grid(), candidates,
              intersections)),
          ownership(ownSceneFluidSurfacePatches(
              surface.definition, state, grid(), candidates,
              intersections, patches)),
          crossings(buildSceneFluidFaceCrossings(
              surface.definition, state, grid(), candidates,
              intersections, patches, ownership)),
          topology(buildSceneFluidFaceTopology(
              surface.definition, state, grid(), candidates,
              intersections, patches, ownership, crossings)),
          graph(buildSceneFluidFaceGraph(
              surface.definition, state, grid(), candidates,
              intersections, patches, ownership, crossings, topology)),
          chains(buildSceneFluidFaceChains(
              surface.definition, state, grid(), candidates,
              intersections, patches, ownership, crossings, topology,
              graph)),
          loops(buildSceneFluidFaceLoops(
              surface.definition, state, grid(), candidates,
              intersections, patches, ownership, crossings, topology,
              graph, chains)) {}
};

SceneFluidFacePartitionSet partitions(const Pipeline& p,
    SceneFluidFacePartitionSettings settings={}, SceneFluidFacePartitionLimits limits={}) {
    return buildSceneFluidFacePartitions(p.surface.definition,p.state,grid(),p.candidates,
        p.intersections,p.patches,p.ownership,p.crossings,p.topology,p.graph,p.chains,p.loops,
        settings,limits);
}

void checkRegionMoments(const SceneFluidFacePartitionSet& partitions,
                        const SceneFluidFacePartition& partition,
                        const SceneFluidActiveFace& face) {
    Vec3 summedMoment;
    for (std::size_t offset = 0;
         offset < partition.regionAreaCount; ++offset) {
        const auto& area = partitions.regionAreas[
            partition.firstRegionArea + offset];
        checkNear(area.firstMomentMeters3.x,
                  area.areaSquareMeters * area.centroidMeters.x,
                  3.0e-15,
                  "face partition: region x first moment matches centroid");
        checkNear(area.firstMomentMeters3.y,
                  area.areaSquareMeters * area.centroidMeters.y,
                  3.0e-15,
                  "face partition: region y first moment matches centroid");
        checkNear(area.firstMomentMeters3.z,
                  area.areaSquareMeters * area.centroidMeters.z,
                  3.0e-15,
                  "face partition: region z first moment matches centroid");
        summedMoment.x += area.firstMomentMeters3.x;
        summedMoment.y += area.firstMomentMeters3.y;
        summedMoment.z += area.firstMomentMeters3.z;
    }
    const auto fluidGrid = grid();
    const auto lower = fluidGrid.lowerMeters();
    const auto spacing = fluidGrid.cellSpacingMeters();
    const Vec3 center{
        face.axis == GridFaceAxis::X
            ? lower.x + static_cast<double>(face.i) * spacing.x
            : lower.x + (static_cast<double>(face.i) + 0.5) * spacing.x,
        face.axis == GridFaceAxis::Y
            ? lower.y + static_cast<double>(face.j) * spacing.y
            : lower.y + (static_cast<double>(face.j) + 0.5) * spacing.y,
        face.axis == GridFaceAxis::Z
            ? lower.z + static_cast<double>(face.k) * spacing.z
            : lower.z + (static_cast<double>(face.k) + 0.5) * spacing.z,
    };
    checkNear(summedMoment.x,
              partition.faceAreaSquareMeters * center.x, 3.0e-15,
              "face partition: region x first moments close the face");
    checkNear(summedMoment.y,
              partition.faceAreaSquareMeters * center.y, 3.0e-15,
              "face partition: region y first moments close the face");
    checkNear(summedMoment.z,
              partition.faceAreaSquareMeters * center.z, 3.0e-15,
              "face partition: region z first moments close the face");
}

void testNestedPartition() {
    Pipeline p;
    auto a=partitions(p), b=partitions(p);
    check(a==b && a.version == sceneFluidFacePartitionVersion
          && a.fingerprint!=0 && a.partitions.size()==1
          && a.loopContainment.size()==2 && a.unresolvedActiveFaceCount==0,
          "face partition: nested loops partition deterministically");
    const auto& part=a.partitions.front();
    check(part.rootExteriorRegionId==1 && part.loopReferenceCount==2
          && part.regionAreaCount==3 && a.segmentPairTestCount==9,
          "face partition: nesting and pair checks remain explicit");
    std::size_t roots=0, children=0;
    for (const auto& c:a.loopContainment) {
        roots += c.parentLoopIndex==noParentFaceLoop;
        children += c.parentLoopIndex!=noParentFaceLoop && c.depth==1;
    }
    check(roots==1 && children==1,"face partition: smallest containing loop is the parent");
    std::map<StableId,double> areas;
    for (std::size_t i=0;i<part.regionAreaCount;++i) {
        const auto& value=a.regionAreas[part.firstRegionArea+i]; areas[value.regionId]=value.areaSquareMeters;
    }
    checkNear(areas[1],0.955,3e-15,"face partition: exterior area closes");
    checkNear(areas[2],0.04,3e-15,"face partition: annular parent region closes");
    checkNear(areas[3],0.005,3e-15,"face partition: nested region area closes");
    checkNear(part.assignedAreaSquareMeters,1.0,3e-15,"face partition: all region areas sum to face area");
    checkNear(part.areaResidualSquareMeters,0.0,3e-15,"face partition: area residual closes");
    checkRegionMoments(a, part, p.topology.activeFaces[part.activeFaceIndex]);
    validateSceneFluidFacePartitions(a,p.surface.definition,p.state,grid(),p.candidates,
        p.intersections,p.patches,p.ownership,p.crossings,p.topology,p.graph,p.chains,p.loops);
}

void testBoundaryOpenChainPartition() {
    Pipeline p(boundarySplitScene());
    const auto first = partitions(p);
    const auto repeated = partitions(p);
    const auto found = std::ranges::find(
        first.partitions,
        SceneFluidFacePartitionKind::BoundaryOpenChain,
        &SceneFluidFacePartition::kind);
    check(first == repeated && found != first.partitions.end()
              && first.partitions.size() == 1
              && first.unresolvedActiveFaceCount + first.partitions.size()
                  == p.topology.activeFaces.size()
              && first.openChainReferences.size() == 1
              && first.segmentPairTestCount > 0,
          "face partition: one boundary chain resolves deterministically");
    if (found == first.partitions.end()) return;
    check(found->rootExteriorRegionId == invalidStableId
              && found->loopReferenceCount == 0
              && found->openChainReferenceCount == 1
              && found->regionAreaCount == 2,
          "face partition: boundary-chain source identity remains explicit");
    const std::size_t chainIndex = first.openChainReferences[
        found->firstOpenChainReference];
    check(chainIndex < p.chains.chains.size()
              && p.chains.chains[chainIndex].kind
                  == SceneFluidFaceChainKind::Open
              && p.chains.chains[chainIndex]
                      .endpointFaceBoundaryMasks[0]
                  != FaceBoundaryNone
              && p.chains.chains[chainIndex]
                      .endpointFaceBoundaryMasks[1]
                  != FaceBoundaryNone,
          "face partition: retained open chain joins two face boundaries");
    std::map<StableId, double> areas;
    for (std::size_t offset = 0;
         offset < found->regionAreaCount; ++offset) {
        const auto& area = first.regionAreas[
            found->firstRegionArea + offset];
        areas.emplace(area.regionId, area.areaSquareMeters);
    }
    checkNear(areas.at(1), 0.4, 3.0e-15,
              "face partition: boundary chain preserves negative-side area");
    checkNear(areas.at(2), 0.6, 3.0e-15,
              "face partition: boundary chain preserves positive-side area");
    checkNear(found->assignedAreaSquareMeters, 1.0, 3.0e-15,
              "face partition: boundary-chain areas close the face");
    checkNear(found->areaResidualSquareMeters, 0.0, 3.0e-15,
              "face partition: boundary-chain residual is exact");
    checkRegionMoments(first, *found,
                       p.topology.activeFaces[found->activeFaceIndex]);
    validateSceneFluidFacePartitions(
        first, p.surface.definition, p.state, grid(), p.candidates,
        p.intersections, p.patches, p.ownership, p.crossings, p.topology,
        p.graph, p.chains, p.loops);

    for (const GridFaceAxis axis : {
             GridFaceAxis::Y, GridFaceAxis::Z}) {
        Pipeline rotated(boundarySplitScene(axis));
        const auto rotatedPartitions = partitions(rotated);
        const auto rotatedFound = std::ranges::find(
            rotatedPartitions.partitions,
            SceneFluidFacePartitionKind::BoundaryOpenChain,
            &SceneFluidFacePartition::kind);
        check(rotatedFound != rotatedPartitions.partitions.end()
                  && rotatedPartitions.partitions.size() == 1,
              "face partition: right-handed Y/Z charts retain boundary closure");
        if (rotatedFound == rotatedPartitions.partitions.end()) continue;
        std::map<StableId, double> rotatedAreas;
        for (std::size_t offset = 0;
             offset < rotatedFound->regionAreaCount; ++offset) {
            const auto& area = rotatedPartitions.regionAreas[
                rotatedFound->firstRegionArea + offset];
            rotatedAreas.emplace(area.regionId, area.areaSquareMeters);
        }
        checkNear(rotatedAreas.at(1), 0.4, 3.0e-15,
                  "face partition: Y/Z chart keeps negative-side area");
        checkNear(rotatedAreas.at(2), 0.6, 3.0e-15,
                   "face partition: Y/Z chart keeps positive-side area");
        checkRegionMoments(
            rotatedPartitions, *rotatedFound,
            rotated.topology.activeFaces[rotatedFound->activeFaceIndex]);
    }

    SceneFluidFacePartitionLimits limits;
    limits.maximumReferences = 2;
    expectLimited(
        [&] { static_cast<void>(partitions(p, {}, limits)); },
        "face partition: boundary-chain references are bounded");
    limits = {};
    limits.maximumSegmentPairTests = first.segmentPairTestCount - 1;
    expectLimited(
        [&] { static_cast<void>(partitions(p, {}, limits)); },
        "face partition: boundary-chain intersection work is bounded");
    auto corrupt = first;
    corrupt.openChainReferences.front() = p.chains.chains.size();
    expectInvalid(
        [&] { validateSceneFluidFacePartitions(
            corrupt, p.surface.definition, p.state, grid(), p.candidates,
            p.intersections, p.patches, p.ownership, p.crossings,
            p.topology, p.graph, p.chains, p.loops); },
        "face partition: corrupt boundary-chain reference is rejected");
}

void testBoundaryJunctionArrangement() {
    Pipeline p(boundaryJunctionScene());
    const auto first = partitions(p);
    const auto repeated = partitions(p);
    const auto found = std::ranges::find(
        first.partitions,
        SceneFluidFacePartitionKind::BoundaryChainArrangement,
        &SceneFluidFacePartition::kind);
    check(first == repeated && found != first.partitions.end(),
          "face partition: three-region boundary junction resolves deterministically");
    if (found == first.partitions.end()) return;
    const auto& activeFace = p.topology.activeFaces[
        found->activeFaceIndex];
    check(activeFace.axis == GridFaceAxis::X
              && activeFace.i == 2
              && activeFace.j == 1
              && activeFace.k == 1
              && found->rootExteriorRegionId == invalidStableId
              && found->loopReferenceCount == 0
              && found->openChainReferenceCount == 3
              && found->regionAreaCount == 3,
          "face partition: junction arrangement retains all source chains");
    std::map<StableId, double> areas;
    for (std::size_t offset = 0;
         offset < found->regionAreaCount; ++offset) {
        const auto& area = first.regionAreas[
            found->firstRegionArea + offset];
        areas.emplace(area.regionId, area.areaSquareMeters);
    }
    checkNear(areas.at(1), 0.5, 3.0e-15,
              "face partition: junction retains lower exterior sector");
    checkNear(areas.at(2), 0.25, 3.0e-15,
              "face partition: junction retains upper-right cell sector");
    checkNear(areas.at(3), 0.25, 3.0e-15,
              "face partition: junction retains upper-left cell sector");
    checkNear(found->assignedAreaSquareMeters, 1.0, 3.0e-15,
              "face partition: junction sectors close the face");
    checkNear(found->areaResidualSquareMeters, 0.0, 3.0e-15,
              "face partition: junction area residual is exact");
    checkRegionMoments(first, *found, activeFace);
    validateSceneFluidFacePartitions(
        first, p.surface.definition, p.state, grid(), p.candidates,
        p.intersections, p.patches, p.ownership, p.crossings, p.topology,
        p.graph, p.chains, p.loops);

    for (const GridFaceAxis axis : {
             GridFaceAxis::Y, GridFaceAxis::Z}) {
        Pipeline rotated(boundaryJunctionScene(axis));
        const auto rotatedPartitions = partitions(rotated);
        const auto rotatedFound = std::ranges::find(
            rotatedPartitions.partitions,
            SceneFluidFacePartitionKind::BoundaryChainArrangement,
            &SceneFluidFacePartition::kind);
        check(rotatedFound != rotatedPartitions.partitions.end(),
              "face partition: right-handed Y/Z charts retain junction sectors");
        if (rotatedFound == rotatedPartitions.partitions.end()) continue;
        std::map<StableId, double> rotatedAreas;
        for (std::size_t offset = 0;
             offset < rotatedFound->regionAreaCount; ++offset) {
            const auto& area = rotatedPartitions.regionAreas[
                rotatedFound->firstRegionArea + offset];
            rotatedAreas.emplace(area.regionId, area.areaSquareMeters);
        }
        checkNear(rotatedAreas.at(1), 0.5, 3.0e-15,
                  "face partition: Y/Z junction keeps exterior sector");
        checkNear(rotatedAreas.at(2), 0.25, 3.0e-15,
                  "face partition: Y/Z junction keeps first cell sector");
        checkNear(rotatedAreas.at(3), 0.25, 3.0e-15,
                   "face partition: Y/Z junction keeps second cell sector");
        checkRegionMoments(
            rotatedPartitions, *rotatedFound,
            rotated.topology.activeFaces[rotatedFound->activeFaceIndex]);
    }

    SceneFluidFacePartitionLimits limits;
    limits.maximumReferences =
        first.openChainReferences.size() + first.loopReferences.size()
        + first.regionAreas.size() - 1;
    expectLimited(
        [&] { static_cast<void>(partitions(p, {}, limits)); },
        "face partition: junction arrangement references are bounded");
    limits = {};
    limits.maximumSegmentPairTests = first.segmentPairTestCount - 1;
    expectLimited(
        [&] { static_cast<void>(partitions(p, {}, limits)); },
        "face partition: junction arrangement intersection work is bounded");

    auto conflictingScene = boundaryJunctionScene();
    std::swap(conflictingScene.triangles[1].negativeSideRegionId,
              conflictingScene.triangles[1].positiveSideRegionId);
    Pipeline conflicting(std::move(conflictingScene));
    expectInvalid(
        [&] { static_cast<void>(partitions(conflicting)); },
        "face partition: junction arrangement rejects conflicting region winding");
}

void testSameRegionSheetDoesNotBlockJunctionArrangement() {
    Pipeline p(boundaryJunctionWithSameRegionSheetScene());
    const auto first = partitions(p);
    const auto repeated = partitions(p);
    const auto found = std::ranges::find_if(
        first.partitions,
        [&](const SceneFluidFacePartition& partition) {
            const auto& face = p.topology.activeFaces[
                partition.activeFaceIndex];
            return face.axis == GridFaceAxis::X
                && face.i == 2 && face.j == 1 && face.k == 1;
        });
    check(first == repeated && found != first.partitions.end()
              && first.ignoredSameRegionChainCount == 1,
          "face partition: same-region sheet is audited without blocking a junction");
    if (found == first.partitions.end()) return;
    check(found->kind
                  == SceneFluidFacePartitionKind::BoundaryChainArrangement
              && found->openChainReferenceCount == 3
              && found->regionAreaCount == 3,
          "face partition: only region-separating junction chains own pressure areas");
    std::map<StableId, double> areas;
    for (std::size_t offset = 0;
         offset < found->regionAreaCount; ++offset) {
        const auto& area = first.regionAreas[
            found->firstRegionArea + offset];
        areas.emplace(area.regionId, area.areaSquareMeters);
    }
    checkNear(areas.at(1), 0.5, 3.0e-15,
              "face partition: ignored same-region sheet preserves exterior area");
    checkNear(areas.at(2), 0.25, 3.0e-15,
              "face partition: ignored same-region sheet preserves first cell area");
    checkNear(areas.at(3), 0.25, 3.0e-15,
              "face partition: ignored same-region sheet preserves second cell area");
    validateSceneFluidFacePartitions(
        first, p.surface.definition, p.state, grid(), p.candidates,
        p.intersections, p.patches, p.ownership, p.crossings, p.topology,
        p.graph, p.chains, p.loops);
}

void testSameRegionSheetOwnsFullRegionArea() {
    Pipeline p(sameRegionSheetScene());
    const auto first = partitions(p);
    const auto repeated = partitions(p);
    check(first == repeated && first.partitions.size() == 1
              && first.unresolvedActiveFaceCount == 0
              && first.ignoredSameRegionChainCount == 1,
          "face partition: isolated same-region sheet resolves deterministically");
    if (first.partitions.empty()) return;
    const auto& partition = first.partitions.front();
    check(partition.kind == SceneFluidFacePartitionKind::SameRegionSheets
              && partition.rootExteriorRegionId == 2
              && partition.loopReferenceCount == 0
              && partition.openChainReferenceCount == 1
              && partition.regionAreaCount == 1,
          "face partition: same-region sheet retains explicit source identity");
    const auto& area = first.regionAreas[partition.firstRegionArea];
    check(area.regionId == 2,
          "face partition: same-region sheet keeps its authored region");
    checkNear(area.areaSquareMeters, 1.0, 3.0e-15,
              "face partition: same-region sheet owns the complete face area");
    checkRegionMoments(
        first, partition,
        p.topology.activeFaces[partition.activeFaceIndex]);
    checkNear(partition.areaResidualSquareMeters, 0.0, 0.0,
              "face partition: same-region sheet closes with zero residual");
    validateSceneFluidFacePartitions(
        first, p.surface.definition, p.state, grid(), p.candidates,
        p.intersections, p.patches, p.ownership, p.crossings, p.topology,
        p.graph, p.chains, p.loops);
}

void testLimitsAndValidation() {
    Pipeline p;
    SceneFluidFacePartitionLimits limits; limits.maximumSegmentPairTests=8;
    expectLimited([&]{(void)partitions(p,{},limits);},"face partition: pair tests are bounded");
    limits={}; limits.maximumPartitions=0;
    expectLimited([&]{(void)partitions(p,{},limits);},"face partition: partition count is bounded");
    limits={}; limits.maximumReferences=4;
    expectLimited([&]{(void)partitions(p,{},limits);},"face partition: result references are bounded");
    limits={};
    limits.maximumPartitionBytes = sizeof(SceneFluidFacePartition)
        + 2 * sizeof(SceneFluidFaceLoopContainment)
        + 2 * sizeof(std::size_t)
        + 3 * sizeof(SceneFluidFaceRegionArea) - 1;
    expectLimited([&]{(void)partitions(p,{},limits);},"face partition: result storage is bounded");
    SceneFluidFacePartitionSettings settings;
    settings.geometryToleranceMeters=1.0e-3;
    expectInvalid([&]{(void)partitions(p,settings);},"face partition: over-broad geometry tolerance is rejected");
    auto accepted=partitions(p); auto corrupt=accepted;
    corrupt.regionAreas.front().areaSquareMeters+=0.01;
    expectInvalid([&]{validateSceneFluidFacePartitions(corrupt,p.surface.definition,p.state,grid(),
        p.candidates,p.intersections,p.patches,p.ownership,p.crossings,p.topology,p.graph,p.chains,p.loops);},
        "face partition: corrupt area is rejected transactionally");
}
}
int main(){testNestedPartition();testBoundaryOpenChainPartition();testBoundaryJunctionArrangement();testSameRegionSheetDoesNotBlockJunctionArrangement();testSameRegionSheetOwnsFullRegionArea();testLimitsAndValidation();if(failures){std::fprintf(stderr,"%d face-partition failures\n",failures);return 1;}std::printf("scene face-partition tests passed\n");}
