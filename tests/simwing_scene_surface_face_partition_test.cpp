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

void testNestedPartition() {
    Pipeline p;
    auto a=partitions(p), b=partitions(p);
    check(a==b && a.fingerprint!=0 && a.partitions.size()==1
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
int main(){testNestedPartition();testBoundaryOpenChainPartition();testLimitsAndValidation();if(failures){std::fprintf(stderr,"%d face-partition failures\n",failures);return 1;}std::printf("scene face-partition tests passed\n");}
