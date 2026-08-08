#include "fluid/scene_surface_face_partition.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>

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
PeriodicCartesianGrid grid(){ return {{4,4,4},{},{4,4,4}}; }

struct Pipeline {
    Scene scene=nestedScene();
    SceneFluidSurfaceAssembly surface=assembleSceneFluidSurface(scene);
    SceneStructureAssembly sa=assembleSceneStructure(scene);
    Structure structure{sa.definition};
    SceneFluidSurfaceState state=captureSceneFluidSurfaceState(surface.definition,sa.mappings,structure);
    SceneFluidGridCandidateSet candidates=buildSceneFluidGridCandidates(surface.definition,state,grid());
    SceneFluidGridIntersectionSet intersections=intersectSceneFluidSurfaceWithGrid(surface.definition,state,grid(),candidates);
    SceneFluidGridPatchSet patches=clipSceneFluidSurfaceToCells(surface.definition,state,grid(),candidates,intersections);
    SceneFluidPatchOwnership ownership=ownSceneFluidSurfacePatches(surface.definition,state,grid(),candidates,intersections,patches);
    SceneFluidFaceCrossingSet crossings=buildSceneFluidFaceCrossings(surface.definition,state,grid(),candidates,intersections,patches,ownership);
    SceneFluidFaceTopology topology=buildSceneFluidFaceTopology(surface.definition,state,grid(),candidates,intersections,patches,ownership,crossings);
    SceneFluidFaceGraph graph=buildSceneFluidFaceGraph(surface.definition,state,grid(),candidates,intersections,patches,ownership,crossings,topology);
    SceneFluidFaceChainSet chains=buildSceneFluidFaceChains(surface.definition,state,grid(),candidates,intersections,patches,ownership,crossings,topology,graph);
    SceneFluidFaceLoopSet loops=buildSceneFluidFaceLoops(surface.definition,state,grid(),candidates,intersections,patches,ownership,crossings,topology,graph,chains);
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
int main(){testNestedPartition();testLimitsAndValidation();if(failures){std::fprintf(stderr,"%d face-partition failures\n",failures);return 1;}std::printf("scene face-partition tests passed\n");}
