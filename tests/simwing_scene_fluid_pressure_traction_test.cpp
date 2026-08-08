#include "scene_fluid_pressure_traction.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {
using namespace simwing::fsi;
using namespace simwing::fsi::fluid;
int failures=0;
void check(bool value,const char* message){if(!value){std::fprintf(stderr,"FAIL: %s\n",message);++failures;}}
void near(double a,double b,double t,const char* message){check(std::isfinite(a)&&std::abs(a-b)<=t,message);}
template<class F> void invalid(F&& f,const char* message){bool caught=false;try{f();}catch(const std::invalid_argument&){caught=true;}check(caught,message);}

Scene scene(bool reverse=false){
    Scene s; s.metadata.designChecksum="sha256:pressure-traction";s.metadata.exporterVersion="test/1";
    s.regions={{1,RegionKind::Outside,"outside"},{2,RegionKind::Cell,"cell"}};
    s.vertices={{10,{1.1,1.1,1.5}},{11,{2.7,1.1,1.5}},{12,{1.1,2.7,1.5}}};
    s.fabricMaterials={{100,"fabric",900,650,220,.015,.041,.02,.0125,2.5e-12}};
    s.triangles={{500,reverse?std::array<StableId,3>{10,12,11}:std::array<StableId,3>{10,11,12},
                  {{{0,0},{1.6,0},{0,1.6}}},1,2,100,900,SurfaceRole::Skin}}; return s;
}
PeriodicCartesianGrid grid(){return {{4,4,4},{},{4,4,4}};}
struct P {
    Scene s; SceneFluidSurfaceAssembly surface; SceneStructureAssembly sa; Structure structure;
    SceneFluidSurfaceTransfer transfer; SceneFluidSurfaceState state; SceneFluidGridCandidateSet candidates;
    SceneFluidGridIntersectionSet intersections; SceneFluidGridPatchSet patches; SceneFluidPatchOwnership ownership;
    SceneFluidQuadratureDefinition quadrature;
    explicit P(bool reverse=false):s(scene(reverse)),surface(assembleSceneFluidSurface(s)),sa(assembleSceneStructure(s)),
      structure(sa.definition),transfer(surface.definition,sa.mappings,structure),
      state(captureSceneFluidSurfaceState(surface.definition,sa.mappings,structure)),
      candidates(buildSceneFluidGridCandidates(surface.definition,state,grid())),
      intersections(intersectSceneFluidSurfaceWithGrid(surface.definition,state,grid(),candidates)),
      patches(clipSceneFluidSurfaceToCells(surface.definition,state,grid(),candidates,intersections)),
      ownership(ownSceneFluidSurfacePatches(surface.definition,state,grid(),candidates,intersections,patches)),
      quadrature(buildSceneFluidQuadrature(surface.definition,state,grid(),candidates,intersections,patches,ownership,transfer)){}
};
std::vector<SceneFluidQuadraturePressure> pressure(const P& p,double negative=100,double positive=20){
    std::vector<SceneFluidQuadraturePressure> v;for(const auto& q:p.quadrature.points)v.push_back({q.stableId,negative,positive});return v;
}
void testPressureDifference(){
    P p; auto samples=pressure(p); auto tractions=buildSceneFluidPressureTractions(p.surface.definition,p.state,p.quadrature,samples);
    check(tractions.size()==3,"pressure traction: one ordered traction per owned patch");
    for(const auto& t:tractions){near(t.tractionPascals.x,0,0,"pressure traction: normal X is zero");near(t.tractionPascals.y,0,0,"pressure traction: normal Y is zero");near(t.tractionPascals.z,80,2e-14,"pressure traction: pressure jump follows oriented normal");}
    auto result=evaluateSceneFluidPressureQuadrature(p.surface.definition,p.state,p.transfer,p.quadrature,samples);
    near(result.diagnostics().integratedSurfaceForceNewtons.z,102.4,3e-13,"pressure traction: pressure force integrates once");
    p.transfer.addLoadsTo(p.structure,result);
    near(p.structure.diagnostics().pendingExternalForceNewtons.z,102.4,3e-13,"pressure traction: pressure load reaches Structure");
    auto equal=pressure(p,50,50);auto zero=evaluateSceneFluidPressureQuadrature(p.surface.definition,p.state,p.transfer,p.quadrature,equal);
    near(zero.diagnostics().integratedSurfaceForceNewtons.z,0,0,"pressure traction: equal one-sided pressure cancels exactly");
}
void testWindingAndValidation(){
    P reversed(true);auto r=evaluateSceneFluidPressureQuadrature(reversed.surface.definition,reversed.state,reversed.transfer,reversed.quadrature,pressure(reversed));
    near(r.diagnostics().integratedSurfaceForceNewtons.z,-102.4,3e-13,"pressure traction: reversed winding reverses force without relabeling sides");
    P p;auto samples=pressure(p);std::swap(samples[0],samples[1]);
    invalid([&]{(void)buildSceneFluidPressureTractions(p.surface.definition,p.state,p.quadrature,samples);},"pressure traction: reordered samples reject");
    samples=pressure(p);samples[0].negativeSidePressurePascals=std::numeric_limits<double>::quiet_NaN();
    invalid([&]{(void)buildSceneFluidPressureTractions(p.surface.definition,p.state,p.quadrature,samples);},"pressure traction: non-finite samples reject");
}
}
int main(){testPressureDifference();testWindingAndValidation();if(failures){std::fprintf(stderr,"%d pressure-traction failures\n",failures);return 1;}std::printf("pressure-traction tests passed\n");}
