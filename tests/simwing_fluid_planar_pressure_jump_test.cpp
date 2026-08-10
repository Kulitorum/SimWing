#include "fluid/planar_pressure_jump.h"
#include "fluid/planar_region_flux.h"
#include "fluid/planar_region_opening_flow.h"
#include "fluid/planar_region_opening_power.h"
#include "fluid/planar_region_fragment.h"
#include "fluid/planar_region_sweep.h"
#include "fluid/projection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

namespace {

using namespace simwing::fsi::fluid;

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
                     "FAIL: %s (actual %.17g, expected %.17g)\n",
                     message, actual, expected);
        ++failures;
    }
}

template<typename Callback>
void expectRejected(Callback&& callback, const char* message) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected, message);
}

PeriodicCartesianGrid grid() {
    return {{4, 2, 2}, {-2.0, -1.0, -1.0}, {2.0, 1.0, 1.0}};
}

ProjectionSettings strictSettings() {
    ProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.01;
    settings.absoluteResidualTolerance = 1.0e-11;
    settings.relativeResidualTolerance = 1.0e-13;
    settings.maximumIterations = 1000;
    return settings;
}

std::vector<PlanarPressureJumpLayerDefinition> pocketLayers() {
    return {
        {
            10, 1, 2,
            {movingPlanarFaceTopologyVersion,
             GridFaceAxis::X, 1, 0},
            -0.8, 70.0,
        },
        {
            20, 2, 1,
            {movingPlanarFaceTopologyVersion,
             GridFaceAxis::X, 2, 0},
            -0.2, -70.0,
        },
    };
}

const GridFacePressureJump* findCrossing(
    const SharpPressureJumpField& field,
    const std::uint64_t surfaceStableId,
    const std::size_t i,
    const std::size_t j = 0,
    const std::size_t k = 0) {
    const auto found = std::ranges::find_if(
        field.faces(),
        [&](const auto& face) {
            return face.surfaceStableId == surfaceStableId
                && face.i == i && face.j == j && face.k == k;
        });
    return found == field.faces().end() ? nullptr : &*found;
}

double maximumVelocity(const MacVelocityField& velocity) {
    double result = 0.0;
    for (std::size_t index = 0; index < velocity.xFaces().size(); ++index) {
        result = std::max({
            result,
            std::abs(velocity.xFaces()[index]),
            std::abs(velocity.yFaces()[index]),
            std::abs(velocity.zFaces()[index]),
        });
    }
    return result;
}

const PlanarPressureRegionSummary* findRegion(
    const StaticPlanarPressureRegionProfile& profile,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find_if(
        profile.regions,
        [&](const auto& region) {
            return region.regionStableId == regionStableId;
        });
    return found == profile.regions.end() ? nullptr : &*found;
}

const PlanarPressureRegionSweepSummary* findSweepRegion(
    const PlanarPressureRegionSweepLedger& ledger,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find_if(
        ledger.regions,
        [&](const auto& region) {
            return region.regionStableId == regionStableId;
        });
    return found == ledger.regions.end() ? nullptr : &*found;
}

const PlanarPressureRegionFluxSummary* findFluxRegion(
    const PlanarPressureRegionFluxCompatibility& compatibility,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find_if(
        compatibility.regions,
        [&](const auto& region) {
            return region.regionStableId == regionStableId;
        });
    return found == compatibility.regions.end() ? nullptr : &*found;
}

const PlanarPressureRegionOpeningBalance* findOpeningRegion(
    const PlanarPressureRegionOpeningFlowAllocation& allocation,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        allocation.regions, regionStableId,
        &PlanarPressureRegionOpeningBalance::regionStableId);
    return found == allocation.regions.end() ? nullptr : &*found;
}

const PlanarPressureRegionFragmentRegionSummary* findFragmentRegion(
    const PlanarPressureRegionFragmentSet& fragments,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        fragments.regions, regionStableId,
        &PlanarPressureRegionFragmentRegionSummary::regionStableId);
    return found == fragments.regions.end() ? nullptr : &*found;
}

void testCanonicalAssemblyAndSameFacePocket() {
    const auto geometry = grid();
    const auto authored = pocketLayers();
    auto reversed = authored;
    std::ranges::reverse(reversed);
    const auto first = makePlanarPressureJumpField(geometry, authored);
    const auto second = makePlanarPressureJumpField(geometry, reversed);
    check(first == second && first.faceCount() == 8,
          "planar jump assembly canonicalizes a complete two-layer pocket");
    check(first.xFaceJumpsPascals()[geometry.cellIndex(1, 0, 0)] == 70.0
              && first.xFaceJumpsPascals()
                     [geometry.cellIndex(2, 0, 0)] == -70.0,
          "planar jump assembly places each separated layer on its owned face");

    const auto translated = translatePlanarPressureJumpLayers(
        geometry, authored, 0.5);
    check(translated.layers.size() == 2
              && translated.rebaseDirections.size() == 2
              && translated.layers[0].surfaceStableId == 10
              && translated.layers[0].topology.faceCoordinate == 2
              && translated.rebaseDirections[0]
                  == PlanarTopologyRebaseDirection::Positive
              && translated.layers[1].surfaceStableId == 20
              && translated.layers[1].topology.faceCoordinate == 2
              && translated.rebaseDirections[1]
                  == PlanarTopologyRebaseDirection::None,
          "rigid layer translation moves one boundary into its adjacent segment");
    checkNear(translated.layers[0].physicalPlaneCoordinateMeters,
              -0.3, 1.0e-15,
              "rigid layer translation retains the first physical position");
    checkNear(translated.layers[1].physicalPlaneCoordinateMeters,
              0.3, 1.0e-15,
              "rigid layer translation retains the second physical position");
    const auto collapsed = makePlanarPressureJumpField(
        geometry, translated.layers);
    const auto* firstCrossing = findCrossing(collapsed, 10, 2);
    const auto* secondCrossing = findCrossing(collapsed, 20, 2);
    check(firstCrossing != nullptr && secondCrossing != nullptr,
          "same-face pocket retains both stable surface layers");
    if (firstCrossing != nullptr && secondCrossing != nullptr) {
        checkNear(firstCrossing->crossingFraction, 0.2, 1.0e-15,
                  "same-face pocket retains its first subcell position");
        checkNear(secondCrossing->crossingFraction, 0.8, 1.0e-15,
                  "same-face pocket retains its second subcell position");
        check(firstCrossing->plusRegionStableId
                      == secondCrossing->minusRegionStableId
                  && firstCrossing->pressureJumpPascals == 70.0
                  && secondCrossing->pressureJumpPascals == -70.0,
              "same-face pocket retains the ordered intermediate region");
    }
    check(collapsed.xFaceJumpsPascals()
                  [geometry.cellIndex(2, 0, 0)] == 0.0,
          "same-face pocket explicitly exposes the zero-net dense-stencil limitation");

    MacVelocityField velocity(geometry);
    CellScalarField pressure(geometry);
    const auto diagnostics = projectVelocityWithPressureJumps(
        geometry, velocity, pressure, collapsed, strictSettings());
    check(diagnostics.converged
              && diagnostics.pressureJumpFaceCount == 8
              && maximumAbsoluteValue(pressure) == 0.0
              && maximumVelocity(velocity) == 0.0,
          "same-face pocket remains metadata-only in the current aggregate projection oracle");
}

void testSeparatedResolutionAndPeriodicWrap() {
    const auto geometry = grid();
    const auto sameFace = translatePlanarPressureJumpLayers(
        geometry, pocketLayers(), 0.5);
    const auto separated = translatePlanarPressureJumpLayers(
        geometry, sameFace.layers, 0.6);
    check(separated.layers[0].topology.faceCoordinate == 2
              && separated.rebaseDirections[0]
                  == PlanarTopologyRebaseDirection::None
              && separated.layers[1].topology.faceCoordinate == 3
              && separated.rebaseDirections[1]
                  == PlanarTopologyRebaseDirection::Positive,
          "rigid layer translation separates the pocket across adjacent faces");
    const auto resolved = makePlanarPressureJumpField(
        geometry, separated.layers);
    MacVelocityField velocity(geometry);
    CellScalarField pressure(geometry);
    const auto diagnostics = projectVelocityWithPressureJumps(
        geometry, velocity, pressure, resolved, strictSettings());
    check(diagnostics.converged && maximumVelocity(velocity) < 1.0e-13,
          "separated pocket projection creates no spurious flow");
    const auto counts = geometry.cellCounts();
    double maximumPressureError = 0.0;
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const double expected = i == 2 ? 52.5 : -17.5;
                maximumPressureError = std::max(
                    maximumPressureError,
                    std::abs(pressure.values()
                                 [geometry.cellIndex(i, j, k)]
                             - expected));
            }
        }
    }
    check(maximumPressureError < 1.0e-12,
          "separated pocket recovers its sharp intermediate pressure region");

    const std::vector<PlanarPressureJumpLayerDefinition> wrapping{
        {
            10, 1, 2,
            {movingPlanarFaceTopologyVersion,
             GridFaceAxis::X, 3, 0},
            1.2, 70.0,
        },
        {
            20, 2, 1,
            {movingPlanarFaceTopologyVersion,
             GridFaceAxis::X, 0, 1},
            1.8, -70.0,
        },
    };
    const auto wrapped = translatePlanarPressureJumpLayers(
        geometry, wrapping, 0.5);
    check(wrapped.layers[0].topology.faceCoordinate == 0
              && wrapped.layers[0].topology.periodicImage == 1
              && wrapped.rebaseDirections[0]
                  == PlanarTopologyRebaseDirection::Positive
              && wrapped.layers[1].topology.faceCoordinate == 0
              && wrapped.layers[1].topology.periodicImage == 1
              && wrapped.rebaseDirections[1]
                  == PlanarTopologyRebaseDirection::None,
          "layer translation preserves a same-face pocket through a periodic wrap");
    const auto wrappedField = makePlanarPressureJumpField(
        geometry, wrapped.layers);
    check(wrappedField.faceCount() == 8
              && wrappedField.xFaceJumpsPascals()
                     [geometry.cellIndex(0, 0, 0)] == 0.0,
          "periodically wrapped pocket retains both layers and zero aggregate jump");

    const std::vector<PlanarPressureJumpLayerDefinition> reverseWrapping{
        {
            10, 1, 2,
            {movingPlanarFaceTopologyVersion,
             GridFaceAxis::X, 3, -1},
            -2.8, 70.0,
        },
        {
            20, 2, 1,
            {movingPlanarFaceTopologyVersion,
             GridFaceAxis::X, 0, 0},
            -2.2, -70.0,
        },
    };
    const auto reverseWrapped = translatePlanarPressureJumpLayers(
        geometry, reverseWrapping, -0.5);
    check(reverseWrapped.layers[0].topology.faceCoordinate == 3
              && reverseWrapped.layers[0].topology.periodicImage == -1
              && reverseWrapped.rebaseDirections[0]
                  == PlanarTopologyRebaseDirection::None
              && reverseWrapped.layers[1].topology.faceCoordinate == 3
              && reverseWrapped.layers[1].topology.periodicImage == -1
              && reverseWrapped.rebaseDirections[1]
                  == PlanarTopologyRebaseDirection::Negative,
          "layer translation preserves a same-face pocket through a negative periodic wrap");
    const auto reverseWrappedField = makePlanarPressureJumpField(
        geometry, reverseWrapped.layers);
    check(reverseWrappedField.faceCount() == 8
              && reverseWrappedField.xFaceJumpsPascals()
                     [geometry.cellIndex(3, 0, 0)] == 0.0,
          "negative periodically wrapped pocket retains both layers and zero aggregate jump");
}

void testStaticRegionalPressureProfile() {
    const auto geometry = grid();
    auto reversed = pocketLayers();
    std::ranges::reverse(reversed);
    const auto first = makeStaticPlanarPressureRegionProfile(
        geometry, pocketLayers());
    const auto second = makeStaticPlanarPressureRegionProfile(
        geometry, reversed);
    check(first == second
              && first.axis == GridFaceAxis::X
              && first.intervals.size() == 2
              && first.regions.size() == 2,
          "static regional profile canonicalizes one complete thin pocket");
    check(first.intervals[0].lowerSurfaceStableId == 10
              && first.intervals[0].upperSurfaceStableId == 20
              && first.intervals[0].regionStableId == 2
              && first.intervals[1].lowerSurfaceStableId == 20
              && first.intervals[1].upperSurfaceStableId == 10
              && first.intervals[1].regionStableId == 1,
          "static regional profile retains the ordered periodic intervals");
    checkNear(first.geometricDomainVolumeCubicMeters, 16.0, 0.0,
              "static regional profile retains the grid domain volume");
    checkNear(first.intervalVolumeCubicMeters, 16.0, 2.0e-15,
              "static regional profile partitions the complete domain volume");
    checkNear(first.volumeClosureResidualCubicMeters, 0.0, 2.0e-15,
              "static regional profile closes its interval volume ledger");
    checkNear(first.intervals[0].volumeCubicMeters, 2.4, 2.0e-15,
              "static regional profile retains the physical thin-pocket volume");
    const auto* exterior = findRegion(first, 1);
    const auto* interior = findRegion(first, 2);
    check(exterior != nullptr && interior != nullptr,
          "static regional profile summarizes both authored regions");
    if (exterior != nullptr && interior != nullptr) {
        checkNear(exterior->volumeCubicMeters, 13.6, 2.0e-15,
                  "static regional profile retains the exterior volume");
        checkNear(interior->volumeCubicMeters, 2.4, 2.0e-15,
                  "static regional profile retains the pocket volume");
        checkNear(exterior->pressurePascals, -10.5, 2.0e-14,
                  "static regional profile gauges the exterior pressure by volume");
        checkNear(interior->pressurePascals, 59.5, 2.0e-14,
                  "static regional profile recovers the intermediate pocket pressure");
        checkNear(interior->pressurePascals - exterior->pressurePascals,
                  70.0, 2.0e-14,
                  "static regional profile preserves the authored pressure jump");
    }
    checkNear(first.achievedVolumeMeanPressurePascals, 0.0, 2.0e-15,
              "static regional profile applies the requested zero pressure gauge");
    const auto shifted = makeStaticPlanarPressureRegionProfile(
        geometry, pocketLayers(), 12.5);
    const auto* shiftedExterior = findRegion(shifted, 1);
    const auto* shiftedInterior = findRegion(shifted, 2);
    check(shiftedExterior != nullptr && shiftedInterior != nullptr,
          "static regional profile retains both regions under a shifted gauge");
    if (shiftedExterior != nullptr && shiftedInterior != nullptr) {
        checkNear(shifted.achievedVolumeMeanPressurePascals,
                  12.5, 4.0e-15,
                  "static regional profile applies a nonzero volume pressure gauge");
        checkNear(shiftedInterior->pressurePascals
                      - shiftedExterior->pressurePascals,
                  70.0, 2.0e-14,
                  "static regional pressure differences are gauge invariant");
    }

    const auto sameFace = translatePlanarPressureJumpLayers(
        geometry, pocketLayers(), 0.5);
    const auto separated = translatePlanarPressureJumpLayers(
        geometry, sameFace.layers, 0.6);
    const auto sameFaceProfile = makeStaticPlanarPressureRegionProfile(
        geometry, sameFace.layers);
    const auto separatedProfile = makeStaticPlanarPressureRegionProfile(
        geometry, separated.layers);
    for (const auto* translated : {&sameFaceProfile, &separatedProfile}) {
        const auto* translatedExterior = findRegion(*translated, 1);
        const auto* translatedInterior = findRegion(*translated, 2);
        check(translatedExterior != nullptr && translatedInterior != nullptr,
              "translated static profile retains both pressure regions");
        if (translatedExterior != nullptr && translatedInterior != nullptr) {
            checkNear(translatedInterior->volumeCubicMeters,
                      2.4, 4.0e-15,
                      "rigid layer translation preserves the thin-pocket volume");
            checkNear(translatedInterior->pressurePascals
                          - translatedExterior->pressurePascals,
                      70.0, 2.0e-14,
                      "rigid layer translation preserves the regional pressure difference");
        }
    }

    auto invalid = pocketLayers();
    invalid[1].pressureJumpPascals = -60.0;
    expectRejected(
        [&] { static_cast<void>(makeStaticPlanarPressureRegionProfile(
            geometry, invalid)); },
        "static regional profile rejects a non-closing pressure cycle");
    const std::vector<PlanarPressureJumpLayerDefinition> repeatedRegion{
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 1, 0},
         -1.4, 10.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 1, 0},
         -0.8, -9.0},
        {30, 1, 3,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 2, 0},
         -0.2, 19.0},
        {40, 3, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 2, 0},
         0.4, -20.0},
    };
    expectRejected(
        [&] { static_cast<void>(makeStaticPlanarPressureRegionProfile(
            geometry, repeatedRegion)); },
        "static regional profile rejects an inconsistent repeated region potential");
    expectRejected(
        [&] { static_cast<void>(makeStaticPlanarPressureRegionProfile(
            geometry, pocketLayers(),
            std::numeric_limits<double>::quiet_NaN())); },
        "static regional profile rejects a non-finite pressure gauge");
}

void testRegionalSweepLedger() {
    const auto geometry = grid();
    const auto previous = pocketLayers();
    const auto rigid = translatePlanarPressureJumpLayers(
        geometry, previous, 0.5);
    const auto first = makePlanarPressureRegionSweepLedger(
        geometry, previous, rigid.layers, 1.0);
    auto reversedPrevious = previous;
    auto reversedCurrent = rigid.layers;
    std::ranges::reverse(reversedPrevious);
    std::ranges::reverse(reversedCurrent);
    const auto second = makePlanarPressureRegionSweepLedger(
        geometry, reversedPrevious, reversedCurrent, 1.0);
    check(first == second
              && first.version == planarPressureRegionSweepVersion
              && first.axis == GridFaceAxis::X
              && first.intervals.size() == 2
              && first.regions.size() == 2
              && first.ownedStorageBytes > 0,
          "regional sweep ledger canonicalizes a complete rigid layer motion");
    for (const auto& interval : first.intervals) {
        checkNear(interval.lowerSurfaceDisplacementMeters,
                  0.5, 1.0e-15,
                  "rigid regional sweep retains each lower-layer displacement");
        checkNear(interval.upperSurfaceDisplacementMeters,
                  0.5, 1.0e-15,
                  "rigid regional sweep retains each upper-layer displacement");
        checkNear(interval.lowerSurfaceVelocityMetersPerSecond,
                  0.5, 1.0e-15,
                  "rigid regional sweep derives the lower-layer velocity");
        checkNear(interval.upperSurfaceVelocityMetersPerSecond,
                  0.5, 1.0e-15,
                  "rigid regional sweep derives the upper-layer velocity");
        checkNear(interval.boundarySweptVolumeCubicMeters,
                  0.0, 4.0e-15,
                  "rigid regional sweep has no interval volume sweep");
        checkNear(interval.surfaceGeometryResidualCubicMeters,
                  0.0, 4.0e-15,
                  "rigid regional sweep closes interval geometry");
    }
    checkNear(first.globalGeometryVolumeChangeCubicMeters,
              0.0, 4.0e-15,
              "rigid regional sweep preserves global volume");
    checkNear(first.globalBoundarySweptVolumeCubicMeters,
              0.0, 4.0e-15,
              "rigid regional sweep has zero global boundary sweep");
    checkNear(first.globalSurfaceGeometryResidualCubicMeters,
              0.0, 4.0e-15,
              "rigid regional sweep closes its global geometry ledger");

    auto breathing = previous;
    breathing[0].physicalPlaneCoordinateMeters -= 0.1;
    breathing[1].physicalPlaneCoordinateMeters += 0.1;
    const auto expanded = makePlanarPressureRegionSweepLedger(
        geometry, previous, breathing, 0.5);
    const auto* pocket = findSweepRegion(expanded, 2);
    const auto* exterior = findSweepRegion(expanded, 1);
    check(pocket != nullptr && exterior != nullptr,
          "breathing regional sweep retains both region ledgers");
    if (pocket != nullptr && exterior != nullptr) {
        checkNear(pocket->previousVolumeCubicMeters,
                  2.4, 2.0e-15,
                  "breathing regional sweep retains the previous pocket volume");
        checkNear(pocket->currentVolumeCubicMeters,
                  3.2, 2.0e-15,
                  "breathing regional sweep retains the expanded pocket volume");
        checkNear(pocket->geometryVolumeChangeCubicMeters,
                  0.8, 4.0e-15,
                  "breathing regional sweep measures pocket volume growth");
        checkNear(pocket->boundarySweptVolumeCubicMeters,
                  0.8, 4.0e-15,
                  "breathing regional sweep independently integrates layer motion");
        checkNear(pocket->surfaceGeometryResidualCubicMeters,
                  0.0, 4.0e-15,
                  "breathing pocket geometry and layer sweep close");
        checkNear(exterior->geometryVolumeChangeCubicMeters,
                  -0.8, 4.0e-15,
                  "breathing regional sweep assigns opposite exterior volume change");
        checkNear(exterior->boundarySweptVolumeCubicMeters,
                  -0.8, 4.0e-15,
                  "breathing regional sweep assigns opposite exterior boundary sweep");
    }
    checkNear(expanded.globalGeometryVolumeChangeCubicMeters,
              0.0, 4.0e-15,
              "breathing regional sweep conserves total periodic volume");
    checkNear(expanded.globalBoundarySweptVolumeCubicMeters,
              0.0, 4.0e-15,
              "breathing regional sweep cancels global boundary motion");
    checkNear(expanded.maximumAbsoluteSurfaceGeometryResidualCubicMeters,
              0.0, 4.0e-15,
              "breathing regional sweep closes every interval geometry ledger");

    const std::vector<PlanarPressureJumpLayerDefinition> wrapping{
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 3, 0},
         1.2, 70.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 1},
         1.8, -70.0},
    };
    const auto wrapped = translatePlanarPressureJumpLayers(
        geometry, wrapping, 0.5);
    const auto wrapLedger = makePlanarPressureRegionSweepLedger(
        geometry, wrapping, wrapped.layers, 1.0);
    check(wrapLedger.currentProfile.intervals[0]
                  .lowerSurfaceStableId == 10
              && wrapped.layers[0].topology.faceCoordinate == 0
              && wrapped.layers[0].topology.periodicImage == 1,
          "regional sweep ledger retains identity through a periodic topology rebase");
    checkNear(wrapLedger.globalSurfaceGeometryResidualCubicMeters,
              0.0, 4.0e-15,
              "periodically rebased regional sweep closes global geometry");

    const std::vector<PlanarPressureJumpLayerDefinition> reverseWrapping{
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 3, -1},
         -2.8, 70.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 0},
         -2.2, -70.0},
    };
    const auto reverseWrapped = translatePlanarPressureJumpLayers(
        geometry, reverseWrapping, -0.5);
    const auto reverseWrapLedger = makePlanarPressureRegionSweepLedger(
        geometry, reverseWrapping, reverseWrapped.layers, 1.0);
    check(reverseWrapped.layers[1].topology.faceCoordinate == 3
              && reverseWrapped.layers[1].topology.periodicImage == -1,
          "regional sweep ledger retains identity through a negative periodic rebase");
    checkNear(reverseWrapLedger.globalSurfaceGeometryResidualCubicMeters,
              0.0, 4.0e-15,
              "negative periodically rebased regional sweep closes global geometry");
}

void testRegionalSweepAllAxesAndRejection() {
    const auto geometry = grid();
    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> previous{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0},
             -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0},
             -0.2, -70.0},
        };
        auto current = previous;
        current[0].physicalPlaneCoordinateMeters -= 0.1;
        current[1].physicalPlaneCoordinateMeters += 0.1;
        const auto ledger = makePlanarPressureRegionSweepLedger(
            geometry, previous, current, 0.5);
        const auto* pocket = findSweepRegion(ledger, 2);
        const double expectedChange = axis == GridFaceAxis::X
            ? 0.8 : 1.6;
        check(pocket != nullptr,
              "regional sweep retains the breathing pocket on each axis");
        if (pocket != nullptr) {
            checkNear(pocket->geometryVolumeChangeCubicMeters,
                      expectedChange, 8.0e-15,
                      "regional sweep uses the physical transverse area on each axis");
            checkNear(pocket->boundarySweptVolumeCubicMeters,
                      expectedChange, 8.0e-15,
                      "regional sweep closes boundary motion on each axis");
        }
    }

    const auto previous = pocketLayers();
    const auto previousCopy = previous;
    auto validCurrent = translatePlanarPressureJumpLayers(
        geometry, previous, 0.1).layers;
    auto invalid = validCurrent;
    invalid[0].surfaceStableId = 30;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, invalid, 1.0)); },
        "regional sweep rejects changed surface identity");
    invalid = validCurrent;
    invalid[0].plusRegionStableId = 3;
    invalid[1].minusRegionStableId = 3;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, invalid, 1.0)); },
        "regional sweep rejects changed region identity");
    invalid = validCurrent;
    invalid[0].pressureJumpPascals = 71.0;
    invalid[1].pressureJumpPascals = -71.0;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, invalid, 1.0)); },
        "regional sweep rejects changed pressure-jump identity");
    invalid = previous;
    invalid[0].physicalPlaneCoordinateMeters = -0.5;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, invalid, 1.0)); },
        "regional sweep rejects exact topology-boundary motion");
    invalid = previous;
    invalid[0].physicalPlaneCoordinateMeters += 1.6;
    invalid[0].topology.faceCoordinate = 3;
    invalid[1].physicalPlaneCoordinateMeters += 1.6;
    invalid[1].topology.faceCoordinate = 3;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, invalid, 1.0)); },
        "regional sweep rejects a skipped topology segment");
    invalid = previous;
    invalid[0].physicalPlaneCoordinateMeters = -0.1;
    invalid[0].topology.faceCoordinate = 2;
    invalid[1].physicalPlaneCoordinateMeters = -0.3;
    invalid[1].topology.faceCoordinate = 2;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, invalid, 1.0)); },
        "regional sweep rejects layer-order crossing");
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, validCurrent, 0.0)); },
        "regional sweep rejects a zero duration");
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, validCurrent, 1.0,
            {1, 2, 1024})); },
        "regional sweep enforces its layer limit");
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, validCurrent, 1.0,
            {2, 1, 1024 * 1024})); },
        "regional sweep enforces its region limit");
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureRegionSweepLedger(
            geometry, previous, validCurrent, 1.0,
            {2, 2, 1})); },
        "regional sweep enforces its byte limit");
    check(previous == previousCopy,
          "rejected regional sweep cannot mutate either endpoint");
}

void testRegionalFluxCompatibility() {
    const auto geometry = grid();
    const auto previous = pocketLayers();
    const auto rigidLayers = translatePlanarPressureJumpLayers(
        geometry, previous, 0.5).layers;
    const auto rigidSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, rigidLayers, 1.0);
    const auto rigid = assessPlanarPressureRegionFluxCompatibility(
        rigidSweep);
    check(rigid.version == planarPressureRegionFluxVersion
              && rigid.fingerprint != 0
              && rigid.sourceSweepVersion
                  == planarPressureRegionSweepVersion
              && rigid.intervals.size() == 2
              && rigid.regions.size() == 2
              && rigid.failedImpermeableIntervalCount == 0
              && rigid.failedContinuityIntervalCount == 0
              && rigid.failedImpermeableRegionCount == 0
              && rigid.failedContinuityRegionCount == 0
              && rigid.allIntervalsImpermeableWithinTolerance
              && rigid.allIntervalsContinuousWithinTolerance
              && rigid.allRegionsImpermeableWithinTolerance
              && rigid.allRegionsContinuousWithinTolerance
              && rigid.ownedStorageBytes > 0,
          "regional flux compatibility accepts an exactly sealed rigid pocket");
    for (const auto& interval : rigid.intervals) {
        checkNear(interval.leastSquaresFluidVelocityMetersPerSecond,
                  0.5, 1.0e-15,
                  "rigid regional flux recovers the common layer velocity");
        checkNear(interval.maximumAbsoluteInterfaceSlipMetersPerSecond,
                  0.0, 1.0e-15,
                  "rigid regional flux has zero material-relative slip");
        checkNear(
            interval.totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond,
            0.0, 1.0e-15,
            "rigid regional flux has zero one-sided leakage");
        checkNear(interval.continuityResidualCubicMeters,
                  0.0, 4.0e-15,
                  "rigid regional flux closes interval continuity");
    }

    auto breathingLayers = previous;
    breathingLayers[0].physicalPlaneCoordinateMeters -= 0.1;
    breathingLayers[1].physicalPlaneCoordinateMeters += 0.1;
    const auto breathingSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, breathingLayers, 0.5);
    const auto breathing = assessPlanarPressureRegionFluxCompatibility(
        breathingSweep);
    check(breathing.failedImpermeableIntervalCount == 2
              && breathing.failedImpermeableRegionCount == 2
              && breathing.failedContinuityIntervalCount == 0
              && breathing.failedContinuityRegionCount == 0
              && !breathing.allIntervalsImpermeableWithinTolerance
              && !breathing.allRegionsImpermeableWithinTolerance
              && breathing.allIntervalsContinuousWithinTolerance
              && breathing.allRegionsContinuousWithinTolerance,
          "regional flux compatibility exposes breathing as necessarily permeable");
    checkNear(breathing.maximumAbsoluteInterfaceSlipMetersPerSecond,
              0.2, 2.0e-15,
              "breathing regional flux reports the minimum interface slip");
    const auto* pocket = findFluxRegion(breathing, 2);
    const auto* exterior = findFluxRegion(breathing, 1);
    check(pocket != nullptr && exterior != nullptr,
          "breathing regional flux retains interior and exterior summaries");
    if (pocket != nullptr && exterior != nullptr) {
        checkNear(pocket->geometryVolumeChangeCubicMeters,
                  0.8, 4.0e-15,
                  "breathing flux retains pocket volume growth");
        checkNear(pocket->outwardRelativeFlowRateCubicMetersPerSecond,
                  -1.6, 8.0e-15,
                  "breathing flux requires inward material-relative pocket flow");
        checkNear(pocket->integratedOutwardRelativeVolumeCubicMeters,
                  -0.8, 4.0e-15,
                  "breathing flux integrates the inward pocket flow");
        checkNear(pocket->continuityResidualCubicMeters,
                  0.0, 4.0e-15,
                  "breathing pocket volume and relative flow close continuity");
        checkNear(
            pocket->totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond,
            1.6, 8.0e-15,
            "breathing flux reports the minimum absolute pocket leakage");
        checkNear(exterior->geometryVolumeChangeCubicMeters,
                  -0.8, 4.0e-15,
                  "breathing flux retains opposite exterior volume change");
        checkNear(exterior->outwardRelativeFlowRateCubicMetersPerSecond,
                  1.6, 8.0e-15,
                  "breathing flux retains opposite exterior relative flow");
    }
    checkNear(breathing.globalGeometryVolumeChangeCubicMeters,
              0.0, 4.0e-15,
              "breathing regional flux preserves global volume");
    checkNear(breathing.globalIntegratedOutwardRelativeVolumeCubicMeters,
              0.0, 4.0e-15,
              "breathing regional flux cancels globally integrated relative flow");
    checkNear(breathing.globalContinuityResidualCubicMeters,
              0.0, 4.0e-15,
              "breathing regional flux closes global continuity");

    const std::vector<PlanarPressureJumpLayerDefinition> wrapping{
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 3, 0},
         1.2, 70.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 1},
         1.8, -70.0},
    };
    const auto wrapped = translatePlanarPressureJumpLayers(
        geometry, wrapping, 0.5);
    const auto wrappedFlux = assessPlanarPressureRegionFluxCompatibility(
        makePlanarPressureRegionSweepLedger(
            geometry, wrapping, wrapped.layers, 1.0));
    check(wrappedFlux.allRegionsImpermeableWithinTolerance
              && wrappedFlux.allRegionsContinuousWithinTolerance,
          "regional flux remains sealed through a periodic rebase");

    const std::vector<PlanarPressureJumpLayerDefinition> reverseWrapping{
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 3, -1},
         -2.8, 70.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 0},
         -2.2, -70.0},
    };
    const auto reverseWrapped = translatePlanarPressureJumpLayers(
        geometry, reverseWrapping, -0.5);
    const auto reverseWrappedFlux =
        assessPlanarPressureRegionFluxCompatibility(
            makePlanarPressureRegionSweepLedger(
                geometry, reverseWrapping,
                reverseWrapped.layers, 1.0));
    check(reverseWrappedFlux.allRegionsImpermeableWithinTolerance
              && reverseWrappedFlux.allRegionsContinuousWithinTolerance,
          "regional flux remains sealed through a negative periodic rebase");
}

void testRegionalFluxAllAxesAndRejection() {
    const auto geometry = grid();
    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> previous{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0},
             -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0},
             -0.2, -70.0},
        };
        auto current = previous;
        current[0].physicalPlaneCoordinateMeters -= 0.1;
        current[1].physicalPlaneCoordinateMeters += 0.1;
        const auto compatibility =
            assessPlanarPressureRegionFluxCompatibility(
                makePlanarPressureRegionSweepLedger(
                    geometry, previous, current, 0.5));
        const auto* pocket = findFluxRegion(compatibility, 2);
        const double expectedAbsoluteFlow = axis == GridFaceAxis::X
            ? 1.6 : 3.2;
        check(pocket != nullptr,
              "regional flux retains the incompatible pocket on each axis");
        if (pocket != nullptr) {
            checkNear(
                pocket->totalAbsoluteInterfaceRelativeFlowRateCubicMetersPerSecond,
                expectedAbsoluteFlow, 1.6e-14,
                "regional flux uses physical interface area on each axis");
            check(!pocket->impermeableWithinTolerance
                      && pocket->continuityWithinTolerance,
                  "regional flux distinguishes impermeability from continuity on each axis");
        }
    }

    const auto previous = pocketLayers();
    const auto current = translatePlanarPressureJumpLayers(
        geometry, previous, 0.1).layers;
    const auto validSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, current, 1.0);
    const auto validCompatibility =
        assessPlanarPressureRegionFluxCompatibility(validSweep);
    const auto repeatedCompatibility =
        assessPlanarPressureRegionFluxCompatibility(validSweep);
    check(validCompatibility == repeatedCompatibility
              && validCompatibility.fingerprint
                  == repeatedCompatibility.fingerprint,
          "regional flux compatibility is deterministic and fingerprinted");
    validatePlanarPressureRegionFluxCompatibility(validCompatibility);

    auto corrupt = validSweep;
    corrupt.intervals[0].lowerSurfaceVelocityMetersPerSecond += 1.0;
    expectRejected(
        [&] { static_cast<void>(
            assessPlanarPressureRegionFluxCompatibility(corrupt)); },
        "regional flux rejects a corrupted source interval ledger");
    corrupt = validSweep;
    corrupt.previousProfile.intervals[0].pressurePascals += 1.0;
    expectRejected(
        [&] { static_cast<void>(
            assessPlanarPressureRegionFluxCompatibility(corrupt)); },
        "regional flux rejects a corrupted source pressure profile");
    auto invalidSettings = PlanarPressureRegionFluxSettings{};
    invalidSettings.absoluteVelocityToleranceMetersPerSecond = 0.0;
    invalidSettings.relativeVelocityTolerance = 0.0;
    expectRejected(
        [&] { static_cast<void>(
            assessPlanarPressureRegionFluxCompatibility(
                validSweep, invalidSettings)); },
        "regional flux rejects a zero velocity-tolerance policy");
    expectRejected(
        [&] { static_cast<void>(
            assessPlanarPressureRegionFluxCompatibility(
                validSweep, {}, {1, 2, 1024 * 1024})); },
        "regional flux enforces its interval limit");
    expectRejected(
        [&] { static_cast<void>(
            assessPlanarPressureRegionFluxCompatibility(
                validSweep, {}, {2, 1, 1024 * 1024})); },
        "regional flux enforces its region limit");
    expectRejected(
        [&] { static_cast<void>(
            assessPlanarPressureRegionFluxCompatibility(
                validSweep, {}, {2, 2, 1})); },
        "regional flux enforces its byte limit");

    auto corruptCompatibility = validCompatibility;
    corruptCompatibility.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation rejects a zero fingerprint");
    corruptCompatibility = validCompatibility;
    corruptCompatibility.intervals[0]
        .leastSquaresFluidVelocityMetersPerSecond += 1.0;
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation recomputes derived interval values");
    corruptCompatibility = validCompatibility;
    corruptCompatibility.intervals[0].previousVolumeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation binds primitive interval volumes");
    corruptCompatibility = validCompatibility;
    corruptCompatibility.intervals[0].lowerSurfaceStableId = 999;
    corruptCompatibility.intervals[1].upperSurfaceStableId = 999;
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation binds foreign surface identities");
    corruptCompatibility = validCompatibility;
    std::swap(corruptCompatibility.regions[0],
              corruptCompatibility.regions[1]);
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation enforces sorted region order");
    corruptCompatibility = validCompatibility;
    corruptCompatibility.failedContinuityRegionCount += 1;
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation recomputes aggregate counts");
    corruptCompatibility = validCompatibility;
    corruptCompatibility.maximumAbsoluteContinuityResidualCubicMeters =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation rejects non-finite aggregates");
    corruptCompatibility = validCompatibility;
    corruptCompatibility.settings.relativeVelocityTolerance *= 2.0;
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            corruptCompatibility); },
        "regional flux validation binds its tolerance policy");
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            validCompatibility, {1, 2, 1024 * 1024}); },
        "regional flux validation enforces its interval limit");
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            validCompatibility, {2, 1, 1024 * 1024}); },
        "regional flux validation enforces its region limit");
    expectRejected(
        [&] { validatePlanarPressureRegionFluxCompatibility(
            validCompatibility, {2, 2, 1}); },
        "regional flux validation enforces its byte limit");
}

void testRegionalOpeningFlowFeasibility() {
    const auto geometry = grid();
    const auto previous = pocketLayers();
    auto current = previous;
    current[0].physicalPlaneCoordinateMeters -= 0.1;
    current[1].physicalPlaneCoordinateMeters += 0.1;
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, current, 0.5);

    const std::vector<PlanarPressureRegionOpeningDefinition> noOpenings;
    const auto sealed = solvePlanarPressureRegionOpeningFlow(
        sweep, noOpenings);
    check(sealed.version == planarPressureRegionOpeningFlowVersion
              && sealed.fingerprint != 0
              && sealed.openings.empty()
              && sealed.regions.size() == 2
              && sealed.components.size() == 2
              && sealed.failedComponentCount == 2
              && sealed.failedRegionCount == 2
              && !sealed.allComponentsFeasible
              && !sealed.allRegionsWithinTolerance,
          "regional opening flow keeps a sealed breathing pocket infeasible");
    checkNear(sealed.globalGeometryVolumeChangeCubicMeters,
              0.0, 4.0e-15,
              "sealed breathing opening flow retains global volume closure");
    checkNear(sealed.globalContinuityResidualCubicMeters,
              0.0, 4.0e-15,
              "sealed breathing reports global cancellation separately from local failure");

    const auto rigidLayers = translatePlanarPressureJumpLayers(
        geometry, previous, 0.5).layers;
    const auto rigid = solvePlanarPressureRegionOpeningFlow(
        makePlanarPressureRegionSweepLedger(
            geometry, previous, rigidLayers, 1.0),
        noOpenings);
    check(rigid.failedComponentCount == 0
              && rigid.failedRegionCount == 0
              && rigid.allComponentsFeasible
              && rigid.allRegionsWithinTolerance,
          "regional opening flow accepts sealed rigid translation");

    const std::vector<PlanarPressureRegionOpeningDefinition> oneOpening{
        {700, 1, 2, 0.5},
    };
    const auto supplied = solvePlanarPressureRegionOpeningFlow(
        sweep, oneOpening);
    check(supplied.openings.size() == 1
              && supplied.components.size() == 1
              && supplied.failedComponentCount == 0
              && supplied.failedRegionCount == 0
              && supplied.allComponentsFeasible
              && supplied.allRegionsWithinTolerance,
          "regional opening flow closes breathing through an authored opening");
    checkNear(
        supplied.openings[0]
            .relativeVolumeFlowRateCubicMetersPerSecond,
        1.6, 8.0e-15,
        "regional opening flow points from exterior into the growing pocket");
    checkNear(supplied.openings[0]
                  .relativeNormalVelocityMetersPerSecond,
              3.2, 1.6e-14,
              "regional opening flow divides volume flow by authored area");
    const auto* exterior = findOpeningRegion(supplied, 1);
    const auto* pocket = findOpeningRegion(supplied, 2);
    check(exterior != nullptr && pocket != nullptr,
          "regional opening flow retains both regional balances");
    if (exterior != nullptr && pocket != nullptr) {
        checkNear(
            exterior->solvedOutwardRelativeFlowRateCubicMetersPerSecond,
            1.6, 8.0e-15,
            "regional opening flow leaves the shrinking exterior region");
        checkNear(
            pocket->solvedOutwardRelativeFlowRateCubicMetersPerSecond,
            -1.6, 8.0e-15,
            "regional opening flow enters the growing pocket region");
        checkNear(exterior->continuityResidualCubicMeters,
                  0.0, 4.0e-15,
                  "regional opening flow closes exterior continuity");
        checkNear(pocket->continuityResidualCubicMeters,
                  0.0, 4.0e-15,
                  "regional opening flow closes pocket continuity");
    }

    const std::vector<PlanarPressureRegionOpeningDefinition> reverseOpening{
        {700, 2, 1, 0.5},
    };
    const auto reversed = solvePlanarPressureRegionOpeningFlow(
        sweep, reverseOpening);
    checkNear(
        reversed.openings[0]
            .relativeVolumeFlowRateCubicMetersPerSecond,
        -1.6, 8.0e-15,
        "regional opening flow changes sign under authored orientation reversal");
    check(reversed.allComponentsFeasible
              && reversed.allRegionsWithinTolerance,
          "orientation reversal preserves regional opening feasibility");

    const std::vector<PlanarPressureRegionOpeningDefinition> parallel{
        {701, 1, 2, 0.375},
        {700, 1, 2, 0.125},
    };
    const auto distributed = solvePlanarPressureRegionOpeningFlow(
        sweep, parallel);
    const std::vector<PlanarPressureRegionOpeningDefinition> sortedParallel{
        parallel[1], parallel[0],
    };
    const auto repeated = solvePlanarPressureRegionOpeningFlow(
        sweep, sortedParallel);
    check(distributed == repeated
              && distributed.fingerprint == repeated.fingerprint
              && distributed.openings[0].openingStableId == 700
              && distributed.openings[1].openingStableId == 701,
          "regional opening flow canonicalizes authored opening order");
    checkNear(
        distributed.openings[0]
            .relativeVolumeFlowRateCubicMetersPerSecond,
        0.4, 4.0e-15,
        "regional opening flow assigns the small aperture by area");
    checkNear(
        distributed.openings[1]
            .relativeVolumeFlowRateCubicMetersPerSecond,
        1.2, 8.0e-15,
        "regional opening flow assigns the large aperture by area");
    checkNear(
        distributed.openings[0]
            .relativeNormalVelocityMetersPerSecond,
        3.2, 1.6e-14,
        "parallel regional openings share one minimum-norm normal velocity");
    checkNear(
        distributed.openings[1]
            .relativeNormalVelocityMetersPerSecond,
        3.2, 1.6e-14,
        "parallel regional opening velocity is area independent");

    const std::vector<PlanarPressureJumpLayerDefinition> threeRegions{
        {30, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 0},
         -1.8, 30.0},
        {40, 2, 3,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 1, 0},
         -0.8, 20.0},
        {50, 3, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 2, 0},
         0.2, -50.0},
    };
    auto movedThreeRegions = threeRegions;
    movedThreeRegions[1].physicalPlaneCoordinateMeters += 0.1;
    const auto threeRegionSweep = makePlanarPressureRegionSweepLedger(
        geometry, threeRegions, movedThreeRegions, 0.5);
    const std::vector<PlanarPressureRegionOpeningDefinition> serialOpenings{
        {800, 3, 1, 0.25},
        {810, 1, 2, 0.5},
    };
    const auto serial = solvePlanarPressureRegionOpeningFlow(
        threeRegionSweep, serialOpenings);
    check(serial.regions.size() == 3
              && serial.components.size() == 1
              && serial.allComponentsFeasible
              && serial.allRegionsWithinTolerance,
          "regional opening flow solves a three-region serial graph");
    checkNear(
        serial.openings[0]
            .relativeVolumeFlowRateCubicMetersPerSecond,
        0.8, 8.0e-15,
        "serial regional opening carries flow from the shrinking region");
    checkNear(
        serial.openings[1]
            .relativeVolumeFlowRateCubicMetersPerSecond,
        0.8, 8.0e-15,
        "serial regional opening carries flow into the growing region");
    const auto* intermediate = findOpeningRegion(serial, 1);
    check(intermediate != nullptr,
          "serial regional opening graph retains its intermediate region");
    if (intermediate != nullptr) {
        checkNear(
            intermediate
                ->solvedOutwardRelativeFlowRateCubicMetersPerSecond,
            0.0, 8.0e-15,
            "serial regional opening graph balances its intermediate region");
    }
}

void testRegionalOpeningFlowAllAxesAndRejection() {
    const auto geometry = grid();
    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> previous{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0},
             -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0},
             -0.2, -70.0},
        };
        auto current = previous;
        current[0].physicalPlaneCoordinateMeters -= 0.1;
        current[1].physicalPlaneCoordinateMeters += 0.1;
        const auto allocation = solvePlanarPressureRegionOpeningFlow(
            makePlanarPressureRegionSweepLedger(
                geometry, previous, current, 0.5),
            std::vector<PlanarPressureRegionOpeningDefinition>{
                {700, 1, 2, 0.5},
            });
        const double expectedFlow = axis == GridFaceAxis::X ? 1.6 : 3.2;
        check(allocation.allComponentsFeasible
                  && allocation.allRegionsWithinTolerance,
              "regional opening flow closes breathing on every axis");
        checkNear(
            allocation.openings[0]
                .relativeVolumeFlowRateCubicMetersPerSecond,
            expectedFlow, 1.6e-14,
            "regional opening flow follows physical swept volume on every axis");
    }

    const auto previous = pocketLayers();
    auto current = previous;
    current[0].physicalPlaneCoordinateMeters -= 0.1;
    current[1].physicalPlaneCoordinateMeters += 0.1;
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, current, 0.5);
    const std::vector<PlanarPressureRegionOpeningDefinition> openings{
        {700, 1, 2, 0.5},
    };
    const auto valid = solvePlanarPressureRegionOpeningFlow(
        sweep, openings);
    validatePlanarPressureRegionOpeningFlow(valid, sweep, openings);

    auto corrupt = valid;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningFlow(
            corrupt, sweep, openings); },
        "regional opening-flow validation rejects fingerprint corruption");
    corrupt = valid;
    corrupt.openings[0]
        .relativeVolumeFlowRateCubicMetersPerSecond += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningFlow(
            corrupt, sweep, openings); },
        "regional opening-flow validation rejects flow corruption");
    corrupt = valid;
    corrupt.regions[0].continuityResidualCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningFlow(
            corrupt, sweep, openings); },
        "regional opening-flow validation rejects region corruption");
    corrupt = valid;
    corrupt.components[0].feasible = false;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningFlow(
            corrupt, sweep, openings); },
        "regional opening-flow validation rejects component corruption");

    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep,
                std::vector<PlanarPressureRegionOpeningDefinition>{
                    {700, 1, 99, 0.5},
                })); },
        "regional opening flow rejects a foreign region");
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep,
                std::vector<PlanarPressureRegionOpeningDefinition>{
                    {700, 1, 1, 0.5},
                })); },
        "regional opening flow rejects a same-region opening");
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep,
                std::vector<PlanarPressureRegionOpeningDefinition>{
                    {700, 1, 2, 0.25},
                    {700, 1, 2, 0.25},
                })); },
        "regional opening flow rejects a duplicate opening identity");
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep,
                std::vector<PlanarPressureRegionOpeningDefinition>{
                    {700, 1, 2, 0.0},
                })); },
        "regional opening flow rejects a zero opening area");

    auto invalidSettings = PlanarPressureRegionOpeningFlowSettings{};
    invalidSettings.relativeCholeskyPivotTolerance = 0.0;
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep, openings, invalidSettings)); },
        "regional opening flow rejects an invalid pivot policy");
    auto limits = PlanarPressureRegionOpeningFlowLimits{};
    limits.maximumIntervals = 1;
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep, openings, {}, limits)); },
        "regional opening flow enforces its interval limit");
    limits = {};
    limits.maximumRegions = 1;
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep, openings, {}, limits)); },
        "regional opening flow enforces its region limit");
    limits = {};
    limits.maximumOpenings = 1;
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep,
                std::vector<PlanarPressureRegionOpeningDefinition>{
                    {700, 1, 2, 0.25},
                    {701, 1, 2, 0.25},
                }, {}, limits)); },
        "regional opening flow enforces its opening limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep, openings, {}, limits)); },
        "regional opening flow enforces its owned byte limit");
    limits = {};
    limits.maximumFactorizationBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                sweep, openings, {}, limits)); },
        "regional opening flow enforces its factorization byte limit");

    const std::vector<PlanarPressureJumpLayerDefinition> threeRegions{
        {30, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 0},
         -1.8, 30.0},
        {40, 2, 3,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 1, 0},
         -0.8, 20.0},
        {50, 3, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 2, 0},
         0.2, -50.0},
    };
    auto movedThreeRegions = threeRegions;
    movedThreeRegions[1].physicalPlaneCoordinateMeters += 0.1;
    const auto threeRegionSweep = makePlanarPressureRegionSweepLedger(
        geometry, threeRegions, movedThreeRegions, 0.5);
    const std::vector<PlanarPressureRegionOpeningDefinition> serialOpenings{
        {800, 3, 1, 0.25},
        {810, 1, 2, 0.5},
    };
    limits = {};
    limits.maximumFactorizationWork = 7;
    expectRejected(
        [&] { static_cast<void>(
            solvePlanarPressureRegionOpeningFlow(
                threeRegionSweep, serialOpenings, {}, limits)); },
        "regional opening flow enforces its factorization work limit");
}

void testRegionalOpeningPressurePower() {
    const auto geometry = grid();
    const auto previous = pocketLayers();
    const std::vector<PlanarPressureRegionOpeningDefinition> noOpenings;
    const auto rigidLayers = translatePlanarPressureJumpLayers(
        geometry, previous, 0.5).layers;
    const auto rigidSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, rigidLayers, 1.0);
    const auto rigidAllocation = solvePlanarPressureRegionOpeningFlow(
        rigidSweep, noOpenings);
    const auto rigid = auditPlanarPressureRegionOpeningPower(
        rigidSweep, noOpenings, rigidAllocation);
    check(rigid.version == planarPressureRegionOpeningPowerVersion
              && rigid.fingerprint != 0
              && rigid.sourceOpeningFlowFingerprint
                  == rigidAllocation.fingerprint
              && rigid.openings.empty()
              && rigid.regions.size() == 2
              && rigid.failedPassiveOpeningCount == 0
              && rigid.allOpeningsPassiveWithinTolerance
              && rigid.pressurePowerClosesWithinTolerance,
          "regional opening power accepts zero-power rigid translation");
    checkNear(rigid.totalOpeningPressurePowerWatts, 0.0, 1.0e-15,
              "rigid regional opening pressure power is zero");
    checkNear(rigid.totalRegionPressureVolumePowerWatts, 0.0, 1.0e-15,
              "rigid regional pressure-volume power is zero");

    auto inflatedLayers = previous;
    inflatedLayers[0].physicalPlaneCoordinateMeters -= 0.1;
    inflatedLayers[1].physicalPlaneCoordinateMeters += 0.1;
    const auto inflationSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, inflatedLayers, 0.5);
    const std::vector<PlanarPressureRegionOpeningDefinition> opening{
        {700, 1, 2, 0.5},
    };
    const auto inflationFlow = solvePlanarPressureRegionOpeningFlow(
        inflationSweep, opening);
    const auto inflation = auditPlanarPressureRegionOpeningPower(
        inflationSweep, opening, inflationFlow);
    check(inflation.openings.size() == 1
              && inflation.failedPassiveOpeningCount == 1
              && !inflation.allOpeningsPassiveWithinTolerance
              && inflation.pressurePowerClosesWithinTolerance,
          "regional opening power exposes uphill inflation flow");
    checkNear(
        inflation.openings[0]
            .negativeToPositivePressureDropPascals,
        -70.0, 1.0e-13,
        "inflation opening retains the static regional pressure rise");
    checkNear(inflation.totalOpeningPressurePowerWatts,
              -112.0, 3.0e-13,
              "inflation opening reports negative pressure power");
    checkNear(inflation.maximumOpeningExternalPowerWatts,
              112.0, 3.0e-13,
              "inflation opening reports its external power requirement");
    checkNear(inflation.minimumNetExternalPowerWatts,
              112.0, 3.0e-13,
              "inflation graph reports its net external power requirement");
    checkNear(inflation.totalRegionPressureVolumePowerWatts,
              112.0, 3.0e-13,
              "inflation regional pressure-volume power is equal and opposite");
    checkNear(inflation.pressurePowerClosureResidualWatts,
              0.0, 6.0e-13,
              "inflation pressure power closes against regional volume work");

    const std::vector<PlanarPressureRegionOpeningDefinition> reverseOpening{
        {700, 2, 1, 0.5},
    };
    const auto reverseFlow = solvePlanarPressureRegionOpeningFlow(
        inflationSweep, reverseOpening);
    const auto reverse = auditPlanarPressureRegionOpeningPower(
        inflationSweep, reverseOpening, reverseFlow);
    checkNear(reverse.totalOpeningPressurePowerWatts,
              -112.0, 3.0e-13,
              "opening pressure power is invariant to authored orientation");
    checkNear(reverse.minimumNetExternalPowerWatts,
              112.0, 3.0e-13,
              "orientation reversal preserves external power demand");

    auto deflatedLayers = previous;
    deflatedLayers[0].physicalPlaneCoordinateMeters += 0.1;
    deflatedLayers[1].physicalPlaneCoordinateMeters -= 0.1;
    const auto deflationSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, deflatedLayers, 0.5);
    const auto deflationFlow = solvePlanarPressureRegionOpeningFlow(
        deflationSweep, opening);
    const auto deflation = auditPlanarPressureRegionOpeningPower(
        deflationSweep, opening, deflationFlow);
    check(deflation.allOpeningsPassiveWithinTolerance
              && deflation.pressurePowerClosesWithinTolerance,
          "regional opening power accepts passive high-to-low deflation");
    checkNear(deflation.totalOpeningPressurePowerWatts,
              112.0, 3.0e-13,
              "deflation releases positive pressure power");
    checkNear(deflation.minimumNetExternalPowerWatts,
              0.0, 1.0e-15,
              "passive deflation needs no external pressure power");
    checkNear(deflation.totalRegionPressureVolumePowerWatts,
              -112.0, 3.0e-13,
              "deflation regional pressure-volume work absorbs released power");

    const std::vector<PlanarPressureRegionOpeningDefinition> parallel{
        {701, 1, 2, 0.375},
        {700, 1, 2, 0.125},
    };
    const auto parallelFlow = solvePlanarPressureRegionOpeningFlow(
        inflationSweep, parallel);
    const auto parallelPower = auditPlanarPressureRegionOpeningPower(
        inflationSweep, parallel, parallelFlow);
    check(parallelPower.failedPassiveOpeningCount == 2
              && parallelPower.openings[0].openingStableId == 700
              && parallelPower.openings[1].openingStableId == 701,
          "parallel regional opening power retains canonical split identities");
    checkNear(parallelPower.openings[0].pressurePowerWatts,
              -28.0, 1.0e-13,
              "small parallel aperture carries its pressure-power share");
    checkNear(parallelPower.openings[1].pressurePowerWatts,
              -84.0, 3.0e-13,
              "large parallel aperture carries its pressure-power share");
    checkNear(parallelPower.summedOpeningExternalPowerDeficitWatts,
              112.0, 3.0e-13,
              "parallel aperture external deficits sum to the graph demand");

    const std::vector<PlanarPressureJumpLayerDefinition> threeRegions{
        {30, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 0},
         -1.8, 30.0},
        {40, 2, 3,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 1, 0},
         -0.8, 20.0},
        {50, 3, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 2, 0},
         0.2, -50.0},
    };
    auto movedThreeRegions = threeRegions;
    movedThreeRegions[1].physicalPlaneCoordinateMeters += 0.1;
    const auto serialSweep = makePlanarPressureRegionSweepLedger(
        geometry, threeRegions, movedThreeRegions, 0.5);
    const std::vector<PlanarPressureRegionOpeningDefinition> serialOpenings{
        {800, 3, 1, 0.25},
        {810, 1, 2, 0.5},
    };
    const auto serialFlow = solvePlanarPressureRegionOpeningFlow(
        serialSweep, serialOpenings);
    const auto serial = auditPlanarPressureRegionOpeningPower(
        serialSweep, serialOpenings, serialFlow);
    check(serial.failedPassiveOpeningCount == 1
              && !serial.allOpeningsPassiveWithinTolerance
              && serial.pressurePowerClosesWithinTolerance,
          "serial opening power distinguishes local uphill and downhill links");
    checkNear(serial.openings[0].pressurePowerWatts,
              40.0, 2.0e-13,
              "serial high-to-low link supplies pressure power");
    checkNear(serial.openings[1].pressurePowerWatts,
              -24.0, 2.0e-13,
              "serial low-to-high link consumes pressure power");
    checkNear(serial.totalOpeningPressurePowerWatts,
              16.0, 3.0e-13,
              "serial opening graph retains its net positive pressure power");
    checkNear(serial.summedOpeningExternalPowerDeficitWatts,
              24.0, 2.0e-13,
              "serial opening graph reports its local uphill deficit");
    checkNear(serial.minimumNetExternalPowerWatts,
              0.0, 1.0e-15,
              "serial opening graph has no net external pressure-power deficit");
    checkNear(serial.totalRegionPressureVolumePowerWatts,
              -16.0, 3.0e-13,
              "serial regional volume work closes net opening power");

    expectRejected(
        [&] { static_cast<void>(
            auditPlanarPressureRegionOpeningPower(
                inflationSweep, noOpenings,
                solvePlanarPressureRegionOpeningFlow(
                    inflationSweep, noOpenings))); },
        "regional opening power rejects kinematically infeasible sources");

    auto corrupt = inflation;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningPower(
            corrupt, inflationSweep, opening, inflationFlow); },
        "regional opening-power validation rejects fingerprint corruption");
    corrupt = inflation;
    corrupt.openings[0].pressurePowerWatts += 1.0;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningPower(
            corrupt, inflationSweep, opening, inflationFlow); },
        "regional opening-power validation rejects opening corruption");
    corrupt = inflation;
    corrupt.regions[0].pressureVolumePowerWatts += 1.0;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningPower(
            corrupt, inflationSweep, opening, inflationFlow); },
        "regional opening-power validation rejects region corruption");
    corrupt = inflation;
    corrupt.minimumNetExternalPowerWatts += 1.0;
    expectRejected(
        [&] { validatePlanarPressureRegionOpeningPower(
            corrupt, inflationSweep, opening, inflationFlow); },
        "regional opening-power validation rejects aggregate corruption");

    auto invalidSettings = PlanarPressureRegionOpeningPowerSettings{};
    invalidSettings.absolutePowerToleranceWatts = 0.0;
    invalidSettings.relativePowerTolerance = 0.0;
    expectRejected(
        [&] { static_cast<void>(
            auditPlanarPressureRegionOpeningPower(
                inflationSweep, opening, inflationFlow,
                invalidSettings)); },
        "regional opening power rejects a zero tolerance policy");
    auto limits = PlanarPressureRegionOpeningPowerLimits{};
    limits.maximumRegions = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditPlanarPressureRegionOpeningPower(
                inflationSweep, opening, inflationFlow, {}, limits)); },
        "regional opening power enforces its region limit");
    limits = {};
    limits.maximumOpenings = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditPlanarPressureRegionOpeningPower(
                inflationSweep, parallel, parallelFlow, {}, limits)); },
        "regional opening power enforces its opening limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditPlanarPressureRegionOpeningPower(
                inflationSweep, opening, inflationFlow, {}, limits)); },
        "regional opening power enforces its owned byte limit");
    limits = {};
    limits.sourceLimits.maximumIntervals = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditPlanarPressureRegionOpeningPower(
                inflationSweep, opening, inflationFlow, {}, limits)); },
        "regional opening power enforces nested source limits");
}

void testPlanarRegionalControlFragments() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto repeated = buildPlanarPressureRegionFragments(
        geometry, sweep);
    check(fragments == repeated
              && fragments.version == planarPressureRegionFragmentVersion
              && fragments.fingerprint != 0
              && fragments.fragments.size() == 24
              && fragments.regions.size() == 2
              && fragments.cells.size() == 16
              && fragments.maximumFragmentsPerCell == 3
              && fragments.ownedStorageBytes > 0,
          "planar regional fragments preserve a deterministic same-cell pocket");
    checkNear(fragments.fragmentVolumeCubicMeters,
              16.0, 4.0e-14,
              "planar regional fragment volume closes the domain");
    checkNear(fragments.domainVolumeClosureResidualCubicMeters,
              0.0, 4.0e-14,
              "planar regional fragments retain zero domain-volume residual");
    checkNear(
        fragments.maximumAbsoluteCellVolumeClosureResidualCubicMeters,
        0.0, 4.0e-14,
        "planar regional fragments close every Cartesian cell volume");
    checkNear(
        fragments
            .maximumAbsoluteCellFirstMomentClosureResidualCubicMetersSquared,
        0.0, 4.0e-14,
        "planar regional fragments close every Cartesian cell centroid moment");

    const auto* exterior = findFragmentRegion(fragments, 1);
    const auto* pocket = findFragmentRegion(fragments, 2);
    check(exterior != nullptr && pocket != nullptr,
          "planar regional fragments retain both pressure regions");
    if (exterior != nullptr && pocket != nullptr) {
        check(exterior->fragmentCount == 20
                  && pocket->fragmentCount == 4,
              "planar regional fragments retain per-region control counts");
        checkNear(exterior->volumeCubicMeters,
                  13.6, 4.0e-14,
                  "planar regional fragments retain exterior volume");
        checkNear(pocket->volumeCubicMeters,
                  2.4, 4.0e-14,
                  "planar regional fragments retain thin-pocket volume");
        checkNear(pocket->pressurePascals - exterior->pressurePascals,
                  70.0, 1.0e-13,
                  "planar regional fragments retain regional pressure difference");
    }

    std::size_t tripleCells = 0;
    for (const auto& cell : fragments.cells) {
        checkNear(cell.fragmentVolumeCubicMeters,
                  1.0, 4.0e-14,
                  "planar regional fragment cell volume is analytic");
        if (cell.fragmentCount == 3) {
            ++tripleCells;
            check(cell.i == 1,
                  "same-cell pocket fragments occupy the authored axial cell");
        }
    }
    check(tripleCells == 4,
          "same-cell pocket creates three controls in every transverse tile");

    std::size_t pocketFragments = 0;
    std::set<std::uint64_t> stableIds;
    for (const auto& fragment : fragments.fragments) {
        stableIds.insert(fragment.stableId);
        if (fragment.regionStableId == 2) {
            ++pocketFragments;
            check(fragment.i == 1
                      && fragment.lowerBoundary.kind
                          == PlanarPressureRegionFragmentBoundaryKind::PressureLayer
                      && fragment.lowerBoundary.surfaceStableId == 10
                      && fragment.upperBoundary.kind
                          == PlanarPressureRegionFragmentBoundaryKind::PressureLayer
                      && fragment.upperBoundary.surfaceStableId == 20,
                  "thin-pocket controls retain both authored layer boundaries");
            checkNear(fragment.volumeCubicMeters,
                      0.6, 4.0e-15,
                      "thin-pocket control owns its physical subcell volume");
            checkNear(fragment.wrappedCentroidMeters.x,
                      -0.5, 4.0e-15,
                      "thin-pocket control owns its physical subcell centroid");
        }
    }
    check(pocketFragments == 4
              && stableIds.size() == fragments.fragments.size(),
          "planar regional fragments have unique stable control identities");
    validatePlanarPressureRegionFragments(
        fragments, geometry, sweep);
}

void testPlanarRegionalFragmentsAllAxesAndRejection() {
    const auto geometry = grid();
    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> layers{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0},
             -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0},
             -0.2, -70.0},
        };
        const auto sweep = makePlanarPressureRegionSweepLedger(
            geometry, layers, layers, 1.0);
        const auto fragments = buildPlanarPressureRegionFragments(
            geometry, sweep);
        const std::size_t expectedCount = axis == GridFaceAxis::X
            ? 24 : 32;
        const double expectedPocketVolume = axis == GridFaceAxis::X
            ? 2.4 : 4.8;
        const auto* pocket = findFragmentRegion(fragments, 2);
        check(fragments.axis == axis
                  && fragments.fragments.size() == expectedCount
                  && fragments.maximumFragmentsPerCell == 3
                  && pocket != nullptr,
              "planar regional fragments preserve same-cell topology on every axis");
        if (pocket != nullptr) {
            checkNear(pocket->volumeCubicMeters,
                      expectedPocketVolume, 8.0e-14,
                      "planar regional fragments use physical volume on every axis");
        }
        checkNear(
            fragments.maximumAbsoluteCellVolumeClosureResidualCubicMeters,
            0.0, 8.0e-14,
            "planar regional fragments close cell volume on every axis");
    }

    const auto layers = pocketLayers();
    const auto translated = translatePlanarPressureJumpLayers(
        geometry, layers, 0.5).layers;
    const auto movedSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, translated, 1.0);
    const auto moved = buildPlanarPressureRegionFragments(
        geometry, movedSweep);
    const auto* movedPocket = findFragmentRegion(moved, 2);
    check(movedPocket != nullptr,
          "translated planar regional fragments retain the pocket region");
    if (movedPocket != nullptr) {
        checkNear(movedPocket->volumeCubicMeters,
                  2.4, 4.0e-14,
                  "translated planar regional fragments preserve pocket volume");
    }
    checkNear(moved.domainVolumeClosureResidualCubicMeters,
              0.0, 4.0e-14,
              "translated planar regional fragments preserve domain volume");

    const std::vector<PlanarPressureJumpLayerDefinition> wrapping{
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 3, 0},
         1.2, 70.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 1},
         1.8, -70.0},
    };
    const auto wrapped = translatePlanarPressureJumpLayers(
        geometry, wrapping, 0.5).layers;
    const auto wrappedFragments = buildPlanarPressureRegionFragments(
        geometry,
        makePlanarPressureRegionSweepLedger(
            geometry, wrapping, wrapped, 1.0));
    const auto* wrappedPocket = findFragmentRegion(wrappedFragments, 2);
    check(wrappedPocket != nullptr,
          "periodically rebased planar fragments retain the pocket region");
    if (wrappedPocket != nullptr) {
        checkNear(wrappedPocket->volumeCubicMeters,
                  2.4, 4.0e-14,
                  "periodically rebased planar fragments preserve pocket volume");
    }
    checkNear(
        wrappedFragments
            .maximumAbsoluteCellVolumeClosureResidualCubicMeters,
        0.0, 4.0e-14,
        "periodically rebased planar fragments close every wrapped cell");

    const std::vector<PlanarPressureJumpLayerDefinition> reverseWrapping{
        {10, 1, 2,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 3, -1},
         -2.8, 70.0},
        {20, 2, 1,
         {movingPlanarFaceTopologyVersion, GridFaceAxis::X, 0, 0},
         -2.2, -70.0},
    };
    const auto reverseWrapped = translatePlanarPressureJumpLayers(
        geometry, reverseWrapping, -0.5).layers;
    const auto reverseWrappedFragments =
        buildPlanarPressureRegionFragments(
            geometry,
            makePlanarPressureRegionSweepLedger(
                geometry, reverseWrapping, reverseWrapped, 1.0));
    checkNear(
        reverseWrappedFragments
            .maximumAbsoluteCellVolumeClosureResidualCubicMeters,
        0.0, 4.0e-14,
        "negative periodically rebased planar fragments close every wrapped cell");

    auto corrupt = moved;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionFragments(
            corrupt, geometry, movedSweep); },
        "planar regional fragment validation rejects fingerprint corruption");
    corrupt = moved;
    corrupt.fragments[0].volumeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragments(
            corrupt, geometry, movedSweep); },
        "planar regional fragment validation rejects control corruption");
    corrupt = moved;
    corrupt.regions[0].volumeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragments(
            corrupt, geometry, movedSweep); },
        "planar regional fragment validation rejects region corruption");
    corrupt = moved;
    corrupt.cells[0].fragmentCount += 1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragments(
            corrupt, geometry, movedSweep); },
        "planar regional fragment validation rejects cell corruption");

    auto limits = PlanarPressureRegionFragmentLimits{};
    limits.maximumIntervals = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragments(
            geometry, movedSweep, limits)); },
        "planar regional fragments enforce the interval limit");
    limits = {};
    limits.maximumRegions = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragments(
            geometry, movedSweep, limits)); },
        "planar regional fragments enforce the region limit");
    limits = {};
    limits.maximumCells = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragments(
            geometry, movedSweep, limits)); },
        "planar regional fragments enforce the cell limit");
    limits = {};
    limits.maximumFragments = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragments(
            geometry, movedSweep, limits)); },
        "planar regional fragments enforce the control limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragments(
            geometry, movedSweep, limits)); },
        "planar regional fragments enforce the byte limit");
}

void testAllAxisAssembly() {
    const auto geometry = grid();
    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> layers{
            {
                10, 1, 2,
                {movingPlanarFaceTopologyVersion,
                 axis, firstFace, 0},
                -0.8, 70.0,
            },
            {
                20, 2, 1,
                {movingPlanarFaceTopologyVersion,
                 axis, secondFace, 0},
                -0.2, -70.0,
            },
        };
        const auto field = makePlanarPressureJumpField(
            geometry, layers);
        const auto counts = geometry.cellCounts();
        const std::size_t transverseTiles = axis == GridFaceAxis::X
            ? counts.y * counts.z
            : (axis == GridFaceAxis::Y
                ? counts.x * counts.z
                : counts.x * counts.y);
        check(field.faceCount() == 2 * transverseTiles,
              "planar jump assembly covers every transverse tile on each axis");
        const auto translated = translatePlanarPressureJumpLayers(
            geometry, layers, 0.5);
        check(translated.layers[0].topology.faceCoordinate == secondFace
                  && translated.layers[1].topology.faceCoordinate
                      == secondFace,
              "planar jump translation forms a same-face pocket on each axis");
        const auto collapsed = makePlanarPressureJumpField(
            geometry, translated.layers);
        const std::size_t index = geometry.cellIndex(
            axis == GridFaceAxis::X ? secondFace : 0,
            axis == GridFaceAxis::Y ? secondFace : 0,
            axis == GridFaceAxis::Z ? secondFace : 0);
        const double aggregate = axis == GridFaceAxis::X
            ? collapsed.xFaceJumpsPascals()[index]
            : (axis == GridFaceAxis::Y
                ? collapsed.yFaceJumpsPascals()[index]
                : collapsed.zFaceJumpsPascals()[index]);
        check(aggregate == 0.0,
              "same-face pocket retains its zero dense jump on each axis");
        const auto profile = makeStaticPlanarPressureRegionProfile(
            geometry, translated.layers);
        const double expectedPocketVolume = axis == GridFaceAxis::X
            ? 2.4 : 4.8;
        const auto* pocket = findRegion(profile, 2);
        check(pocket != nullptr,
              "static regional profile retains the pocket on each axis");
        if (pocket != nullptr) {
            checkNear(pocket->volumeCubicMeters,
                      expectedPocketVolume, 4.0e-15,
                      "static regional profile uses the physical transverse area on each axis");
        }
    }
}

void testTransactionalRejection() {
    const auto geometry = grid();
    const auto valid = pocketLayers();
    auto invalid = valid;
    invalid[1].minusRegionStableId = 3;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureJumpField(
            geometry, invalid)); },
        "planar jump assembly rejects a broken ordered region chain");
    invalid = valid;
    invalid[1].plusRegionStableId = 3;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureJumpField(
            geometry, invalid)); },
        "planar jump assembly rejects a nonperiodic outer region chain");
    invalid = valid;
    invalid[1].surfaceStableId = invalid[0].surfaceStableId;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureJumpField(
            geometry, invalid)); },
        "planar jump assembly rejects duplicate surface identity");
    invalid = valid;
    invalid[1].physicalPlaneCoordinateMeters =
        invalid[0].physicalPlaneCoordinateMeters;
    invalid[1].topology = invalid[0].topology;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureJumpField(
            geometry, invalid)); },
        "planar jump assembly rejects coincident physical layers");
    invalid = valid;
    invalid[1].physicalPlaneCoordinateMeters = 3.2;
    invalid[1].topology.periodicImage = 1;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureJumpField(
            geometry, invalid)); },
        "planar jump assembly rejects layers spanning a complete period");
    invalid = valid;
    invalid[1].topology.axis = GridFaceAxis::Y;
    expectRejected(
        [&] { static_cast<void>(makePlanarPressureJumpField(
            geometry, invalid)); },
        "planar jump assembly rejects mixed layer axes");

    expectRejected(
        [&] { static_cast<void>(translatePlanarPressureJumpLayers(
            geometry, valid, 0.3)); },
        "layer translation rejects exact topology boundaries");
    expectRejected(
        [&] { static_cast<void>(translatePlanarPressureJumpLayers(
            geometry, valid, 1.6)); },
        "layer translation rejects a skipped topology segment");
    expectRejected(
        [&] { static_cast<void>(translatePlanarPressureJumpLayers(
            geometry, valid,
            std::numeric_limits<double>::quiet_NaN())); },
        "layer translation rejects non-finite displacement");
    check(valid == pocketLayers(),
          "rejected layer translation cannot mutate the caller topology");

    auto reversed = valid;
    std::ranges::reverse(reversed);
    const auto first = translatePlanarPressureJumpLayers(
        geometry, valid, 0.5);
    const auto second = translatePlanarPressureJumpLayers(
        geometry, reversed, 0.5);
    check(first == second,
          "layer translation canonicalizes authored order deterministically");
}

} // namespace

int main() {
    testCanonicalAssemblyAndSameFacePocket();
    testSeparatedResolutionAndPeriodicWrap();
    testStaticRegionalPressureProfile();
    testRegionalSweepLedger();
    testRegionalSweepAllAxesAndRejection();
    testRegionalFluxCompatibility();
    testRegionalFluxAllAxesAndRejection();
    testRegionalOpeningFlowFeasibility();
    testRegionalOpeningFlowAllAxesAndRejection();
    testRegionalOpeningPressurePower();
    testPlanarRegionalControlFragments();
    testPlanarRegionalFragmentsAllAxesAndRejection();
    testAllAxisAssembly();
    testTransactionalRejection();
    if (failures != 0) {
        std::fprintf(stderr,
                     "%d planar pressure-jump check(s) failed\n",
                     failures);
        return 1;
    }
    std::puts("all planar pressure-jump checks passed");
    return 0;
}
