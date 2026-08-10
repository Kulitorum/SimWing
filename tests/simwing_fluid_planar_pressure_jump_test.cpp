#include "fluid/planar_pressure_jump.h"
#include "fluid/planar_region_sweep.h"
#include "fluid/projection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
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
