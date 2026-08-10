#include "fluid/planar_pressure_jump.h"
#include "fluid/planar_region_flux.h"
#include "fluid/planar_region_opening_flow.h"
#include "fluid/planar_region_opening_power.h"
#include "fluid/planar_region_fragment.h"
#include "fluid/planar_region_fragment_pressure_operator.h"
#include "fluid/planar_region_fragment_pressure_projection.h"
#include "fluid/planar_region_fragment_projection_energy.h"
#include "fluid/planar_region_fragment_pressure_solve.h"
#include "fluid/planar_region_fragment_topology.h"
#include "fluid/planar_region_fragment_velocity_metric.h"
#include "fluid/planar_region_fragment_velocity_state.h"
#include "fluid/planar_region_fragment_volume_rate.h"
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

const PlanarPressureRegionFragmentComponent* findFragmentComponent(
    const PlanarPressureRegionFragmentTopology& topology,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        topology.components, regionStableId,
        &PlanarPressureRegionFragmentComponent::regionStableId);
    return found == topology.components.end() ? nullptr : &*found;
}

const PlanarPressureRegionFragmentPressureOperatorComponent*
findFragmentOperatorComponent(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        pressureOperator.components, regionStableId,
        &PlanarPressureRegionFragmentPressureOperatorComponent::
            regionStableId);
    return found == pressureOperator.components.end() ? nullptr : &*found;
}

const PlanarPressureRegionFragmentRegionVolumeRate*
findFragmentVolumeRateRegion(
    const PlanarPressureRegionFragmentVolumeRateSet& rates,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        rates.regions, regionStableId,
        &PlanarPressureRegionFragmentRegionVolumeRate::regionStableId);
    return found == rates.regions.end() ? nullptr : &*found;
}

const PlanarPressureRegionFragmentComponentVolumeRate*
findFragmentVolumeRateComponent(
    const PlanarPressureRegionFragmentVolumeRateSet& rates,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        rates.components, regionStableId,
        &PlanarPressureRegionFragmentComponentVolumeRate::regionStableId);
    return found == rates.components.end() ? nullptr : &*found;
}

const PlanarPressureRegionFragmentVelocityMetricComponent*
findFragmentVelocityMetricComponent(
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        metric.components, regionStableId,
        &PlanarPressureRegionFragmentVelocityMetricComponent::regionStableId);
    return found == metric.components.end() ? nullptr : &*found;
}

const PlanarPressureRegionFragmentVelocityStateComponent*
findFragmentVelocityStateComponent(
    const PlanarPressureRegionFragmentVelocityState& state,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        state.components, regionStableId,
        &PlanarPressureRegionFragmentVelocityStateComponent::regionStableId);
    return found == state.components.end() ? nullptr : &*found;
}

std::vector<double> regionalMetricVelocityFromTopologyLinks(
    const PlanarPressureRegionFragmentVelocityMetric& metric,
    const std::vector<double>& linkVelocityMetersPerSecond,
    const PlanarPressureRegionFragmentVolumeRateSet* volumeRates = nullptr) {
    std::vector<double> result(metric.dofs.size(), 0.0);
    for (const auto& dof : metric.dofs) {
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                SharedRegionGrid) {
            result[dof.dofIndex] =
                linkVelocityMetersPerSecond[dof.sourceFaceLinkIndex];
        } else if (volumeRates != nullptr) {
            const auto& rate =
                volumeRates->fragments[dof.ownerFragmentIndex];
            result[dof.dofIndex] = dof.kind
                    == PlanarPressureRegionFragmentVelocityDofKind::
                        PressureLayerMinusTrace
                ? rate.upperBoundaryVelocityMetersPerSecond
                : rate.lowerBoundaryVelocityMetersPerSecond;
        }
    }
    return result;
}

double dotProduct(const std::vector<double>& first,
                  const std::vector<double>& second) {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result += first[index] * second[index];
    }
    return result;
}

double fragmentCorrectionVolumeMean(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentSet& fragments,
    const PlanarPressureRegionFragmentPressureOperatorComponent& component,
    const std::vector<double>& values) {
    double moment = 0.0;
    for (std::size_t offset = 0;
         offset < component.fragmentCount; ++offset) {
        const std::size_t index =
            pressureOperator.componentFragmentIndices[
                component.firstFragmentMember + offset];
        moment += fragments.fragments[index].volumeCubicMeters
            * values[index];
    }
    return moment / component.totalVolumeCubicMeters;
}

void subtractFragmentCorrectionVolumeMeans(
    const PlanarPressureRegionFragmentPressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentSet& fragments,
    std::vector<double>& values) {
    for (const auto& component : pressureOperator.components) {
        const double mean = fragmentCorrectionVolumeMean(
            pressureOperator, fragments, component, values);
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            values[pressureOperator.componentFragmentIndices[
                component.firstFragmentMember + offset]] -= mean;
        }
    }
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

void testPlanarRegionalFragmentTopology() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto repeated = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    check(topology == repeated
              && topology.version
                  == planarPressureRegionFragmentTopologyVersion
              && topology.fingerprint != 0
              && topology.sourceFragmentFingerprint
                  == fragments.fingerprint
              && topology.links.size() == 72
              && topology.fragments.size() == 24
              && topology.components.size() == 2
              && topology.sameRegionGridLinkCount == 64
              && topology.pressureLayerWallLinkCount == 8
              && topology.periodicGridLinkCount == 28
              && topology.ownedStorageBytes > 0,
          "planar regional fragment topology pairs the complete canonical graph");
    checkNear(topology.totalUniqueFaceAreaSquareMeters,
              56.0, 2.0e-13,
              "planar regional fragment topology retains unique face area");
    checkNear(topology.totalIncidentFaceAreaSquareMeters,
              112.0, 2.0e-13,
              "planar regional fragment topology retains two-sided face area");
    checkNear(topology.totalExpectedFragmentBoundaryAreaSquareMeters,
              112.0, 2.0e-13,
              "planar regional fragment topology closes analytic prism area");
    checkNear(
        topology
            .maximumAbsoluteFragmentBoundaryAreaClosureResidualSquareMeters,
        0.0, 2.0e-13,
        "planar regional fragment topology closes every prism boundary");

    const auto* exterior = findFragmentComponent(topology, 1);
    const auto* pocket = findFragmentComponent(topology, 2);
    check(exterior != nullptr && pocket != nullptr,
          "planar regional fragment topology retains both components");
    if (exterior != nullptr && pocket != nullptr) {
        check(exterior->fragmentCount == 20
                  && exterior->sameRegionGridLinkCount == 56
                  && pocket->fragmentCount == 4
                  && pocket->sameRegionGridLinkCount == 8,
              "planar regional fragment topology retains component graph counts");
        checkNear(exterior->volumeCubicMeters,
                  13.6, 4.0e-14,
                  "planar regional fragment topology retains exterior volume");
        checkNear(pocket->volumeCubicMeters,
                  2.4, 4.0e-14,
                  "planar regional fragment topology retains pocket volume");
    }

    std::size_t surface10Count = 0;
    std::size_t surface20Count = 0;
    std::set<std::uint64_t> stableIds;
    for (const auto& link : topology.links) {
        stableIds.insert(link.stableId);
        check(link.linkIndex < topology.links.size()
                  && link.minusFragmentIndex < fragments.fragments.size()
                  && link.plusFragmentIndex < fragments.fragments.size()
                  && link.centerDistanceMeters > 0.0,
              "planar regional fragment link retains bounded source identity");
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            check(link.surfaceStableId == 0
                      && link.minusRegionStableId
                          == link.plusRegionStableId
                      && link.minusComponentIndex
                          == link.plusComponentIndex
                      && link.pressureJumpPascals == 0.0
                      && link.sameRegionGeometryWeightMeters > 0.0,
                  "planar regional grid link connects only one pressure component");
        } else {
            check(link.sameRegionGeometryWeightMeters == 0.0
                      && link.minusComponentIndex
                          != link.plusComponentIndex
                      && link.i == 1
                      && link.areaSquareMeters == 1.0
                      && !link.crossesPeriodicBoundary,
                  "planar pressure-layer wall remains nonconductive and two-sided");
            checkNear(link.centerDistanceMeters,
                      0.4, 4.0e-15,
                      "planar pressure-layer wall retains control-center distance");
            if (link.surfaceStableId == 10) {
                ++surface10Count;
                checkNear(link.pressureJumpPascals,
                          70.0, 1.0e-13,
                          "first pressure-layer wall retains positive jump");
                checkNear(link.wrappedCentroidMeters.x,
                          -0.8, 4.0e-15,
                          "first pressure-layer wall retains physical position");
            } else if (link.surfaceStableId == 20) {
                ++surface20Count;
                checkNear(link.pressureJumpPascals,
                          -70.0, 1.0e-13,
                          "second pressure-layer wall retains negative jump");
                checkNear(link.wrappedCentroidMeters.x,
                          -0.2, 4.0e-15,
                          "second pressure-layer wall retains physical position");
            }
        }
    }
    check(surface10Count == 4 && surface20Count == 4
              && stableIds.size() == topology.links.size(),
          "planar regional fragment topology retains every unique layer tile");
    for (const auto& summary : topology.fragments) {
        check(summary.incidentFaceCount == 6
                  && summary.sameRegionGridFaceCount
                          + summary.pressureLayerWallFaceCount
                      == 6,
              "planar regional fragment topology closes six faces per prism");
        if (summary.regionStableId == 2) {
            check(summary.sameRegionGridFaceCount == 4
                      && summary.pressureLayerWallFaceCount == 2,
                  "thin-pocket controls retain four grid faces and two walls");
        }
    }
    validatePlanarPressureRegionFragmentTopology(
        topology, geometry, sweep, fragments);

    const auto translatedLayers = translatePlanarPressureJumpLayers(
        geometry, layers, 0.1).layers;
    const auto translatedSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, translatedLayers, 1.0);
    const auto translatedFragments = buildPlanarPressureRegionFragments(
        geometry, translatedSweep);
    const auto translatedTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, translatedSweep, translatedFragments);
    std::set<std::uint64_t> translatedLinkIds;
    for (const auto& link : translatedTopology.links) {
        translatedLinkIds.insert(link.stableId);
    }
    check(translatedLinkIds == stableIds
              && translatedTopology.components[0].stableId
                  == topology.components[0].stableId
              && translatedTopology.components[1].stableId
                  == topology.components[1].stableId,
          "within-segment layer motion preserves graph and component identity");
}

void testPlanarRegionalFragmentTopologyAxesAndRejection() {
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
        const auto topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        const auto counts = geometry.cellCounts();
        const std::size_t transverseTiles = axis == GridFaceAxis::X
            ? counts.y * counts.z
            : (axis == GridFaceAxis::Y
                ? counts.x * counts.z
                : counts.x * counts.y);
        check(topology.profileAxis == axis
                  && topology.links.size() == 3 * fragments.fragments.size()
                  && topology.pressureLayerWallLinkCount
                      == 2 * transverseTiles
                  && topology.components.size() == 2
                  && topology.periodicGridLinkCount > 0,
              "planar regional fragment topology closes every axis");
        checkNear(
            topology
                .maximumAbsoluteFragmentBoundaryAreaClosureResidualSquareMeters,
            0.0, 3.0e-13,
            "planar regional fragment topology closes area on every axis");
        validatePlanarPressureRegionFragmentTopology(
            topology, geometry, sweep, fragments);
    }

    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);

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
    const auto wrappedSweep = makePlanarPressureRegionSweepLedger(
        geometry, wrapping, wrapped, 1.0);
    const auto wrappedFragments = buildPlanarPressureRegionFragments(
        geometry, wrappedSweep);
    const auto wrappedTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, wrappedSweep, wrappedFragments);
    check(wrappedTopology.links.size()
                  == 3 * wrappedFragments.fragments.size()
              && wrappedTopology.components.size() == 2,
          "positive periodic rebase preserves the regional fragment graph");

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
    const auto reverseSweep = makePlanarPressureRegionSweepLedger(
        geometry, reverseWrapping, reverseWrapped, 1.0);
    const auto reverseFragments = buildPlanarPressureRegionFragments(
        geometry, reverseSweep);
    const auto reverseTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, reverseSweep, reverseFragments);
    check(reverseTopology.links.size()
                  == 3 * reverseFragments.fragments.size()
              && reverseTopology.components.size() == 2,
          "negative periodic rebase preserves the regional fragment graph");

    auto corrupt = topology;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentTopology(
            corrupt, geometry, sweep, fragments); },
        "regional fragment topology rejects fingerprint corruption");
    corrupt = topology;
    corrupt.links[0].areaSquareMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentTopology(
            corrupt, geometry, sweep, fragments); },
        "regional fragment topology rejects face-link corruption");
    corrupt = topology;
    corrupt.fragments[0].incidentFaceCount += 1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentTopology(
            corrupt, geometry, sweep, fragments); },
        "regional fragment topology rejects incidence corruption");
    corrupt = topology;
    corrupt.components[0].volumeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentTopology(
            corrupt, geometry, sweep, fragments); },
        "regional fragment topology rejects component corruption");

    auto limits = PlanarPressureRegionFragmentTopologyLimits{};
    limits.maximumLinks = topology.links.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentTopology(
                geometry, sweep, fragments, limits)); },
        "regional fragment topology enforces the link limit");
    limits = {};
    limits.maximumComponents = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentTopology(
                geometry, sweep, fragments, limits)); },
        "regional fragment topology enforces the component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentTopology(
                geometry, sweep, fragments, limits)); },
        "regional fragment topology enforces the byte limit");
    limits = {};
    limits.fragmentLimits.maximumFragments = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentTopology(
                geometry, sweep, fragments, limits)); },
        "regional fragment topology enforces nested source limits");
    auto corruptSource = fragments;
    corruptSource.fragments[0].stableId ^= 1U;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentTopology(
                geometry, sweep, corruptSource)); },
        "regional fragment topology rejects mutated source geometry");
}

void testPlanarRegionalFragmentVolumeRates() {
    const auto geometry = grid();
    const auto previousLayers = pocketLayers();
    auto breathingLayers = previousLayers;
    breathingLayers[0].physicalPlaneCoordinateMeters -= 0.1;
    breathingLayers[1].physicalPlaneCoordinateMeters += 0.1;
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, previousLayers, breathingLayers, 0.5);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto rates = buildPlanarPressureRegionFragmentVolumeRates(
        geometry, sweep, fragments, topology);
    const auto repeated = buildPlanarPressureRegionFragmentVolumeRates(
        geometry, sweep, fragments, topology);
    check(rates == repeated
              && rates.version
                  == planarPressureRegionFragmentVolumeRateVersion
              && rates.fingerprint != 0
              && rates.sourceFragmentFingerprint == fragments.fingerprint
              && rates.sourceTopologyFingerprint == topology.fingerprint
              && rates.sourceSweepVersion == sweep.version
              && rates.axis == GridFaceAxis::X
              && rates.durationSeconds == 0.5
              && rates.topologyStable
              && rates.fragments.size() == 24
              && rates.cells.size() == 16
              && rates.regions.size() == 2
              && rates.components.size() == 2
              && rates.ownedStorageBytes > 0,
          "regional fragment volume rates publish a deterministic local ledger");
    checkNear(rates.maximumAbsoluteFragmentVolumeChangeCubicMeters,
              0.2, 3.0e-16,
              "breathing volume rates retain the largest local volume change");
    checkNear(rates.maximumAbsoluteFragmentVolumeRateCubicMetersPerSecond,
              0.4, 6.0e-16,
              "breathing volume rates retain the largest local dV/dt");
    checkNear(rates.maximumAbsoluteCellClosureResidualCubicMeters,
              0.0, 3.0e-16,
              "breathing fragment rates close every fixed Cartesian cell");
    checkNear(rates.maximumAbsoluteRegionClosureResidualCubicMeters,
              0.0, 3.0e-15,
              "breathing fragment rates close the source region sweep");
    checkNear(
        rates.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond,
        1.6, 3.0e-15,
        "breathing fragment rates retain component volume demand");
    checkNear(rates.globalGeometryVolumeChangeCubicMeters,
              0.0, 3.0e-15,
              "breathing fragment rates conserve the periodic domain");
    checkNear(rates.globalGeometryVolumeChangeRateCubicMetersPerSecond,
              0.0, 6.0e-15,
              "breathing fragment dV/dt cancels globally");

    const auto* pocket = findFragmentVolumeRateRegion(rates, 2);
    const auto* exterior = findFragmentVolumeRateRegion(rates, 1);
    const auto* pocketComponent = findFragmentVolumeRateComponent(rates, 2);
    const auto* exteriorComponent = findFragmentVolumeRateComponent(rates, 1);
    check(pocket != nullptr && exterior != nullptr
              && pocketComponent != nullptr && exteriorComponent != nullptr,
          "regional fragment volume rates retain both pressure components");
    if (pocket != nullptr && exterior != nullptr
        && pocketComponent != nullptr && exteriorComponent != nullptr) {
        checkNear(pocket->previousVolumeCubicMeters,
                  2.4, 4.0e-14,
                  "breathing pocket rates reconstruct previous volume");
        checkNear(pocket->currentVolumeCubicMeters,
                  3.2, 4.0e-14,
                  "breathing pocket rates retain current volume");
        checkNear(pocket->geometryVolumeChangeCubicMeters,
                  0.8, 4.0e-15,
                  "breathing pocket fragment changes sum exactly");
        checkNear(pocket->geometryVolumeChangeRateCubicMetersPerSecond,
                  1.6, 8.0e-15,
                  "breathing pocket fragment rates sum exactly");
        checkNear(exterior->geometryVolumeChangeCubicMeters,
                  -0.8, 4.0e-15,
                  "breathing exterior fragment changes oppose the pocket");
        checkNear(
            exterior->geometryVolumeChangeRateCubicMetersPerSecond,
            -1.6, 8.0e-15,
            "breathing exterior fragment rates oppose the pocket");
        checkNear(pocketComponent->geometryVolumeChangeCubicMeters,
                  pocket->geometryVolumeChangeCubicMeters, 0.0,
                  "pocket component and region share one volume ledger");
        checkNear(exteriorComponent->geometryVolumeChangeCubicMeters,
                  exterior->geometryVolumeChangeCubicMeters, 0.0,
                  "exterior component and region share one volume ledger");
    }
    std::size_t pocketFragmentCount = 0;
    for (const auto& rate : rates.fragments) {
        checkNear(rate.currentVolumeCubicMeters
                      - rate.previousVolumeCubicMeters,
                  rate.geometryVolumeChangeCubicMeters, 0.0,
                  "local fragment volume rate retains exact endpoint identity");
        checkNear(rate.geometryVolumeChangeRateCubicMetersPerSecond,
                  rate.geometryVolumeChangeCubicMeters / 0.5, 0.0,
                  "local fragment dV/dt uses the exact epoch duration");
        if (rate.regionStableId == 2) {
            ++pocketFragmentCount;
            checkNear(rate.lowerBoundaryDisplacementMeters,
                      -0.1, 3.0e-16,
                      "pocket control retains lower-wall displacement");
            checkNear(rate.upperBoundaryDisplacementMeters,
                      0.1, 3.0e-16,
                      "pocket control retains upper-wall displacement");
            checkNear(rate.lowerBoundaryVelocityMetersPerSecond,
                      -0.2, 6.0e-16,
                      "pocket control retains lower-wall velocity");
            checkNear(rate.upperBoundaryVelocityMetersPerSecond,
                      0.2, 6.0e-16,
                      "pocket control retains upper-wall velocity");
            checkNear(rate.previousVolumeCubicMeters,
                      0.6, 2.0e-15,
                      "pocket control reconstructs its previous local volume");
            checkNear(rate.currentVolumeCubicMeters,
                      0.8, 2.0e-15,
                      "pocket control retains its expanded local volume");
        }
    }
    check(pocketFragmentCount == 4,
          "breathing volume rates retain all four pocket controls");
    validatePlanarPressureRegionFragmentVolumeRates(
        rates, geometry, sweep, fragments, topology);

    const auto staticSweep = makePlanarPressureRegionSweepLedger(
        geometry, previousLayers, previousLayers, 1.0);
    const auto staticFragments = buildPlanarPressureRegionFragments(
        geometry, staticSweep);
    const auto staticTopology = buildPlanarPressureRegionFragmentTopology(
        geometry, staticSweep, staticFragments);
    const auto staticRates = buildPlanarPressureRegionFragmentVolumeRates(
        geometry, staticSweep, staticFragments, staticTopology);
    check(staticRates.maximumAbsoluteFragmentVolumeChangeCubicMeters == 0.0
              && staticRates
                      .maximumAbsoluteFragmentVolumeRateCubicMetersPerSecond
                  == 0.0
              && staticRates
                      .maximumAbsoluteComponentVolumeRateCubicMetersPerSecond
                  == 0.0,
          "static regional fragments publish exact zero geometry rates");
    for (const auto& rate : staticRates.fragments) {
        check(rate.previousVolumeCubicMeters == rate.currentVolumeCubicMeters,
              "static regional fragments preserve exact endpoint volume");
    }

    const auto translatedLayers = translatePlanarPressureJumpLayers(
        geometry, previousLayers, 0.1).layers;
    const auto translatedSweep = makePlanarPressureRegionSweepLedger(
        geometry, previousLayers, translatedLayers, 1.0);
    const auto translatedFragments = buildPlanarPressureRegionFragments(
        geometry, translatedSweep);
    const auto translatedTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, translatedSweep, translatedFragments);
    const auto translatedRates =
        buildPlanarPressureRegionFragmentVolumeRates(
            geometry, translatedSweep, translatedFragments,
            translatedTopology);
    checkNear(
        translatedRates.maximumAbsoluteFragmentVolumeChangeCubicMeters,
        0.1, 3.0e-16,
        "rigid translation retains local exchange across fixed grid faces");
    checkNear(
        translatedRates.maximumAbsoluteComponentVolumeRateCubicMetersPerSecond,
        0.0, 3.0e-15,
        "rigid translation preserves every component volume");
    checkNear(translatedRates.maximumAbsoluteCellClosureResidualCubicMeters,
              0.0, 3.0e-16,
              "rigid translation closes every fixed cell");
}

void testPlanarRegionalFragmentVolumeRatesAxesAndRejection() {
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
        const auto sweep = makePlanarPressureRegionSweepLedger(
            geometry, previous, current, 0.5);
        const auto fragments = buildPlanarPressureRegionFragments(
            geometry, sweep);
        const auto topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        const auto rates = buildPlanarPressureRegionFragmentVolumeRates(
            geometry, sweep, fragments, topology);
        const double crossSection = axis == GridFaceAxis::X ? 4.0 : 8.0;
        const auto* pocket = findFragmentVolumeRateRegion(rates, 2);
        check(rates.axis == axis && rates.topologyStable
                  && pocket != nullptr,
              "regional fragment volume rates retain every profile axis");
        if (pocket != nullptr) {
            checkNear(pocket->geometryVolumeChangeCubicMeters,
                      0.2 * crossSection, 8.0e-15,
                      "all-axis pocket volume change scales by cross-section");
            checkNear(
                pocket->geometryVolumeChangeRateCubicMetersPerSecond,
                0.4 * crossSection, 1.6e-14,
                "all-axis pocket dV/dt scales by cross-section");
        }
        checkNear(rates.maximumAbsoluteCellClosureResidualCubicMeters,
                  0.0, 4.0e-16,
                  "all-axis fragment rates close each fixed cell");
        checkNear(rates.globalGeometryVolumeChangeCubicMeters,
                  0.0, 8.0e-15,
                  "all-axis fragment rates conserve domain volume");
    }

    const auto previous = pocketLayers();
    const auto stableCurrent = translatePlanarPressureJumpLayers(
        geometry, previous, 0.1).layers;
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, stableCurrent, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto rates = buildPlanarPressureRegionFragmentVolumeRates(
        geometry, sweep, fragments, topology);

    const auto rebasedCurrent = translatePlanarPressureJumpLayers(
        geometry, previous, 0.5).layers;
    const auto rebasedSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, rebasedCurrent, 1.0);
    const auto rebasedFragments = buildPlanarPressureRegionFragments(
        geometry, rebasedSweep);
    const auto rebasedTopology = buildPlanarPressureRegionFragmentTopology(
        geometry, rebasedSweep, rebasedFragments);
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, rebasedSweep, rebasedFragments,
                rebasedTopology)); },
        "regional fragment volume rates reject a Cartesian topology rebase");

    auto corrupt = rates;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVolumeRates(
            corrupt, geometry, sweep, fragments, topology); },
        "regional fragment volume rates reject fingerprint corruption");
    corrupt = rates;
    corrupt.fragments[0].previousVolumeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVolumeRates(
            corrupt, geometry, sweep, fragments, topology); },
        "regional fragment volume rates reject control corruption");
    corrupt = rates;
    corrupt.cells[0].fragmentCount += 1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVolumeRates(
            corrupt, geometry, sweep, fragments, topology); },
        "regional fragment volume rates reject cell corruption");
    corrupt = rates;
    corrupt.regions[0].geometryVolumeChangeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVolumeRates(
            corrupt, geometry, sweep, fragments, topology); },
        "regional fragment volume rates reject region corruption");
    corrupt = rates;
    corrupt.components[0].geometryVolumeChangeRateCubicMetersPerSecond += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVolumeRates(
            corrupt, geometry, sweep, fragments, topology); },
        "regional fragment volume rates reject component corruption");

    auto limits = PlanarPressureRegionFragmentVolumeRateLimits{};
    limits.maximumFragments = rates.fragments.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology, limits)); },
        "regional fragment volume rates enforce the fragment limit");
    limits = {};
    limits.maximumCells = rates.cells.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology, limits)); },
        "regional fragment volume rates enforce the cell limit");
    limits = {};
    limits.maximumRegions = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology, limits)); },
        "regional fragment volume rates enforce the region limit");
    limits = {};
    limits.maximumComponents = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology, limits)); },
        "regional fragment volume rates enforce the component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology, limits)); },
        "regional fragment volume rates enforce the byte limit");
    limits = {};
    limits.topologyLimits.maximumLinks = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology, limits)); },
        "regional fragment volume rates enforce nested topology limits");
    auto corruptTopology = topology;
    corruptTopology.links[0].stableId ^= 1U;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, corruptTopology)); },
        "regional fragment volume rates reject mutated source topology");
}

void testPlanarRegionalFragmentVelocityMetric() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    const auto repeated = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    check(metric == repeated
              && metric.version
                  == planarPressureRegionFragmentVelocityMetricVersion
              && metric.fingerprint != 0
              && metric.sourceFragmentFingerprint == fragments.fingerprint
              && metric.sourceTopologyFingerprint == topology.fingerprint
              && metric.profileAxis == GridFaceAxis::X
              && metric.dofs.size() == 80
              && metric.fragments.size() == 24
              && metric.components.size() == 2
              && metric.sharedRegionGridDofCount == 64
              && metric.pressureLayerTraceDofCount == 16
              && metric.ownedStorageBytes > 0,
          "regional velocity metric publishes the complete canonical basis");
    checkNear(metric.sharedRegionGridDualVolumeCubicMeters,
              44.8, 2.0e-13,
              "regional velocity metric retains shared grid dual volume");
    checkNear(metric.pressureLayerTraceDualVolumeCubicMeters,
              3.2, 2.0e-14,
              "regional velocity metric retains both wall-side half volumes");
    checkNear(metric.totalDualVolumeCubicMeters,
              48.0, 2.0e-13,
              "regional velocity metric closes three component volumes");
    for (const double value : {
             metric.dualVolumeByAxisCubicMeters.x,
             metric.dualVolumeByAxisCubicMeters.y,
             metric.dualVolumeByAxisCubicMeters.z}) {
        checkNear(value, 16.0, 8.0e-14,
                  "regional velocity metric closes domain volume per axis");
    }
    checkNear(
        std::max({
            std::abs(
                metric.domainVolumeClosureResidualByAxisCubicMeters.x),
            std::abs(
                metric.domainVolumeClosureResidualByAxisCubicMeters.y),
            std::abs(
                metric.domainVolumeClosureResidualByAxisCubicMeters.z)}),
        0.0, 8.0e-14,
        "regional velocity metric diagnoses bounded domain closure");
    checkNear(
        metric.maximumAbsoluteFragmentVolumeClosureResidualCubicMeters,
        0.0, 4.0e-15,
        "regional velocity metric closes every fragment on all axes");
    checkNear(
        metric.maximumAbsoluteComponentVolumeClosureResidualCubicMeters,
        0.0, 8.0e-14,
        "regional velocity metric closes every component on all axes");

    const auto* exterior = findFragmentVelocityMetricComponent(metric, 1);
    const auto* pocket = findFragmentVelocityMetricComponent(metric, 2);
    check(exterior != nullptr && pocket != nullptr,
          "regional velocity metric retains both pressure components");
    if (exterior != nullptr && pocket != nullptr) {
        check(exterior->sharedGridDofCount == 56
                  && exterior->pressureLayerTraceDofCount == 8
                  && pocket->sharedGridDofCount == 8
                  && pocket->pressureLayerTraceDofCount == 8,
              "regional velocity metric splits grid and wall-trace ownership");
        for (const double value : {
                 exterior->dualVolumeByAxisCubicMeters.x,
                 exterior->dualVolumeByAxisCubicMeters.y,
                 exterior->dualVolumeByAxisCubicMeters.z}) {
            checkNear(value, 13.6, 6.0e-14,
                      "exterior velocity metric closes each axis");
        }
        for (const double value : {
                 pocket->dualVolumeByAxisCubicMeters.x,
                 pocket->dualVolumeByAxisCubicMeters.y,
                 pocket->dualVolumeByAxisCubicMeters.z}) {
            checkNear(value, 2.4, 2.0e-14,
                      "pocket velocity metric closes each axis");
        }
    }

    std::set<std::uint64_t> stableIds;
    std::size_t surface10TraceCount = 0;
    std::size_t surface20TraceCount = 0;
    for (const auto& dof : metric.dofs) {
        stableIds.insert(dof.stableId);
        check(dof.dofIndex < metric.dofs.size()
                  && dof.sourceFaceLinkIndex < topology.links.size()
                  && dof.ownerFragmentIndex < fragments.fragments.size()
                  && dof.oppositeFragmentIndex < fragments.fragments.size()
                  && dof.areaSquareMeters > 0.0
                  && dof.ownerHalfDistanceMeters > 0.0
                  && dof.ownerDualVolumeCubicMeters > 0.0
                  && dof.dualVolumeCubicMeters > 0.0,
              "regional velocity DOF retains bounded geometric ownership");
        const auto& source = topology.links[dof.sourceFaceLinkIndex];
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::SharedRegionGrid) {
            check(dof.surfaceStableId == 0
                      && dof.oppositeHalfDistanceMeters > 0.0
                      && dof.oppositeDualVolumeCubicMeters > 0.0,
                  "shared grid velocity owns both adjacent half volumes");
            checkNear(dof.dualVolumeCubicMeters,
                      source.areaSquareMeters
                          * source.centerDistanceMeters,
                      4.0e-15,
                      "shared velocity dual is area times center distance");
        } else {
            check(dof.surfaceStableId != 0
                      && dof.oppositeHalfDistanceMeters == 0.0
                      && dof.oppositeDualVolumeCubicMeters == 0.0,
                  "wall trace velocity owns only one fluid-side half volume");
            checkNear(dof.dualVolumeCubicMeters,
                      dof.areaSquareMeters
                          * dof.ownerHalfDistanceMeters,
                      0.0,
                      "wall trace dual is its one-sided area times distance");
            if (dof.surfaceStableId == 10) ++surface10TraceCount;
            if (dof.surfaceStableId == 20) ++surface20TraceCount;
        }
    }
    check(stableIds.size() == metric.dofs.size()
              && surface10TraceCount == 8
              && surface20TraceCount == 8,
          "regional velocity metric retains unique shared and wall trace IDs");
    for (const auto& fragment : metric.fragments) {
        check(fragment.velocityDofIncidenceCount == 6,
              "regional velocity metric assigns six incidences per fragment");
        for (const double value : {
                 fragment.dualVolumeByAxisCubicMeters.x,
                 fragment.dualVolumeByAxisCubicMeters.y,
                 fragment.dualVolumeByAxisCubicMeters.z}) {
            checkNear(value, fragment.sourceVolumeCubicMeters, 4.0e-15,
                      "fragment velocity metric recovers its volume per axis");
        }
    }
    validatePlanarPressureRegionFragmentVelocityMetric(
        metric, geometry, sweep, fragments, topology);

    const auto translatedLayers = translatePlanarPressureJumpLayers(
        geometry, layers, 0.1).layers;
    const auto translatedSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, translatedLayers, 1.0);
    const auto translatedFragments = buildPlanarPressureRegionFragments(
        geometry, translatedSweep);
    const auto translatedTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, translatedSweep, translatedFragments);
    const auto translatedMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, translatedSweep, translatedFragments,
            translatedTopology);
    std::set<std::uint64_t> translatedIds;
    for (const auto& dof : translatedMetric.dofs) {
        translatedIds.insert(dof.stableId);
    }
    check(translatedIds == stableIds
              && translatedMetric.fingerprint != metric.fingerprint,
          "within-segment motion preserves velocity identity and updates metric");
    checkNear(translatedMetric.totalDualVolumeCubicMeters,
              48.0, 2.0e-13,
              "translated velocity metric preserves total dual volume");

    auto breathingLayers = layers;
    breathingLayers[0].physicalPlaneCoordinateMeters -= 0.1;
    breathingLayers[1].physicalPlaneCoordinateMeters += 0.1;
    const auto breathingSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, breathingLayers, 0.5);
    const auto breathingFragments = buildPlanarPressureRegionFragments(
        geometry, breathingSweep);
    const auto breathingTopology = buildPlanarPressureRegionFragmentTopology(
        geometry, breathingSweep, breathingFragments);
    const auto breathingMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, breathingSweep, breathingFragments,
            breathingTopology);
    checkNear(breathingMetric.pressureLayerTraceDualVolumeCubicMeters,
              3.6, 3.0e-14,
              "breathing geometry updates its one-sided wall inertia");
    checkNear(breathingMetric.totalDualVolumeCubicMeters,
              48.0, 2.0e-13,
              "breathing velocity metric preserves total dual volume");
}

void testPlanarRegionalFragmentVelocityMetricAxesAndRejection() {
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
        const auto topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        const auto metric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology);
        check(metric.profileAxis == axis
                  && metric.sharedRegionGridDofCount
                      == topology.sameRegionGridLinkCount
                  && metric.pressureLayerTraceDofCount
                      == 2 * topology.pressureLayerWallLinkCount,
              "regional velocity metric retains all-axis DOF ownership");
        checkNear(metric.totalDualVolumeCubicMeters,
                  48.0, 3.0e-13,
                  "all-axis regional velocity metric closes total inertia");
        for (const double value : {
                 metric.dualVolumeByAxisCubicMeters.x,
                 metric.dualVolumeByAxisCubicMeters.y,
                 metric.dualVolumeByAxisCubicMeters.z}) {
            checkNear(value, 16.0, 1.2e-13,
                      "all-axis regional velocity metric closes domain volume");
        }
        for (const auto& component : metric.components) {
            for (const double value : {
                     component.dualVolumeByAxisCubicMeters.x,
                     component.dualVolumeByAxisCubicMeters.y,
                     component.dualVolumeByAxisCubicMeters.z}) {
                checkNear(value, component.sourceVolumeCubicMeters, 8.0e-14,
                          "all-axis component velocity metric closes volume");
            }
        }
    }

    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    auto corrupt = metric;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityMetric(
            corrupt, geometry, sweep, fragments, topology); },
        "regional velocity metric rejects fingerprint corruption");
    corrupt = metric;
    corrupt.dofs[0].dualVolumeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityMetric(
            corrupt, geometry, sweep, fragments, topology); },
        "regional velocity metric rejects DOF corruption");
    corrupt = metric;
    corrupt.fragments[0].velocityDofIncidenceCount += 1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityMetric(
            corrupt, geometry, sweep, fragments, topology); },
        "regional velocity metric rejects fragment corruption");
    corrupt = metric;
    corrupt.components[0].dualVolumeByAxisCubicMeters.x += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityMetric(
            corrupt, geometry, sweep, fragments, topology); },
        "regional velocity metric rejects component corruption");

    auto limits = PlanarPressureRegionFragmentVelocityMetricLimits{};
    limits.maximumDofs = metric.dofs.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology, limits)); },
        "regional velocity metric enforces the DOF limit");
    limits = {};
    limits.maximumFragments = metric.fragments.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology, limits)); },
        "regional velocity metric enforces the fragment limit");
    limits = {};
    limits.maximumComponents = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology, limits)); },
        "regional velocity metric enforces the component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology, limits)); },
        "regional velocity metric enforces the byte limit");
    limits = {};
    limits.topologyLimits.maximumLinks = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology, limits)); },
        "regional velocity metric enforces nested topology limits");
    auto corruptTopology = topology;
    corruptTopology.links[0].areaSquareMeters += 0.1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, corruptTopology)); },
        "regional velocity metric rejects mutated source topology");
}

void testPlanarRegionalFragmentVelocityState() {
    constexpr double density = 1.25;
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    std::vector<double> velocity(metric.dofs.size(), 0.0);
    for (const auto& dof : metric.dofs) {
        switch (dof.axis) {
        case GridFaceAxis::X:
            velocity[dof.dofIndex] = 2.0;
            break;
        case GridFaceAxis::Y:
            velocity[dof.dofIndex] = -0.5;
            break;
        case GridFaceAxis::Z:
            velocity[dof.dofIndex] = 0.25;
            break;
        }
    }
    const auto state = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, velocity, density);
    const auto repeated = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, velocity, density);
    check(state == repeated
              && state.version
                  == planarPressureRegionFragmentVelocityStateVersion
              && state.fingerprint != 0
              && state.sourceMetricFingerprint == metric.fingerprint
              && state.sourceFragmentFingerprint == fragments.fingerprint
              && state.sourceTopologyFingerprint == topology.fingerprint
              && state.profileAxis == GridFaceAxis::X
              && state.densityKgPerCubicMeter == density
              && state.samples.size() == metric.dofs.size()
              && state.fragments.size() == fragments.fragments.size()
              && state.components.size() == topology.components.size()
              && state.sharedRegionGridSampleCount == 64
              && state.pressureLayerTraceSampleCount == 16
              && state.maximumAbsoluteVelocityMetersPerSecond == 2.0
              && state.ownedStorageBytes > 0,
          "regional velocity state publishes its complete diagonal inertia state");
    checkNear(state.massByAxisKilograms.x, 20.0, 2.0e-13,
              "regional velocity state closes X mass");
    checkNear(state.massByAxisKilograms.y, 20.0, 2.0e-13,
              "regional velocity state closes Y mass");
    checkNear(state.massByAxisKilograms.z, 20.0, 2.0e-13,
              "regional velocity state closes Z mass");
    checkNear(state.momentumKilogramMetersPerSecond.x, 40.0, 4.0e-13,
              "regional velocity state closes X momentum");
    checkNear(state.momentumKilogramMetersPerSecond.y, -10.0, 2.0e-13,
              "regional velocity state closes Y momentum");
    checkNear(state.momentumKilogramMetersPerSecond.z, 5.0, 1.0e-13,
              "regional velocity state closes Z momentum");
    checkNear(state.kineticEnergyByAxisJoules.x, 40.0, 4.0e-13,
              "regional velocity state closes X kinetic energy");
    checkNear(state.kineticEnergyByAxisJoules.y, 2.5, 5.0e-14,
              "regional velocity state closes Y kinetic energy");
    checkNear(state.kineticEnergyByAxisJoules.z, 0.625, 2.0e-14,
              "regional velocity state closes Z kinetic energy");
    checkNear(state.kineticEnergyJoules, 43.125, 5.0e-13,
              "regional velocity state closes total kinetic energy");
    checkNear(
        std::max({
            std::abs(state.domainMassClosureResidualByAxisKilograms.x),
            std::abs(state.domainMassClosureResidualByAxisKilograms.y),
            std::abs(state.domainMassClosureResidualByAxisKilograms.z),
            state.maximumAbsoluteFragmentMassClosureResidualKilograms,
            state.maximumAbsoluteComponentMassClosureResidualKilograms}),
        0.0, 2.0e-13,
        "regional velocity state closes every diagonal mass ledger");

    const auto* exterior = findFragmentVelocityStateComponent(state, 1);
    const auto* pocket = findFragmentVelocityStateComponent(state, 2);
    check(exterior != nullptr && pocket != nullptr,
          "regional velocity state retains both pressure components");
    if (exterior != nullptr && pocket != nullptr) {
        for (const double mass : {
                 exterior->massByAxisKilograms.x,
                 exterior->massByAxisKilograms.y,
                 exterior->massByAxisKilograms.z}) {
            checkNear(mass, 17.0, 2.0e-13,
                      "exterior regional velocity state closes mass");
        }
        for (const double mass : {
                 pocket->massByAxisKilograms.x,
                 pocket->massByAxisKilograms.y,
                 pocket->massByAxisKilograms.z}) {
            checkNear(mass, 3.0, 3.0e-14,
                      "pocket regional velocity state closes mass");
        }
        checkNear(exterior->momentumKilogramMetersPerSecond.x,
                  34.0, 4.0e-13,
                  "exterior regional velocity state closes momentum");
        checkNear(pocket->momentumKilogramMetersPerSecond.x,
                  6.0, 6.0e-14,
                  "pocket regional velocity state closes momentum");
        checkNear(exterior->kineticEnergyJoules,
                  36.65625, 5.0e-13,
                  "exterior regional velocity state closes energy");
        checkNear(pocket->kineticEnergyJoules,
                  6.46875, 8.0e-14,
                  "pocket regional velocity state closes energy");
    }

    for (std::size_t index = 0; index < state.samples.size(); ++index) {
        const auto& sample = state.samples[index];
        const auto& dof = metric.dofs[index];
        check(sample.dofIndex == index
                  && sample.stableId == dof.stableId
                  && sample.kind == dof.kind
                  && sample.sourceFaceLinkStableId
                      == dof.sourceFaceLinkStableId
                  && sample.axis == dof.axis
                  && sample.surfaceStableId == dof.surfaceStableId
                  && sample.componentIndex == dof.componentIndex
                  && sample.regionStableId == dof.regionStableId
                  && sample.normalVelocityMetersPerSecond == velocity[index],
              "regional velocity sample remains bound to its metric DOF");
        const double expectedMass =
            density * dof.dualVolumeCubicMeters;
        checkNear(sample.normalMomentumKilogramMetersPerSecond,
                  expectedMass * velocity[index], 2.0e-14,
                  "regional velocity sample owns scalar normal momentum");
        checkNear(sample.kineticEnergyJoules,
                  0.5 * expectedMass * velocity[index] * velocity[index],
                  2.0e-14,
                  "regional velocity sample owns diagonal kinetic energy");
    }
    for (std::size_t index = 0; index < state.fragments.size(); ++index) {
        const double expectedMass = density
            * fragments.fragments[index].volumeCubicMeters;
        for (const double mass : {
                 state.fragments[index].massByAxisKilograms.x,
                 state.fragments[index].massByAxisKilograms.y,
                 state.fragments[index].massByAxisKilograms.z}) {
            checkNear(mass, expectedMass, 2.0e-14,
                      "regional velocity fragment closes physical mass");
        }
    }
    validatePlanarPressureRegionFragmentVelocityState(
        state, geometry, sweep, fragments, topology, metric);

    const std::vector<double> zeroVelocity(metric.dofs.size(), 0.0);
    const auto zero = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, zeroVelocity, density);
    check(zero.momentumKilogramMetersPerSecond == Vector3{}
              && zero.kineticEnergyByAxisJoules == Vector3{}
              && zero.kineticEnergyJoules == 0.0
              && zero.maximumAbsoluteVelocityMetersPerSecond == 0.0,
          "zero regional velocity has exact zero momentum and energy");

    std::vector<double> independentWallVelocity(
        metric.dofs.size(), 0.0);
    for (const auto& dof : metric.dofs) {
        if (dof.surfaceStableId != 10) continue;
        if (dof.kind
            == PlanarPressureRegionFragmentVelocityDofKind::
                PressureLayerMinusTrace) {
            independentWallVelocity[dof.dofIndex] = 1.0;
        } else if (dof.kind
                   == PlanarPressureRegionFragmentVelocityDofKind::
                       PressureLayerPlusTrace) {
            independentWallVelocity[dof.dofIndex] = -2.0;
        }
    }
    const auto wall = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric,
        independentWallVelocity, density);
    const auto* wallExterior = findFragmentVelocityStateComponent(wall, 1);
    const auto* wallPocket = findFragmentVelocityStateComponent(wall, 2);
    check(wallExterior != nullptr && wallPocket != nullptr,
          "independent wall traces retain regional ownership");
    if (wallExterior != nullptr && wallPocket != nullptr) {
        checkNear(wallExterior->momentumKilogramMetersPerSecond.x,
                  0.5, 2.0e-14,
                  "minus wall trace contributes only exterior momentum");
        checkNear(wallPocket->momentumKilogramMetersPerSecond.x,
                  -3.0, 3.0e-14,
                  "plus wall trace contributes only pocket momentum");
        checkNear(wallExterior->kineticEnergyJoules,
                  0.25, 2.0e-14,
                  "minus wall trace contributes only exterior energy");
        checkNear(wallPocket->kineticEnergyJoules,
                  3.0, 3.0e-14,
                  "plus wall trace contributes only pocket energy");
    }
    checkNear(wall.momentumKilogramMetersPerSecond.x,
              -2.5, 4.0e-14,
              "independent wall traces close aggregate momentum");
    checkNear(wall.kineticEnergyJoules, 3.25, 4.0e-14,
              "independent wall traces close aggregate energy");
}

void testPlanarRegionalFragmentVelocityStateMotionAndRejection() {
    constexpr double density = 1.25;
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    std::vector<double> velocity(metric.dofs.size(), 0.0);
    for (const auto& dof : metric.dofs) {
        velocity[dof.dofIndex] = dof.axis == GridFaceAxis::X ? 1.0
            : dof.axis == GridFaceAxis::Y ? -0.5 : 0.25;
    }
    const auto state = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, velocity, density);

    auto breathingLayers = layers;
    breathingLayers[0].physicalPlaneCoordinateMeters -= 0.1;
    breathingLayers[1].physicalPlaneCoordinateMeters += 0.1;
    const auto breathingSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, breathingLayers, 0.5);
    const auto breathingFragments = buildPlanarPressureRegionFragments(
        geometry, breathingSweep);
    const auto breathingTopology = buildPlanarPressureRegionFragmentTopology(
        geometry, breathingSweep, breathingFragments);
    const auto breathingMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, breathingSweep, breathingFragments,
            breathingTopology);
    std::vector<double> breathingVelocity(
        breathingMetric.dofs.size(), 0.0);
    bool identitiesStable = breathingMetric.dofs.size() == metric.dofs.size();
    for (std::size_t index = 0;
         index < breathingMetric.dofs.size(); ++index) {
        const auto& dof = breathingMetric.dofs[index];
        breathingVelocity[index] = dof.axis == GridFaceAxis::X ? 1.0
            : dof.axis == GridFaceAxis::Y ? -0.5 : 0.25;
        identitiesStable = identitiesStable
            && dof.stableId == metric.dofs[index].stableId;
    }
    const auto breathing = buildPlanarPressureRegionFragmentVelocityState(
        geometry, breathingSweep, breathingFragments, breathingTopology,
        breathingMetric, breathingVelocity, density);
    check(identitiesStable && breathing.fingerprint != state.fingerprint,
          "topology-stable motion preserves velocity identity and updates state");
    checkNear(breathing.momentumKilogramMetersPerSecond.x,
              state.momentumKilogramMetersPerSecond.x, 2.0e-13,
              "breathing preserves uniform X momentum globally");
    checkNear(breathing.momentumKilogramMetersPerSecond.y,
              state.momentumKilogramMetersPerSecond.y, 2.0e-13,
              "breathing preserves uniform Y momentum globally");
    checkNear(breathing.momentumKilogramMetersPerSecond.z,
              state.momentumKilogramMetersPerSecond.z, 2.0e-13,
              "breathing preserves uniform Z momentum globally");
    checkNear(breathing.kineticEnergyJoules,
              state.kineticEnergyJoules, 3.0e-13,
              "breathing preserves uniform kinetic energy globally");
    const auto* breathingExterior =
        findFragmentVelocityStateComponent(breathing, 1);
    const auto* breathingPocket =
        findFragmentVelocityStateComponent(breathing, 2);
    check(breathingExterior != nullptr && breathingPocket != nullptr,
          "breathing velocity state retains both components");
    if (breathingExterior != nullptr && breathingPocket != nullptr) {
        checkNear(breathingExterior->massByAxisKilograms.x,
                  16.0, 2.0e-13,
                  "breathing exterior state follows current volume");
        checkNear(breathingPocket->massByAxisKilograms.x,
                  4.0, 5.0e-14,
                  "breathing pocket state follows current volume");
    }

    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> axisLayers{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0},
             -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0},
             -0.2, -70.0},
        };
        const auto axisSweep = makePlanarPressureRegionSweepLedger(
            geometry, axisLayers, axisLayers, 1.0);
        const auto axisFragments = buildPlanarPressureRegionFragments(
            geometry, axisSweep);
        const auto axisTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, axisSweep, axisFragments);
        const auto axisMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, axisSweep, axisFragments, axisTopology);
        std::vector<double> axisVelocity(axisMetric.dofs.size(), 0.0);
        for (const auto& dof : axisMetric.dofs) {
            axisVelocity[dof.dofIndex] =
                dof.axis == GridFaceAxis::X ? 1.0
                : dof.axis == GridFaceAxis::Y ? -0.5 : 0.25;
        }
        const auto axisState =
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, axisSweep, axisFragments, axisTopology,
                axisMetric, axisVelocity, density);
        check(axisState.profileAxis == axis,
              "all-axis regional velocity state retains profile axis");
        checkNear(axisState.momentumKilogramMetersPerSecond.x,
                  20.0, 3.0e-13,
                  "all-axis regional velocity state closes X momentum");
        checkNear(axisState.momentumKilogramMetersPerSecond.y,
                  -10.0, 2.0e-13,
                  "all-axis regional velocity state closes Y momentum");
        checkNear(axisState.momentumKilogramMetersPerSecond.z,
                  5.0, 1.0e-13,
                  "all-axis regional velocity state closes Z momentum");
        checkNear(axisState.kineticEnergyJoules,
                  13.125, 3.0e-13,
                  "all-axis regional velocity state closes kinetic energy");
    }

    auto corrupt = state;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityState(
            corrupt, geometry, sweep, fragments, topology, metric); },
        "regional velocity state rejects fingerprint corruption");
    corrupt = state;
    corrupt.samples[0].normalMomentumKilogramMetersPerSecond += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityState(
            corrupt, geometry, sweep, fragments, topology, metric); },
        "regional velocity state rejects sample corruption");
    corrupt = state;
    corrupt.fragments[0].massByAxisKilograms.x += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityState(
            corrupt, geometry, sweep, fragments, topology, metric); },
        "regional velocity state rejects fragment corruption");
    corrupt = state;
    corrupt.components[0].kineticEnergyJoules += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityState(
            corrupt, geometry, sweep, fragments, topology, metric); },
        "regional velocity state rejects component corruption");

    auto wrongSize = velocity;
    wrongSize.pop_back();
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, wrongSize,
                density)); },
        "regional velocity state rejects a wrong-sized field");
    auto nonFinite = velocity;
    nonFinite[0] = std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, nonFinite,
                density)); },
        "regional velocity state rejects a non-finite field");
    auto overflowing = velocity;
    overflowing[0] = std::numeric_limits<double>::max();
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, overflowing,
                density)); },
        "regional velocity state rejects overflowing energy");
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, velocity,
                0.0)); },
        "regional velocity state rejects zero density");

    auto limits = PlanarPressureRegionFragmentVelocityStateLimits{};
    limits.maximumSamples = state.samples.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, velocity,
                density, limits)); },
        "regional velocity state enforces the sample limit");
    limits = {};
    limits.maximumFragments = state.fragments.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, velocity,
                density, limits)); },
        "regional velocity state enforces the fragment limit");
    limits = {};
    limits.maximumComponents = state.components.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, velocity,
                density, limits)); },
        "regional velocity state enforces the component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, velocity,
                density, limits)); },
        "regional velocity state enforces the byte limit");
    limits = {};
    limits.maximumWorkingBytes = 1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentVelocityState(
            state, geometry, sweep, fragments, topology, metric, limits); },
        "regional velocity state bounds validation storage");
    limits = {};
    limits.metricLimits.maximumDofs = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric, velocity,
                density, limits)); },
        "regional velocity state enforces nested metric limits");
    auto corruptMetric = metric;
    corruptMetric.dofs[0].dualVolumeCubicMeters += 0.1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, corruptMetric,
                velocity, density)); },
        "regional velocity state rejects a mutated source metric");
}

void testPlanarRegionalFragmentPressureOperator() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
    const auto repeated =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
    check(pressureOperator == repeated
              && pressureOperator.version
                  == planarPressureRegionFragmentPressureOperatorVersion
              && pressureOperator.fingerprint != 0
              && pressureOperator.sourceFragmentFingerprint
                  == fragments.fingerprint
              && pressureOperator.sourceTopologyFingerprint
                  == topology.fingerprint
              && pressureOperator.rows.size() == 24
              && pressureOperator.entries.size() == 128
              && pressureOperator.components.size() == 2
              && pressureOperator.componentFragmentIndices.size() == 24
              && pressureOperator.includedSameRegionGridLinkCount == 64
              && pressureOperator.excludedPressureLayerWallLinkCount == 8
              && pressureOperator.ownedStorageBytes > 0,
          "regional fragment pressure operator retains the complete conductive graph");
    checkNear(pressureOperator.totalGeometryWeightMeters,
              160.0 / 3.0, 2.0e-13,
              "regional fragment pressure operator retains unique graph weight");
    checkNear(pressureOperator.totalDiagonalGeometryWeightMeters,
              320.0 / 3.0, 4.0e-13,
              "regional fragment pressure operator closes two-sided diagonal weight");

    const auto* exterior = findFragmentOperatorComponent(
        pressureOperator, 1);
    const auto* pocket = findFragmentOperatorComponent(
        pressureOperator, 2);
    check(exterior != nullptr && pocket != nullptr,
          "regional fragment pressure operator retains both gauges");
    if (exterior != nullptr && pocket != nullptr) {
        check(exterior->fragmentCount == 20
                  && pocket->fragmentCount == 4
                  && pressureOperator.rows[
                      exterior->gaugeFragmentIndex].isGauge
                  && pressureOperator.rows[
                      pocket->gaugeFragmentIndex].isGauge,
              "regional fragment pressure operator retains one gauge per component");
        checkNear(exterior->totalVolumeCubicMeters,
                  13.6, 4.0e-14,
                  "regional fragment pressure operator retains exterior volume");
        checkNear(pocket->totalVolumeCubicMeters,
                  2.4, 4.0e-14,
                  "regional fragment pressure operator retains pocket volume");
        checkNear(exterior->totalGeometryWeightMeters,
                  728.0 / 15.0, 2.0e-13,
                  "regional fragment pressure operator retains exterior graph weight");
        checkNear(pocket->totalGeometryWeightMeters,
                  4.8, 4.0e-14,
                  "regional fragment pressure operator retains pocket graph weight");
    }

    std::vector<double> regionalPressure;
    regionalPressure.reserve(fragments.fragments.size());
    for (const auto& fragment : fragments.fragments) {
        regionalPressure.push_back(fragment.pressurePascals);
    }
    const auto regionalAction =
        applyPlanarPressureRegionFragmentPressureOperator(
            pressureOperator, regionalPressure);
    for (const double value : regionalAction) {
        checkNear(value, 0.0, 0.0,
                  "regional static pressure is an exact component-null mode");
    }

    std::vector<double> first(pressureOperator.rows.size());
    std::vector<double> second(pressureOperator.rows.size());
    for (std::size_t index = 0; index < first.size(); ++index) {
        const double sample = static_cast<double>(index + 1);
        first[index] = std::sin(0.37 * sample) + 0.05 * sample;
        second[index] = std::cos(0.23 * sample) - 0.03 * sample;
    }
    const auto firstAction =
        applyPlanarPressureRegionFragmentPressureOperator(
            pressureOperator, first);
    const auto secondAction =
        applyPlanarPressureRegionFragmentPressureOperator(
            pressureOperator, second);
    checkNear(dotProduct(first, secondAction),
              dotProduct(second, firstAction), 3.0e-13,
              "regional fragment pressure operator is symmetric");
    double expectedEnergy = 0.0;
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        const double difference = first[link.minusFragmentIndex]
            - first[link.plusFragmentIndex];
        expectedEnergy += link.sameRegionGeometryWeightMeters
            * difference * difference;
    }
    check(expectedEnergy > 0.0,
          "nonconstant regional fragment pressure has positive energy");
    checkNear(dotProduct(first, firstAction), expectedEnergy, 4.0e-13,
              "regional fragment pressure energy equals its link sum");
    for (const auto& component : pressureOperator.components) {
        double actionSum = 0.0;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            actionSum += firstAction[
                pressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset]];
        }
        checkNear(actionSum, 0.0, 4.0e-14,
                  "regional fragment pressure component conserves row sum");
    }

    std::set<std::uint64_t> rowIds;
    std::set<std::uint64_t> entryIds;
    for (const auto& row : pressureOperator.rows) {
        rowIds.insert(row.stableId);
        check(row.fragmentStableId
                  == fragments.fragments[row.fragmentIndex].stableId
                  && row.diagonalGeometryWeightMeters > 0.0,
              "regional pressure row remains bound to its fragment");
    }
    for (const auto& entry : pressureOperator.entries) {
        entryIds.insert(entry.stableId);
        const auto& link = topology.links[entry.sourceFaceLinkIndex];
        check(link.kind
                  == PlanarPressureRegionFragmentFaceKind::SameRegionGrid
                  && entry.sourceFaceLinkStableId == link.stableId
                  && entry.geometryWeightMeters
                      == link.sameRegionGeometryWeightMeters,
              "regional pressure entry excludes every fabric wall");
    }
    check(rowIds.size() == pressureOperator.rows.size()
              && entryIds.size() == pressureOperator.entries.size(),
          "regional pressure rows and directed entries have unique identity");
    validatePlanarPressureRegionFragmentPressureOperator(
        pressureOperator, geometry, sweep, fragments, topology);

    const auto translatedLayers = translatePlanarPressureJumpLayers(
        geometry, layers, 0.1).layers;
    const auto translatedSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, translatedLayers, 1.0);
    const auto translatedFragments = buildPlanarPressureRegionFragments(
        geometry, translatedSweep);
    const auto translatedTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, translatedSweep, translatedFragments);
    const auto translatedOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, translatedSweep, translatedFragments,
            translatedTopology);
    std::set<std::uint64_t> translatedRowIds;
    std::set<std::uint64_t> translatedEntryIds;
    for (const auto& row : translatedOperator.rows) {
        translatedRowIds.insert(row.stableId);
    }
    for (const auto& entry : translatedOperator.entries) {
        translatedEntryIds.insert(entry.stableId);
    }
    check(translatedRowIds == rowIds
              && translatedEntryIds == entryIds
              && translatedOperator.components[0].stableId
                  == pressureOperator.components[0].stableId
              && translatedOperator.components[1].stableId
                  == pressureOperator.components[1].stableId
              && translatedOperator.fingerprint
                  != pressureOperator.fingerprint,
          "within-segment motion preserves operator identity while updating metrics");
}

void testPlanarRegionalFragmentPressureOperatorAxesAndRejection() {
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
        const auto topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        const auto pressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology);
        std::vector<double> regionalPressure;
        for (const auto& fragment : fragments.fragments) {
            regionalPressure.push_back(fragment.pressurePascals);
        }
        const auto action =
            applyPlanarPressureRegionFragmentPressureOperator(
                pressureOperator, regionalPressure);
        check(pressureOperator.rows.size() == fragments.fragments.size()
                  && pressureOperator.entries.size()
                      == 2 * topology.sameRegionGridLinkCount
                  && pressureOperator.components.size() == 2
                  && pressureOperator.excludedPressureLayerWallLinkCount
                      == topology.pressureLayerWallLinkCount
                  && pressureOperator.totalGeometryWeightMeters > 0.0,
              "regional fragment pressure operator closes on every axis");
        for (const double value : action) {
            checkNear(value, 0.0, 0.0,
                      "regional pressure jump remains a null mode on every axis");
        }
    }

    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);

    auto corrupt = pressureOperator;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentPressureOperator(
            corrupt, geometry, sweep, fragments, topology); },
        "regional pressure operator rejects fingerprint corruption");
    corrupt = pressureOperator;
    corrupt.rows[0].diagonalGeometryWeightMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentPressureOperator(
            corrupt, geometry, sweep, fragments, topology); },
        "regional pressure operator rejects row corruption");
    corrupt = pressureOperator;
    corrupt.entries[0].geometryWeightMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentPressureOperator(
            corrupt, geometry, sweep, fragments, topology); },
        "regional pressure operator rejects entry corruption");
    corrupt = pressureOperator;
    corrupt.components[0].totalVolumeCubicMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentPressureOperator(
            corrupt, geometry, sweep, fragments, topology); },
        "regional pressure operator rejects component corruption");
    corrupt = pressureOperator;
    corrupt.componentFragmentIndices[0] =
        corrupt.componentFragmentIndices[1];
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentPressureOperator(
            corrupt, geometry, sweep, fragments, topology); },
        "regional pressure operator rejects member corruption");
    expectRejected(
        [&] { static_cast<void>(
            applyPlanarPressureRegionFragmentPressureOperator(
                corrupt, std::vector<double>(fragments.fragments.size()))); },
        "regional pressure application rejects a corrupted operator");
    expectRejected(
        [&] { static_cast<void>(
            applyPlanarPressureRegionFragmentPressureOperator(
                pressureOperator, std::vector<double>(1))); },
        "regional pressure application rejects a wrong-sized field");
    auto nonfinite = std::vector<double>(fragments.fragments.size(), 0.0);
    nonfinite[0] = std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(
            applyPlanarPressureRegionFragmentPressureOperator(
                pressureOperator, nonfinite)); },
        "regional pressure application rejects non-finite values");

    auto limits =
        PlanarPressureRegionFragmentPressureOperatorLimits{};
    limits.maximumRows = pressureOperator.rows.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology, limits)); },
        "regional pressure operator enforces the row limit");
    limits = {};
    limits.maximumEntries = pressureOperator.entries.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology, limits)); },
        "regional pressure operator enforces the entry limit");
    limits = {};
    limits.maximumComponents = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology, limits)); },
        "regional pressure operator enforces the component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology, limits)); },
        "regional pressure operator enforces the byte limit");
    limits = {};
    limits.topologyLimits.maximumLinks = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology, limits)); },
        "regional pressure operator enforces nested topology limits");
    auto corruptTopology = topology;
    corruptTopology.links[0].stableId ^= 1U;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, corruptTopology)); },
        "regional pressure operator rejects mutated topology");
}

void testPlanarRegionalFragmentPressureCorrectionSolve() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);

    std::vector<double> manufactured(pressureOperator.rows.size(), 0.0);
    for (std::size_t index = 0; index < manufactured.size(); ++index) {
        const double sample = static_cast<double>(index + 1);
        manufactured[index] =
            std::sin(0.31 * sample) + 0.07 * std::cos(0.17 * sample);
    }
    subtractFragmentCorrectionVolumeMeans(
        pressureOperator, fragments, manufactured);
    const auto rightHandSide =
        applyPlanarPressureRegionFragmentPressureOperator(
            pressureOperator, manufactured);
    std::vector<double> firstCorrection(
        pressureOperator.rows.size(), 0.0);
    std::vector<double> secondCorrection = firstCorrection;
    PlanarPressureRegionFragmentPressureSolveSettings settings;
    settings.absoluteResidualTolerancePascalsMeters = 1.0e-13;
    settings.relativeResidualTolerance = 1.0e-12;
    settings.maximumIterations = 200;
    const auto first =
        solvePlanarPressureRegionFragmentPressureCorrection(
            pressureOperator, geometry, sweep, fragments, topology,
            rightHandSide, firstCorrection, settings);
    const auto second =
        solvePlanarPressureRegionFragmentPressureCorrection(
            pressureOperator, geometry, sweep, fragments, topology,
            rightHandSide, secondCorrection, settings);
    check(first == second && firstCorrection == secondCorrection
              && first.compatible && first.converged && first.finite
              && first.pressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && first.fragmentFingerprint == fragments.fingerprint
              && first.rowCount == pressureOperator.rows.size()
              && first.componentCount == 2
              && first.components.size() == 2
              && first.iterationCount > 0
              && first.iterationCount <= settings.maximumIterations,
          "regional pressure correction solve converges deterministically");
    double maximumRecoveryError = 0.0;
    for (std::size_t index = 0; index < manufactured.size(); ++index) {
        maximumRecoveryError = std::max(
            maximumRecoveryError,
            std::abs(firstCorrection[index] - manufactured[index]));
    }
    check(maximumRecoveryError < 3.0e-11,
          "regional pressure correction solve recovers the manufactured field");
    const auto recoveredRightHandSide =
        applyPlanarPressureRegionFragmentPressureOperator(
            pressureOperator, firstCorrection);
    double maximumResidual = 0.0;
    for (std::size_t index = 0; index < rightHandSide.size(); ++index) {
        maximumResidual = std::max(
            maximumResidual,
            std::abs(recoveredRightHandSide[index]
                     - rightHandSide[index]));
    }
    check(maximumResidual < 2.0e-12
              && first.finalResidualL2PascalsMeters < 1.0e-13
              && first.finalResidualMaximumPascalsMeters < 2.0e-12,
          "regional pressure correction solve closes its recomputed residual");
    for (const auto& component : pressureOperator.components) {
        checkNear(
            fragmentCorrectionVolumeMean(
                pressureOperator, fragments, component, firstCorrection),
            0.0, 3.0e-16,
            "regional pressure correction has roundoff-zero volume mean");
        checkNear(first.components[component.componentIndex]
                      .rightHandSideSumPascalsMeters,
                  0.0, 2.0e-14,
                  "regional pressure correction RHS is component-compatible");
    }
    check(first.maximumAbsoluteCorrectionVolumeMeanPascals < 3.0e-16,
          "regional pressure correction diagnoses its bounded gauge residual");

    std::vector<double> constantCorrection(
        pressureOperator.rows.size(), 0.0);
    for (const auto& component : pressureOperator.components) {
        const double constant = component.regionStableId == 1 ? 12.0 : -7.0;
        for (std::size_t offset = 0;
             offset < component.fragmentCount; ++offset) {
            constantCorrection[
                pressureOperator.componentFragmentIndices[
                    component.firstFragmentMember + offset]] = constant;
        }
    }
    const std::vector<double> zeroRightHandSide(
        pressureOperator.rows.size(), 0.0);
    const auto zero =
        solvePlanarPressureRegionFragmentPressureCorrection(
            pressureOperator, geometry, sweep, fragments, topology,
            zeroRightHandSide, constantCorrection, settings);
    check(zero.compatible && zero.converged && zero.finite
              && zero.iterationCount == 0,
          "zero regional RHS accepts without an artificial pressure mode");
    for (const double value : constantCorrection) {
        checkNear(value, 0.0, 2.0e-29,
                  "zero regional RHS removes only correction gauges");
    }
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            continue;
        }
        const double totalMinus = fragments.fragments[
            link.minusFragmentIndex].pressurePascals
            + constantCorrection[link.minusFragmentIndex];
        const double totalPlus = fragments.fragments[
            link.plusFragmentIndex].pressurePascals
            + constantCorrection[link.plusFragmentIndex];
        checkNear(totalPlus - totalMinus,
                  link.pressureJumpPascals, 0.0,
                  "zero correction preserves the authored static wall jump");
    }
}

void testPlanarRegionalFragmentPressureCorrectionRollback() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
    std::vector<double> manufactured(pressureOperator.rows.size(), 0.0);
    for (std::size_t index = 0; index < manufactured.size(); ++index) {
        manufactured[index] =
            std::sin(0.23 * static_cast<double>(index + 1));
    }
    subtractFragmentCorrectionVolumeMeans(
        pressureOperator, fragments, manufactured);
    const auto compatibleRightHandSide =
        applyPlanarPressureRegionFragmentPressureOperator(
            pressureOperator, manufactured);

    auto incompatibleRightHandSide = compatibleRightHandSide;
    incompatibleRightHandSide[
        pressureOperator.components[0].gaugeFragmentIndex] += 1.0e-4;
    std::vector<double> warmStart(
        pressureOperator.rows.size(), 0.125);
    const auto originalWarmStart = warmStart;
    const auto incompatible =
        solvePlanarPressureRegionFragmentPressureCorrection(
            pressureOperator, geometry, sweep, fragments, topology,
            incompatibleRightHandSide, warmStart);
    check(!incompatible.compatible && !incompatible.converged
              && incompatible.finite && warmStart == originalWarmStart
              && incompatible
                      .maximumAbsoluteComponentCompatibilityPascalsMeters
                  > 9.9e-5,
          "incompatible regional pressure RHS rolls back transactionally");

    PlanarPressureRegionFragmentPressureSolveSettings truncatedSettings;
    truncatedSettings.absoluteResidualTolerancePascalsMeters = 1.0e-16;
    truncatedSettings.relativeResidualTolerance = 0.0;
    truncatedSettings.maximumIterations = 1;
    warmStart = originalWarmStart;
    const auto truncated =
        solvePlanarPressureRegionFragmentPressureCorrection(
            pressureOperator, geometry, sweep, fragments, topology,
            compatibleRightHandSide, warmStart, truncatedSettings);
    check(truncated.compatible && !truncated.converged
              && truncated.iterationCount == 1
              && warmStart == originalWarmStart,
          "truncated regional pressure solve preserves its warm start");

    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> axisLayers{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0},
             -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0},
             -0.2, -70.0},
        };
        const auto axisSweep = makePlanarPressureRegionSweepLedger(
            geometry, axisLayers, axisLayers, 1.0);
        const auto axisFragments = buildPlanarPressureRegionFragments(
            geometry, axisSweep);
        const auto axisTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, axisSweep, axisFragments);
        const auto axisOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, axisSweep, axisFragments, axisTopology);
        std::vector<double> correction(axisOperator.rows.size(), 3.0);
        const auto diagnostics =
            solvePlanarPressureRegionFragmentPressureCorrection(
                axisOperator, geometry, axisSweep, axisFragments,
                axisTopology,
                std::vector<double>(axisOperator.rows.size(), 0.0),
                correction);
        check(diagnostics.compatible && diagnostics.converged
                  && diagnostics.finite
                  && diagnostics.iterationCount == 0,
              "zero regional pressure correction closes on every axis");
        for (const double value : correction) {
            checkNear(value, 0.0, 2.0e-29,
                      "all-axis regional correction gauge is roundoff-zero");
        }
    }

    expectRejected(
        [&] { auto correction = originalWarmStart;
            static_cast<void>(
                solvePlanarPressureRegionFragmentPressureCorrection(
                    pressureOperator, geometry, sweep, fragments,
                    topology, std::vector<double>(1), correction)); },
        "regional pressure solve rejects a wrong-sized RHS");
    expectRejected(
        [&] { auto correction = std::vector<double>(1);
            static_cast<void>(
                solvePlanarPressureRegionFragmentPressureCorrection(
                    pressureOperator, geometry, sweep, fragments,
                    topology, compatibleRightHandSide, correction)); },
        "regional pressure solve rejects a wrong-sized correction");
    auto nonfinite = compatibleRightHandSide;
    nonfinite[0] = std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { auto correction = originalWarmStart;
            static_cast<void>(
                solvePlanarPressureRegionFragmentPressureCorrection(
                    pressureOperator, geometry, sweep, fragments,
                    topology, nonfinite, correction)); },
        "regional pressure solve rejects non-finite fields");
    auto invalidSettings =
        PlanarPressureRegionFragmentPressureSolveSettings{};
    invalidSettings.maximumIterations = 0;
    expectRejected(
        [&] { auto correction = originalWarmStart;
            static_cast<void>(
                solvePlanarPressureRegionFragmentPressureCorrection(
                    pressureOperator, geometry, sweep, fragments,
                    topology, compatibleRightHandSide, correction,
                    invalidSettings)); },
        "regional pressure solve rejects a zero iteration bound");
    invalidSettings = {};
    invalidSettings.absoluteResidualTolerancePascalsMeters =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { auto correction = originalWarmStart;
            static_cast<void>(
                solvePlanarPressureRegionFragmentPressureCorrection(
                    pressureOperator, geometry, sweep, fragments,
                    topology, compatibleRightHandSide, correction,
                    invalidSettings)); },
        "regional pressure solve rejects invalid tolerances");
    auto corruptOperator = pressureOperator;
    corruptOperator.entries[0].geometryWeightMeters += 0.1;
    expectRejected(
        [&] { auto correction = originalWarmStart;
            static_cast<void>(
                solvePlanarPressureRegionFragmentPressureCorrection(
                    corruptOperator, geometry, sweep, fragments,
                    topology, compatibleRightHandSide, correction)); },
        "regional pressure solve rejects a mutated operator");
}

void testPlanarRegionalFragmentPressureProjection() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);

    PlanarPressureRegionFragmentPressureProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.01;
    settings.absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
    settings.relativeContinuityTolerance = 1.0e-10;
    settings.pressureSolve.absoluteResidualTolerancePascalsMeters = 1.0e-13;
    settings.pressureSolve.relativeResidualTolerance = 1.0e-12;
    settings.pressureSolve.maximumIterations = 200;

    std::vector<double> manufactured(pressureOperator.rows.size(), 0.0);
    for (std::size_t index = 0; index < manufactured.size(); ++index) {
        const double sample = static_cast<double>(index + 1);
        manufactured[index] =
            0.8 * std::sin(0.29 * sample)
            + 0.13 * std::cos(0.11 * sample);
    }
    subtractFragmentCorrectionVolumeMeans(
        pressureOperator, fragments, manufactured);

    std::vector<double> predictedVelocity(topology.links.size(), 0.0);
    for (const auto& link : topology.links) {
        if (link.kind
            != PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            continue;
        }
        predictedVelocity[link.linkIndex] =
            -settings.timeStepSeconds / settings.densityKgPerCubicMeter
            * (manufactured[link.minusFragmentIndex]
               - manufactured[link.plusFragmentIndex])
            / link.centerDistanceMeters;
    }
    const auto originalPredictedVelocity = predictedVelocity;
    std::vector<double> pressureCorrection(
        pressureOperator.rows.size(), 0.0);
    auto repeatedVelocity = predictedVelocity;
    auto repeatedPressure = pressureCorrection;
    const auto diagnostics =
        projectStaticPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            predictedVelocity, pressureCorrection, settings);
    const auto repeated =
        projectStaticPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            repeatedVelocity, repeatedPressure, settings);
    check(diagnostics == repeated
              && predictedVelocity == repeatedVelocity
              && pressureCorrection == repeatedPressure
              && diagnostics.accepted && diagnostics.finite
              && diagnostics.staticGeometry
              && !diagnostics.usesMovingVolumeRates
              && diagnostics.volumeRateFingerprint == 0
              && diagnostics.pressureSolve.compatible
              && diagnostics.pressureSolve.converged
              && diagnostics.pressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && diagnostics.topologyFingerprint == topology.fingerprint
              && diagnostics.fragmentFingerprint == fragments.fingerprint
              && diagnostics.fragmentCount == fragments.fragments.size()
              && diagnostics.linkCount == topology.links.size()
              && diagnostics.projectedSameRegionGridLinkCount == 64
              && diagnostics.sealedPressureLayerWallLinkCount == 8
              && diagnostics.workingStorageBytes
                  == 9 * sizeof(double) * fragments.fragments.size()
                     + sizeof(double) * topology.links.size()
                     + sizeof(
                         PlanarPressureRegionFragmentPressureSolveComponentDiagnostics)
                         * pressureOperator.components.size(),
          "static regional face projection is deterministic and source-bound");
    double maximumPressureRecoveryError = 0.0;
    double maximumCorrectedVelocity = 0.0;
    for (std::size_t index = 0; index < pressureCorrection.size(); ++index) {
        maximumPressureRecoveryError = std::max(
            maximumPressureRecoveryError,
            std::abs(pressureCorrection[index] - manufactured[index]));
    }
    for (const double velocity : predictedVelocity) {
        maximumCorrectedVelocity = std::max(
            maximumCorrectedVelocity, std::abs(velocity));
    }
    check(maximumPressureRecoveryError < 3.0e-11
              && maximumCorrectedVelocity < 3.0e-14
              && diagnostics
                      .predictedNetOutwardFlowMaximumCubicMetersPerSecond
                  > 1.0e-4
              && diagnostics
                      .correctedNetOutwardFlowMaximumCubicMetersPerSecond
                  < 3.0e-14
              && diagnostics
                      .correctedContinuityResidualMaximumCubicMetersPerSecond
                  == diagnostics
                      .correctedNetOutwardFlowMaximumCubicMetersPerSecond
              && diagnostics
                      .maximumAbsoluteCorrectedComponentBalanceCubicMetersPerSecond
                  < 3.0e-14,
          "regional pressure projection cancels manufactured divergence");
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            checkNear(predictedVelocity[link.linkIndex], 0.0, 0.0,
                      "regional pressure projection keeps layer walls sealed");
        }
    }

    std::vector<double> uniformVelocity(topology.links.size(), 0.0);
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            continue;
        }
        switch (link.axis) {
        case GridFaceAxis::X:
            uniformVelocity[link.linkIndex] = 0.0;
            break;
        case GridFaceAxis::Y:
            uniformVelocity[link.linkIndex] = -0.11;
            break;
        case GridFaceAxis::Z:
            uniformVelocity[link.linkIndex] = 0.08;
            break;
        }
    }
    const auto originalUniformVelocity = uniformVelocity;
    std::vector<double> uniformPressure(
        pressureOperator.rows.size(), 5.0);
    const auto uniformDiagnostics =
        projectStaticPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            uniformVelocity, uniformPressure, settings);
    double maximumUniformVelocityChange = 0.0;
    for (std::size_t index = 0; index < uniformVelocity.size(); ++index) {
        maximumUniformVelocityChange = std::max(
            maximumUniformVelocityChange,
            std::abs(uniformVelocity[index]
                     - originalUniformVelocity[index]));
    }
    check(uniformDiagnostics.accepted
              && uniformDiagnostics.pressureSolve.iterationCount == 0
              && uniformDiagnostics
                      .predictedNetOutwardFlowMaximumCubicMetersPerSecond
                  == 0.0
              && uniformDiagnostics
                      .correctedNetOutwardFlowMaximumCubicMetersPerSecond
                  == 0.0
              && maximumUniformVelocityChange == 0.0,
          "regional face projection preserves uniform wall-tangential flow");
    for (const double pressure : uniformPressure) {
        checkNear(pressure, 0.0, 2.0e-29,
                  "uniform flow removes only regional correction gauges");
    }

    check(originalPredictedVelocity != predictedVelocity,
          "accepted regional projection publishes corrected link velocity");
}

void testPlanarRegionalFragmentProjectionEnergyAudit() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);

    PlanarPressureRegionFragmentPressureProjectionSettings projectionSettings;
    projectionSettings.densityKgPerCubicMeter = 1.2;
    projectionSettings.timeStepSeconds = 0.01;
    projectionSettings.absoluteContinuityToleranceCubicMetersPerSecond =
        1.0e-12;
    projectionSettings.relativeContinuityTolerance = 1.0e-10;
    projectionSettings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-13;
    projectionSettings.pressureSolve.relativeResidualTolerance = 1.0e-12;
    projectionSettings.pressureSolve.maximumIterations = 200;
    PlanarPressureRegionFragmentProjectionEnergySettings auditSettings;
    auditSettings.densityKgPerCubicMeter =
        projectionSettings.densityKgPerCubicMeter;
    auditSettings.timeStepSeconds = projectionSettings.timeStepSeconds;
    auditSettings.absoluteContinuityToleranceCubicMetersPerSecond =
        projectionSettings
            .absoluteContinuityToleranceCubicMetersPerSecond;
    auditSettings.relativeContinuityTolerance =
        projectionSettings.relativeContinuityTolerance;

    std::vector<double> manufactured(pressureOperator.rows.size(), 0.0);
    for (std::size_t index = 0; index < manufactured.size(); ++index) {
        const double sample = static_cast<double>(index + 1);
        manufactured[index] =
            0.8 * std::sin(0.29 * sample)
            + 0.13 * std::cos(0.11 * sample);
    }
    subtractFragmentCorrectionVolumeMeans(
        pressureOperator, fragments, manufactured);
    std::vector<double> linkVelocity(topology.links.size(), 0.0);
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            linkVelocity[link.linkIndex] =
                -projectionSettings.timeStepSeconds
                / projectionSettings.densityKgPerCubicMeter
                * (manufactured[link.minusFragmentIndex]
                   - manufactured[link.plusFragmentIndex])
                / link.centerDistanceMeters;
        }
    }
    const auto beforeLinkVelocity = linkVelocity;
    std::vector<double> pressure(pressureOperator.rows.size(), 0.0);
    const auto projection =
        projectStaticPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            linkVelocity, pressure, projectionSettings);
    check(projection.accepted,
          "projection-energy fixture has an accepted pressure correction");
    const auto beforeVelocity = regionalMetricVelocityFromTopologyLinks(
        metric, beforeLinkVelocity);
    const auto afterVelocity = regionalMetricVelocityFromTopologyLinks(
        metric, linkVelocity);
    const auto before = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, beforeVelocity,
        projectionSettings.densityKgPerCubicMeter);
    const auto after = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, afterVelocity,
        projectionSettings.densityKgPerCubicMeter);
    const auto audit =
        auditStaticPlanarPressureRegionFragmentProjectionEnergy(
            geometry, sweep, fragments, topology, metric, before, after,
            pressure, auditSettings);
    const auto repeated =
        auditStaticPlanarPressureRegionFragmentProjectionEnergy(
            geometry, sweep, fragments, topology, metric, before, after,
            pressure, auditSettings);
    check(audit == repeated
              && audit.version
                  == planarPressureRegionFragmentProjectionEnergyVersion
              && audit.fingerprint != 0
              && audit.accepted
              && audit.nonIncreasingKineticEnergy
              && audit.sourceMetricFingerprint == metric.fingerprint
              && audit.sourceTopologyFingerprint == topology.fingerprint
              && audit.sourceFragmentFingerprint == fragments.fingerprint
              && audit.volumeRateFingerprint == 0
              && audit.beforeVelocityStateFingerprint == before.fingerprint
              && audit.afterVelocityStateFingerprint == after.fingerprint
              && audit.staticGeometry
              && !audit.usesMovingVolumeRates
              && audit.pressureCorrectionPascals == pressure
              && audit.corrections.size() == 64
              && audit.components.size() == 2
              && audit.pressureLayerTraceCount == 16
              && audit
                      .maximumAbsoluteWallTraceVelocityResidualMetersPerSecond
                  == 0.0
              && audit
                      .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond
                  == 0.0
              && audit.ownedStorageBytes > 0
              && audit.workingStorageBytes
                  == 2 * sizeof(double) * fragments.fragments.size()
                     + sizeof(double) * metric.components.size(),
          "static projection energy audit is deterministic and source-bound");
    check(audit.predictedContinuityResidualMaximumCubicMetersPerSecond
              > 1.0e-4
              && audit
                      .correctedContinuityResidualMaximumCubicMetersPerSecond
                  < 3.0e-14
              && audit.maximumAbsoluteVelocityChangeResidualMetersPerSecond
                  < 3.0e-16
              && audit
                      .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond
                  < 3.0e-17
              && audit.maximumAbsolutePressureGaugePascals < 3.0e-16,
          "static projection energy audit closes continuity, impulse, and gauge");
    check(audit.kineticEnergyBeforeJoules > 0.0
              && audit.kineticEnergyAfterJoules < 1.0e-26
              && audit.kineticEnergyChangeJoules < 0.0
              && audit.midpointPressureWorkJoules < 0.0
              && audit.kineticEnergyRemovedJoules > 0.0,
          "manufactured gradient projection removes its kinetic energy");
    checkNear(audit.workEnergyResidualJoules, 0.0, 3.0e-18,
              "midpoint pressure work closes total kinetic-energy change");
    checkNear(audit.kineticEnergyRemovedJoules,
              audit.kineticEnergyBeforeJoules
                  - audit.kineticEnergyAfterJoules,
              0.0,
              "projection energy audit retains exact removed energy");
    checkNear(audit.correctionKineticEnergyJoules,
              audit.kineticEnergyRemovedJoules, 3.0e-18,
              "static correction kinetic energy equals removed energy");
    checkNear(audit.finalPressureWorkJoules, 0.0, 3.0e-18,
              "static final divergence has zero correction-pressure work");
    checkNear(audit.geometryPressureWorkJoules, 0.0, 0.0,
              "static geometry contributes exact zero pressure work");
    checkNear(audit.affineEnergyResidualJoules, 0.0, 3.0e-18,
              "static affine identity reduces to orthogonal projection");
    checkNear(
        std::max({
            std::abs(
                audit.momentumImpulseResidualKilogramMetersPerSecond.x),
            std::abs(
                audit.momentumImpulseResidualKilogramMetersPerSecond.y),
            std::abs(
                audit.momentumImpulseResidualKilogramMetersPerSecond.z)}),
        0.0, 3.0e-16,
        "pressure impulse closes aggregate diagonal momentum change");
    for (const auto& correction : audit.corrections) {
        check(correction.correctionIndex < audit.corrections.size()
                  && correction.dofIndex < metric.dofs.size()
                  && correction.dofStableId
                      == metric.dofs[correction.dofIndex].stableId
                  && correction.diagonalMassKilograms > 0.0,
              "projection energy correction retains metric identity and mass");
        checkNear(correction.velocityChangeResidualMetersPerSecond,
                  0.0, 3.0e-16,
                  "projection correction matches pressure gradient");
        checkNear(
            correction
                .momentumImpulseResidualKilogramMetersPerSecond,
            0.0, 3.0e-17,
            "projection correction matches pressure impulse");
        checkNear(correction.workEnergyResidualJoules,
                  0.0, 3.0e-18,
                  "projection correction closes midpoint work and energy");
    }
    const auto exterior = std::ranges::find(
        audit.components, 1,
        &PlanarPressureRegionFragmentProjectionEnergyComponent::
            regionStableId);
    const auto pocket = std::ranges::find(
        audit.components, 2,
        &PlanarPressureRegionFragmentProjectionEnergyComponent::
            regionStableId);
    check(exterior != audit.components.end()
              && pocket != audit.components.end()
              && exterior->correctionCount == 56
              && pocket->correctionCount == 8,
          "projection energy audit retains separate component corrections");
    validateStaticPlanarPressureRegionFragmentProjectionEnergyAudit(
        audit, geometry, sweep, fragments, topology, metric, before, after);

    std::vector<double> uniformLinks(topology.links.size(), 0.0);
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            continue;
        }
        if (link.axis == GridFaceAxis::Y) {
            uniformLinks[link.linkIndex] = -0.11;
        } else if (link.axis == GridFaceAxis::Z) {
            uniformLinks[link.linkIndex] = 0.08;
        }
    }
    const auto uniformBeforeLinks = uniformLinks;
    std::vector<double> uniformPressure(
        pressureOperator.rows.size(), 5.0);
    const auto uniformProjection =
        projectStaticPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            uniformLinks, uniformPressure, projectionSettings);
    const auto uniformBeforeVelocity =
        regionalMetricVelocityFromTopologyLinks(
            metric, uniformBeforeLinks);
    const auto uniformAfterVelocity =
        regionalMetricVelocityFromTopologyLinks(metric, uniformLinks);
    const auto uniformBefore =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            uniformBeforeVelocity,
            projectionSettings.densityKgPerCubicMeter);
    const auto uniformAfter =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            uniformAfterVelocity,
            projectionSettings.densityKgPerCubicMeter);
    const auto uniformAudit =
        auditStaticPlanarPressureRegionFragmentProjectionEnergy(
            geometry, sweep, fragments, topology, metric, uniformBefore,
            uniformAfter, uniformPressure, auditSettings);
    check(uniformProjection.accepted && uniformAudit.accepted
              && uniformAudit.kineticEnergyBeforeJoules > 0.0
              && uniformAudit.kineticEnergyAfterJoules
                  == uniformAudit.kineticEnergyBeforeJoules
              && uniformAudit.kineticEnergyChangeJoules == 0.0
              && uniformAudit.midpointPressureWorkJoules == 0.0
              && uniformAudit.pressureImpulseKilogramMetersPerSecond
                  == Vector3{},
          "projection energy audit preserves the tangential null field exactly");
}

void testPlanarRegionalFragmentProjectionEnergyAxesAndRejection() {
    const auto geometry = grid();
    PlanarPressureRegionFragmentPressureProjectionSettings projectionSettings;
    projectionSettings.densityKgPerCubicMeter = 1.2;
    projectionSettings.timeStepSeconds = 0.01;
    projectionSettings.absoluteContinuityToleranceCubicMetersPerSecond =
        1.0e-12;
    projectionSettings.relativeContinuityTolerance = 1.0e-10;
    projectionSettings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-13;
    projectionSettings.pressureSolve.relativeResidualTolerance = 1.0e-12;
    projectionSettings.pressureSolve.maximumIterations = 200;
    PlanarPressureRegionFragmentProjectionEnergySettings auditSettings;
    auditSettings.densityKgPerCubicMeter =
        projectionSettings.densityKgPerCubicMeter;
    auditSettings.timeStepSeconds = projectionSettings.timeStepSeconds;

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
        const auto topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        const auto pressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology);
        const auto metric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology);
        std::vector<double> manufactured(
            pressureOperator.rows.size(), 0.0);
        for (std::size_t index = 0; index < manufactured.size(); ++index) {
            manufactured[index] =
                std::sin(0.19 * static_cast<double>(index + 1));
        }
        subtractFragmentCorrectionVolumeMeans(
            pressureOperator, fragments, manufactured);
        std::vector<double> linkVelocity(topology.links.size(), 0.0);
        for (const auto& link : topology.links) {
            if (link.kind
                == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
                linkVelocity[link.linkIndex] =
                    -projectionSettings.timeStepSeconds
                    / projectionSettings.densityKgPerCubicMeter
                    * (manufactured[link.minusFragmentIndex]
                       - manufactured[link.plusFragmentIndex])
                    / link.centerDistanceMeters;
            }
        }
        const auto beforeLinks = linkVelocity;
        std::vector<double> pressure(pressureOperator.rows.size(), 0.0);
        const auto projection =
            projectStaticPlanarPressureRegionFragmentFaceVelocities(
                pressureOperator, geometry, sweep, fragments, topology,
                linkVelocity, pressure, projectionSettings);
        const auto beforeVelocity =
            regionalMetricVelocityFromTopologyLinks(metric, beforeLinks);
        const auto afterVelocity =
            regionalMetricVelocityFromTopologyLinks(metric, linkVelocity);
        const auto before =
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric,
                beforeVelocity,
                projectionSettings.densityKgPerCubicMeter);
        const auto after =
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric,
                afterVelocity,
                projectionSettings.densityKgPerCubicMeter);
        const auto audit =
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, before, after,
                pressure, auditSettings);
        check(projection.accepted && audit.accepted
                  && audit.corrections.size()
                      == metric.sharedRegionGridDofCount
                  && audit.pressureLayerTraceCount
                      == metric.pressureLayerTraceDofCount
                  && audit.kineticEnergyRemovedJoules > 0.0
                  && audit
                          .correctedContinuityResidualMaximumCubicMetersPerSecond
                      < 3.0e-14,
              "projection energy audit closes a nontrivial field on every axis");
    }

    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    const std::vector<double> zeroVelocity(metric.dofs.size(), 0.0);
    const auto zero = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, zeroVelocity,
        projectionSettings.densityKgPerCubicMeter);
    const std::vector<double> zeroPressure(fragments.fragments.size(), 0.0);
    const auto audit =
        auditStaticPlanarPressureRegionFragmentProjectionEnergy(
            geometry, sweep, fragments, topology, metric, zero, zero,
            zeroPressure, auditSettings);
    check(audit.accepted && audit.kineticEnergyRemovedJoules == 0.0,
          "zero projection has an exact zero energy certificate");

    auto corruptAudit = audit;
    corruptAudit.fingerprint = 0;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentProjectionEnergyAudit(
                corruptAudit, geometry, sweep, fragments, topology, metric,
                zero, zero);
        },
        "projection energy audit rejects fingerprint corruption");
    corruptAudit = audit;
    corruptAudit.corrections[0].pressureImpulseKilogramMetersPerSecond +=
        0.1;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentProjectionEnergyAudit(
                corruptAudit, geometry, sweep, fragments, topology, metric,
                zero, zero);
        },
        "projection energy audit rejects correction corruption");
    corruptAudit = audit;
    corruptAudit.components[0].workEnergyResidualJoules += 0.1;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentProjectionEnergyAudit(
                corruptAudit, geometry, sweep, fragments, topology, metric,
                zero, zero);
        },
        "projection energy audit rejects component corruption");

    auto changedVelocity = zeroVelocity;
    const auto sharedDof = std::ranges::find_if(
        metric.dofs,
        [](const auto& dof) {
            return dof.kind
                == PlanarPressureRegionFragmentVelocityDofKind::
                    SharedRegionGrid;
        });
    changedVelocity[sharedDof->dofIndex] = 0.01;
    const auto inconsistentAfter =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            changedVelocity,
            projectionSettings.densityKgPerCubicMeter);
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero,
                inconsistentAfter, zeroPressure, auditSettings)); },
        "projection energy audit rejects a velocity change without pressure impulse");

    auto wallVelocity = zeroVelocity;
    const auto wallDof = std::ranges::find_if(
        metric.dofs,
        [](const auto& dof) {
            return dof.kind
                != PlanarPressureRegionFragmentVelocityDofKind::
                    SharedRegionGrid;
        });
    wallVelocity[wallDof->dofIndex] = 0.01;
    const auto unsealed = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, wallVelocity,
        projectionSettings.densityKgPerCubicMeter);
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, unsealed,
                unsealed, zeroPressure, auditSettings)); },
        "projection energy audit rejects an unsealed one-sided wall trace");

    auto shiftedPressure = zeroPressure;
    std::ranges::fill(shiftedPressure, 0.01);
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                shiftedPressure, auditSettings)); },
        "projection energy audit rejects a nonzero correction gauge");
    auto wrongPressure = zeroPressure;
    wrongPressure.pop_back();
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                wrongPressure, auditSettings)); },
        "projection energy audit rejects a wrong-sized pressure field");
    auto nonFinitePressure = zeroPressure;
    nonFinitePressure[0] = std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                nonFinitePressure, auditSettings)); },
        "projection energy audit rejects non-finite pressure");

    auto corruptState = zero;
    corruptState.samples[0].kineticEnergyJoules += 0.1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, corruptState,
                zero, zeroPressure, auditSettings)); },
        "projection energy audit rejects a corrupted velocity state");
    auto invalidSettings = auditSettings;
    invalidSettings.densityKgPerCubicMeter = 1.3;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, invalidSettings)); },
        "projection energy audit rejects density foreign to its states");
    invalidSettings = auditSettings;
    invalidSettings.timeStepSeconds = 0.0;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, invalidSettings)); },
        "projection energy audit rejects a nonpositive time step");
    invalidSettings = auditSettings;
    invalidSettings.absoluteEnergyResidualToleranceJoules = 0.0;
    invalidSettings.relativeEnergyResidualTolerance = 0.0;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, invalidSettings)); },
        "projection energy audit rejects an empty energy tolerance");

    const auto movingLayers = translatePlanarPressureJumpLayers(
        geometry, layers, 0.1).layers;
    const auto movingSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, movingLayers, 1.0);
    const auto movingFragments = buildPlanarPressureRegionFragments(
        geometry, movingSweep);
    const auto movingTopology = buildPlanarPressureRegionFragmentTopology(
        geometry, movingSweep, movingFragments);
    const auto movingMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, movingSweep, movingFragments, movingTopology);
    const std::vector<double> movingVelocity(
        movingMetric.dofs.size(), 0.0);
    const auto movingState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, movingSweep, movingFragments, movingTopology,
            movingMetric, movingVelocity,
            projectionSettings.densityKgPerCubicMeter);
    const std::vector<double> movingPressure(
        movingFragments.fragments.size(), 0.0);
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, movingSweep, movingFragments, movingTopology,
                movingMetric, movingState, movingState, movingPressure,
                auditSettings)); },
        "static projection energy audit rejects moving geometry");

    auto limits = PlanarPressureRegionFragmentProjectionEnergyLimits{};
    limits.maximumCorrections = metric.sharedRegionGridDofCount - 1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, auditSettings, limits)); },
        "projection energy audit enforces the correction limit");
    limits = {};
    limits.maximumPressureSamples = zeroPressure.size() - 1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, auditSettings, limits)); },
        "projection energy audit enforces the pressure-sample limit");
    limits = {};
    limits.maximumComponents = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, auditSettings, limits)); },
        "projection energy audit enforces the component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, auditSettings, limits)); },
        "projection energy audit enforces the owned-byte limit");
    limits = {};
    limits.maximumWorkingBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, auditSettings, limits)); },
        "projection energy audit enforces the working-byte limit");
    limits = {};
    limits.velocityStateLimits.maximumSamples = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, metric, zero, zero,
                zeroPressure, auditSettings, limits)); },
        "projection energy audit enforces nested state limits");
    auto corruptMetric = metric;
    corruptMetric.dofs[0].dualVolumeCubicMeters += 0.1;
    expectRejected(
        [&] { static_cast<void>(
            auditStaticPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, corruptMetric, zero,
                zero, zeroPressure, auditSettings)); },
        "projection energy audit rejects a mutated source metric");
}

void testPlanarRegionalFragmentPressureProjectionAxesAndRollback() {
    const auto geometry = grid();
    PlanarPressureRegionFragmentPressureProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.01;
    settings.pressureSolve.absoluteResidualTolerancePascalsMeters = 1.0e-13;
    settings.pressureSolve.relativeResidualTolerance = 1.0e-12;
    settings.pressureSolve.maximumIterations = 200;

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
        const auto topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        const auto pressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology);
        std::vector<double> manufactured(
            pressureOperator.rows.size(), 0.0);
        for (std::size_t index = 0; index < manufactured.size(); ++index) {
            manufactured[index] =
                std::sin(0.19 * static_cast<double>(index + 1));
        }
        subtractFragmentCorrectionVolumeMeans(
            pressureOperator, fragments, manufactured);
        std::vector<double> velocity(topology.links.size(), 0.0);
        for (const auto& link : topology.links) {
            if (link.kind
                == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
                velocity[link.linkIndex] =
                    -settings.timeStepSeconds
                    / settings.densityKgPerCubicMeter
                    * (manufactured[link.minusFragmentIndex]
                       - manufactured[link.plusFragmentIndex])
                    / link.centerDistanceMeters;
            }
        }
        std::vector<double> pressure(pressureOperator.rows.size(), 0.0);
        const auto diagnostics =
            projectStaticPlanarPressureRegionFragmentFaceVelocities(
                pressureOperator, geometry, sweep, fragments, topology,
                velocity, pressure, settings);
        check(diagnostics.accepted
                  && diagnostics
                          .correctedNetOutwardFlowMaximumCubicMetersPerSecond
                      < 3.0e-14,
              "regional face projection closes manufactured flow on every axis");
    }

    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
    std::vector<double> manufactured(pressureOperator.rows.size(), 0.0);
    for (std::size_t index = 0; index < manufactured.size(); ++index) {
        manufactured[index] =
            std::sin(0.31 * static_cast<double>(index + 1));
    }
    subtractFragmentCorrectionVolumeMeans(
        pressureOperator, fragments, manufactured);
    std::vector<double> velocity(topology.links.size(), 0.0);
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::SameRegionGrid) {
            velocity[link.linkIndex] =
                -settings.timeStepSeconds / settings.densityKgPerCubicMeter
                * (manufactured[link.minusFragmentIndex]
                   - manufactured[link.plusFragmentIndex])
                / link.centerDistanceMeters;
        }
    }
    std::vector<double> pressure(pressureOperator.rows.size(), 0.25);
    const auto originalVelocity = velocity;
    const auto originalPressure = pressure;
    auto truncatedSettings = settings;
    truncatedSettings.pressureSolve.absoluteResidualTolerancePascalsMeters =
        1.0e-16;
    truncatedSettings.pressureSolve.relativeResidualTolerance = 0.0;
    truncatedSettings.pressureSolve.maximumIterations = 1;
    const auto truncated =
        projectStaticPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            velocity, pressure, truncatedSettings);
    check(!truncated.accepted && truncated.pressureSolve.compatible
              && !truncated.pressureSolve.converged
              && velocity == originalVelocity
              && pressure == originalPressure,
          "truncated regional face projection rolls back both state vectors");

    auto overstrictContinuitySettings = settings;
    overstrictContinuitySettings
        .absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-30;
    overstrictContinuitySettings.relativeContinuityTolerance = 0.0;
    velocity = originalVelocity;
    pressure = originalPressure;
    const auto overstrictContinuity =
        projectStaticPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            velocity, pressure, overstrictContinuitySettings);
    check(!overstrictContinuity.accepted
              && overstrictContinuity.pressureSolve.converged
              && overstrictContinuity
                      .correctedNetOutwardFlowMaximumCubicMetersPerSecond
                  > overstrictContinuity
                      .continuityToleranceCubicMetersPerSecond
              && velocity == originalVelocity
              && pressure == originalPressure,
          "failed regional continuity check rolls back both state vectors");

    auto movingLayers = translatePlanarPressureJumpLayers(
        geometry, layers, 0.1).layers;
    const auto movingSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, movingLayers, 1.0);
    const auto movingFragments = buildPlanarPressureRegionFragments(
        geometry, movingSweep);
    const auto movingTopology = buildPlanarPressureRegionFragmentTopology(
        geometry, movingSweep, movingFragments);
    const auto movingOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, movingSweep, movingFragments, movingTopology);
    expectRejected(
        [&] {
            auto movingVelocity =
                std::vector<double>(movingTopology.links.size(), 0.0);
            auto movingPressure =
                std::vector<double>(movingOperator.rows.size(), 0.0);
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    movingOperator, geometry, movingSweep, movingFragments,
                    movingTopology, movingVelocity, movingPressure));
        },
        "static regional face projection rejects moving layer geometry");

    const auto wall = std::ranges::find_if(
        topology.links,
        [](const auto& link) {
            return link.kind
                == PlanarPressureRegionFragmentFaceKind::PressureLayerWall;
        });
    expectRejected(
        [&] {
            auto invalidVelocity = originalVelocity;
            invalidVelocity[wall->linkIndex] = 1.0e-15;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    invalidVelocity, candidatePressure, settings));
        },
        "regional face projection rejects nonzero layer-wall flow");
    expectRejected(
        [&] {
            auto invalidVelocity = std::vector<double>(1, 0.0);
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    invalidVelocity, candidatePressure, settings));
        },
        "regional face projection rejects wrong-sized link velocity");
    expectRejected(
        [&] {
            auto invalidVelocity = originalVelocity;
            invalidVelocity[0] = std::numeric_limits<double>::quiet_NaN();
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    invalidVelocity, candidatePressure, settings));
        },
        "regional face projection rejects non-finite link velocity");
    auto invalidSettings = settings;
    invalidSettings.densityKgPerCubicMeter = 0.0;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    candidateVelocity, candidatePressure, invalidSettings));
        },
        "regional face projection rejects nonpositive density");
    invalidSettings = settings;
    invalidSettings.timeStepSeconds =
        std::numeric_limits<double>::quiet_NaN();
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    candidateVelocity, candidatePressure, invalidSettings));
        },
        "regional face projection rejects non-finite time step");
    invalidSettings = settings;
    invalidSettings.absoluteContinuityToleranceCubicMetersPerSecond = 0.0;
    invalidSettings.relativeContinuityTolerance = 0.0;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    candidateVelocity, candidatePressure, invalidSettings));
        },
        "regional face projection rejects an empty continuity tolerance");
    auto corruptOperator = pressureOperator;
    corruptOperator.rows[0].diagonalGeometryWeightMeters += 0.1;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    corruptOperator, geometry, sweep, fragments, topology,
                    candidateVelocity, candidatePressure, settings));
        },
        "regional face projection rejects a mutated operator");
    auto limits =
        PlanarPressureRegionFragmentPressureProjectionLimits{};
    limits.maximumWorkingBytes = 1;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    candidateVelocity, candidatePressure, settings, limits));
        },
        "regional face projection enforces its working-storage limit");
    limits = {};
    limits.pressureOperatorLimits.maximumRows = 1;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectStaticPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    candidateVelocity, candidatePressure, settings, limits));
        },
        "regional face projection enforces nested operator limits");
}

void testPlanarRegionalMovingFragmentPressureProjection() {
    const auto geometry = grid();
    const auto previousLayers = pocketLayers();
    const auto currentLayers = translatePlanarPressureJumpLayers(
        geometry, previousLayers, 0.1).layers;
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, previousLayers, currentLayers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
    const auto volumeRates = buildPlanarPressureRegionFragmentVolumeRates(
        geometry, sweep, fragments, topology);

    PlanarPressureRegionFragmentPressureProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 1.0;
    settings.absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-12;
    settings.relativeContinuityTolerance = 1.0e-10;
    settings.pressureSolve.absoluteResidualTolerancePascalsMeters = 1.0e-13;
    settings.pressureSolve.relativeResidualTolerance = 1.0e-12;
    settings.pressureSolve.maximumIterations = 200;

    std::vector<double> velocity(topology.links.size(), 0.0);
    const auto beforeProjectionVelocity = velocity;
    std::vector<double> pressure(pressureOperator.rows.size(), 0.0);
    auto repeatedVelocity = velocity;
    auto repeatedPressure = pressure;
    const auto diagnostics =
        projectMovingPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            volumeRates, velocity, pressure, settings);
    const auto repeated =
        projectMovingPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            volumeRates, repeatedVelocity, repeatedPressure, settings);
    check(diagnostics == repeated
              && velocity == repeatedVelocity
              && pressure == repeatedPressure
              && diagnostics.accepted && diagnostics.finite
              && !diagnostics.staticGeometry
              && diagnostics.usesMovingVolumeRates
              && diagnostics.volumeRateFingerprint
                  == volumeRates.fingerprint
              && diagnostics.pressureSolve.compatible
              && diagnostics.pressureSolve.converged,
          "moving regional projection is deterministic and rate-bound");
    checkNear(
        diagnostics.maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond,
        0.1, 4.0e-16,
        "moving regional projection retains local geometry dV/dt");
    checkNear(
        diagnostics.predictedContinuityResidualMaximumCubicMetersPerSecond,
        0.1, 4.0e-16,
        "moving regional projection starts from geometry continuity demand");
    check(diagnostics.maximumAbsoluteVelocityCorrectionMetersPerSecond
              > 1.0e-3
              && diagnostics
                      .correctedContinuityResidualMaximumCubicMetersPerSecond
                  < 1.0e-11
              && diagnostics
                      .maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond
                  < 1.0e-11,
          "moving regional projection redistributes rigid-motion volume");
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            checkNear(velocity[link.linkIndex], 0.0, 0.0,
                      "moving regional projection keeps wall-relative flow zero");
        }
    }

    const auto beforeMetricVelocity =
        regionalMetricVelocityFromTopologyLinks(
            metric, beforeProjectionVelocity, &volumeRates);
    const auto afterMetricVelocity =
        regionalMetricVelocityFromTopologyLinks(
            metric, velocity, &volumeRates);
    const auto beforeState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            beforeMetricVelocity, settings.densityKgPerCubicMeter);
    const auto afterState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            afterMetricVelocity, settings.densityKgPerCubicMeter);
    PlanarPressureRegionFragmentProjectionEnergySettings energySettings;
    energySettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    energySettings.timeStepSeconds = settings.timeStepSeconds;
    energySettings.absoluteContinuityToleranceCubicMetersPerSecond =
        settings.absoluteContinuityToleranceCubicMetersPerSecond;
    energySettings.relativeContinuityTolerance =
        settings.relativeContinuityTolerance;
    const auto energy =
        auditMovingPlanarPressureRegionFragmentProjectionEnergy(
            geometry, sweep, fragments, topology, volumeRates, metric,
            beforeState, afterState, pressure, energySettings);
    const auto repeatedEnergy =
        auditMovingPlanarPressureRegionFragmentProjectionEnergy(
            geometry, sweep, fragments, topology, volumeRates, metric,
            beforeState, afterState, pressure, energySettings);
    check(energy == repeatedEnergy && energy.accepted
              && !energy.staticGeometry
              && energy.usesMovingVolumeRates
              && energy.volumeRateFingerprint == volumeRates.fingerprint
              && energy.pressureLayerTraceCount == 16
              && energy
                      .maximumAbsoluteWallTraceVelocityResidualMetersPerSecond
                  == 0.0
              && energy
                      .maximumAbsoluteGeometryVolumeRateCubicMetersPerSecond
                  == 0.1
              && energy.kineticEnergyChangeJoules > 0.0
              && !energy.nonIncreasingKineticEnergy
              && energy.correctionKineticEnergyJoules > 0.0
              && energy.geometryPressureWorkJoules > 0.0,
          "moving projection energy audit binds material walls and affine volume work");
    checkNear(energy.kineticEnergyChangeJoules,
              energy.correctionKineticEnergyJoules, 2.0e-13,
              "zero-flow moving projection adds exactly correction kinetic energy");
    checkNear(energy.geometryPressureWorkJoules,
              2.0 * energy.correctionKineticEnergyJoules, 4.0e-13,
              "moving geometry supplies twice the correction energy at the final state");
    checkNear(energy.finalPressureWorkJoules,
              energy.geometryPressureWorkJoules, 4.0e-13,
              "final pressure work closes against geometry volume work");
    checkNear(energy.midpointPressureWorkJoules,
              energy.kineticEnergyChangeJoules, 2.0e-13,
              "moving midpoint pressure work closes kinetic-energy change");
    checkNear(energy.affineEnergyResidualJoules, 0.0, 4.0e-13,
              "moving affine pressure-work identity closes");
    checkNear(energy.finalGeometryWorkResidualJoules, 0.0, 4.0e-13,
              "moving final pressure work closes geometry work");
    validateMovingPlanarPressureRegionFragmentProjectionEnergyAudit(
        energy, geometry, sweep, fragments, topology, volumeRates, metric,
        beforeState, afterState);
    auto corruptEnergy = energy;
    corruptEnergy.volumeRateFingerprint = 0;
    expectRejected(
        [&] {
            validateMovingPlanarPressureRegionFragmentProjectionEnergyAudit(
                corruptEnergy, geometry, sweep, fragments, topology,
                volumeRates, metric, beforeState, afterState);
        },
        "moving projection energy audit rejects volume-rate binding corruption");
    auto wrongWallVelocity = beforeMetricVelocity;
    const auto firstWallDof = std::ranges::find_if(
        metric.dofs,
        [](const auto& dof) {
            return dof.kind
                != PlanarPressureRegionFragmentVelocityDofKind::
                    SharedRegionGrid;
        });
    wrongWallVelocity[firstWallDof->dofIndex] = 0.0;
    const auto wrongWallState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            wrongWallVelocity, settings.densityKgPerCubicMeter);
    expectRejected(
        [&] { static_cast<void>(
            auditMovingPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, volumeRates, metric,
                wrongWallState, afterState, pressure, energySettings)); },
        "moving projection energy audit rejects a trace foreign to wall motion");
    auto mismatchedEnergySettings = energySettings;
    mismatchedEnergySettings.timeStepSeconds = 0.5;
    expectRejected(
        [&] { static_cast<void>(
            auditMovingPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, volumeRates, metric,
                beforeState, afterState, pressure,
                mismatchedEnergySettings)); },
        "moving projection energy audit rejects a mismatched duration");
    auto corruptEnergyRates = volumeRates;
    corruptEnergyRates.fragments[0]
        .geometryVolumeChangeRateCubicMetersPerSecond += 0.1;
    expectRejected(
        [&] { static_cast<void>(
            auditMovingPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, corruptEnergyRates,
                metric, beforeState, afterState, pressure,
                energySettings)); },
        "moving projection energy audit rejects corrupted volume rates");
    auto energyLimits =
        PlanarPressureRegionFragmentProjectionEnergyLimits{};
    energyLimits.volumeRateLimits.maximumFragments = 1;
    expectRejected(
        [&] { static_cast<void>(
            auditMovingPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, volumeRates, metric,
                beforeState, afterState, pressure, energySettings,
                energyLimits)); },
        "moving projection energy audit enforces nested volume-rate limits");

    const auto balancedVelocity = velocity;
    pressure.assign(pressure.size(), 0.0);
    const auto balanced =
        projectMovingPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            volumeRates, velocity, pressure, settings);
    double maximumBalancedVelocityChange = 0.0;
    for (std::size_t index = 0; index < velocity.size(); ++index) {
        maximumBalancedVelocityChange = std::max(
            maximumBalancedVelocityChange,
            std::abs(velocity[index] - balancedVelocity[index]));
    }
    check(balanced.accepted
              && balanced
                      .predictedContinuityResidualMaximumCubicMetersPerSecond
                  < 1.0e-11
              && maximumBalancedVelocityChange < 1.0e-10,
          "moving regional projection preserves a prebalanced flow field");

    auto breathingLayers = previousLayers;
    breathingLayers[0].physicalPlaneCoordinateMeters -= 0.1;
    breathingLayers[1].physicalPlaneCoordinateMeters += 0.1;
    const auto breathingSweep = makePlanarPressureRegionSweepLedger(
        geometry, previousLayers, breathingLayers, 0.5);
    const auto breathingFragments = buildPlanarPressureRegionFragments(
        geometry, breathingSweep);
    const auto breathingTopology = buildPlanarPressureRegionFragmentTopology(
        geometry, breathingSweep, breathingFragments);
    const auto breathingOperator =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, breathingSweep, breathingFragments,
            breathingTopology);
    const auto breathingRates =
        buildPlanarPressureRegionFragmentVolumeRates(
            geometry, breathingSweep, breathingFragments,
            breathingTopology);
    auto breathingSettings = settings;
    breathingSettings.timeStepSeconds = 0.5;
    std::vector<double> breathingVelocity(
        breathingTopology.links.size(), 0.0);
    std::vector<double> breathingPressure(
        breathingOperator.rows.size(), 0.25);
    const auto originalBreathingVelocity = breathingVelocity;
    const auto originalBreathingPressure = breathingPressure;
    const auto breathing =
        projectMovingPlanarPressureRegionFragmentFaceVelocities(
            breathingOperator, geometry, breathingSweep,
            breathingFragments, breathingTopology, breathingRates,
            breathingVelocity, breathingPressure, breathingSettings);
    check(!breathing.accepted && breathing.finite
              && !breathing.pressureSolve.compatible
              && !breathing.pressureSolve.converged
              && breathingVelocity == originalBreathingVelocity
              && breathingPressure == originalBreathingPressure,
          "sealed breathing projection rejects incompatible component volume");
    checkNear(
        breathing.maximumAbsolutePredictedComponentContinuityResidualCubicMetersPerSecond,
        1.6, 4.0e-15,
        "sealed breathing exposes its component continuity deficit");
    checkNear(
        breathing.pressureSolve
            .maximumAbsoluteComponentCompatibilityPascalsMeters,
        3.84, 1.0e-14,
        "sealed breathing maps volume deficit into the physical pressure RHS");
    const auto breathingMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, breathingSweep, breathingFragments,
            breathingTopology);
    const auto breathingMetricVelocity =
        regionalMetricVelocityFromTopologyLinks(
            breathingMetric, originalBreathingVelocity, &breathingRates);
    const auto breathingState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, breathingSweep, breathingFragments,
            breathingTopology, breathingMetric, breathingMetricVelocity,
            breathingSettings.densityKgPerCubicMeter);
    const std::vector<double> zeroBreathingPressure(
        breathingFragments.fragments.size(), 0.0);
    auto breathingEnergySettings = energySettings;
    breathingEnergySettings.timeStepSeconds =
        breathingSettings.timeStepSeconds;
    expectRejected(
        [&] { static_cast<void>(
            auditMovingPlanarPressureRegionFragmentProjectionEnergy(
                geometry, breathingSweep, breathingFragments,
                breathingTopology, breathingRates, breathingMetric,
                breathingState, breathingState, zeroBreathingPressure,
                breathingEnergySettings)); },
        "moving projection energy audit rejects incompatible breathing continuity");

    auto truncatedSettings = settings;
    truncatedSettings.pressureSolve.absoluteResidualTolerancePascalsMeters =
        1.0e-16;
    truncatedSettings.pressureSolve.relativeResidualTolerance = 0.0;
    truncatedSettings.pressureSolve.maximumIterations = 1;
    velocity.assign(velocity.size(), 0.0);
    pressure.assign(pressure.size(), 0.25);
    const auto originalVelocity = velocity;
    const auto originalPressure = pressure;
    const auto truncated =
        projectMovingPlanarPressureRegionFragmentFaceVelocities(
            pressureOperator, geometry, sweep, fragments, topology,
            volumeRates, velocity, pressure, truncatedSettings);
    check(!truncated.accepted && truncated.pressureSolve.compatible
              && !truncated.pressureSolve.converged
              && velocity == originalVelocity
              && pressure == originalPressure,
          "truncated moving regional projection rolls back both fields");

    auto mismatchedSettings = settings;
    mismatchedSettings.timeStepSeconds = 0.5;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectMovingPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    volumeRates, candidateVelocity, candidatePressure,
                    mismatchedSettings));
        },
        "moving regional projection rejects a mismatched epoch duration");
    auto corruptRates = volumeRates;
    corruptRates.fragments[0]
        .geometryVolumeChangeRateCubicMetersPerSecond += 0.1;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectMovingPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    corruptRates, candidateVelocity, candidatePressure,
                    settings));
        },
        "moving regional projection rejects corrupted volume rates");
    auto limits =
        PlanarPressureRegionFragmentPressureProjectionLimits{};
    limits.volumeRateLimits.maximumFragments = 1;
    expectRejected(
        [&] {
            auto candidateVelocity = originalVelocity;
            auto candidatePressure = originalPressure;
            static_cast<void>(
                projectMovingPlanarPressureRegionFragmentFaceVelocities(
                    pressureOperator, geometry, sweep, fragments, topology,
                    volumeRates, candidateVelocity, candidatePressure,
                    settings, limits));
        },
        "moving regional projection enforces nested volume-rate limits");
}

void testPlanarRegionalMovingFragmentPressureProjectionAllAxes() {
    const auto geometry = grid();
    PlanarPressureRegionFragmentPressureProjectionSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 1.0;
    settings.pressureSolve.absoluteResidualTolerancePascalsMeters = 1.0e-13;
    settings.pressureSolve.relativeResidualTolerance = 1.0e-12;
    settings.pressureSolve.maximumIterations = 200;
    PlanarPressureRegionFragmentProjectionEnergySettings energySettings;
    energySettings.densityKgPerCubicMeter = settings.densityKgPerCubicMeter;
    energySettings.timeStepSeconds = settings.timeStepSeconds;
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
        const auto current = translatePlanarPressureJumpLayers(
            geometry, previous, 0.1).layers;
        const auto sweep = makePlanarPressureRegionSweepLedger(
            geometry, previous, current, 1.0);
        const auto fragments = buildPlanarPressureRegionFragments(
            geometry, sweep);
        const auto topology = buildPlanarPressureRegionFragmentTopology(
            geometry, sweep, fragments);
        const auto metric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology);
        const auto pressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology);
        const auto volumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology);
        std::vector<double> velocity(topology.links.size(), 0.0);
        const auto beforeLinks = velocity;
        std::vector<double> pressure(pressureOperator.rows.size(), 0.0);
        const auto diagnostics =
            projectMovingPlanarPressureRegionFragmentFaceVelocities(
                pressureOperator, geometry, sweep, fragments, topology,
                volumeRates, velocity, pressure, settings);
        const auto beforeVelocity =
            regionalMetricVelocityFromTopologyLinks(
                metric, beforeLinks, &volumeRates);
        const auto afterVelocity =
            regionalMetricVelocityFromTopologyLinks(
                metric, velocity, &volumeRates);
        const auto before =
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric,
                beforeVelocity, settings.densityKgPerCubicMeter);
        const auto after =
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, sweep, fragments, topology, metric,
                afterVelocity, settings.densityKgPerCubicMeter);
        const auto energy =
            auditMovingPlanarPressureRegionFragmentProjectionEnergy(
                geometry, sweep, fragments, topology, volumeRates, metric,
                before, after, pressure, energySettings);
        check(diagnostics.accepted
                  && diagnostics.usesMovingVolumeRates
                  && diagnostics
                          .correctedContinuityResidualMaximumCubicMetersPerSecond
                      < 1.0e-11
                  && energy.accepted
                  && energy.usesMovingVolumeRates
                  && energy
                          .maximumAbsoluteWallTraceVelocityResidualMetersPerSecond
                      == 0.0
                  && std::abs(energy.affineEnergyResidualJoules)
                      < 4.0e-13
                  && std::abs(energy.finalGeometryWorkResidualJoules)
                      < 4.0e-13,
              "moving regional projection and affine energy close on every axis");
    }
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
    testPlanarRegionalFragmentTopology();
    testPlanarRegionalFragmentTopologyAxesAndRejection();
    testPlanarRegionalFragmentVolumeRates();
    testPlanarRegionalFragmentVolumeRatesAxesAndRejection();
    testPlanarRegionalFragmentVelocityMetric();
    testPlanarRegionalFragmentVelocityMetricAxesAndRejection();
    testPlanarRegionalFragmentVelocityState();
    testPlanarRegionalFragmentVelocityStateMotionAndRejection();
    testPlanarRegionalFragmentPressureOperator();
    testPlanarRegionalFragmentPressureOperatorAxesAndRejection();
    testPlanarRegionalFragmentPressureCorrectionSolve();
    testPlanarRegionalFragmentPressureCorrectionRollback();
    testPlanarRegionalFragmentPressureProjection();
    testPlanarRegionalFragmentProjectionEnergyAudit();
    testPlanarRegionalFragmentProjectionEnergyAxesAndRejection();
    testPlanarRegionalFragmentPressureProjectionAxesAndRollback();
    testPlanarRegionalMovingFragmentPressureProjection();
    testPlanarRegionalMovingFragmentPressureProjectionAllAxes();
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
