#include "fluid/planar_pressure_jump.h"
#include "fluid/planar_region_flux.h"
#include "fluid/planar_region_opening_flow.h"
#include "fluid/planar_region_opening_power.h"
#include "fluid/planar_region_fragment.h"
#include "fluid/planar_region_fragment_accepted_state.h"
#include "fluid/planar_region_fragment_opening_accepted_state.h"
#include "fluid/planar_region_fragment_opening_accepted_state_persistence.h"
#include "fluid/planar_region_fragment_opening_continuation.h"
#include "fluid/planar_region_fragment_opening_continuation_momentum_audit.h"
#include "fluid/planar_region_fragment_opening_velocity_metric.h"
#include "fluid/planar_region_fragment_opening_velocity_state.h"
#include "fluid/planar_region_fragment_opening_load_state.h"
#include "fluid/planar_region_fragment_opening_momentum_transport.h"
#include "fluid/planar_region_fragment_opening_momentum_prediction.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle_state.h"
#include "fluid/planar_region_fragment_opening_momentum_cycle_state_persistence.h"
#include "fluid/planar_region_fragment_opening_momentum_pressure_epoch.h"
#include "fluid/planar_region_fragment_opening_pressure_state.h"
#include "fluid/planar_region_fragment_opening.h"
#include "fluid/planar_region_fragment_opening_flux.h"
#include "fluid/planar_region_fragment_opening_pressure_operator.h"
#include "fluid/planar_region_fragment_opening_pressure_epoch.h"
#include "fluid/planar_region_fragment_opening_pressure_projection.h"
#include "fluid/planar_region_fragment_opening_pressure_step.h"
#include "fluid/planar_region_fragment_opening_resistance.h"
#include "fluid/planar_region_fragment_opening_surface_load.h"
#include "fluid/planar_region_fragment_pressure_operator.h"
#include "fluid/planar_region_fragment_pressure_jump_energy.h"
#include "fluid/planar_region_fragment_pressure_projection.h"
#include "fluid/planar_region_fragment_projection_energy.h"
#include "fluid/planar_region_fragment_pressure_solve.h"
#include "fluid/planar_region_fragment_pressure_state.h"
#include "fluid/planar_region_fragment_surface_load.h"
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
#include <numeric>
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

void refreshOpeningAcceptedStatePersistenceChecksum(
    std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 24;
    constexpr std::size_t checksumOffset = 16;
    constexpr std::uint64_t offsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t value = offsetBasis;
    for (std::size_t index = envelopeBytes; index < bytes.size(); ++index) {
        value ^= bytes[index];
        value *= prime;
    }
    for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
        bytes[checksumOffset + byte] = static_cast<std::uint8_t>(
            value >> (8U * byte));
    }
}

void refreshOpeningMomentumCycleStatePersistenceChecksum(
    std::vector<std::uint8_t>& bytes) {
    constexpr std::size_t envelopeBytes = 24;
    constexpr std::size_t checksumOffset = 16;
    constexpr std::uint64_t offsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t value = offsetBasis;
    for (std::size_t index = envelopeBytes; index < bytes.size(); ++index) {
        value ^= bytes[index];
        value *= prime;
    }
    for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
        bytes[checksumOffset + byte] = static_cast<std::uint8_t>(
            value >> (8U * byte));
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

const PlanarPressureRegionFragmentPressureJumpEnergyComponent*
findFragmentPressureJumpEnergyComponent(
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& audit,
    const std::uint64_t regionStableId) {
    const auto found = std::ranges::find(
        audit.components, regionStableId,
        &PlanarPressureRegionFragmentPressureJumpEnergyComponent::
            regionStableId);
    return found == audit.components.end() ? nullptr : &*found;
}

double pressureJumpSurfaceWork(
    const PlanarPressureRegionFragmentPressureJumpEnergyAudit& audit,
    const std::uint64_t surfaceStableId) {
    double result = 0.0;
    for (const auto& layer : audit.layers) {
        if (layer.surfaceStableId == surfaceStableId) {
            result += layer.pressureJumpWorkToFluidJoules;
        }
    }
    return result;
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

template<typename PressureOperator, typename Component>
double fragmentCorrectionVolumeMean(
    const PressureOperator& pressureOperator,
    const PlanarPressureRegionFragmentSet& fragments,
    const Component& component,
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

template<typename PressureOperator>
void subtractFragmentCorrectionVolumeMeans(
    const PressureOperator& pressureOperator,
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

void testPlanarRegionalFragmentPressureJumpEnergyAudit() {
    constexpr double density = 1.2;
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 0.5);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    const auto state = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric,
        std::vector<double>(metric.dofs.size(), 0.0), density);
    PlanarPressureRegionFragmentPressureJumpEnergySettings settings;
    settings.timeStepSeconds = 0.5;

    const auto audit =
        auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, sweep, fragments, topology, metric, state, settings);
    const auto repeated =
        auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, sweep, fragments, topology, metric, state, settings);
    check(audit == repeated
              && audit.version
                  == planarPressureRegionFragmentPressureJumpEnergyVersion
              && audit.fingerprint != 0
              && audit.sourceMetricFingerprint == metric.fingerprint
              && audit.sourceVelocityStateFingerprint == state.fingerprint
              && audit.sourceTopologyFingerprint == topology.fingerprint
              && audit.sourceFragmentFingerprint == fragments.fingerprint
              && audit.volumeRateFingerprint == 0
              && audit.staticGeometry
              && !audit.usesMovingVolumeRates
              && audit.accepted
              && audit.layers.size() == 8
              && audit.components.size() == 2
              && audit.pressureLayerTraceCount == 16
              && audit.ownedStorageBytes > 0
              && audit.workingStorageBytes
                  == 2 * topology.links.size() * sizeof(std::size_t),
          "static regional pressure-jump energy audit is deterministic and source-bound");
    for (const auto& layer : audit.layers) {
        const double sign = layer.surfaceStableId == 10 ? 1.0 : -1.0;
        check(layer.minusTraceDofIndex != layer.plusTraceDofIndex
                  && layer.areaSquareMeters == 1.0
                  && layer.axis == GridFaceAxis::X
                  && layer.authoredPressureJumpPascals == sign * 70.0
                  && layer.reconstructedPressureJumpPascals
                      == layer.authoredPressureJumpPascals
                  && layer.pressureJumpResidualPascals == 0.0
                  && layer.materialWallVelocityMetersPerSecond == 0.0
                  && layer.minusTraceVelocityMetersPerSecond == 0.0
                  && layer.plusTraceVelocityMetersPerSecond == 0.0
                  && layer.maximumAbsoluteWallVelocityResidualMetersPerSecond
                      == 0.0,
              "static pressure-jump layer retains exact wall sources");
        checkNear(
            layer.authoredPressureJumpForceOnFluidNewtons.x,
            sign * 70.0, 0.0,
            "static pressure-jump layer reconstructs signed fluid force");
        checkNear(
            layer.pressureForceOnSheetNewtons.x,
            -sign * 70.0, 0.0,
            "static pressure-jump layer publishes opposite sheet force");
        checkNear(
            layer.pressureJumpImpulseOnFluidNewtonSeconds.x,
            sign * 35.0, 0.0,
            "static pressure-jump layer integrates fluid impulse");
        check(layer.resolvedPressureForceOnFluidNewtons
                      == layer.authoredPressureJumpForceOnFluidNewtons
                  && layer.pressureForceClosureResidualNewtons == Vector3{}
                  && layer.actionReactionForceResidualNewtons == Vector3{}
                  && layer.actionReactionImpulseResidualNewtonSeconds
                      == Vector3{}
                  && layer.resolvedPressurePowerToFluidWatts == 0.0
                  && layer.authoredPressureJumpPowerToFluidWatts == 0.0
                  && layer.pressureJumpWorkToFluidJoules == 0.0
                  && layer.pressureWorkToSheetJoules == 0.0,
              "static pressure-jump layer closes action, reaction, and zero work");
    }
    const auto* exterior = findFragmentPressureJumpEnergyComponent(audit, 1);
    const auto* pocket = findFragmentPressureJumpEnergyComponent(audit, 2);
    check(exterior != nullptr && pocket != nullptr,
          "pressure-jump audit retains both pressure components");
    if (exterior != nullptr && pocket != nullptr) {
        check(exterior->pressureLayerSideCount == 8
                  && pocket->pressureLayerSideCount == 8
                  && exterior->resolvedPressureForceOnFluidNewtons
                      == Vector3{}
                  && pocket->resolvedPressureForceOnFluidNewtons
                      == Vector3{}
                  && exterior->pressureWorkToFluidJoules == 0.0
                  && pocket->pressureWorkToFluidJoules == 0.0
                  && exterior->workGeometryResidualJoules == 0.0
                  && pocket->workGeometryResidualJoules == 0.0,
              "static pressure components close boundary force and work independently");
    }
    check(audit.resolvedPressureForceOnFluidNewtons == Vector3{}
              && audit.authoredPressureJumpForceOnFluidNewtons == Vector3{}
              && audit.pressureForceOnSheetNewtons == Vector3{}
              && audit.actionReactionForceResidualNewtons == Vector3{}
              && audit.pressureJumpImpulseOnFluidNewtonSeconds == Vector3{}
              && audit.pressureJumpWorkToFluidJoules == 0.0
              && audit.pressureWorkToSheetJoules == 0.0
              && audit.geometryPressureWorkToFluidJoules == 0.0
              && audit.workGeometryResidualJoules == 0.0,
          "static authored pressure jumps are globally internal and do no work");
    validateStaticPlanarPressureRegionFragmentPressureJumpEnergyAudit(
        audit, geometry, sweep, fragments, topology, metric, state);

    auto corrupt = audit;
    corrupt.fingerprint = 0;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentPressureJumpEnergyAudit(
                corrupt, geometry, sweep, fragments, topology, metric,
                state);
        },
        "pressure-jump energy audit rejects fingerprint corruption");
    corrupt = audit;
    corrupt.layers[0].pressureJumpWorkToFluidJoules += 0.1;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentPressureJumpEnergyAudit(
                corrupt, geometry, sweep, fragments, topology, metric,
                state);
        },
        "pressure-jump energy audit rejects layer corruption");
    corrupt = audit;
    corrupt.components[0].pressureWorkToFluidJoules += 0.1;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentPressureJumpEnergyAudit(
                corrupt, geometry, sweep, fragments, topology, metric,
                state);
        },
        "pressure-jump energy audit rejects component corruption");

    auto movingWallVelocity = std::vector<double>(metric.dofs.size(), 0.0);
    const auto wallDof = std::ranges::find_if(
        metric.dofs,
        [](const auto& dof) {
            return dof.kind
                != PlanarPressureRegionFragmentVelocityDofKind::
                    SharedRegionGrid;
        });
    movingWallVelocity[wallDof->dofIndex] = 0.01;
    const auto inconsistentWallState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, sweep, fragments, topology, metric,
            movingWallVelocity, density);
    expectRejected(
        [&] {
            static_cast<void>(
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric,
                    inconsistentWallState, settings));
        },
        "static pressure-jump energy audit rejects a moving wall trace");

    auto invalidSettings = settings;
    invalidSettings.timeStepSeconds = 0.0;
    expectRejected(
        [&] {
            static_cast<void>(
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric, state,
                    invalidSettings));
        },
        "pressure-jump energy audit rejects an invalid time step");
    auto limits = PlanarPressureRegionFragmentPressureJumpEnergyLimits{};
    limits.maximumLayers = audit.layers.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric, state,
                    settings, limits));
        },
        "pressure-jump energy audit enforces the layer limit");
    limits = {};
    limits.maximumComponents = audit.components.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric, state,
                    settings, limits));
        },
        "pressure-jump energy audit enforces the component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric, state,
                    settings, limits));
        },
        "pressure-jump energy audit enforces the owned byte limit");
    limits = {};
    limits.maximumWorkingBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric, state,
                    settings, limits));
        },
        "pressure-jump energy audit enforces the working byte limit");
    limits = {};
    limits.velocityStateLimits.maximumSamples = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, metric, state,
                    settings, limits));
        },
        "pressure-jump energy audit enforces nested velocity-state limits");
}

void testPlanarRegionalMovingPressureJumpEnergyAudit() {
    constexpr double density = 1.2;
    const auto geometry = grid();
    const auto previous = pocketLayers();
    const auto current = translatePlanarPressureJumpLayers(
        geometry, previous, 0.1).layers;
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, current, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto rates = buildPlanarPressureRegionFragmentVolumeRates(
        geometry, sweep, fragments, topology);
    const auto metric = buildPlanarPressureRegionFragmentVelocityMetric(
        geometry, sweep, fragments, topology);
    const auto wallVelocity = regionalMetricVelocityFromTopologyLinks(
        metric, std::vector<double>(topology.links.size(), 0.0), &rates);
    const auto state = buildPlanarPressureRegionFragmentVelocityState(
        geometry, sweep, fragments, topology, metric, wallVelocity, density);
    PlanarPressureRegionFragmentPressureJumpEnergySettings settings;
    settings.timeStepSeconds = 1.0;
    const auto audit =
        auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, sweep, fragments, topology, rates, metric, state,
            settings);
    const auto repeated =
        auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, sweep, fragments, topology, rates, metric, state,
            settings);
    check(audit == repeated && audit.accepted && !audit.staticGeometry
              && audit.usesMovingVolumeRates
              && audit.volumeRateFingerprint == rates.fingerprint
              && audit.maximumAbsolutePressureJumpResidualPascals == 0.0
              && audit.maximumAbsoluteWallVelocityResidualMetersPerSecond
                  == 0.0,
          "moving pressure-jump energy audit is deterministic and rate-bound");
    checkNear(pressureJumpSurfaceWork(audit, 10), 28.0, 3.0e-14,
              "rigid translation receives positive work at the first layer");
    checkNear(pressureJumpSurfaceWork(audit, 20), -28.0, 3.0e-14,
              "rigid translation returns work at the second layer");
    checkNear(audit.pressureJumpWorkToFluidJoules, 0.0, 6.0e-14,
              "rigid translation has zero net authored pressure work");
    checkNear(audit.geometryPressureWorkToFluidJoules, 0.0, 6.0e-14,
              "rigid translation has zero net geometry pressure work");
    checkNear(audit.workGeometryResidualJoules, 0.0, 6.0e-14,
              "rigid translation closes authored and geometry work");
    check(audit.authoredPressureJumpForceOnFluidNewtons == Vector3{}
              && audit.pressureForceOnSheetNewtons == Vector3{}
              && audit.actionReactionForceResidualNewtons == Vector3{},
          "rigid pressure-jump forces remain globally internal");
    validateMovingPlanarPressureRegionFragmentPressureJumpEnergyAudit(
        audit, geometry, sweep, fragments, topology, rates, metric, state);

    auto breathing = previous;
    breathing[0].physicalPlaneCoordinateMeters -= 0.1;
    breathing[1].physicalPlaneCoordinateMeters += 0.1;
    const auto breathingSweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, breathing, 0.5);
    const auto breathingFragments = buildPlanarPressureRegionFragments(
        geometry, breathingSweep);
    const auto breathingTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, breathingSweep, breathingFragments);
    const auto breathingRates =
        buildPlanarPressureRegionFragmentVolumeRates(
            geometry, breathingSweep, breathingFragments,
            breathingTopology);
    const auto breathingMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, breathingSweep, breathingFragments,
            breathingTopology);
    const auto breathingVelocity = regionalMetricVelocityFromTopologyLinks(
        breathingMetric,
        std::vector<double>(breathingTopology.links.size(), 0.0),
        &breathingRates);
    const auto breathingState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, breathingSweep, breathingFragments,
            breathingTopology, breathingMetric, breathingVelocity,
            density);
    PlanarPressureRegionFragmentPressureJumpEnergySettings breathingSettings;
    breathingSettings.timeStepSeconds = 0.5;
    const auto breathingAudit =
        auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, breathingSweep, breathingFragments,
            breathingTopology, breathingRates, breathingMetric,
            breathingState, breathingSettings);
    check(breathingAudit.accepted
              && breathingAudit.authoredPressureJumpForceOnFluidNewtons
                  == Vector3{}
              && breathingAudit.pressureForceOnSheetNewtons == Vector3{},
          "sealed breathing retains zero aggregate pressure force");
    checkNear(pressureJumpSurfaceWork(breathingAudit, 10), -28.0, 3.0e-14,
              "breathing transfers work through the first moving layer");
    checkNear(pressureJumpSurfaceWork(breathingAudit, 20), -28.0, 3.0e-14,
              "breathing transfers work through the second moving layer");
    checkNear(breathingAudit.pressureJumpWorkToFluidJoules,
              -56.0, 8.0e-14,
              "breathing pressure jump removes the analytic fluid work");
    checkNear(breathingAudit.pressureWorkToSheetJoules,
              56.0, 8.0e-14,
              "breathing pressure jump gives opposite work to the sheet");
    checkNear(breathingAudit.geometryPressureWorkToFluidJoules,
              -56.0, 8.0e-14,
              "breathing geometry reconstructs analytic pressure work");
    checkNear(breathingAudit.workGeometryResidualJoules,
              0.0, 8.0e-14,
              "breathing closes authored pressure work to geometry");
    const auto* breathingExterior =
        findFragmentPressureJumpEnergyComponent(breathingAudit, 1);
    const auto* breathingPocket =
        findFragmentPressureJumpEnergyComponent(breathingAudit, 2);
    check(breathingExterior != nullptr && breathingPocket != nullptr,
          "breathing pressure-jump audit retains both components");
    if (breathingExterior != nullptr && breathingPocket != nullptr) {
        checkNear(breathingExterior->geometryPressureWorkToFluidJoules,
                  -11.2, 4.0e-14,
                  "breathing exterior closes its gauge-dependent work");
        checkNear(breathingPocket->geometryPressureWorkToFluidJoules,
                  -44.8, 6.0e-14,
                  "breathing pocket closes its pressure work");
        checkNear(breathingExterior->workGeometryResidualJoules,
                  0.0, 4.0e-14,
                  "breathing exterior wall work matches volume work");
        checkNear(breathingPocket->workGeometryResidualJoules,
                  0.0, 6.0e-14,
                  "breathing pocket wall work matches volume work");
    }

    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> axisPrevious{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0},
             -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0},
             -0.2, -70.0},
        };
        const auto axisCurrent = translatePlanarPressureJumpLayers(
            geometry, axisPrevious, 0.1).layers;
        const auto axisSweep = makePlanarPressureRegionSweepLedger(
            geometry, axisPrevious, axisCurrent, 1.0);
        const auto axisFragments = buildPlanarPressureRegionFragments(
            geometry, axisSweep);
        const auto axisTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, axisSweep, axisFragments);
        const auto axisRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, axisSweep, axisFragments, axisTopology);
        const auto axisMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, axisSweep, axisFragments, axisTopology);
        const auto axisVelocity = regionalMetricVelocityFromTopologyLinks(
            axisMetric,
            std::vector<double>(axisTopology.links.size(), 0.0),
            &axisRates);
        const auto axisState =
            buildPlanarPressureRegionFragmentVelocityState(
                geometry, axisSweep, axisFragments, axisTopology,
                axisMetric, axisVelocity, density);
        const auto axisAudit =
            auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                geometry, axisSweep, axisFragments, axisTopology,
                axisRates, axisMetric, axisState, settings);
        check(axisAudit.accepted && axisAudit.layers.size()
                  == axisTopology.pressureLayerWallLinkCount
                  && axisAudit.authoredPressureJumpForceOnFluidNewtons
                      == Vector3{},
              "all-axis rigid pressure-jump audit retains internal force closure");
        checkNear(pressureJumpSurfaceWork(axisAudit, 10),
                  axis == GridFaceAxis::X ? 28.0 : 56.0,
                  8.0e-14,
                  "all-axis rigid pressure-jump audit closes first-layer work");
        checkNear(axisAudit.pressureJumpWorkToFluidJoules,
                  0.0, 2.0e-13,
                  "all-axis rigid pressure-jump audit closes net work");
    }

    auto corrupt = breathingAudit;
    corrupt.layers[0].authoredPressureJumpPascals += 1.0;
    expectRejected(
        [&] {
            validateMovingPlanarPressureRegionFragmentPressureJumpEnergyAudit(
                corrupt, geometry, breathingSweep, breathingFragments,
                breathingTopology, breathingRates, breathingMetric,
                breathingState);
        },
        "moving pressure-jump energy audit rejects layer corruption");
    auto wrongBreathingVelocity = breathingVelocity;
    const auto breathingWallDof = std::ranges::find_if(
        breathingMetric.dofs,
        [](const auto& dof) {
            return dof.kind
                != PlanarPressureRegionFragmentVelocityDofKind::
                    SharedRegionGrid;
        });
    wrongBreathingVelocity[breathingWallDof->dofIndex] += 0.01;
    const auto wrongBreathingState =
        buildPlanarPressureRegionFragmentVelocityState(
            geometry, breathingSweep, breathingFragments,
            breathingTopology, breathingMetric, wrongBreathingVelocity,
            density);
    expectRejected(
        [&] {
            static_cast<void>(
                auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology, breathingRates, breathingMetric,
                    wrongBreathingState, breathingSettings));
        },
        "moving pressure-jump energy audit rejects a wrong wall trace");
    auto wrongDuration = breathingSettings;
    wrongDuration.timeStepSeconds = 1.0;
    expectRejected(
        [&] {
            static_cast<void>(
                auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology, breathingRates, breathingMetric,
                    breathingState, wrongDuration));
        },
        "moving pressure-jump energy audit rejects a duration mismatch");
    auto corruptRates = breathingRates;
    corruptRates.fingerprint = 0;
    expectRejected(
        [&] {
            static_cast<void>(
                auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology, corruptRates, breathingMetric,
                    breathingState, breathingSettings));
        },
        "moving pressure-jump energy audit rejects corrupted volume rates");
    auto limits = PlanarPressureRegionFragmentPressureJumpEnergyLimits{};
    limits.volumeRateLimits.maximumFragments = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology, breathingRates, breathingMetric,
                    breathingState, breathingSettings, limits));
        },
        "moving pressure-jump energy audit enforces nested volume-rate limits");
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

void testPlanarRegionalFragmentOpeningPressureOperator() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto base = buildPlanarPressureRegionFragmentPressureOperator(
        geometry, sweep, fragments, topology);
    const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions{{100, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.5}};
    const auto openings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, definitions);
    const auto op = buildPlanarPressureRegionFragmentOpeningPressureOperator(
        base, geometry, sweep, fragments, topology, definitions, openings);
    const auto repeated = buildPlanarPressureRegionFragmentOpeningPressureOperator(
        base, geometry, sweep, fragments, topology, definitions, openings);
    check(op == repeated && op.fingerprint != 0
              && op.sourceBaseOperatorFingerprint == base.fingerprint
              && op.sourceOpeningFingerprint == openings.fingerprint
              && op.rows.size() == 24 && op.entries.size() == 130
              && op.components.size() == 1
              && op.componentFragmentIndices.size() == 24
              && op.includedSameRegionGridLinkCount == 64
              && op.includedOpeningPatchCount == 1,
          "opening pressure operator is deterministic and joins sealed components");
    checkNear(op.totalPressureLayerWallAreaSquareMeters, 8.0, 2.0e-15,
              "opening pressure operator retains total wall area");
    checkNear(op.totalOpeningAreaSquareMeters, 0.5, 2.0e-15,
              "opening pressure operator includes exact aperture area");
    checkNear(op.totalSolidPressureLayerWallAreaSquareMeters, 7.5, 2.0e-15,
              "opening pressure operator excludes remaining solid wall area");
    checkNear(op.wallAreaPartitionResidualSquareMeters, 0.0, 2.0e-15,
              "opening pressure operator closes wall area");
    checkNear(op.openingGeometryWeightMeters, 1.25, 2.0e-15,
              "opening pressure operator uses aperture area over center distance");
    checkNear(op.totalGeometryWeightMeters, 655.0 / 12.0, 4.0e-14,
              "opening pressure operator adds aperture and grid weights");
    checkNear(op.totalDiagonalGeometryWeightMeters, 655.0 / 6.0, 8.0e-14,
              "opening pressure operator closes two-sided diagonal weight");

    std::vector<double> constant(op.rows.size(), 3.0);
    const auto constantAction =
        applyPlanarPressureRegionFragmentOpeningPressureOperator(op, constant);
    check(std::ranges::all_of(
              constantAction,
              [](const double value) { return std::abs(value) < 3.0e-14; }),
          "opening pressure operator retains one connected constant null mode");
    std::vector<double> regional;
    regional.reserve(fragments.fragments.size());
    for (const auto& fragment : fragments.fragments)
        regional.push_back(fragment.pressurePascals);
    const auto regionalAction =
        applyPlanarPressureRegionFragmentOpeningPressureOperator(op, regional);
    checkNear(dotProduct(regional, regionalAction), 6125.0, 2.0e-12,
              "opening pressure operator exposes authored pressure across aperture");
    checkNear(std::accumulate(regionalAction.begin(), regionalAction.end(), 0.0),
              0.0, 2.0e-13,
              "opening pressure operator has zero integrated connected action");
    std::vector<double> first(op.rows.size());
    std::vector<double> second(op.rows.size());
    for (std::size_t index = 0; index < first.size(); ++index) {
        first[index] = 0.25 * static_cast<double>(index) - 1.0;
        second[index] = std::sin(0.3 * static_cast<double>(index));
    }
    const auto firstAction =
        applyPlanarPressureRegionFragmentOpeningPressureOperator(op, first);
    const auto secondAction =
        applyPlanarPressureRegionFragmentOpeningPressureOperator(op, second);
    checkNear(dotProduct(first, secondAction), dotProduct(second, firstAction),
              2.0e-12,
              "opening pressure operator is symmetric");
    check(dotProduct(first, firstAction) >= 0.0,
          "opening pressure operator has nonnegative graph energy");
    validatePlanarPressureRegionFragmentOpeningPressureOperator(
        op, base, geometry, sweep, fragments, topology, definitions, openings);

    const auto emptyOpenings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, {});
    const auto empty = buildPlanarPressureRegionFragmentOpeningPressureOperator(
        base, geometry, sweep, fragments, topology, {}, emptyOpenings);
    const auto baseAction = applyPlanarPressureRegionFragmentPressureOperator(
        base, regional);
    bool emptyMatchesBase = empty.entries.size() == base.entries.size();
    const auto emptyAction =
        applyPlanarPressureRegionFragmentOpeningPressureOperator(
            empty, regional);
    for (std::size_t index = 0; index < baseAction.size(); ++index) {
        emptyMatchesBase = emptyMatchesBase
            && std::abs(emptyAction[index] - baseAction[index]) < 3.0e-13;
    }
    check(emptyMatchesBase
              && empty.components.size() == base.components.size()
              && empty.openingGeometryWeightMeters == 0.0,
          "empty opening operator preserves the sealed base action within roundoff");

    const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition> shared{
        {102, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.2},
        {101, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.3}};
    auto reversedShared = shared;
    std::ranges::reverse(reversedShared);
    const auto sharedOpenings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, shared);
    const auto reversedOpenings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, reversedShared);
    const auto sharedOp = buildPlanarPressureRegionFragmentOpeningPressureOperator(
        base, geometry, sweep, fragments, topology, shared, sharedOpenings);
    const auto reversedOp = buildPlanarPressureRegionFragmentOpeningPressureOperator(
        base, geometry, sweep, fragments, topology, reversedShared,
        reversedOpenings);
    check(sharedOp == reversedOp && sharedOp.entries.size() == 132,
          "opening pressure operator canonicalizes parallel patches on one tile");
    checkNear(sharedOp.openingGeometryWeightMeters, 1.25, 2.0e-15,
              "parallel opening patches add their geometric weights");

    for (const GridFaceAxis axis
         : {GridFaceAxis::X, GridFaceAxis::Y, GridFaceAxis::Z}) {
        const std::size_t firstFace = axis == GridFaceAxis::X ? 1 : 0;
        const std::size_t secondFace = axis == GridFaceAxis::X ? 2 : 1;
        const std::vector<PlanarPressureJumpLayerDefinition> axisLayers{
            {10, 1, 2,
             {movingPlanarFaceTopologyVersion, axis, firstFace, 0}, -0.8, 70.0},
            {20, 2, 1,
             {movingPlanarFaceTopologyVersion, axis, secondFace, 0}, -0.2, -70.0}};
        const auto axisSweep = makePlanarPressureRegionSweepLedger(
            geometry, axisLayers, axisLayers, 1.0);
        const auto axisFragments = buildPlanarPressureRegionFragments(
            geometry, axisSweep);
        const auto axisTopology = buildPlanarPressureRegionFragmentTopology(
            geometry, axisSweep, axisFragments);
        const auto axisBase = buildPlanarPressureRegionFragmentPressureOperator(
            geometry, axisSweep, axisFragments, axisTopology);
        const auto wall = std::ranges::find_if(
            axisTopology.links, [](const auto& link) {
                return link.kind
                    == PlanarPressureRegionFragmentFaceKind::PressureLayerWall;
            });
        check(wall != axisTopology.links.end(),
              "opening pressure operator finds a wall on every axis");
        if (wall == axisTopology.links.end()) continue;
        const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
            axisDefinitions{{100, 1000, wall->surfaceStableId, wall->axis,
                wall->i, wall->j, wall->k, wall->minusRegionStableId,
                wall->plusRegionStableId, 0.5 * wall->areaSquareMeters}};
        const auto axisOpenings = buildPlanarPressureRegionFragmentOpenings(
            geometry, axisSweep, axisFragments, axisTopology, axisDefinitions);
        const auto axisOp =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                axisBase, geometry, axisSweep, axisFragments, axisTopology,
                axisDefinitions, axisOpenings);
        check(axisOp.components.size() == 1
                  && axisOp.entries.size() == axisBase.entries.size() + 2,
              "opening pressure operator joins components on every axis");
        checkNear(axisOp.openingGeometryWeightMeters,
                  0.5 * wall->areaSquareMeters / wall->centerDistanceMeters,
                  4.0e-15,
                  "opening pressure operator uses physical geometry on every axis");
    }

    auto corrupt = op;
    corrupt.entries.back().geometryWeightMeters += 0.1;
    expectRejected([&] {
        validatePlanarPressureRegionFragmentOpeningPressureOperator(
            corrupt, base, geometry, sweep, fragments, topology,
            definitions, openings);
    }, "opening pressure operator rejects entry corruption");
    auto limits = PlanarPressureRegionFragmentOpeningPressureOperatorLimits{};
    limits.maximumEntries = 129;
    expectRejected([&] {
        static_cast<void>(buildPlanarPressureRegionFragmentOpeningPressureOperator(
            base, geometry, sweep, fragments, topology, definitions,
            openings, limits));
    }, "opening pressure operator enforces entry limits");
    limits = {};
    limits.maximumWorkingBytes = 1;
    expectRejected([&] {
        static_cast<void>(buildPlanarPressureRegionFragmentOpeningPressureOperator(
            base, geometry, sweep, fragments, topology, definitions,
            openings, limits));
    }, "opening pressure operator enforces working storage limits");
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
              && !first.usesOpeningPressureOperator
              && first.pressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && first.basePressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && first.openingFingerprint == 0
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

void testPlanarRegionalFragmentOpeningPressureCorrectionSolve() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto base = buildPlanarPressureRegionFragmentPressureOperator(
        geometry, sweep, fragments, topology);
    const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions{{100, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.5}};
    const auto openings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, definitions);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            base, geometry, sweep, fragments, topology, definitions,
            openings);

    std::vector<double> manufactured;
    manufactured.reserve(fragments.fragments.size());
    for (const auto& fragment : fragments.fragments) {
        manufactured.push_back(fragment.pressurePascals);
    }
    subtractFragmentCorrectionVolumeMeans(
        pressureOperator, fragments, manufactured);
    const auto rightHandSide =
        applyPlanarPressureRegionFragmentOpeningPressureOperator(
            pressureOperator, manufactured);

    std::vector<double> sealedCorrection(base.rows.size(), 0.0);
    const auto sealed = solvePlanarPressureRegionFragmentPressureCorrection(
        base, geometry, sweep, fragments, topology, rightHandSide,
        sealedCorrection);
    check(!sealed.compatible && !sealed.converged && sealed.finite,
          "sealed regional pressure solve rejects cross-wall aperture transfer");

    PlanarPressureRegionFragmentPressureSolveSettings settings;
    settings.absoluteResidualTolerancePascalsMeters = 1.0e-13;
    settings.relativeResidualTolerance = 0.0;
    settings.maximumIterations = 300;
    std::vector<double> firstCorrection(
        pressureOperator.rows.size(), 0.0);
    auto secondCorrection = firstCorrection;
    const auto first =
        solvePlanarPressureRegionFragmentOpeningPressureCorrection(
            pressureOperator, base, geometry, sweep, fragments, topology,
            definitions, openings, rightHandSide, firstCorrection, settings);
    const auto second =
        solvePlanarPressureRegionFragmentOpeningPressureCorrection(
            pressureOperator, base, geometry, sweep, fragments, topology,
            definitions, openings, rightHandSide, secondCorrection, settings);
    check(first == second && firstCorrection == secondCorrection
              && first.compatible && first.converged && first.finite
              && first.usesOpeningPressureOperator
              && first.pressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && first.basePressureOperatorFingerprint == base.fingerprint
              && first.openingFingerprint == openings.fingerprint
              && first.fragmentFingerprint == fragments.fingerprint
              && first.componentCount == 1
              && first.components.size() == 1
              && first.iterationCount > 0
              && first.iterationCount <= settings.maximumIterations,
          "opening regional pressure solve accepts cross-wall transfer deterministically");
    double maximumRecoveryError = 0.0;
    for (std::size_t index = 0; index < manufactured.size(); ++index) {
        maximumRecoveryError = std::max(
            maximumRecoveryError,
            std::abs(firstCorrection[index] - manufactured[index]));
    }
    const auto recoveredRightHandSide =
        applyPlanarPressureRegionFragmentOpeningPressureOperator(
            pressureOperator, firstCorrection);
    double maximumResidual = 0.0;
    for (std::size_t index = 0; index < rightHandSide.size(); ++index) {
        maximumResidual = std::max(
            maximumResidual,
            std::abs(recoveredRightHandSide[index]
                     - rightHandSide[index]));
    }
    check(maximumRecoveryError < 3.0e-10,
          "opening regional pressure solve recovers its aperture-coupled field");
    check(maximumResidual < 2.0e-11,
          "opening regional pressure solve closes its recomputed residual");
    check(first.finalResidualL2PascalsMeters < 1.0e-12,
          "opening regional pressure solve reports a small RMS residual");
    check(first.finalResidualMaximumPascalsMeters < 2.0e-11,
          "opening regional pressure solve reports a small maximum residual");
    checkNear(
        fragmentCorrectionVolumeMean(
            pressureOperator, fragments, pressureOperator.components[0],
            firstCorrection),
        0.0, 2.0e-15,
        "opening regional pressure correction has one volume-weighted gauge");
    checkNear(first.components[0].rightHandSideSumPascalsMeters,
              0.0, 2.0e-13,
              "opening regional pressure RHS is connected-component compatible");

    const auto emptyOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, {});
    const auto emptyOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            base, geometry, sweep, fragments, topology, {}, emptyOpenings);
    std::vector<double> sealedManufactured(base.rows.size(), 0.0);
    for (std::size_t index = 0; index < sealedManufactured.size(); ++index) {
        const double sample = static_cast<double>(index + 1);
        sealedManufactured[index] =
            std::sin(0.23 * sample) + 0.04 * std::cos(0.37 * sample);
    }
    subtractFragmentCorrectionVolumeMeans(
        base, fragments, sealedManufactured);
    const auto sealedRightHandSide =
        applyPlanarPressureRegionFragmentPressureOperator(
            base, sealedManufactured);
    std::vector<double> baseRecovered(base.rows.size(), 0.0);
    std::vector<double> emptyRecovered(base.rows.size(), 0.0);
    const auto baseDiagnostics =
        solvePlanarPressureRegionFragmentPressureCorrection(
            base, geometry, sweep, fragments, topology,
            sealedRightHandSide, baseRecovered, settings);
    const auto emptyDiagnostics =
        solvePlanarPressureRegionFragmentOpeningPressureCorrection(
            emptyOperator, base, geometry, sweep, fragments, topology,
            {}, emptyOpenings, sealedRightHandSide, emptyRecovered, settings);
    double maximumEmptyDifference = 0.0;
    for (std::size_t index = 0; index < baseRecovered.size(); ++index) {
        maximumEmptyDifference = std::max(
            maximumEmptyDifference,
            std::abs(baseRecovered[index] - emptyRecovered[index]));
    }
    check(baseDiagnostics.compatible && baseDiagnostics.converged
              && emptyDiagnostics.compatible && emptyDiagnostics.converged
              && emptyDiagnostics.componentCount
                  == baseDiagnostics.componentCount
              && maximumEmptyDifference < 3.0e-11,
          "empty opening pressure solve preserves the sealed solution");

    auto incompatibleRightHandSide = rightHandSide;
    incompatibleRightHandSide[
        pressureOperator.components[0].gaugeFragmentIndex] += 1.0e-4;
    std::vector<double> warmStart(pressureOperator.rows.size(), 0.125);
    const auto originalWarmStart = warmStart;
    const auto incompatible =
        solvePlanarPressureRegionFragmentOpeningPressureCorrection(
            pressureOperator, base, geometry, sweep, fragments, topology,
            definitions, openings, incompatibleRightHandSide, warmStart,
            settings);
    check(!incompatible.compatible && !incompatible.converged
              && incompatible.finite && warmStart == originalWarmStart,
          "incompatible opening pressure RHS rolls back transactionally");

    auto truncatedSettings = settings;
    truncatedSettings.absoluteResidualTolerancePascalsMeters = 1.0e-16;
    truncatedSettings.relativeResidualTolerance = 0.0;
    truncatedSettings.maximumIterations = 1;
    warmStart = originalWarmStart;
    const auto truncated =
        solvePlanarPressureRegionFragmentOpeningPressureCorrection(
            pressureOperator, base, geometry, sweep, fragments, topology,
            definitions, openings, rightHandSide, warmStart,
            truncatedSettings);
    check(truncated.compatible && !truncated.converged
              && truncated.iterationCount == 1
              && warmStart == originalWarmStart,
          "truncated opening pressure solve preserves its warm start");

    auto corrupt = pressureOperator;
    corrupt.entries.back().geometryWeightMeters += 0.1;
    expectRejected([&] {
        auto correction = originalWarmStart;
        static_cast<void>(
            solvePlanarPressureRegionFragmentOpeningPressureCorrection(
                corrupt, base, geometry, sweep, fragments, topology,
                definitions, openings, rightHandSide, correction, settings));
    }, "opening pressure solve rejects a corrupted augmented operator");
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

    PlanarPressureRegionFragmentPressureJumpEnergySettings jumpSettings;
    jumpSettings.timeStepSeconds = projectionSettings.timeStepSeconds;
    const auto jumpAudit =
        auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, sweep, fragments, topology, metric, after,
            jumpSettings);
    const auto pressureState =
        composeStaticPlanarPressureRegionFragmentPressureState(
            geometry, sweep, fragments, topology, metric, before, after,
            audit, jumpAudit);
    const auto repeatedPressureState =
        composeStaticPlanarPressureRegionFragmentPressureState(
            geometry, sweep, fragments, topology, metric, before, after,
            audit, jumpAudit);
    check(pressureState == repeatedPressureState
              && pressureState.version
                  == planarPressureRegionFragmentPressureStateVersion
              && pressureState.fingerprint != 0
              && pressureState.accepted
              && pressureState.staticGeometry
              && !pressureState.usesMovingVolumeRates
              && pressureState.sourceFragmentFingerprint
                  == fragments.fingerprint
              && pressureState.sourceTopologyFingerprint
                  == topology.fingerprint
              && pressureState.sourceMetricFingerprint
                  == metric.fingerprint
              && pressureState.sourceProjectionEnergyFingerprint
                  == audit.fingerprint
              && pressureState.sourcePressureJumpEnergyFingerprint
                  == jumpAudit.fingerprint
              && pressureState.volumeRateFingerprint == 0
              && pressureState.timeStepSeconds
                  == projectionSettings.timeStepSeconds
              && pressureState.controls.size()
                  == fragments.fragments.size()
              && pressureState.walls.size() == 8
              && pressureState.components.size() == 2
              && pressureState.ownedStorageBytes > 0
              && pressureState.workingStorageBytes
                  == topology.links.size() * sizeof(std::size_t),
          "static regional pressure state composes both accepted pressure ledgers");
    double maximumCorrectionWallJump = 0.0;
    for (std::size_t index = 0;
         index < pressureState.controls.size(); ++index) {
        const auto& control = pressureState.controls[index];
        check(control.fragmentIndex == index
                  && control.fragmentStableId
                      == fragments.fragments[index].stableId
                  && control.authoredPressurePascals
                      == fragments.fragments[index].pressurePascals
                  && control.correctionPressurePascals == pressure[index],
              "composed pressure control retains authored and correction identity");
        checkNear(control.totalPressurePascals,
                  control.authoredPressurePascals
                      + control.correctionPressurePascals,
                  2.0e-16,
                  "composed pressure control closes total pressure");
    }
    for (const auto& wall : pressureState.walls) {
        maximumCorrectionWallJump = std::max(
            maximumCorrectionWallJump,
            std::abs(wall.correctionPressureJumpPascals));
        checkNear(wall.totalPressureJumpPascals,
                  wall.authoredPressureJumpPascals
                      + wall.correctionPressureJumpPascals,
                  2.0e-14,
                  "composed wall pressure jump retains both pressure owners");
        checkNear(
            std::max({
                std::abs(wall.pressureForceSplitResidualNewtons.x),
                std::abs(wall.pressureForceSplitResidualNewtons.y),
                std::abs(wall.pressureForceSplitResidualNewtons.z)}),
            0.0, 3.0e-14,
            "composed wall sheet force closes its pressure split");
        check(wall.materialWallVelocityMetersPerSecond == 0.0
                  && wall.totalPressureWorkToFluidJoules == 0.0
                  && wall.totalPressureWorkToSheetJoules == 0.0,
              "static composed wall retains exact zero pressure work");
    }
    check(maximumCorrectionWallJump > 1.0e-3
              && pressureState.maximumAbsoluteCorrectionGaugePascals
                  < 3.0e-16
              && pressureState.totalPressureWorkToFluidJoules == 0.0
              && pressureState.totalGeometryPressureWorkJoules == 0.0
              && pressureState.wallGeometryWorkResidualJoules == 0.0,
          "static composed pressure keeps transient wall jump and zero work distinct");
    validateStaticPlanarPressureRegionFragmentPressureState(
        pressureState, geometry, sweep, fragments, topology, metric,
        before, after, audit, jumpAudit);

    const auto surfaceLoads =
        capturePlanarPressureRegionFragmentSurfaceLoads(pressureState);
    const auto repeatedSurfaceLoads =
        capturePlanarPressureRegionFragmentSurfaceLoads(pressureState);
    check(surfaceLoads == repeatedSurfaceLoads
              && surfaceLoads.version
                  == planarPressureRegionFragmentSurfaceLoadVersion
              && surfaceLoads.fingerprint != 0
              && surfaceLoads.accepted
              && surfaceLoads.sourcePressureStateFingerprint
                  == pressureState.fingerprint
              && surfaceLoads.sourceTopologyFingerprint
                  == topology.fingerprint
              && surfaceLoads.staticGeometry
              && !surfaceLoads.usesMovingVolumeRates
              && surfaceLoads.timeStepSeconds
                  == projectionSettings.timeStepSeconds
              && surfaceLoads.tiles.size() == 8
              && surfaceLoads.surfaces.size() == 2
              && surfaceLoads.totalAreaSquareMeters == 8.0
              && surfaceLoads.ownedStorageBytes > 0
              && surfaceLoads.workingStorageBytes
                  == 2 * sizeof(std::uint64_t),
          "static regional surface loads capture a deterministic authored handoff");
    check(surfaceLoads.surfaces[0].surfaceStableId == 10
              && surfaceLoads.surfaces[1].surfaceStableId == 20,
          "regional surface loads sort stable authored surface identity");
    for (std::size_t index = 0;
         index < surfaceLoads.surfaces.size(); ++index) {
        const auto& surface = surfaceLoads.surfaces[index];
        const double sign = index == 0 ? -1.0 : 1.0;
        const double coordinate = index == 0 ? -0.8 : -0.2;
        check(surface.tileCount == 4
                  && surface.areaSquareMeters == 4.0
                  && surface.axis == GridFaceAxis::X,
              "regional surface load aggregates every transverse tile");
        checkNear(surface.areaWeightedCentroidMeters.x,
                  coordinate, 2.0e-16,
                  "regional surface load retains the authored plane centroid");
        checkNear(surface.areaWeightedCentroidMeters.y, 0.0, 0.0,
                  "regional surface load closes transverse Y centroid");
        checkNear(surface.areaWeightedCentroidMeters.z, 0.0, 0.0,
                  "regional surface load closes transverse Z centroid");
        checkNear(surface.authoredPressureForceOnSheetNewtons.x,
                  sign * 280.0, 0.0,
                  "regional surface load sums authored sheet force");
        check(surface.totalPressureWorkToSheetJoules == 0.0,
              "static regional surface load does exact zero work");
    }
    for (const auto& tile : surfaceLoads.tiles) {
        checkNear(
            tile.totalPressureTractionOnSheetPascals.x
                * tile.areaSquareMeters,
            tile.totalPressureForceOnSheetNewtons.x, 3.0e-14,
            "regional surface-load tile reconstructs force from traction");
        checkNear(
            tile.totalPressureImpulseOnSheetNewtonSeconds.x,
            tile.totalPressureForceOnSheetNewtons.x
                * surfaceLoads.timeStepSeconds,
            3.0e-16,
            "regional surface-load tile integrates pressure impulse");
        check(tile.totalPressureWorkToSheetJoules == 0.0,
              "static regional surface-load tile retains zero work");
    }
    check(surfaceLoads.authoredPressureForceOnSheetNewtons
                  == pressureState.authoredPressureForceOnSheetNewtons
              && surfaceLoads.correctionPressureForceOnSheetNewtons
                  == pressureState.correctionPressureForceOnSheetNewtons
              && surfaceLoads.totalPressureForceOnSheetNewtons
                  == pressureState.totalPressureForceOnSheetNewtons
              && surfaceLoads.totalPressureWorkToSheetJoules
                  == pressureState.totalPressureWorkToSheetJoules
              && surfaceLoads.sourceWorkResidualJoules == 0.0,
          "regional surface-load handoff closes exactly to its pressure state");
    validatePlanarPressureRegionFragmentSurfaceLoadLedgerIntegrity(
        surfaceLoads);
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressureState);

    const std::vector<
        PlanarPressureRegionFragmentOpeningPatchDefinition> loadOpening{
        {100, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.5},
    };
    const auto loadOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, loadOpening);
    const auto openingAdjustedLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressureState, geometry, sweep, fragments,
            topology, loadOpening, loadOpenings);
    const auto repeatedOpeningAdjustedLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressureState, geometry, sweep, fragments,
            topology, loadOpening, loadOpenings);
    check(openingAdjustedLoads == repeatedOpeningAdjustedLoads
              && openingAdjustedLoads.version
                  == planarPressureRegionFragmentOpeningSurfaceLoadVersion
              && openingAdjustedLoads.fingerprint != 0
              && openingAdjustedLoads.accepted
              && openingAdjustedLoads.sourceSurfaceLoadFingerprint
                  == surfaceLoads.fingerprint
              && openingAdjustedLoads.sourcePressureStateFingerprint
                  == pressureState.fingerprint
              && openingAdjustedLoads.sourceOpeningFingerprint
                  == loadOpenings.fingerprint
              && openingAdjustedLoads.sourceTopologyFingerprint
                  == topology.fingerprint
              && openingAdjustedLoads.staticGeometry
              && !openingAdjustedLoads.usesMovingVolumeRates
              && openingAdjustedLoads.tiles.size()
                  == surfaceLoads.tiles.size()
              && openingAdjustedLoads.surfaces.size()
                  == surfaceLoads.surfaces.size()
              && openingAdjustedLoads.totalWallAreaSquareMeters == 8.0
              && openingAdjustedLoads.totalOpeningAreaSquareMeters == 0.5
              && openingAdjustedLoads.totalSolidAreaSquareMeters == 7.5
              && openingAdjustedLoads.wallAreaPartitionResidualSquareMeters
                  == 0.0
              && openingAdjustedLoads.ownedStorageBytes > 0
              && openingAdjustedLoads.workingStorageBytes
                  == surfaceLoads.tiles.size()
                     * (sizeof(std::pair<std::size_t, std::size_t>)
                        + sizeof(std::size_t)),
          "opening-aware surface load is deterministic and source-bound");
    const auto adjustedTile = std::ranges::find(
        openingAdjustedLoads.tiles,
        loadOpenings.patches[0].sourceFaceLinkStableId,
        &PlanarPressureRegionFragmentOpeningSurfaceLoadTile::
            sourceFaceLinkStableId);
    const auto sourceTile = std::ranges::find(
        surfaceLoads.tiles,
        loadOpenings.patches[0].sourceFaceLinkStableId,
        &PlanarPressureRegionFragmentSurfaceLoadTile::
            sourceFaceLinkStableId);
    check(adjustedTile != openingAdjustedLoads.tiles.end()
              && sourceTile != surfaceLoads.tiles.end(),
          "opening-aware surface load retains the exact source wall tile");
    if (adjustedTile != openingAdjustedLoads.tiles.end()
        && sourceTile != surfaceLoads.tiles.end()) {
        check(adjustedTile->touchedByOpening
                  && adjustedTile->sourceOpeningPartitionIndex == 0
                  && adjustedTile->openingPatchCount == 1
                  && adjustedTile->wallAreaSquareMeters == 1.0
                  && adjustedTile->openingAreaSquareMeters == 0.5
                  && adjustedTile->solidAreaSquareMeters == 0.5
                  && adjustedTile->openingAreaFraction == 0.5,
              "opening-aware surface load partitions the touched tile exactly");
        check(adjustedTile->authoredPressureTractionOnSheetPascals
                      == sourceTile
                             ->authoredPressureTractionOnSheetPascals
                  && adjustedTile
                         ->correctionPressureTractionOnSheetPascals
                      == sourceTile
                             ->correctionPressureTractionOnSheetPascals
                  && adjustedTile->totalPressureTractionOnSheetPascals
                      == sourceTile->totalPressureTractionOnSheetPascals,
              "opening-aware load preserves pressure-jump traction without an absolute-pressure gauge path");
        for (const auto& component : {
                 std::pair{
                     adjustedTile
                         ->openingRemovedTotalPressureForceOnSheetNewtons.x,
                     sourceTile->totalPressureForceOnSheetNewtons.x},
                 std::pair{
                     adjustedTile->solidTotalPressureForceOnSheetNewtons.x,
                     sourceTile->totalPressureForceOnSheetNewtons.x},
                 std::pair{
                     adjustedTile
                         ->openingRemovedTotalPressureImpulseOnSheetNewtonSeconds.x,
                     sourceTile
                         ->totalPressureImpulseOnSheetNewtonSeconds.x},
                 std::pair{
                     adjustedTile
                         ->solidTotalPressureImpulseOnSheetNewtonSeconds.x,
                     sourceTile
                         ->totalPressureImpulseOnSheetNewtonSeconds.x}}) {
            checkNear(component.first, 0.5 * component.second, 3.0e-14,
                      "half-open tile halves retained and removed load ownership");
        }
        const auto cross = [](const Vector3& first, const Vector3& second) {
            return Vector3{
                first.y * second.z - first.z * second.y,
                first.z * second.x - first.x * second.z,
                first.x * second.y - first.y * second.x,
            };
        };
        const Vector3 sourceMoment = cross(
            sourceTile->wrappedCentroidMeters,
            sourceTile->totalPressureForceOnSheetNewtons);
        const Vector3 partitionedMoment{
            adjustedTile
                    ->openingRemovedTotalPressureMomentOnSheetNewtonMeters.x
                + adjustedTile
                      ->solidTotalPressureMomentOnSheetNewtonMeters.x,
            adjustedTile
                    ->openingRemovedTotalPressureMomentOnSheetNewtonMeters.y
                + adjustedTile
                      ->solidTotalPressureMomentOnSheetNewtonMeters.y,
            adjustedTile
                    ->openingRemovedTotalPressureMomentOnSheetNewtonMeters.z
                + adjustedTile
                      ->solidTotalPressureMomentOnSheetNewtonMeters.z,
        };
        checkNear(partitionedMoment.x, sourceMoment.x, 3.0e-14,
                  "opening-aware tile closes pressure moment X");
        checkNear(partitionedMoment.y, sourceMoment.y, 3.0e-14,
                  "opening-aware tile closes pressure moment Y");
        checkNear(partitionedMoment.z, sourceMoment.z, 3.0e-14,
                  "opening-aware tile closes pressure moment Z");
    }
    check(openingAdjustedLoads.surfaces[0].surfaceStableId == 10
              && openingAdjustedLoads.surfaces[0].tileCount == 4
              && openingAdjustedLoads.surfaces[0].openingTouchedTileCount
                  == 1
              && openingAdjustedLoads.surfaces[0].fullyOpenTileCount == 0
              && openingAdjustedLoads.surfaces[0].wallAreaSquareMeters
                  == 4.0
              && openingAdjustedLoads.surfaces[0].openingAreaSquareMeters
                  == 0.5
              && openingAdjustedLoads.surfaces[0].solidAreaSquareMeters
                  == 3.5
              && openingAdjustedLoads.surfaces[1].openingAreaSquareMeters
                  == 0.0
              && openingAdjustedLoads.surfaces[1].solidAreaSquareMeters
                  == 4.0,
          "opening-aware surface summaries retain touched and untouched area");
    check(openingAdjustedLoads.maximumAbsoluteAreaPartitionResidualSquareMeters
                  < 1.0e-15
              && openingAdjustedLoads
                         .maximumAbsoluteForcePartitionResidualNewtons
                  < 1.0e-12
              && openingAdjustedLoads
                         .maximumAbsoluteForceSplitResidualNewtons
                  < 1.0e-12
              && openingAdjustedLoads
                         .maximumAbsoluteImpulsePartitionResidualNewtonSeconds
                  < 1.0e-12
              && openingAdjustedLoads
                         .maximumAbsoluteMomentPartitionResidualNewtonMeters
                  < 1.0e-12
              && openingAdjustedLoads
                         .maximumAbsoluteSurfaceAggregationResidualNewtons
                  < 1.0e-12
              && openingAdjustedLoads.workPartitionResidualJoules == 0.0,
          "opening-aware load closes area, force, impulse, moment, and work");
    validatePlanarPressureRegionFragmentOpeningSurfaceLoadLedgerIntegrity(
        openingAdjustedLoads);
    validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
        openingAdjustedLoads, surfaceLoads, pressureState, geometry, sweep,
        fragments, topology, loadOpening, loadOpenings);

    auto centroidLoadOpening = loadOpening;
    centroidLoadOpening[0].authoredWrappedCentroidMeters =
        Vector3{-0.8, -0.75, -0.5};
    const auto centroidLoadOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, centroidLoadOpening);
    const auto centroidAdjustedLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressureState, geometry, sweep, fragments,
            topology, centroidLoadOpening, centroidLoadOpenings);
    const auto centroidAdjustedTile = std::ranges::find(
        centroidAdjustedLoads.tiles,
        centroidLoadOpenings.patches[0].sourceFaceLinkStableId,
        &PlanarPressureRegionFragmentOpeningSurfaceLoadTile::
            sourceFaceLinkStableId);
    check(centroidAdjustedTile != centroidAdjustedLoads.tiles.end()
              && centroidAdjustedTile->hasExactSubtileCentroids,
          "opening-aware surface load retains exact sub-tile centroid provenance");
    if (centroidAdjustedTile != centroidAdjustedLoads.tiles.end()) {
        checkNear(
            centroidAdjustedTile->openingAreaWeightedCentroidMeters.y,
            -0.75, 0.0,
            "opening-aware surface load retains the aperture centroid");
        checkNear(
            centroidAdjustedTile->solidAreaWeightedCentroidMeters.y,
            -0.25, 0.0,
            "opening-aware surface load derives the retained centroid");
        const auto centroidCross = [](const Vector3& first,
                                      const Vector3& second) {
            return Vector3{
                first.y * second.z - first.z * second.y,
                first.z * second.x - first.x * second.z,
                first.x * second.y - first.y * second.x,
            };
        };
        const Vector3 expectedOpeningMoment = centroidCross(
            centroidAdjustedTile->openingAreaWeightedCentroidMeters,
            centroidAdjustedTile
                ->openingRemovedTotalPressureForceOnSheetNewtons);
        const Vector3 expectedSolidMoment = centroidCross(
            centroidAdjustedTile->solidAreaWeightedCentroidMeters,
            centroidAdjustedTile->solidTotalPressureForceOnSheetNewtons);
        check(centroidAdjustedTile
                      ->openingRemovedTotalPressureMomentOnSheetNewtonMeters
                  == expectedOpeningMoment
                  && centroidAdjustedTile
                         ->solidTotalPressureMomentOnSheetNewtonMeters
                      == expectedSolidMoment,
              "opening-aware surface load uses distinct exact moment arms");
    }
    checkNear(
        centroidAdjustedLoads.maximumAbsoluteMomentPartitionResidualNewtonMeters,
        0.0, 3.0e-14,
        "opening-aware exact centroid moments close to the source wall");
    validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
        centroidAdjustedLoads, surfaceLoads, pressureState, geometry, sweep,
        fragments, topology, centroidLoadOpening, centroidLoadOpenings);

    const std::vector<
        PlanarPressureRegionFragmentOpeningPatchDefinition> fullLoadOpening{
        {100, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 1.0},
    };
    const auto fullLoadOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, fullLoadOpening);
    const auto fullyOpenLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressureState, geometry, sweep, fragments,
            topology, fullLoadOpening, fullLoadOpenings);
    const auto fullyOpenTile = std::ranges::find(
        fullyOpenLoads.tiles,
        fullLoadOpenings.patches[0].sourceFaceLinkStableId,
        &PlanarPressureRegionFragmentOpeningSurfaceLoadTile::
            sourceFaceLinkStableId);
    check(fullyOpenTile != fullyOpenLoads.tiles.end()
              && fullyOpenTile->solidAreaSquareMeters == 0.0
              && fullyOpenTile->solidAuthoredPressureForceOnSheetNewtons
                  == Vector3{}
              && fullyOpenTile->solidCorrectionPressureForceOnSheetNewtons
                  == Vector3{}
              && fullyOpenTile->solidTotalPressureForceOnSheetNewtons
                  == Vector3{}
              && fullyOpenTile
                         ->solidTotalPressureImpulseOnSheetNewtonSeconds
                  == Vector3{}
              && fullyOpenTile
                         ->solidTotalPressureMomentOnSheetNewtonMeters
                  == Vector3{}
              && fullyOpenTile->solidTotalPressureWorkToSheetJoules == 0.0
              && fullyOpenLoads.surfaces[0].fullyOpenTileCount == 1,
          "a fully open tile carries exactly zero solid-wall pressure load");

    std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        sharedLoadOpening{
            {111, 2000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.2},
            {110, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.3},
        };
    const auto sharedLoadOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, sharedLoadOpening);
    const auto sharedOpeningLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressureState, geometry, sweep, fragments,
            topology, sharedLoadOpening, sharedLoadOpenings);
    auto reversedSharedLoadOpening = sharedLoadOpening;
    std::ranges::reverse(reversedSharedLoadOpening);
    const auto reversedSharedLoadOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology,
            reversedSharedLoadOpening);
    const auto reversedSharedOpeningLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressureState, geometry, sweep, fragments,
            topology, reversedSharedLoadOpening,
            reversedSharedLoadOpenings);
    check(sharedLoadOpenings == reversedSharedLoadOpenings
              && sharedOpeningLoads == reversedSharedOpeningLoads
              && sharedOpeningLoads.tiles[0].openingPatchCount == 2
              && sharedOpeningLoads.tiles[0].openingAreaSquareMeters == 0.5,
          "opening-aware load canonicalizes multiple patches on one tile");

    const auto emptyLoadOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, {});
    const auto emptyOpeningLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            surfaceLoads, pressureState, geometry, sweep, fragments,
            topology, {}, emptyLoadOpenings);
    bool emptyPreservesSource =
        emptyOpeningLoads.totalOpeningAreaSquareMeters == 0.0
        && emptyOpeningLoads.totalSolidAreaSquareMeters
            == surfaceLoads.totalAreaSquareMeters
        && emptyOpeningLoads.solidAuthoredPressureForceOnSheetNewtons
            == surfaceLoads.authoredPressureForceOnSheetNewtons
        && emptyOpeningLoads.solidCorrectionPressureForceOnSheetNewtons
            == surfaceLoads.correctionPressureForceOnSheetNewtons
        && emptyOpeningLoads.solidTotalPressureForceOnSheetNewtons
            == surfaceLoads.totalPressureForceOnSheetNewtons
        && emptyOpeningLoads.solidTotalPressureImpulseOnSheetNewtonSeconds
            == surfaceLoads.totalPressureImpulseOnSheetNewtonSeconds
        && emptyOpeningLoads.solidTotalPressureWorkToSheetJoules
            == surfaceLoads.totalPressureWorkToSheetJoules;
    for (std::size_t index = 0;
         index < surfaceLoads.tiles.size(); ++index) {
        emptyPreservesSource = emptyPreservesSource
            && !emptyOpeningLoads.tiles[index].touchedByOpening
            && emptyOpeningLoads.tiles[index].openingAreaSquareMeters == 0.0
            && emptyOpeningLoads.tiles[index].solidAreaSquareMeters
                == surfaceLoads.tiles[index].areaSquareMeters
            && emptyOpeningLoads.tiles[index]
                   .solidTotalPressureForceOnSheetNewtons
                == surfaceLoads.tiles[index]
                       .totalPressureForceOnSheetNewtons;
    }
    check(emptyPreservesSource,
          "an empty opening overlay preserves every sealed load bit-exactly");

    auto corruptOpeningAdjustedLoads = openingAdjustedLoads;
    corruptOpeningAdjustedLoads.tiles[0]
        .solidTotalPressureForceOnSheetNewtons.x += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningSurfaceLoadLedgerIntegrity(
                corruptOpeningAdjustedLoads);
        },
        "opening-aware surface-load integrity rejects tile corruption");
    auto corruptLoadOpenings = loadOpenings;
    corruptLoadOpenings.partitions[0].solidAreaSquareMeters += 0.1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                    surfaceLoads, pressureState, geometry, sweep, fragments,
                    topology, loadOpening, corruptLoadOpenings));
        },
        "opening-aware surface load rejects opening partition corruption");
    auto openingLoadLimits =
        PlanarPressureRegionFragmentOpeningSurfaceLoadLimits{};
    openingLoadLimits.maximumTiles = surfaceLoads.tiles.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                    surfaceLoads, pressureState, geometry, sweep, fragments,
                    topology, loadOpening, loadOpenings,
                    openingLoadLimits));
        },
        "opening-aware surface load enforces its tile limit");
    openingLoadLimits = {};
    openingLoadLimits.maximumOwnedBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                    surfaceLoads, pressureState, geometry, sweep, fragments,
                    topology, loadOpening, loadOpenings,
                    openingLoadLimits));
        },
        "opening-aware surface load enforces its owned-byte limit");
    openingLoadLimits = {};
    openingLoadLimits.maximumWorkingBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                    surfaceLoads, pressureState, geometry, sweep, fragments,
                    topology, loadOpening, loadOpenings,
                    openingLoadLimits));
        },
        "opening-aware surface load enforces its working-byte limit");
    openingLoadLimits = {};
    openingLoadLimits.openingLimits.maximumPatches = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                    surfaceLoads, pressureState, geometry, sweep, fragments,
                    topology, sharedLoadOpening, sharedLoadOpenings,
                    openingLoadLimits));
        },
        "opening-aware surface load enforces nested opening limits");

    const auto acceptedState =
        capturePlanarPressureRegionFragmentAcceptedState(
            geometry, sweep, fragments, topology, metric, after,
            pressureState, surfaceLoads);
    const auto repeatedAcceptedState =
        capturePlanarPressureRegionFragmentAcceptedState(
            geometry, sweep, fragments, topology, metric, after,
            pressureState, surfaceLoads);
    check(acceptedState == repeatedAcceptedState
              && acceptedState.version
                  == planarPressureRegionFragmentAcceptedStateVersion
              && acceptedState.fingerprint != 0
              && acceptedState.accepted
              && acceptedState.staticGeometry
              && !acceptedState.usesMovingVolumeRates
              && acceptedState.sourceVelocityStateFingerprint
                  == after.fingerprint
              && acceptedState.sourcePressureStateFingerprint
                  == pressureState.fingerprint
              && acceptedState.sourceSurfaceLoadFingerprint
                  == surfaceLoads.fingerprint
              && acceptedState.velocity == after
              && acceptedState.pressure == pressureState
              && acceptedState.surfaceLoads == surfaceLoads
              && acceptedState.fluidMomentumKilogramMetersPerSecond
                  == after.momentumKilogramMetersPerSecond
              && acceptedState.fluidKineticEnergyJoules
                  == after.kineticEnergyJoules
              && acceptedState.pressureForceOnSheetNewtons
                  == surfaceLoads.totalPressureForceOnSheetNewtons
              && acceptedState.pressureImpulseOnSheetNewtonSeconds
                  == surfaceLoads.totalPressureImpulseOnSheetNewtonSeconds
              && acceptedState.pressureWorkToSheetJoules == 0.0
              && acceptedState.ownedStorageBytes
                  == after.ownedStorageBytes
                     + pressureState.ownedStorageBytes
                     + surfaceLoads.ownedStorageBytes,
          "static regional accepted state atomically retains all continuation products");
    validatePlanarPressureRegionFragmentAcceptedStateIntegrity(
        acceptedState);
    validatePlanarPressureRegionFragmentAcceptedState(
        acceptedState, geometry, sweep, fragments, topology, metric);
    auto corruptAcceptedState = acceptedState;
    corruptAcceptedState.velocity.samples[0]
        .normalVelocityMetersPerSecond += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentAcceptedStateIntegrity(
                corruptAcceptedState);
        },
        "regional accepted state rejects nested velocity corruption");
    corruptAcceptedState = acceptedState;
    corruptAcceptedState.pressure.controls[0].totalPressurePascals += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentAcceptedStateIntegrity(
                corruptAcceptedState);
        },
        "regional accepted state rejects nested pressure corruption");
    corruptAcceptedState = acceptedState;
    corruptAcceptedState.surfaceLoads.tiles[0]
        .totalPressureWorkToSheetJoules += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentAcceptedStateIntegrity(
                corruptAcceptedState);
        },
        "regional accepted state rejects nested load corruption");
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentAcceptedState(
                    geometry, sweep, fragments, topology, metric, before,
                    pressureState, surfaceLoads));
        },
        "regional accepted state rejects a foreign velocity endpoint");
    auto acceptedLimits =
        PlanarPressureRegionFragmentAcceptedStateLimits{};
    acceptedLimits.maximumOwnedBytes = acceptedState.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentAcceptedState(
                    geometry, sweep, fragments, topology, metric, after,
                    pressureState, surfaceLoads, acceptedLimits));
        },
        "regional accepted state enforces its aggregate byte limit");
    acceptedLimits = {};
    acceptedLimits.velocityStateLimits.maximumSamples = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentAcceptedState(
                    geometry, sweep, fragments, topology, metric, after,
                    pressureState, surfaceLoads, acceptedLimits));
        },
        "regional accepted state enforces nested velocity limits");
    acceptedLimits = {};
    acceptedLimits.surfaceLoadLimits.maximumTiles = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentAcceptedState(
                    geometry, sweep, fragments, topology, metric, after,
                    pressureState, surfaceLoads, acceptedLimits));
        },
        "regional accepted state enforces nested surface-load limits");
    auto corruptSurfaceLoads = surfaceLoads;
    corruptSurfaceLoads.tiles[0].totalPressureForceOnSheetNewtons.x += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentSurfaceLoadLedgerIntegrity(
                corruptSurfaceLoads);
        },
        "regional surface-load ledger rejects tile corruption");
    auto surfaceLoadLimits =
        PlanarPressureRegionFragmentSurfaceLoadLimits{};
    surfaceLoadLimits.maximumTiles = surfaceLoads.tiles.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentSurfaceLoads(
                    pressureState, surfaceLoadLimits));
        },
        "regional surface-load ledger enforces the tile limit");
    surfaceLoadLimits = {};
    surfaceLoadLimits.maximumSurfaces = surfaceLoads.surfaces.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentSurfaceLoads(
                    pressureState, surfaceLoadLimits));
        },
        "regional surface-load ledger enforces the surface limit");
    surfaceLoadLimits = {};
    surfaceLoadLimits.maximumOwnedBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentSurfaceLoads(
                    pressureState, surfaceLoadLimits));
        },
        "regional surface-load ledger enforces the owned byte limit");
    surfaceLoadLimits = {};
    surfaceLoadLimits.maximumWorkingBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentSurfaceLoads(
                    pressureState, surfaceLoadLimits));
        },
        "regional surface-load ledger enforces the working byte limit");

    auto corruptPressureState = pressureState;
    corruptPressureState.fingerprint = 0;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentPressureState(
                corruptPressureState, geometry, sweep, fragments, topology,
                metric, before, after, audit, jumpAudit);
        },
        "composed pressure state rejects fingerprint corruption");
    corruptPressureState = pressureState;
    corruptPressureState.walls[0].totalPressureJumpPascals += 0.1;
    expectRejected(
        [&] {
            validateStaticPlanarPressureRegionFragmentPressureState(
                corruptPressureState, geometry, sweep, fragments, topology,
                metric, before, after, audit, jumpAudit);
        },
        "composed pressure state rejects wall corruption");
    auto mismatchedJumpSettings = jumpSettings;
    mismatchedJumpSettings.timeStepSeconds *= 2.0;
    const auto mismatchedJumpAudit =
        auditStaticPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, sweep, fragments, topology, metric, after,
            mismatchedJumpSettings);
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, mismatchedJumpAudit));
        },
        "composed pressure state rejects mismatched audit durations");
    auto pressureStateLimits =
        PlanarPressureRegionFragmentPressureStateLimits{};
    pressureStateLimits.maximumControls =
        pressureState.controls.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, jumpAudit, pressureStateLimits));
        },
        "composed pressure state enforces the control limit");
    pressureStateLimits = {};
    pressureStateLimits.maximumWalls = pressureState.walls.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, jumpAudit, pressureStateLimits));
        },
        "composed pressure state enforces the wall limit");
    pressureStateLimits = {};
    pressureStateLimits.maximumComponents =
        pressureState.components.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, jumpAudit, pressureStateLimits));
        },
        "composed pressure state enforces the component limit");
    pressureStateLimits = {};
    pressureStateLimits.maximumOwnedBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, jumpAudit, pressureStateLimits));
        },
        "composed pressure state enforces the owned byte limit");
    pressureStateLimits = {};
    pressureStateLimits.maximumWorkingBytes = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, jumpAudit, pressureStateLimits));
        },
        "composed pressure state enforces the working byte limit");
    pressureStateLimits = {};
    pressureStateLimits.projectionEnergyLimits.maximumCorrections = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, jumpAudit, pressureStateLimits));
        },
        "composed pressure state enforces nested projection-audit limits");
    pressureStateLimits = {};
    pressureStateLimits.pressureJumpLimits.maximumLayers = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeStaticPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, metric, before,
                    after, audit, jumpAudit, pressureStateLimits));
        },
        "composed pressure state enforces nested jump-audit limits");

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
    PlanarPressureRegionFragmentPressureJumpEnergySettings jumpSettings;
    jumpSettings.timeStepSeconds = settings.timeStepSeconds;
    const auto jumpEnergy =
        auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
            geometry, sweep, fragments, topology, volumeRates, metric,
            afterState, jumpSettings);
    const auto pressureState =
        composeMovingPlanarPressureRegionFragmentPressureState(
            geometry, sweep, fragments, topology, volumeRates, metric,
            beforeState, afterState, energy, jumpEnergy);
    const auto repeatedPressureState =
        composeMovingPlanarPressureRegionFragmentPressureState(
            geometry, sweep, fragments, topology, volumeRates, metric,
            beforeState, afterState, energy, jumpEnergy);
    check(pressureState == repeatedPressureState
              && pressureState.accepted
              && !pressureState.staticGeometry
              && pressureState.usesMovingVolumeRates
              && pressureState.volumeRateFingerprint
                  == volumeRates.fingerprint
              && pressureState.sourceProjectionEnergyFingerprint
                  == energy.fingerprint
              && pressureState.sourcePressureJumpEnergyFingerprint
                  == jumpEnergy.fingerprint
              && pressureState.maximumAbsoluteCorrectionGaugePascals
                  < 3.0e-13,
          "moving composed pressure state is deterministic and source-bound");
    checkNear(pressureState.authoredPressureWorkToFluidJoules,
              0.0, 6.0e-14,
              "rigid composed pressure retains zero net authored work");
    checkNear(pressureState.correctionPressureWorkToFluidJoules,
              energy.geometryPressureWorkJoules, 4.0e-13,
              "rigid composed pressure exposes correction wall work");
    checkNear(pressureState.totalPressureWorkToFluidJoules,
              energy.geometryPressureWorkJoules, 4.0e-13,
              "rigid composed pressure sums authored and correction work");
    checkNear(pressureState.totalGeometryPressureWorkJoules,
              energy.geometryPressureWorkJoules, 4.0e-13,
              "rigid composed pressure retains total geometry work");
    checkNear(pressureState.wallGeometryWorkResidualJoules,
              0.0, 4.0e-13,
              "rigid composed pressure closes total wall and geometry work");
    checkNear(pressureState.totalPressureWorkToSheetJoules,
              -pressureState.totalPressureWorkToFluidJoules, 0.0,
              "rigid composed pressure publishes opposite sheet work");
    for (const auto& component : pressureState.components) {
        checkNear(component.wallGeometryWorkResidualJoules,
                  0.0, 4.0e-13,
                  "moving composed pressure closes component wall work");
    }
    validateMovingPlanarPressureRegionFragmentPressureState(
        pressureState, geometry, sweep, fragments, topology, volumeRates,
        metric, beforeState, afterState, energy, jumpEnergy);
    const auto surfaceLoads =
        capturePlanarPressureRegionFragmentSurfaceLoads(pressureState);
    check(surfaceLoads.accepted && !surfaceLoads.staticGeometry
              && surfaceLoads.usesMovingVolumeRates
              && surfaceLoads.sourcePressureStateFingerprint
                  == pressureState.fingerprint
              && surfaceLoads.surfaces.size() == 2
              && surfaceLoads.tiles.size() == 8
              && surfaceLoads.totalAreaSquareMeters == 8.0,
          "moving regional surface loads retain pressure-state provenance");
    checkNear(surfaceLoads.totalPressureWorkToSheetJoules,
              pressureState.totalPressureWorkToSheetJoules, 0.0,
              "moving regional surface loads retain exact sheet work");
    checkNear(
        surfaceLoads.totalPressureImpulseOnSheetNewtonSeconds.x,
        surfaceLoads.totalPressureForceOnSheetNewtons.x
            * surfaceLoads.timeStepSeconds,
        3.0e-14,
        "moving regional surface loads integrate total sheet impulse");
    checkNear(surfaceLoads.sourceWorkResidualJoules, 0.0, 0.0,
              "moving regional surface loads close source work exactly");
    validatePlanarPressureRegionFragmentSurfaceLoads(
        surfaceLoads, pressureState);
    const auto acceptedState =
        capturePlanarPressureRegionFragmentAcceptedState(
            geometry, sweep, fragments, topology, metric, afterState,
            pressureState, surfaceLoads);
    check(acceptedState.accepted && !acceptedState.staticGeometry
              && acceptedState.usesMovingVolumeRates
              && acceptedState.sourceVelocityStateFingerprint
                  == afterState.fingerprint
              && acceptedState.pressureWorkToSheetJoules
                  == surfaceLoads.totalPressureWorkToSheetJoules
              && acceptedState.pressureImpulseOnSheetNewtonSeconds
                  == surfaceLoads.totalPressureImpulseOnSheetNewtonSeconds,
          "moving regional accepted state retains the projected endpoint and sheet reaction");
    validatePlanarPressureRegionFragmentAcceptedState(
        acceptedState, geometry, sweep, fragments, topology, metric);
    auto corruptPressureState = pressureState;
    corruptPressureState.controls[0].totalPressurePascals += 0.1;
    expectRejected(
        [&] {
            validateMovingPlanarPressureRegionFragmentPressureState(
                corruptPressureState, geometry, sweep, fragments, topology,
                volumeRates, metric, beforeState, afterState, energy,
                jumpEnergy);
        },
        "moving composed pressure state rejects control corruption");
    auto mismatchedJumpSettings = jumpSettings;
    mismatchedJumpSettings.timeStepSeconds = 0.5;
    expectRejected(
        [&] {
            const auto mismatchedJump =
                auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                    geometry, sweep, fragments, topology, volumeRates,
                    metric, afterState, mismatchedJumpSettings);
            static_cast<void>(
                composeMovingPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, volumeRates,
                    metric, beforeState, afterState, energy,
                    mismatchedJump));
        },
        "moving composed pressure state rejects mismatched audit duration");
    auto pressureStateLimits =
        PlanarPressureRegionFragmentPressureStateLimits{};
    pressureStateLimits.projectionEnergyLimits.volumeRateLimits
        .maximumFragments = 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composeMovingPlanarPressureRegionFragmentPressureState(
                    geometry, sweep, fragments, topology, volumeRates,
                    metric, beforeState, afterState, energy, jumpEnergy,
                    pressureStateLimits));
        },
        "moving composed pressure state enforces nested volume-rate limits");
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
    PlanarPressureRegionFragmentPressureJumpEnergySettings jumpSettings;
    jumpSettings.timeStepSeconds = settings.timeStepSeconds;
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
        const auto jumpEnergy =
            auditMovingPlanarPressureRegionFragmentPressureJumpEnergy(
                geometry, sweep, fragments, topology, volumeRates, metric,
                after, jumpSettings);
        const auto pressureState =
            composeMovingPlanarPressureRegionFragmentPressureState(
                geometry, sweep, fragments, topology, volumeRates, metric,
                before, after, energy, jumpEnergy);
        const auto surfaceLoads =
            capturePlanarPressureRegionFragmentSurfaceLoads(pressureState);
        const auto wall = std::ranges::find_if(
            topology.links,
            [](const auto& link) {
                return link.kind
                        == PlanarPressureRegionFragmentFaceKind::
                            PressureLayerWall
                    && link.surfaceStableId == 10;
            });
        check(wall != topology.links.end(),
              "moving opening-aware load finds its wall on every axis");
        std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
            openingDefinitions;
        if (wall != topology.links.end()) {
            openingDefinitions.push_back({
                100,
                1000,
                wall->surfaceStableId,
                wall->axis,
                wall->i,
                wall->j,
                wall->k,
                wall->minusRegionStableId,
                wall->plusRegionStableId,
                0.5 * wall->areaSquareMeters,
            });
        }
        const auto openings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, sweep, fragments, topology,
                openingDefinitions);
        const auto openingAdjustedLoads =
            capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                surfaceLoads, pressureState, geometry, sweep, fragments,
                topology, openingDefinitions, openings);
        const auto acceptedState =
            capturePlanarPressureRegionFragmentAcceptedState(
                geometry, sweep, fragments, topology, metric, after,
                pressureState, surfaceLoads);
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
                      < 4.0e-13
                  && pressureState.accepted
                  && pressureState.usesMovingVolumeRates
                  && std::abs(pressureState.wallGeometryWorkResidualJoules)
                      < 5.0e-13
                  && surfaceLoads.accepted
                  && surfaceLoads.surfaces.size() == 2
                  && surfaceLoads.sourceWorkResidualJoules == 0.0
                  && openingAdjustedLoads.accepted
                  && openingAdjustedLoads.usesMovingVolumeRates
                  && wall != topology.links.end()
                  && openingAdjustedLoads.totalOpeningAreaSquareMeters
                      == 0.5 * wall->areaSquareMeters
                  && openingAdjustedLoads.wallAreaPartitionResidualSquareMeters
                      == 0.0
                  && acceptedState.accepted
                  && acceptedState.usesMovingVolumeRates,
              "moving regional projection and accepted endpoint close on every axis");
        if (wall != topology.links.end()) {
            const auto sourceTile = std::ranges::find(
                surfaceLoads.tiles, wall->stableId,
                &PlanarPressureRegionFragmentSurfaceLoadTile::
                    sourceFaceLinkStableId);
            const auto adjustedTile = std::ranges::find(
                openingAdjustedLoads.tiles, wall->stableId,
                &PlanarPressureRegionFragmentOpeningSurfaceLoadTile::
                    sourceFaceLinkStableId);
            check(sourceTile != surfaceLoads.tiles.end()
                      && adjustedTile != openingAdjustedLoads.tiles.end(),
                  "moving opening-aware load retains source identity on every axis");
            if (sourceTile != surfaceLoads.tiles.end()
                && adjustedTile != openingAdjustedLoads.tiles.end()) {
                checkNear(
                    adjustedTile->solidTotalPressureForceOnSheetNewtons.x,
                    0.5 * sourceTile->totalPressureForceOnSheetNewtons.x,
                    5.0e-13,
                    "moving opening-aware solid force closes X");
                checkNear(
                    adjustedTile->solidTotalPressureForceOnSheetNewtons.y,
                    0.5 * sourceTile->totalPressureForceOnSheetNewtons.y,
                    5.0e-13,
                    "moving opening-aware solid force closes Y");
                checkNear(
                    adjustedTile->solidTotalPressureForceOnSheetNewtons.z,
                    0.5 * sourceTile->totalPressureForceOnSheetNewtons.z,
                    5.0e-13,
                    "moving opening-aware solid force closes Z");
                checkNear(
                    adjustedTile->solidTotalPressureWorkToSheetJoules,
                    0.5 * sourceTile->totalPressureWorkToSheetJoules,
                    5.0e-13,
                    "moving opening-aware solid work follows retained area");
                checkNear(
                    adjustedTile
                        ->openingRemovedTotalPressureWorkToSheetJoules,
                    0.5 * sourceTile->totalPressureWorkToSheetJoules,
                    5.0e-13,
                    "moving opening-aware removed work follows aperture area");
                const double axialSourceForce =
                    axis == GridFaceAxis::X
                    ? sourceTile->totalPressureForceOnSheetNewtons.x
                    : axis == GridFaceAxis::Y
                    ? sourceTile->totalPressureForceOnSheetNewtons.y
                    : sourceTile->totalPressureForceOnSheetNewtons.z;
                check(std::abs(axialSourceForce) > 1.0
                          && sourceTile->totalPressureWorkToSheetJoules
                              != 0.0,
                      "moving opening-aware load exercises oriented force and nonzero work");
            }
        }
    }
}

void testPlanarRegionalFragmentOpeningTopology() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);

    const auto empty = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, {});
    check(empty.patches.empty() && empty.partitions.empty()
              && empty.openings.empty()
              && empty.baseComponents.size() == 2
              && empty.connectedComponents.size() == 2
              && empty.connectedComponents[0].baseComponentCount == 1
              && empty.connectedComponents[1].baseComponentCount == 1
              && empty.totalOpeningAreaSquareMeters == 0.0
              && empty.totalTouchedWallAreaSquareMeters == 0.0
              && empty.totalSolidAreaOnTouchedWallsSquareMeters == 0.0
              && empty.wallAreaPartitionResidualSquareMeters == 0.0,
          "regional fragment openings preserve the sealed base components when empty");

    const std::vector<
        PlanarPressureRegionFragmentOpeningPatchDefinition> single{
        {100, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.5},
    };
    const auto opened = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, single);
    const auto repeated = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, single);
    check(opened == repeated
              && opened.version
                  == planarPressureRegionFragmentOpeningVersion
              && opened.fingerprint != 0
              && opened.sourceFragmentFingerprint == fragments.fingerprint
              && opened.sourceTopologyFingerprint == topology.fingerprint
              && opened.profileAxis == GridFaceAxis::X
              && opened.patches.size() == 1
              && opened.partitions.size() == 1
              && opened.openings.size() == 1
              && opened.baseComponents.size() == 2
              && opened.connectedComponents.size() == 1
              && opened.ownedStorageBytes > 0
              && opened.workingStorageBytes > 0,
          "regional fragment opening capture is deterministic and source-bound");
    const auto& patch = opened.patches[0];
    check(patch.patchStableId == 100
              && patch.openingStableId == 1000
              && patch.surfaceStableId == 10
              && patch.axis == GridFaceAxis::X
              && patch.i == 1 && patch.j == 0 && patch.k == 0
              && patch.negativeSideRegionStableId == 1
              && patch.positiveSideRegionStableId == 2
              && patch.areaSquareMeters == 0.5
              && patch.sourceWallAreaSquareMeters == 1.0
              && patch.sourceWallAreaFraction == 0.5
              && patch.minusBaseComponentIndex
                  != patch.plusBaseComponentIndex,
          "regional fragment opening retains its exact wall tile and orientation");
    const auto& partition = opened.partitions[0];
    check(partition.sourceFaceLinkStableId
                  == patch.sourceFaceLinkStableId
              && partition.openingPatchCount == 1
              && partition.wallAreaSquareMeters == 1.0
              && partition.openingAreaSquareMeters == 0.5
              && partition.solidAreaSquareMeters == 0.5
              && partition.openingAreaFraction == 0.5,
          "regional fragment opening partitions one wall without area repair");
    check(opened.openings[0].openingStableId == 1000
              && opened.openings[0].patchCount == 1
              && opened.openings[0].areaSquareMeters == 0.5
              && opened.connectedComponents[0].baseComponentCount == 2
              && opened.connectedComponents[0].fragmentCount
                  == fragments.fragments.size()
              && opened.connectedComponents[0].openingPatchCount == 1,
          "regional fragment opening merges both pressure components exactly");
    checkNear(opened.connectedComponents[0].volumeCubicMeters,
              16.0, 2.0e-15,
              "regional fragment opening component retains total domain volume");
    check(opened.baseComponents[0].connectedComponentIndex == 0
              && opened.baseComponents[1].connectedComponentIndex == 0
              && opened.baseComponents[0].connectedComponentStableId
                  == opened.connectedComponents[0].stableId
              && opened.baseComponents[1].connectedComponentStableId
                  == opened.connectedComponents[0].stableId,
          "regional fragment opening publishes the complete base-component union");
    check(opened.totalOpeningAreaSquareMeters == 0.5
              && opened.totalTouchedWallAreaSquareMeters == 1.0
              && opened.totalSolidAreaOnTouchedWallsSquareMeters == 0.5
              && opened.wallAreaPartitionResidualSquareMeters == 0.0,
          "regional fragment opening closes its global wall-area partition");
    validatePlanarPressureRegionFragmentOpenings(
        opened, geometry, sweep, fragments, topology, single);

    auto centroidAuthored = single;
    centroidAuthored[0].authoredWrappedCentroidMeters =
        Vector3{-0.8, -0.75, -0.5};
    const auto centroidOpened = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, centroidAuthored);
    check(centroidOpened.patches[0].usesAuthoredCentroid
              && centroidOpened.patches[0].wrappedCentroidMeters
                  == *centroidAuthored[0].authoredWrappedCentroidMeters
              && centroidOpened.partitions[0].hasExactSubtileCentroids,
          "regional fragment opening retains authored sub-tile centroid provenance");
    checkNear(
        centroidOpened.partitions[0]
            .openingAreaWeightedCentroidMeters.y,
        -0.75, 0.0,
        "regional fragment opening retains exact aperture first moment");
    checkNear(
        centroidOpened.partitions[0]
            .solidAreaWeightedCentroidMeters.y,
        -0.25, 0.0,
        "regional fragment opening derives exact retained-solid first moment");
    checkNear(
        centroidOpened.partitions[0]
            .solidAreaWeightedCentroidMeters.z,
        -0.5, 0.0,
        "regional fragment opening preserves the unaffected centroid coordinate");
    validatePlanarPressureRegionFragmentOpenings(
        centroidOpened, geometry, sweep, fragments, topology,
        centroidAuthored);

    std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition> multiple{
        {103, 2000, 20, GridFaceAxis::X, 1, 0, 1, 2, 1, 0.5},
        {102, 1000, 10, GridFaceAxis::X, 1, 1, 0, 1, 2, 0.25},
        {101, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.25},
    };
    auto reversed = multiple;
    std::ranges::reverse(reversed);
    const auto multi = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, multiple);
    const auto reversedMulti = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, reversed);
    check(multi == reversedMulti && multi.patches.size() == 3
              && multi.partitions.size() == 3
              && multi.openings.size() == 2
              && multi.openings[0].openingStableId == 1000
              && multi.openings[0].patchCount == 2
              && multi.openings[0].areaSquareMeters == 0.5
              && multi.openings[1].openingStableId == 2000
              && multi.openings[1].patchCount == 1
              && multi.openings[1].areaSquareMeters == 0.5
              && multi.totalOpeningAreaSquareMeters == 1.0
              && multi.totalTouchedWallAreaSquareMeters == 3.0
              && multi.totalSolidAreaOnTouchedWallsSquareMeters == 2.0,
          "regional fragment openings canonicalize multi-patch authored order");

    const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        sharedTile{
            {110, 3000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.3},
            {111, 4000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.2},
        };
    const auto shared = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, sharedTile);
    check(shared.partitions.size() == 1
              && shared.partitions[0].openingPatchCount == 2
              && shared.partitions[0].openingAreaSquareMeters == 0.5
              && shared.partitions[0].solidAreaSquareMeters == 0.5,
          "regional fragment openings aggregate multiple apertures on one wall tile");

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
        const auto wall = std::ranges::find_if(
            axisTopology.links,
            [](const auto& link) {
                return link.kind
                    == PlanarPressureRegionFragmentFaceKind::
                        PressureLayerWall;
            });
        check(wall != axisTopology.links.end(),
              "regional fragment opening finds a wall on every profile axis");
        if (wall == axisTopology.links.end()) continue;
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition> definition{
            {100, 1000, wall->surfaceStableId, wall->axis,
             wall->i, wall->j, wall->k,
             wall->minusRegionStableId, wall->plusRegionStableId,
             0.5 * wall->areaSquareMeters},
        };
        const auto axisOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, axisSweep, axisFragments, axisTopology,
                definition);
        check(axisOpenings.profileAxis == axis
                  && axisOpenings.patches.size() == 1
                  && axisOpenings.patches[0].sourceFaceLinkStableId
                      == wall->stableId
                  && axisOpenings.partitions[0].openingAreaFraction == 0.5
                  && axisOpenings.connectedComponents.size() == 1,
              "regional fragment opening ownership closes on every profile axis");
    }

    auto corrupt = opened;
    corrupt.patches[0].areaSquareMeters += 0.1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentOpenings(
            corrupt, geometry, sweep, fragments, topology, single); },
        "regional fragment-opening validation rejects patch corruption");
    corrupt = opened;
    corrupt.baseComponents[0].connectedComponentStableId += 1;
    expectRejected(
        [&] { validatePlanarPressureRegionFragmentOpenings(
            corrupt, geometry, sweep, fragments, topology, single); },
        "regional fragment-opening validation rejects component corruption");

    auto invalid = single;
    invalid.push_back(
        {100, 2000, 20, GridFaceAxis::X, 1, 0, 0, 2, 1, 0.25});
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject duplicate patch identity");
    invalid = single;
    invalid[0].surfaceStableId = 999;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject a missing wall tile");
    invalid = single;
    invalid[0].openingStableId = 0;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject invalid stable identity");
    invalid = single;
    invalid[0].axis = GridFaceAxis::Y;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject a foreign wall orientation");
    invalid = single;
    std::swap(invalid[0].negativeSideRegionStableId,
              invalid[0].positiveSideRegionStableId);
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject reversed side regions");
    invalid = single;
    invalid[0].areaSquareMeters = 1.01;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject a patch larger than its wall");
    invalid = single;
    invalid[0].areaSquareMeters = 0.0;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject nonpositive patch area");
    invalid = sharedTile;
    invalid[1].areaSquareMeters = 0.8;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject aggregate overfill on one wall");
    invalid = centroidAuthored;
    invalid[0].authoredWrappedCentroidMeters =
        Vector3{-0.7, -0.75, -0.5};
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject a centroid off its wall plane");
    invalid = centroidAuthored;
    invalid[0].authoredWrappedCentroidMeters =
        Vector3{-0.8, 0.1, -0.5};
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject a centroid outside its wall tile");
    invalid = centroidAuthored;
    invalid[0].areaSquareMeters = 0.9;
    invalid[0].authoredWrappedCentroidMeters =
        Vector3{-0.8, -0.9, -0.5};
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject an impossible retained-solid centroid");
    invalid = centroidAuthored;
    invalid[0].areaSquareMeters = 1.0;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject a full-tile first-moment mismatch");
    invalid = multiple;
    invalid[2].openingStableId = 2000;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, invalid)); },
        "regional fragment openings reject one opening spread across surfaces");

    auto limits = PlanarPressureRegionFragmentOpeningLimits{};
    limits.maximumPatches = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, multiple, limits)); },
        "regional fragment openings enforce the patch limit");
    limits = {};
    limits.maximumPartitions = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, multiple, limits)); },
        "regional fragment openings enforce the partition limit");
    limits = {};
    limits.maximumOpenings = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, multiple, limits)); },
        "regional fragment openings enforce the opening limit");
    limits = {};
    limits.maximumConnectedComponents = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, {}, limits)); },
        "regional fragment openings enforce the connected-component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, single, limits)); },
        "regional fragment openings enforce the owned-byte limit");
    limits = {};
    limits.maximumWorkingBytes = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, single, limits)); },
        "regional fragment openings enforce the working-byte limit");
    limits = {};
    limits.topologyLimits.maximumLinks = 1;
    expectRejected(
        [&] { static_cast<void>(buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, single, limits)); },
        "regional fragment openings enforce nested topology limits");
}

void testPlanarRegionalFragmentOpeningFluxState() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 1.0);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);

    const std::vector<
        PlanarPressureRegionFragmentOpeningPatchDefinition> definitions{
        {100, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.5},
    };
    const auto openings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, definitions);
    const std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        samples{{100, 3.2}};
    const auto state = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, definitions, openings,
        samples);
    const auto repeated =
        buildPlanarPressureRegionFragmentOpeningFluxState(
            geometry, sweep, fragments, topology, definitions, openings,
            samples);
    check(state == repeated
              && state.version
                  == planarPressureRegionFragmentOpeningFluxVersion
              && state.fingerprint != 0
              && state.sourceOpeningFingerprint == openings.fingerprint
              && state.sourceFragmentFingerprint == fragments.fingerprint
              && state.sourceTopologyFingerprint == topology.fingerprint
              && state.profileAxis == GridFaceAxis::X
              && state.patches.size() == 1
              && state.openings.size() == 1
              && state.fragments.size() == fragments.fragments.size()
              && state.baseComponents.size() == topology.components.size()
              && state.connectedComponents.size() == 1
              && state.ownedStorageBytes > 0
              && state.workingStorageBytes > 0,
          "regional fragment opening flux is deterministic and source-bound");
    const auto& patch = state.patches[0];
    check(patch.patchStableId == 100
              && patch.openingStableId == 1000
              && patch.sourceFaceLinkStableId
                  == openings.patches[0].sourceFaceLinkStableId
              && patch.minusFragmentIndex
                  == openings.patches[0].minusFragmentIndex
              && patch.plusFragmentIndex
                  == openings.patches[0].plusFragmentIndex
              && patch.areaSquareMeters == 0.5
              && patch.relativeNormalVelocityMetersPerSecond == 3.2,
          "regional fragment opening flux retains exact patch provenance");
    checkNear(patch.relativeVolumeFlowRateCubicMetersPerSecond,
              1.6, 8.0e-15,
              "regional fragment opening flux uses area times relative velocity");
    checkNear(
        state.fragments[patch.minusFragmentIndex]
            .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
        1.6, 8.0e-15,
        "regional fragment opening flux leaves the minus fragment");
    checkNear(
        state.fragments[patch.plusFragmentIndex]
            .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
        -1.6, 8.0e-15,
        "regional fragment opening flux enters the plus fragment");
    checkNear(
        state.baseComponents[patch.minusBaseComponentIndex]
            .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
        1.6, 8.0e-15,
        "regional fragment opening flux leaves the minus base component");
    checkNear(
        state.baseComponents[patch.plusBaseComponentIndex]
            .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
        -1.6, 8.0e-15,
        "regional fragment opening flux enters the plus base component");
    check(state.openings[0].patchCount == 1
              && state.openings[0].areaSquareMeters == 0.5
              && state.openings[0]
                     .relativeVolumeFlowRateCubicMetersPerSecond == 1.6
              && state.openings[0]
                     .areaWeightedRelativeNormalVelocityMetersPerSecond
                  == 3.2,
          "regional fragment opening flux publishes its opening aggregate");
    checkNear(
        state.connectedComponents[0]
            .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
        0.0, state.conservationToleranceCubicMetersPerSecond,
        "regional fragment opening flux closes its connected component");
    checkNear(state.globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond,
              0.0, state.conservationToleranceCubicMetersPerSecond,
              "regional fragment opening flux closes globally");
    validatePlanarPressureRegionFragmentOpeningFluxState(
        state, geometry, sweep, fragments, topology, definitions, openings,
        samples);
    validatePlanarPressureRegionFragmentOpeningFluxStateIntegrity(state);

    const auto emptyOpenings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, {});
    const auto empty = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, {}, emptyOpenings, {});
    check(empty.patches.empty() && empty.openings.empty()
              && empty.fragments.size() == fragments.fragments.size()
              && empty.baseComponents.size() == 2
              && empty.connectedComponents.size() == 2
              && empty.maximumAbsoluteRelativeNormalVelocityMetersPerSecond
                  == 0.0
              && empty.totalAbsolutePatchVolumeFlowRateCubicMetersPerSecond
                  == 0.0
              && empty.globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond
                  == 0.0
              && empty.conservationToleranceCubicMetersPerSecond == 0.0,
          "regional fragment opening flux preserves a sealed zero-flow overlay");

    const std::vector<
        PlanarPressureRegionFragmentOpeningPatchDefinition> multiple{
        {103, 2000, 20, GridFaceAxis::X, 1, 0, 1, 2, 1, 0.5},
        {102, 1000, 10, GridFaceAxis::X, 1, 1, 0, 1, 2, 0.25},
        {101, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.25},
    };
    const auto multipleOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, multiple);
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        multipleSamples{{103, -1.0}, {101, 4.0}, {102, 2.0}};
    auto reversedSamples = multipleSamples;
    std::ranges::reverse(reversedSamples);
    const auto multi = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, multiple, multipleOpenings,
        multipleSamples);
    const auto canonicalMulti =
        buildPlanarPressureRegionFragmentOpeningFluxState(
            geometry, sweep, fragments, topology, multiple,
            multipleOpenings, reversedSamples);
    check(multi == canonicalMulti && multi.patches.size() == 3
              && multi.patches[0].patchStableId == 101
              && multi.patches[1].patchStableId == 102
              && multi.patches[2].patchStableId == 103
              && multi.openings[0]
                     .relativeVolumeFlowRateCubicMetersPerSecond == 1.5
              && multi.openings[0]
                     .areaWeightedRelativeNormalVelocityMetersPerSecond
                  == 3.0
              && multi.openings[1]
                     .relativeVolumeFlowRateCubicMetersPerSecond == -0.5
              && multi.openings[1]
                     .areaWeightedRelativeNormalVelocityMetersPerSecond
                  == -1.0,
          "regional fragment opening flux canonicalizes samples and aggregates patches");

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
        const auto wall = std::ranges::find_if(
            axisTopology.links,
            [](const auto& link) {
                return link.kind
                    == PlanarPressureRegionFragmentFaceKind::
                        PressureLayerWall;
            });
        check(wall != axisTopology.links.end(),
              "regional fragment opening flux finds a wall on every axis");
        if (wall == axisTopology.links.end()) continue;
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition>
            axisDefinitions{
                {100, 1000, wall->surfaceStableId, wall->axis,
                 wall->i, wall->j, wall->k,
                 wall->minusRegionStableId, wall->plusRegionStableId,
                 0.5 * wall->areaSquareMeters},
            };
        const auto axisOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, axisSweep, axisFragments, axisTopology,
                axisDefinitions);
        const std::vector<
            PlanarPressureRegionFragmentOpeningVelocitySample>
            axisSamples{{100, 1.25}};
        const auto axisState =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, axisSweep, axisFragments, axisTopology,
                axisDefinitions, axisOpenings, axisSamples);
        check(axisState.profileAxis == axis
                  && axisState.patches[0]
                         .relativeVolumeFlowRateCubicMetersPerSecond
                      == 0.625 * wall->areaSquareMeters
                  && axisState.globalOutwardRelativeVolumeFlowRateCubicMetersPerSecond
                      == 0.0,
              "regional fragment opening flux uses physical patch area on every axis");
    }

    auto current = layers;
    current[0].physicalPlaneCoordinateMeters -= 0.1;
    current[1].physicalPlaneCoordinateMeters += 0.1;
    const auto breathingSweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, current, 0.5);
    const auto breathingFragments = buildPlanarPressureRegionFragments(
        geometry, breathingSweep);
    const auto breathingTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, breathingSweep, breathingFragments);
    const auto breathingWall = std::ranges::find_if(
        breathingTopology.links,
        [](const auto& link) {
            return link.kind
                    == PlanarPressureRegionFragmentFaceKind::PressureLayerWall
                && link.surfaceStableId == 10 && link.j == 0 && link.k == 0;
        });
    check(breathingWall != breathingTopology.links.end(),
          "regional fragment opening flux finds the breathing intake tile");
    if (breathingWall != breathingTopology.links.end()) {
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition>
            breathingDefinitions{
                {100, 1000, breathingWall->surfaceStableId,
                 breathingWall->axis, breathingWall->i, breathingWall->j,
                 breathingWall->k, breathingWall->minusRegionStableId,
                 breathingWall->plusRegionStableId, 0.5},
            };
        const auto breathingOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, breathingSweep, breathingFragments,
                breathingTopology, breathingDefinitions);
        const auto breathingState =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, breathingSweep, breathingFragments,
                breathingTopology, breathingDefinitions,
                breathingOpenings, samples);
        const auto volumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, breathingSweep, breathingFragments,
                breathingTopology);
        for (const auto& component : volumeRates.components) {
            checkNear(
                component.geometryVolumeChangeRateCubicMetersPerSecond
                    + breathingState.baseComponents[component.componentIndex]
                          .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
                0.0, 1.6e-14,
                "prescribed regional opening flux makes breathing component-compatible");
        }
    }

    auto corrupt = state;
    corrupt.patches[0].relativeVolumeFlowRateCubicMetersPerSecond += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningFluxStateIntegrity(
                corrupt);
        },
        "regional fragment-opening flux integrity rejects corruption");
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample> invalid;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, invalid)); },
        "regional fragment opening flux rejects a missing sample");
    invalid = {{999, 3.2}};
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, invalid)); },
        "regional fragment opening flux rejects a foreign sample");
    invalid = {{100, std::numeric_limits<double>::quiet_NaN()}};
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, invalid)); },
        "regional fragment opening flux rejects non-finite velocity");
    invalid = {{101, 1.0}, {101, 2.0}, {103, 3.0}};
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, multiple,
                multipleOpenings, invalid)); },
        "regional fragment opening flux rejects duplicate sample identity");

    auto limits = PlanarPressureRegionFragmentOpeningFluxLimits{};
    limits.maximumPatches = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, multiple,
                multipleOpenings, multipleSamples, limits)); },
        "regional fragment opening flux enforces the patch limit");
    limits = {};
    limits.maximumOpenings = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, multiple,
                multipleOpenings, multipleSamples, limits)); },
        "regional fragment opening flux enforces the opening limit");
    limits = {};
    limits.maximumFragments = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, samples, limits)); },
        "regional fragment opening flux enforces the fragment limit");
    limits = {};
    limits.maximumBaseComponents = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, samples, limits)); },
        "regional fragment opening flux enforces the base-component limit");
    limits = {};
    limits.maximumConnectedComponents = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, {}, emptyOpenings,
                {}, limits)); },
        "regional fragment opening flux enforces the connected-component limit");
    limits = {};
    limits.maximumOwnedBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, samples, limits)); },
        "regional fragment opening flux enforces the owned-byte limit");
    limits = {};
    limits.maximumWorkingBytes = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, samples, limits)); },
        "regional fragment opening flux enforces the working-byte limit");
    limits = {};
    limits.openingLimits.maximumPatches = 1;
    expectRejected(
        [&] { static_cast<void>(
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, multiple,
                multipleOpenings, multipleSamples, limits)); },
        "regional fragment opening flux enforces nested opening limits");
}

void testPlanarRegionalOpeningFluxPressureProjection() {
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
        const auto pressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, sweep, fragments, topology);
        const auto volumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology);
        const auto wall = std::ranges::find_if(
            topology.links,
            [](const auto& link) {
                return link.kind
                        == PlanarPressureRegionFragmentFaceKind::
                            PressureLayerWall
                    && link.surfaceStableId == 10;
            });
        check(wall != topology.links.end(),
              "opening-flux moving projection finds an intake wall on every axis");
        if (wall == topology.links.end()) continue;
        const double patchArea = 0.5 * wall->areaSquareMeters;
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition> definitions{
            {100, 1000, wall->surfaceStableId, wall->axis,
             wall->i, wall->j, wall->k,
             wall->minusRegionStableId, wall->plusRegionStableId,
             patchArea},
        };
        const auto openings = buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, definitions);
        const double requiredFlow =
            -volumeRates.components[wall->minusComponentIndex]
                 .geometryVolumeChangeRateCubicMetersPerSecond;
        const std::vector<
            PlanarPressureRegionFragmentOpeningVelocitySample> samples{
            {100, requiredFlow / patchArea},
        };
        const auto openingFlux =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, samples);
        checkNear(
            volumeRates.components[wall->minusComponentIndex]
                    .geometryVolumeChangeRateCubicMetersPerSecond
                + openingFlux.baseComponents[wall->minusComponentIndex]
                      .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
            0.0, 1.6e-14,
            "opening flux balances the shrinking projection component");
        checkNear(
            volumeRates.components[wall->plusComponentIndex]
                    .geometryVolumeChangeRateCubicMetersPerSecond
                + openingFlux.baseComponents[wall->plusComponentIndex]
                      .outwardRelativeVolumeFlowRateCubicMetersPerSecond,
            0.0, 1.6e-14,
            "opening flux balances the expanding projection component");

        PlanarPressureRegionFragmentPressureProjectionSettings settings;
        settings.densityKgPerCubicMeter = 1.2;
        settings.timeStepSeconds = 0.5;
        settings.absoluteContinuityToleranceCubicMetersPerSecond = 1.0e-11;
        settings.relativeContinuityTolerance = 1.0e-10;
        settings.pressureSolve.absoluteResidualTolerancePascalsMeters =
            1.0e-13;
        settings.pressureSolve.relativeResidualTolerance = 1.0e-12;
        settings.pressureSolve.maximumIterations = 300;
        std::vector<double> velocity(topology.links.size(), 0.0);
        std::vector<double> pressure(pressureOperator.rows.size(), 0.0);
        auto repeatedVelocity = velocity;
        auto repeatedPressure = pressure;
        const auto diagnostics =
            projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                pressureOperator, geometry, sweep, fragments, topology,
                volumeRates, definitions, openings, openingFlux, samples,
                velocity, pressure, settings);
        const auto repeated =
            projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                pressureOperator, geometry, sweep, fragments, topology,
                volumeRates, definitions, openings, openingFlux, samples,
                repeatedVelocity, repeatedPressure, settings);
        check(diagnostics == repeated && velocity == repeatedVelocity
                  && pressure == repeatedPressure
                  && diagnostics.accepted && diagnostics.finite
                  && diagnostics.usesMovingVolumeRates
                  && diagnostics.usesOpeningFlux
                  && diagnostics.volumeRateFingerprint
                      == volumeRates.fingerprint
                  && diagnostics.openingFingerprint == openings.fingerprint
                  && diagnostics.openingFluxFingerprint
                      == openingFlux.fingerprint
                  && diagnostics.openingPatchCount == 1
                  && diagnostics.pressureSolve.compatible
                  && diagnostics.pressureSolve.converged,
              "prescribed opening flux makes breathing projection deterministic and acceptable");
        checkNear(
            diagnostics.maximumAbsoluteOpeningFragmentOutwardFlowRateCubicMetersPerSecond,
            std::abs(requiredFlow), 1.6e-14,
            "moving projection reports the prescribed fragment opening flow");
        check(
            diagnostics.correctedContinuityResidualMaximumCubicMetersPerSecond
                    <= diagnostics.continuityToleranceCubicMetersPerSecond
                && diagnostics
                       .maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond
                    <= diagnostics.continuityToleranceCubicMetersPerSecond,
            "opening-flux moving projection closes local and component continuity");
        check(std::ranges::any_of(
                  velocity, [](const double value) { return value != 0.0; })
                  && std::ranges::any_of(
                      pressure,
                      [](const double value) { return value != 0.0; }),
              "opening-flux moving projection redistributes the localized intake flow");
        for (const auto& link : topology.links) {
            if (link.kind
                == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
                check(velocity[link.linkIndex] == 0.0,
                      "opening flux remains separate from sealed topology-link velocity");
            }
        }

        std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
            reversedSamples{{100, -samples[0]
                                      .relativeNormalVelocityMetersPerSecond}};
        const auto reversedFlux =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, reversedSamples);
        std::vector<double> rejectedVelocity(topology.links.size(), 0.0);
        std::vector<double> rejectedPressure(
            pressureOperator.rows.size(), 0.25);
        const auto originalRejectedVelocity = rejectedVelocity;
        const auto originalRejectedPressure = rejectedPressure;
        const auto reversed =
            projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                pressureOperator, geometry, sweep, fragments, topology,
                volumeRates, definitions, openings, reversedFlux,
                reversedSamples, rejectedVelocity, rejectedPressure,
                settings);
        check(!reversed.accepted && reversed.finite
                  && !reversed.pressureSolve.compatible
                  && !reversed.pressureSolve.converged
                  && rejectedVelocity == originalRejectedVelocity
                  && rejectedPressure == originalRejectedPressure,
              "reversed prescribed opening flux rejects and rolls back projection");
        checkNear(
            reversed.maximumAbsolutePredictedComponentContinuityResidualCubicMetersPerSecond,
            2.0 * std::abs(requiredFlow), 3.2e-14,
            "reversed opening flux doubles the breathing component deficit");

        if (axis == GridFaceAxis::X) {
            auto corruptFlux = openingFlux;
            corruptFlux.patches[0]
                .relativeVolumeFlowRateCubicMetersPerSecond += 0.1;
            expectRejected(
                [&] {
                    auto candidateVelocity = originalRejectedVelocity;
                    auto candidatePressure = originalRejectedPressure;
                    static_cast<void>(
                        projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                            pressureOperator, geometry, sweep, fragments,
                            topology, volumeRates, definitions, openings,
                            corruptFlux, samples, candidateVelocity,
                            candidatePressure, settings));
                },
                "opening-flux moving projection rejects corrupted flux state");
            expectRejected(
                [&] {
                    auto candidateVelocity = originalRejectedVelocity;
                    auto candidatePressure = originalRejectedPressure;
                    static_cast<void>(
                        projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                            pressureOperator, geometry, sweep, fragments,
                            topology, volumeRates, definitions, openings,
                            openingFlux, reversedSamples, candidateVelocity,
                            candidatePressure, settings));
                },
                "opening-flux moving projection rejects foreign velocity samples");
            auto limits =
                PlanarPressureRegionFragmentPressureProjectionLimits{};
            limits.openingFluxLimits.maximumFragments = 1;
            expectRejected(
                [&] {
                    auto candidateVelocity = originalRejectedVelocity;
                    auto candidatePressure = originalRejectedPressure;
                    static_cast<void>(
                        projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                            pressureOperator, geometry, sweep, fragments,
                            topology, volumeRates, definitions, openings,
                            openingFlux, samples, candidateVelocity,
                            candidatePressure, settings, limits));
                },
                "opening-flux moving projection enforces nested flux limits");
        }
    }
}

void testPlanarRegionalPressureDrivenOpeningProjection() {
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
        const auto base = buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
        const auto volumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology);
        const auto wall = std::ranges::find_if(
            topology.links,
            [](const auto& link) {
                return link.kind
                        == PlanarPressureRegionFragmentFaceKind::
                            PressureLayerWall
                    && link.surfaceStableId == 10;
            });
        check(wall != topology.links.end(),
              "pressure-driven opening projection finds an intake wall on every axis");
        if (wall == topology.links.end()) continue;
        const double patchArea = 0.5 * wall->areaSquareMeters;
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition> definitions{
            {100, 1000, wall->surfaceStableId, wall->axis,
             wall->i, wall->j, wall->k,
             wall->minusRegionStableId, wall->plusRegionStableId,
             patchArea},
        };
        const auto openings = buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, definitions);
        const auto baseVelocityMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology);
        const auto openingVelocityMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, sweep, fragments, topology,
                baseVelocityMetric, definitions, openings);
        const auto apertureDof = std::ranges::find_if(
            openingVelocityMetric.dofs,
            [](const auto& dof) {
                return dof.kind
                    == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                        OpeningPatch;
            });
        check(apertureDof != openingVelocityMetric.dofs.end()
                  && openingVelocityMetric.profileAxis == axis
                  && apertureDof->axis == axis
                  && apertureDof->sourceOpeningPatchStableId == 100
                  && openingVelocityMetric.openingPatchDofCount == 1
                  && openingVelocityMetric.solidWallTraceDofCount
                      == 2 * topology.pressureLayerWallLinkCount,
              "opening velocity metric retains the same aperture inertia contract on every axis");
        checkNear(
            std::max({
                std::abs(openingVelocityMetric
                             .domainVolumeClosureResidualByAxisCubicMeters.x),
                std::abs(openingVelocityMetric
                             .domainVolumeClosureResidualByAxisCubicMeters.y),
                std::abs(openingVelocityMetric
                             .domainVolumeClosureResidualByAxisCubicMeters.z),
                openingVelocityMetric
                    .maximumAbsoluteFragmentVolumeClosureResidualCubicMeters,
                openingVelocityMetric
                    .maximumAbsoluteComponentVolumeClosureResidualCubicMeters}),
            0.0, 2.0e-13,
            "opening velocity metric closes all-axis fragment inertia");
        validatePlanarPressureRegionFragmentOpeningVelocityMetric(
            openingVelocityMetric, geometry, sweep, fragments, topology,
            baseVelocityMetric, definitions, openings);
        const Vector3 uniformVelocity{0.3, -0.2, 0.1};
        std::vector<double> uniformNormal(
            openingVelocityMetric.dofs.size(), 0.0);
        std::vector<double> uniformMaterial(
            openingVelocityMetric.dofs.size(), 0.0);
        std::vector<double> uniformRelative(
            openingVelocityMetric.dofs.size(), 0.0);
        for (const auto& dof : openingVelocityMetric.dofs) {
            const double component = dof.axis == GridFaceAxis::X
                ? uniformVelocity.x
                : dof.axis == GridFaceAxis::Y
                ? uniformVelocity.y : uniformVelocity.z;
            uniformNormal[dof.dofIndex] = component;
            if (dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SharedRegionGrid) {
                uniformRelative[dof.dofIndex] = component;
            } else {
                uniformMaterial[dof.dofIndex] = component;
            }
        }
        const auto uniformState =
            buildPlanarPressureRegionFragmentOpeningVelocityState(
                openingVelocityMetric, uniformNormal, uniformMaterial,
                uniformRelative, 1.2);
        const double uniformMass = 1.2
            * geometry.cellVolumeCubicMeters()
            * static_cast<double>(geometry.cellCount());
        checkNear(uniformState.momentumKilogramMetersPerSecond.x,
                  uniformMass * uniformVelocity.x, 2.0e-13,
                  "all-axis opening metric preserves uniform X momentum");
        checkNear(uniformState.momentumKilogramMetersPerSecond.y,
                  uniformMass * uniformVelocity.y, 2.0e-13,
                  "all-axis opening metric preserves uniform Y momentum");
        checkNear(uniformState.momentumKilogramMetersPerSecond.z,
                  uniformMass * uniformVelocity.z, 2.0e-13,
                  "all-axis opening metric preserves uniform Z momentum");
        checkNear(uniformState.staggeringKineticEnergyJoules,
                  0.0, 3.0e-13,
                  "all-axis uniform opening flow has zero staggering energy");
        const auto pressureOperator =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                base, geometry, sweep, fragments, topology, definitions,
                openings);
        std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
            samples{{100, 0.0}};
        auto openingFlux =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, definitions,
                openings, samples);
        const auto predictedOpeningFlux = openingFlux;
        std::vector<double> velocity(topology.links.size(), 0.0);
        std::vector<double> pressure(pressureOperator.rows.size(), 0.0);
        auto repeatedSamples = samples;
        auto repeatedOpeningFlux = openingFlux;
        auto repeatedVelocity = velocity;
        auto repeatedPressure = pressure;

        PlanarPressureRegionFragmentOpeningPressureProjectionSettings settings;
        settings.densityKgPerCubicMeter = 1.2;
        settings.timeStepSeconds = 0.5;
        settings.absoluteContinuityToleranceCubicMetersPerSecond = 2.0e-11;
        settings.relativeContinuityTolerance = 1.0e-10;
        settings.absoluteMomentumResidualToleranceKilogramMetersPerSecond =
            2.0e-12;
        settings.relativeMomentumResidualTolerance = 1.0e-10;
        settings.absoluteEnergyResidualToleranceJoules = 2.0e-11;
        settings.relativeEnergyResidualTolerance = 1.0e-10;
        settings.pressureSolve.absoluteResidualTolerancePascalsMeters =
            1.0e-13;
        settings.pressureSolve.relativeResidualTolerance = 0.0;
        settings.pressureSolve.maximumIterations = 300;
        const auto diagnostics =
            projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                pressureOperator, base, geometry, sweep, fragments, topology,
                volumeRates, definitions, openings, velocity, samples,
                openingFlux, pressure, settings);
        const auto repeated =
            projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                pressureOperator, base, geometry, sweep, fragments, topology,
                volumeRates, definitions, openings, repeatedVelocity,
                repeatedSamples, repeatedOpeningFlux, repeatedPressure,
                settings);
        check(diagnostics == repeated && velocity == repeatedVelocity
                  && samples == repeatedSamples
                  && openingFlux == repeatedOpeningFlux
                  && pressure == repeatedPressure
                  && diagnostics.accepted && diagnostics.finite
                  && diagnostics.energyAccepted
                  && diagnostics.pressureSolve.compatible
                  && diagnostics.pressureSolve.converged
                  && diagnostics.pressureSolve.usesOpeningPressureOperator
                  && diagnostics.pressureOperatorFingerprint
                      == pressureOperator.fingerprint
                  && diagnostics.basePressureOperatorFingerprint
                      == base.fingerprint
                  && diagnostics.volumeRateFingerprint
                      == volumeRates.fingerprint
                  && diagnostics.openingFingerprint == openings.fingerprint
                  && diagnostics.predictedOpeningFluxFingerprint
                      == predictedOpeningFlux.fingerprint
                  && diagnostics.correctedOpeningFluxFingerprint
                      == openingFlux.fingerprint,
              "pressure-driven opening projection is deterministic and source-bound");
        check(diagnostics.predictedContinuityResidualMaximumCubicMetersPerSecond
                  > 0.0
                  && diagnostics
                         .correctedContinuityResidualMaximumCubicMetersPerSecond
                      <= diagnostics.continuityToleranceCubicMetersPerSecond
                  && diagnostics
                         .maximumAbsoluteCorrectedComponentContinuityResidualCubicMetersPerSecond
                      <= diagnostics.continuityToleranceCubicMetersPerSecond,
              "pressure-driven opening projection closes local and connected continuity");
        const double requiredFlow =
            -volumeRates.components[wall->minusComponentIndex]
                 .geometryVolumeChangeRateCubicMetersPerSecond;
        check(samples[0].relativeNormalVelocityMetersPerSecond > 0.0,
              "pressure-driven intake accelerates from rest toward the growing pocket");
        checkNear(
            openingFlux.patches[0]
                .relativeVolumeFlowRateCubicMetersPerSecond,
            requiredFlow, 3.0e-11,
            "pressure-driven intake supplies the complete breathing flow");
        checkNear(
            openingFlux.baseComponents[wall->minusComponentIndex]
                    .outwardRelativeVolumeFlowRateCubicMetersPerSecond
                + volumeRates.components[wall->minusComponentIndex]
                      .geometryVolumeChangeRateCubicMetersPerSecond,
            0.0, 3.0e-11,
            "pressure-driven intake balances the shrinking base component");
        checkNear(
            openingFlux.baseComponents[wall->plusComponentIndex]
                    .outwardRelativeVolumeFlowRateCubicMetersPerSecond
                + volumeRates.components[wall->plusComponentIndex]
                      .geometryVolumeChangeRateCubicMetersPerSecond,
            0.0, 3.0e-11,
            "pressure-driven intake balances the growing base component");
        check(diagnostics.maximumAbsoluteOpeningVelocityCorrectionMetersPerSecond
                      > 0.0
                  && diagnostics.maximumAbsoluteGridVelocityCorrectionMetersPerSecond
                      > 0.0
                  && diagnostics.correctionKineticEnergyJoules > 0.0
                  && std::abs(diagnostics.workEnergyResidualJoules)
                      <= diagnostics.energyResidualToleranceJoules
                  && diagnostics.maximumAbsoluteWorkEnergyResidualJoules
                      <= diagnostics.energyResidualToleranceJoules
                  && std::abs(diagnostics.affineEnergyResidualJoules)
                      <= diagnostics.energyResidualToleranceJoules
                  && diagnostics
                         .maximumAbsoluteMomentumImpulseResidualKilogramMetersPerSecond
                      <= diagnostics
                             .momentumResidualToleranceKilogramMetersPerSecond,
              "pressure-driven opening projection closes impulse and moving energy");
        for (const auto& link : topology.links) {
            if (link.kind
                == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
                check(velocity[link.linkIndex] == 0.0,
                      "pressure-driven openings do not leak through remaining wall links");
            }
        }

        if (axis == GridFaceAxis::X) {
            auto truncatedSettings = settings;
            truncatedSettings.pressureSolve
                .absoluteResidualTolerancePascalsMeters = 1.0e-16;
            truncatedSettings.pressureSolve.relativeResidualTolerance = 0.0;
            truncatedSettings.pressureSolve.maximumIterations = 1;
            std::vector<double> rejectedVelocity(topology.links.size(), 0.0);
            std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
                rejectedSamples{{100, 0.0}};
            auto rejectedFlux =
                buildPlanarPressureRegionFragmentOpeningFluxState(
                    geometry, sweep, fragments, topology, definitions,
                    openings, rejectedSamples);
            std::vector<double> rejectedPressure(
                pressureOperator.rows.size(), 0.125);
            const auto originalVelocity = rejectedVelocity;
            const auto originalSamples = rejectedSamples;
            const auto originalFlux = rejectedFlux;
            const auto originalPressure = rejectedPressure;
            const auto truncated =
                projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                    pressureOperator, base, geometry, sweep, fragments,
                    topology, volumeRates, definitions, openings,
                    rejectedVelocity, rejectedSamples, rejectedFlux,
                    rejectedPressure, truncatedSettings);
            check(!truncated.accepted && truncated.finite
                      && truncated.pressureSolve.compatible
                      && !truncated.pressureSolve.converged
                      && rejectedVelocity == originalVelocity
                      && rejectedSamples == originalSamples
                      && rejectedFlux == originalFlux
                      && rejectedPressure == originalPressure,
                  "truncated pressure-driven opening projection rolls back every field");

            const auto emptyOpenings =
                buildPlanarPressureRegionFragmentOpenings(
                    geometry, sweep, fragments, topology, {});
            const auto emptyOperator =
                buildPlanarPressureRegionFragmentOpeningPressureOperator(
                    base, geometry, sweep, fragments, topology, {},
                    emptyOpenings);
            std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
                emptySamples;
            auto emptyFlux =
                buildPlanarPressureRegionFragmentOpeningFluxState(
                    geometry, sweep, fragments, topology, {}, emptyOpenings,
                    emptySamples);
            std::vector<double> emptyVelocity(topology.links.size(), 0.0);
            std::vector<double> emptyPressure(emptyOperator.rows.size(), 0.25);
            const auto originalEmptyPressure = emptyPressure;
            const auto sealed =
                projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                    emptyOperator, base, geometry, sweep, fragments,
                    topology, volumeRates, {}, emptyOpenings, emptyVelocity,
                    emptySamples, emptyFlux, emptyPressure, settings);
            check(!sealed.accepted && sealed.finite
                      && !sealed.pressureSolve.compatible
                      && emptyPressure == originalEmptyPressure,
                  "empty opening overlay preserves sealed breathing incompatibility");

            auto corruptFlux = predictedOpeningFlux;
            corruptFlux.patches[0]
                .relativeVolumeFlowRateCubicMetersPerSecond += 0.1;
            expectRejected([&] {
                auto candidateVelocity = originalVelocity;
                auto candidateSamples = originalSamples;
                auto candidatePressure = originalPressure;
                static_cast<void>(
                    projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                        pressureOperator, base, geometry, sweep, fragments,
                        topology, volumeRates, definitions, openings,
                        candidateVelocity, candidateSamples, corruptFlux,
                        candidatePressure, settings));
            }, "pressure-driven opening projection rejects corrupted flux state");
            auto limits =
                PlanarPressureRegionFragmentOpeningPressureProjectionLimits{};
            limits.maximumWorkingBytes = 1;
            expectRejected([&] {
                auto candidateVelocity = originalVelocity;
                auto candidateSamples = originalSamples;
                auto candidateFlux = originalFlux;
                auto candidatePressure = originalPressure;
                static_cast<void>(
                    projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                        pressureOperator, base, geometry, sweep, fragments,
                        topology, volumeRates, definitions, openings,
                        candidateVelocity, candidateSamples, candidateFlux,
                        candidatePressure, settings, limits));
            }, "pressure-driven opening projection enforces its working limit");

            expectRejected([&] {
                auto candidateVelocity = originalVelocity;
                candidateVelocity[wall->linkIndex] = 0.1;
                auto candidateSamples = originalSamples;
                auto candidateFlux = originalFlux;
                auto candidatePressure = originalPressure;
                static_cast<void>(
                    projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                        pressureOperator, base, geometry, sweep, fragments,
                        topology, volumeRates, definitions, openings,
                        candidateVelocity, candidateSamples, candidateFlux,
                        candidatePressure, settings));
            }, "pressure-driven opening projection rejects topology-wall leakage");
            auto mismatchedSettings = settings;
            mismatchedSettings.timeStepSeconds = 0.25;
            expectRejected([&] {
                auto candidateVelocity = originalVelocity;
                auto candidateSamples = originalSamples;
                auto candidateFlux = originalFlux;
                auto candidatePressure = originalPressure;
                static_cast<void>(
                    projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                        pressureOperator, base, geometry, sweep, fragments,
                        topology, volumeRates, definitions, openings,
                        candidateVelocity, candidateSamples, candidateFlux,
                        candidatePressure, mismatchedSettings));
            }, "pressure-driven opening projection rejects duration mismatch");
            auto invalidSettings = settings;
            invalidSettings.absoluteEnergyResidualToleranceJoules = 0.0;
            invalidSettings.relativeEnergyResidualTolerance = 0.0;
            expectRejected([&] {
                auto candidateVelocity = originalVelocity;
                auto candidateSamples = originalSamples;
                auto candidateFlux = originalFlux;
                auto candidatePressure = originalPressure;
                static_cast<void>(
                    projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
                        pressureOperator, base, geometry, sweep, fragments,
                        topology, volumeRates, definitions, openings,
                        candidateVelocity, candidateSamples, candidateFlux,
                        candidatePressure, invalidSettings));
            }, "pressure-driven opening projection rejects empty energy tolerances");
        }
    }
}

void testPlanarRegionalOpeningResistance() {
    const auto geometry = grid();
    const auto layers = pocketLayers();
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, layers, layers, 0.02);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions{{100, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.5}};
    const auto openings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, definitions);
    const std::vector<
        PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistance{{100, {8.0, 3.0}}};
    PlanarPressureRegionFragmentOpeningResistanceSettings settings;
    settings.densityKgPerCubicMeter = 1.2;
    settings.timeStepSeconds = 0.02;

    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        samples{{100, 3.2}};
    auto openingFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, definitions, openings, samples);
    const auto sourceFlux = openingFlux;
    auto repeatedSamples = samples;
    auto repeatedFlux = openingFlux;
    const auto diagnostics =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, definitions, openings,
            resistance, samples, openingFlux, settings);
    const auto repeated =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, definitions, openings,
            resistance, repeatedSamples, repeatedFlux, settings);
    check(diagnostics == repeated && samples == repeatedSamples
              && openingFlux == repeatedFlux
              && diagnostics.accepted && diagnostics.finite
              && diagnostics.nonIncreasingKineticEnergy
              && !diagnostics.usesAuthoredPressureDrive
              && diagnostics.authoredPressureWorkJoules == 0.0
              && diagnostics.authoredPressureForceOnOpeningFluidNewtons
                  == Vector3{}
              && diagnostics.sourceOpeningFingerprint == openings.fingerprint
              && diagnostics.sourceOpeningFluxFingerprint
                  == sourceFlux.fingerprint
              && diagnostics.resultOpeningFluxFingerprint
                  == openingFlux.fingerprint
              && diagnostics.resistanceDefinitionFingerprint != 0
              && diagnostics.patches.size() == 1
              && diagnostics.zeroResistancePatchCount == 0,
          "opening resistance advances deterministically with source provenance");
    check(samples[0].relativeNormalVelocityMetersPerSecond > 0.0
              && samples[0].relativeNormalVelocityMetersPerSecond < 3.2
              && openingFlux.patches[0]
                     .relativeVolumeFlowRateCubicMetersPerSecond
                  > 0.0
              && openingFlux.patches[0]
                     .relativeVolumeFlowRateCubicMetersPerSecond
                  < sourceFlux.patches[0]
                        .relativeVolumeFlowRateCubicMetersPerSecond,
          "opening resistance passively reduces positive aperture flow");
    check(diagnostics.patches[0].plugFlow.midpointPressureDropPascals > 0.0
              && diagnostics.patches[0].plugFlow.endpointPressureDropPascals
                  > 0.0
              && diagnostics.dissipatedEnergyJoules > 0.0
              && diagnostics.kineticEnergyChangeJoules < 0.0
              && std::abs(diagnostics.energyResidualJoules)
                  <= diagnostics.energyToleranceJoules
              && diagnostics
                     .maximumAbsoluteMomentumResidualKilogramMetersPerSecond
                  < 1.0e-14,
          "opening resistance closes pressure impulse and passive energy loss");

    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        reverseSamples{{100, -3.2}};
    auto reverseFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, definitions, openings,
        reverseSamples);
    const auto reverse =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, definitions, openings,
            resistance, reverseSamples, reverseFlux, settings);
    check(reverse.accepted
              && reverseSamples[0].relativeNormalVelocityMetersPerSecond < 0.0,
          "opening resistance preserves reverse-flow orientation");
    checkNear(
        reverseSamples[0].relativeNormalVelocityMetersPerSecond,
        -samples[0].relativeNormalVelocityMetersPerSecond, 2.0e-15,
        "opening resistance is odd in aperture velocity");
    checkNear(reverse.dissipatedEnergyJoules,
              diagnostics.dissipatedEnergyJoules, 2.0e-15,
              "opening resistance dissipates equally in either orientation");
    check(reverse.patches[0].plugFlow.midpointPressureDropPascals < 0.0,
          "opening resistance pressure drop follows reverse flow");

    const std::vector<
        PlanarPressureRegionFragmentOpeningResistanceDefinition>
        zeroResistance{{100, {0.0, 0.0}}};
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        identitySamples{{100, 3.2}};
    auto identityFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, definitions, openings,
        identitySamples);
    const auto originalIdentitySamples = identitySamples;
    const auto originalIdentityFlux = identityFlux;
    const auto identity =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, definitions, openings,
            zeroResistance, identitySamples, identityFlux, settings);
    check(identity.accepted && identity.zeroResistancePatchCount == 1
              && identity.patches[0].zeroResistanceIdentity
              && identitySamples == originalIdentitySamples
              && identityFlux == originalIdentityFlux
              && identity.dissipatedEnergyJoules == 0.0
              && identity.energyResidualJoules == 0.0,
          "zero opening resistance is a bit-exact inviscid identity");

    auto drivenSettings = settings;
    drivenSettings.useAuthoredPressureDrive = true;
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        drivenSamples{{100, 0.0}};
    auto drivenFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, definitions, openings,
        drivenSamples);
    auto repeatedDrivenSamples = drivenSamples;
    auto repeatedDrivenFlux = drivenFlux;
    const auto driven =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, definitions, openings,
            resistance, drivenSamples, drivenFlux, drivenSettings);
    const auto repeatedDriven =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, definitions, openings,
            resistance, repeatedDrivenSamples, repeatedDrivenFlux,
            drivenSettings);
    check(driven == repeatedDriven
              && drivenSamples == repeatedDrivenSamples
              && drivenFlux == repeatedDrivenFlux
              && driven.accepted && driven.finite
              && driven.usesAuthoredPressureDrive
              && !driven.nonIncreasingKineticEnergy
              && driven.patches.size() == 1
              && !driven.patches[0].zeroResistanceIdentity
              && driven.patches[0].authoredPressureJumpPascals == 70.0
              && driven.patches[0].drivingPressureRisePascals == -70.0
              && driven.maximumAbsoluteAuthoredPressureJumpPascals == 70.0
              && driven.maximumAbsoluteDrivingPressureRisePascals == 70.0,
          "authored opening pressure drive is deterministic and explicit");
    check(drivenSamples[0].relativeNormalVelocityMetersPerSecond < 0.0
              && driven.patches[0]
                         .authoredPressureForceOnOpeningFluidNewtons.x
                  == -35.0
              && driven.patches[0]
                         .authoredPressureForceOnOpeningFluidNewtons.y
                  == 0.0
              && driven.patches[0]
                         .authoredPressureForceOnOpeningFluidNewtons.z
                  == 0.0
              && driven.authoredPressureForceOnOpeningFluidNewtons
                  == driven.patches[0]
                         .authoredPressureForceOnOpeningFluidNewtons
              && driven.authoredPressureImpulseOnOpeningFluidNewtonSeconds
                  == driven.patches[0]
                         .authoredPressureImpulseOnOpeningFluidNewtonSeconds
              && driven.authoredPressureWorkJoules > 0.0
              && driven.dissipatedEnergyJoules > 0.0
              && driven.kineticEnergyChangeJoules > 0.0
              && std::abs(driven.energyResidualJoules)
                  <= driven.energyToleranceJoules,
          "authored pressure jump accelerates aperture fluid with closed impulse and work");
    checkNear(
        driven.authoredPressureImpulseOnOpeningFluidNewtonSeconds.x,
        -0.7, 2.0e-16,
        "authored opening pressure impulse uses area times rise times duration");

    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        inviscidDrivenSamples{{100, 0.0}};
    auto inviscidDrivenFlux =
        buildPlanarPressureRegionFragmentOpeningFluxState(
            geometry, sweep, fragments, topology, definitions, openings,
            inviscidDrivenSamples);
    const auto inviscidDriven =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, definitions, openings,
            zeroResistance, inviscidDrivenSamples, inviscidDrivenFlux,
            drivenSettings);
    check(inviscidDriven.accepted
              && inviscidDriven.zeroResistancePatchCount == 1
              && !inviscidDriven.patches[0].zeroResistanceIdentity
              && inviscidDrivenSamples[0]
                     .relativeNormalVelocityMetersPerSecond
                  < 0.0
              && inviscidDriven.dissipatedEnergyJoules == 0.0
              && inviscidDriven.authoredPressureWorkJoules > 0.0
              && std::abs(inviscidDriven.energyResidualJoules)
                  <= inviscidDriven.energyToleranceJoules,
          "zero resistance remains an active inviscid authored-pressure drive");

    const auto oppositeWall = std::ranges::find_if(
        topology.links,
        [](const auto& link) {
            return link.kind
                    == PlanarPressureRegionFragmentFaceKind::PressureLayerWall
                && link.surfaceStableId == 20;
        });
    check(oppositeWall != topology.links.end(),
          "authored opening drive finds the opposite wall");
    if (oppositeWall != topology.links.end()) {
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition>
            oppositeDefinitions{{
                200, 2000, oppositeWall->surfaceStableId,
                oppositeWall->axis, oppositeWall->i, oppositeWall->j,
                oppositeWall->k, oppositeWall->minusRegionStableId,
                oppositeWall->plusRegionStableId, 0.5,
            }};
        const auto oppositeOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, sweep, fragments, topology,
                oppositeDefinitions);
        const std::vector<
            PlanarPressureRegionFragmentOpeningResistanceDefinition>
            oppositeResistance{{200, {8.0, 3.0}}};
        std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
            oppositeSamples{{200, 0.0}};
        auto oppositeFlux =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology,
                oppositeDefinitions, oppositeOpenings, oppositeSamples);
        const auto opposite =
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology,
                oppositeDefinitions, oppositeOpenings,
                oppositeResistance, oppositeSamples, oppositeFlux,
                drivenSettings);
        check(opposite.accepted
                  && opposite.patches[0].authoredPressureJumpPascals
                      == -70.0
                  && opposite.patches[0].drivingPressureRisePascals
                      == 70.0
                  && oppositeSamples[0]
                         .relativeNormalVelocityMetersPerSecond
                      > 0.0,
              "reversed authored wall jump reverses aperture acceleration");
        checkNear(
            oppositeSamples[0].relativeNormalVelocityMetersPerSecond,
            -drivenSamples[0].relativeNormalVelocityMetersPerSecond,
            2.0e-15,
            "opposite authored opening drive is sign symmetric");
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
            geometry, axisLayers, axisLayers, settings.timeStepSeconds);
        const auto axisFragments = buildPlanarPressureRegionFragments(
            geometry, axisSweep);
        const auto axisTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, axisSweep, axisFragments);
        const auto axisWall = std::ranges::find_if(
            axisTopology.links,
            [](const auto& link) {
                return link.kind
                        == PlanarPressureRegionFragmentFaceKind::
                            PressureLayerWall
                    && link.surfaceStableId == 10;
            });
        check(axisWall != axisTopology.links.end(),
              "authored opening drive finds its wall on every axis");
        if (axisWall == axisTopology.links.end()) continue;
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition>
            axisDefinitions{{
                100, 1000, axisWall->surfaceStableId, axisWall->axis,
                axisWall->i, axisWall->j, axisWall->k,
                axisWall->minusRegionStableId,
                axisWall->plusRegionStableId,
                0.5 * axisWall->areaSquareMeters,
            }};
        const auto axisOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, axisSweep, axisFragments, axisTopology,
                axisDefinitions);
        std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
            axisSamples{{100, 0.0}};
        auto axisFlux =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, axisSweep, axisFragments, axisTopology,
                axisDefinitions, axisOpenings, axisSamples);
        const auto axisDriven =
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, axisSweep, axisFragments, axisTopology,
                axisDefinitions, axisOpenings, resistance, axisSamples,
                axisFlux, drivenSettings);
        const Vector3 expectedForce = axis == GridFaceAxis::X
            ? Vector3{-35.0, 0.0, 0.0}
            : axis == GridFaceAxis::Y
            ? Vector3{0.0, -35.0, 0.0}
            : Vector3{0.0, 0.0, -35.0};
        check(axisDriven.accepted
                  && axisDriven.authoredPressureForceOnOpeningFluidNewtons
                      == expectedForce
                  && axisDriven.patches[0]
                         .authoredPressureForceOnOpeningFluidNewtons
                      == expectedForce
                  && axisSamples[0]
                         .relativeNormalVelocityMetersPerSecond
                      < 0.0
                  && std::abs(axisDriven.energyResidualJoules)
                      <= axisDriven.energyToleranceJoules,
              "authored aperture force and work close on every axis");
    }

    const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        parallelDefinitions{
            {102, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.2},
            {101, 1000, 10, GridFaceAxis::X, 1, 0, 0, 1, 2, 0.3},
        };
    const auto parallelOpenings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, parallelDefinitions);
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        parallelSamples{{102, 3.2}, {101, 3.2}};
    auto parallelFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, parallelDefinitions,
        parallelOpenings, parallelSamples);
    std::vector<PlanarPressureRegionFragmentOpeningResistanceDefinition>
        parallelResistance{{102, {8.0, 3.0}}, {101, {8.0, 3.0}}};
    auto reversedResistance = parallelResistance;
    std::ranges::reverse(reversedResistance);
    auto reorderedSamples = parallelSamples;
    auto reorderedFlux = parallelFlux;
    const auto parallel =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, parallelDefinitions,
            parallelOpenings, parallelResistance, parallelSamples,
            parallelFlux, settings);
    const auto reordered =
        advancePlanarPressureRegionFragmentOpeningResistance(
            geometry, sweep, fragments, topology, parallelDefinitions,
            parallelOpenings, reversedResistance, reorderedSamples,
            reorderedFlux, settings);
    check(parallel == reordered && parallelSamples == reorderedSamples
              && parallelFlux == reorderedFlux
              && parallel.patches.size() == 2,
          "opening resistance canonicalizes multi-patch coefficient order");
    checkNear(parallelSamples[0].relativeNormalVelocityMetersPerSecond,
              parallelSamples[1].relativeNormalVelocityMetersPerSecond,
              2.0e-15,
              "parallel aperture patches retain one resistance decay speed");
    checkNear(
        parallelFlux.patches[0]
                .relativeVolumeFlowRateCubicMetersPerSecond
            / parallelFlux.patches[0].areaSquareMeters,
        parallelFlux.patches[1]
                .relativeVolumeFlowRateCubicMetersPerSecond
            / parallelFlux.patches[1].areaSquareMeters,
        2.0e-15,
        "parallel aperture resistance remains area-consistent");

    const auto originalRejectedSamples = originalIdentitySamples;
    const auto originalRejectedFlux = originalIdentityFlux;
    expectRejected([&] {
        auto candidateSamples = originalRejectedSamples;
        auto candidateFlux = originalRejectedFlux;
        static_cast<void>(
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology, definitions, openings,
                {}, candidateSamples, candidateFlux, settings));
    }, "opening resistance rejects a missing patch coefficient");
    auto invalidResistance = resistance;
    invalidResistance[0].resistance.linearPascalSecondsPerMeter = -1.0;
    expectRejected([&] {
        auto candidateSamples = originalRejectedSamples;
        auto candidateFlux = originalRejectedFlux;
        static_cast<void>(
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology, definitions, openings,
                invalidResistance, candidateSamples, candidateFlux,
                settings));
    }, "opening resistance rejects a negative coefficient");
    std::vector<PlanarPressureRegionFragmentOpeningResistanceDefinition>
        foreignResistance{{999, {8.0, 3.0}}};
    expectRejected([&] {
        auto candidateSamples = originalRejectedSamples;
        auto candidateFlux = originalRejectedFlux;
        static_cast<void>(
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology, definitions, openings,
                foreignResistance, candidateSamples, candidateFlux,
                settings));
    }, "opening resistance rejects a foreign patch identity");
    const std::vector<
        PlanarPressureRegionFragmentOpeningResistanceDefinition>
        duplicateResistance{{101, {8.0, 3.0}}, {101, {8.0, 3.0}}};
    expectRejected([&] {
        auto candidateSamples = std::vector<
            PlanarPressureRegionFragmentOpeningVelocitySample>{
            {101, 3.2}, {102, 3.2}};
        auto candidateFlux =
            buildPlanarPressureRegionFragmentOpeningFluxState(
                geometry, sweep, fragments, topology, parallelDefinitions,
                parallelOpenings, candidateSamples);
        static_cast<void>(
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology, parallelDefinitions,
                parallelOpenings, duplicateResistance, candidateSamples,
                candidateFlux, settings));
    }, "opening resistance rejects a duplicate patch identity");
    auto corruptFlux = originalRejectedFlux;
    corruptFlux.patches[0].relativeVolumeFlowRateCubicMetersPerSecond += 0.1;
    expectRejected([&] {
        auto candidateSamples = originalRejectedSamples;
        static_cast<void>(
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology, definitions, openings,
                resistance, candidateSamples, corruptFlux, settings));
    }, "opening resistance rejects a corrupted flux source");
    auto limits = PlanarPressureRegionFragmentOpeningResistanceLimits{};
    limits.maximumWorkingBytes = 1;
    expectRejected([&] {
        auto candidateSamples = originalRejectedSamples;
        auto candidateFlux = originalRejectedFlux;
        static_cast<void>(
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology, definitions, openings,
                resistance, candidateSamples, candidateFlux, settings,
                limits));
    }, "opening resistance enforces its working-storage limit");
    auto invalidSettings = settings;
    invalidSettings.timeStepSeconds = 0.0;
    auto retainedSamples = originalRejectedSamples;
    auto retainedFlux = originalRejectedFlux;
    expectRejected([&] {
        static_cast<void>(
            advancePlanarPressureRegionFragmentOpeningResistance(
                geometry, sweep, fragments, topology, definitions, openings,
                resistance, retainedSamples, retainedFlux,
                invalidSettings));
    }, "opening resistance rejects a zero time step");
    check(retainedSamples == originalRejectedSamples
              && retainedFlux == originalRejectedFlux,
          "opening resistance rejection preserves samples and flux");
}

void testPlanarRegionalOpeningMomentumTransport() {
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
        current[0].physicalPlaneCoordinateMeters += 0.05;
        current[1].physicalPlaneCoordinateMeters += 0.05;
        const auto sourceSweep = makePlanarPressureRegionSweepLedger(
            geometry, previous, previous, 0.5);
        const auto sourceFragments = buildPlanarPressureRegionFragments(
            geometry, sourceSweep);
        const auto sourceTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, sourceSweep, sourceFragments);
        const auto targetSweep = makePlanarPressureRegionSweepLedger(
            geometry, previous, current, 0.5);
        const auto targetFragments = buildPlanarPressureRegionFragments(
            geometry, targetSweep);
        const auto targetTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, targetSweep, targetFragments);
        const auto targetVolumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, targetSweep, targetFragments, targetTopology);
        const auto targetWall = std::ranges::find_if(
            targetTopology.links,
            [](const auto& link) {
                return link.kind
                        == PlanarPressureRegionFragmentFaceKind::
                            PressureLayerWall
                    && link.surfaceStableId == 10;
            });
        check(targetWall != targetTopology.links.end(),
              "opening momentum transport finds its moving wall");
        if (targetWall == targetTopology.links.end()) continue;
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition>
            definitions{{
                100, 1000, targetWall->surfaceStableId, targetWall->axis,
                targetWall->i, targetWall->j, targetWall->k,
                targetWall->minusRegionStableId,
                targetWall->plusRegionStableId,
                0.5 * targetWall->areaSquareMeters,
            }};
        const auto sourceOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, sourceSweep, sourceFragments, sourceTopology,
                definitions);
        const auto targetOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, targetSweep, targetFragments, targetTopology,
                definitions);
        const auto sourceBaseMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sourceSweep, sourceFragments, sourceTopology);
        const auto targetBaseMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, targetSweep, targetFragments, targetTopology);
        const auto sourceMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, sourceSweep, sourceFragments, sourceTopology,
                sourceBaseMetric, definitions, sourceOpenings);
        const auto targetMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, targetSweep, targetFragments, targetTopology,
                targetBaseMetric, definitions, targetOpenings);

        Vector3 uniformVelocity{0.03, -0.04, 0.05};
        if (axis == GridFaceAxis::X) uniformVelocity.x = 0.1;
        if (axis == GridFaceAxis::Y) uniformVelocity.y = 0.1;
        if (axis == GridFaceAxis::Z) uniformVelocity.z = 0.1;
        const auto makeUniformState = [&](const auto& metric) {
            std::vector<double> normal(metric.dofs.size(), 0.0);
            std::vector<double> material(metric.dofs.size(), 0.0);
            std::vector<double> relative(metric.dofs.size(), 0.0);
            for (const auto& dof : metric.dofs) {
                const double component = dof.axis == GridFaceAxis::X
                    ? uniformVelocity.x
                    : dof.axis == GridFaceAxis::Y
                    ? uniformVelocity.y : uniformVelocity.z;
                normal[dof.dofIndex] = component;
                if (dof.kind
                    == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                        SharedRegionGrid) {
                    relative[dof.dofIndex] = component;
                } else {
                    material[dof.dofIndex] = component;
                }
            }
            return buildPlanarPressureRegionFragmentOpeningVelocityState(
                metric, normal, material, relative, 1.2);
        };
        const auto sourceState = makeUniformState(sourceMetric);
        const auto targetFlowState = makeUniformState(targetMetric);
        const auto transport =
            advancePlanarPressureRegionFragmentOpeningMomentum(
                sourceState, sourceMetric, targetFlowState, targetMetric,
                geometry, targetSweep, targetFragments, targetTopology,
                targetVolumeRates);
        const auto repeatedTransport =
            advancePlanarPressureRegionFragmentOpeningMomentum(
                sourceState, sourceMetric, targetFlowState, targetMetric,
                geometry, targetSweep, targetFragments, targetTopology,
                targetVolumeRates);
        check(transport == repeatedTransport && transport.diagnostics.accepted
                  && transport.diagnostics.finite
                  && transport.sourceStateFingerprint
                      == sourceState.fingerprint
                  && transport.sourceTransportFingerprint == 0
                  && transport.diagnostics.failureStage
                      == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                          None
                  && transport.controls.size()
                      == targetFragments.fragments.size()
                  && transport.diagnostics.openingDofCount == 1
                  && transport.diagnostics
                         .maximumContinuityResidualCubicMetersPerSecond
                      <= transport.diagnostics
                             .continuityToleranceCubicMetersPerSecond
                  && transport.diagnostics
                         .momentumResidualNormKilogramMetersPerSecond
                      < 2.0e-13,
              "opening momentum transport deterministically accepts moving uniform flow on every axis");
        for (const auto& control : transport.controls) {
            checkNear(control.velocityMetersPerSecond.x,
                      uniformVelocity.x, 2.0e-14,
                      "opening momentum transport preserves uniform X velocity");
            checkNear(control.velocityMetersPerSecond.y,
                      uniformVelocity.y, 2.0e-14,
                      "opening momentum transport preserves uniform Y velocity");
            checkNear(control.velocityMetersPerSecond.z,
                      uniformVelocity.z, 2.0e-14,
                      "opening momentum transport preserves uniform Z velocity");
        }
        checkNear(transport.diagnostics.maximumVelocityChangeMetersPerSecond,
                  0.0, 3.0e-14,
                  "opening momentum transport satisfies the moving free-stream GCL");
        checkNear(transport.diagnostics.advectiveKineticEnergyLossJoules,
                  0.0, 4.0e-13,
                  "uniform opening momentum transport has zero mixing loss");
        validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
            transport);
        validatePlanarPressureRegionFragmentOpeningMomentumTransport(
            transport, sourceState, sourceMetric, targetFlowState,
            targetMetric, geometry, targetSweep, targetFragments,
            targetTopology, targetVolumeRates);

        auto predictionCurrent = current;
        predictionCurrent[0].physicalPlaneCoordinateMeters += 0.05;
        predictionCurrent[1].physicalPlaneCoordinateMeters += 0.05;
        const auto predictionSweep = makePlanarPressureRegionSweepLedger(
            geometry, current, predictionCurrent, 0.5);
        const auto predictionFragments =
            buildPlanarPressureRegionFragments(geometry, predictionSweep);
        const auto predictionTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, predictionSweep, predictionFragments);
        const auto predictionVolumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, predictionSweep, predictionFragments,
                predictionTopology);
        const auto predictionOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, predictionSweep, predictionFragments,
                predictionTopology, definitions);
        const auto predictionBaseMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, predictionSweep, predictionFragments,
                predictionTopology);
        const auto predictionMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, predictionSweep, predictionFragments,
                predictionTopology, predictionBaseMetric, definitions,
                predictionOpenings);
        const auto prediction =
            predictPlanarPressureRegionFragmentOpeningMomentum(
                transport, targetMetric, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                predictionBaseMetric, predictionMetric);
        const auto repeatedPrediction =
            predictPlanarPressureRegionFragmentOpeningMomentum(
                transport, targetMetric, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                predictionBaseMetric, predictionMetric);
        check(prediction == repeatedPrediction
                  && prediction.diagnostics.finite
                  && prediction.diagnostics.fragmentCount
                      == predictionFragments.fragments.size()
                  && prediction.diagnostics.dofCount
                      == predictionMetric.dofs.size()
                  && prediction.diagnostics.openingDofCount == 1
                  && prediction.diagnostics
                         .maximumAbsoluteVolumeChangeCubicMeters
                      > 0.0
                  && prediction.diagnostics
                         .maximumEndpointNormalVelocityJumpMetersPerSecond
                      < 5.0e-14
                  && prediction.diagnostics
                         .maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond
                      < 5.0e-14
                  && prediction.predictedVelocityState
                         .staggeringKineticEnergyJoules
                      < 5.0e-13,
              "transported momentum predicts the consecutive moving uniform state on every axis");
        for (const auto& sample
             : prediction.predictedVelocityState.samples) {
            const double expected = sample.axis == GridFaceAxis::X
                ? uniformVelocity.x
                : sample.axis == GridFaceAxis::Y
                ? uniformVelocity.y : uniformVelocity.z;
            checkNear(sample.normalVelocityMetersPerSecond, expected,
                      3.0e-14,
                      "momentum prediction preserves uniform absolute velocity");
            if (sample.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SharedRegionGrid) {
                check(sample.materialNormalVelocityMetersPerSecond == 0.0,
                      "momentum prediction keeps fixed-grid material velocity zero");
            } else if (sample.kind
                       != PlanarPressureRegionFragmentOpeningVelocityDofKind::
                           OpeningPatch) {
                check(sample.relativeNormalVelocityMetersPerSecond == 0.0,
                      "momentum prediction keeps solid relative velocity zero");
            }
        }
        checkNear(
            prediction.diagnostics
                .geometricMomentumChangeKilogramMetersPerSecond.x,
            0.0, 3.0e-13,
            "uniform momentum prediction has no global geometric X impulse");
        checkNear(
            prediction.diagnostics
                .geometricMomentumChangeKilogramMetersPerSecond.y,
            0.0, 3.0e-13,
            "uniform momentum prediction has no global geometric Y impulse");
        checkNear(
            prediction.diagnostics
                .geometricMomentumChangeKilogramMetersPerSecond.z,
            0.0, 3.0e-13,
            "uniform momentum prediction has no global geometric Z impulse");
        checkNear(
            prediction.diagnostics
                .reconstructionMomentumChangeKilogramMetersPerSecond.x,
            0.0, 3.0e-13,
            "uniform momentum prediction has no reconstruction X impulse");
        checkNear(
            prediction.diagnostics
                .reconstructionMomentumChangeKilogramMetersPerSecond.y,
            0.0, 3.0e-13,
            "uniform momentum prediction has no reconstruction Y impulse");
        checkNear(
            prediction.diagnostics
                .reconstructionMomentumChangeKilogramMetersPerSecond.z,
            0.0, 3.0e-13,
            "uniform momentum prediction has no reconstruction Z impulse");
        validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
            prediction, transport, targetMetric, geometry,
            predictionSweep, predictionFragments, predictionTopology,
            predictionVolumeRates, definitions, predictionOpenings,
            predictionBaseMetric, predictionMetric);

        const auto predictionBasePressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, predictionSweep, predictionFragments,
                predictionTopology);
        const auto predictionPressureOperator =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology, definitions,
                predictionOpenings);
        const std::vector<
            PlanarPressureRegionFragmentOpeningResistanceDefinition>
            zeroResistance{{100, {0.0, 0.0}}};
        PlanarPressureRegionFragmentOpeningPressureStepSettings
            momentumPressureSettings;
        momentumPressureSettings.projection.densityKgPerCubicMeter = 1.2;
        momentumPressureSettings.projection.timeStepSeconds = 0.5;
        momentumPressureSettings.projection.pressureSolve
            .absoluteResidualTolerancePascalsMeters = 1.0e-13;
        momentumPressureSettings.projection.pressureSolve
            .relativeResidualTolerance = 1.0e-12;
        momentumPressureSettings.projection.pressureSolve
            .maximumIterations = 300;
        const auto momentumPressureEpoch =
            acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                prediction, predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, predictionBaseMetric, predictionMetric,
                momentumPressureSettings);
        const auto repeatedMomentumPressureEpoch =
            acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                prediction, predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, predictionBaseMetric, predictionMetric,
                momentumPressureSettings);
        check(momentumPressureEpoch == repeatedMomentumPressureEpoch
                  && momentumPressureEpoch.diagnostics.accepted
                  && momentumPressureEpoch.diagnostics
                         .usedTransportedPrediction
                  && momentumPressureEpoch.diagnostics
                         .usedColdPressureStart
                  && !momentumPressureEpoch.diagnostics
                          .usedWarmPressureStart
                  && momentumPressureEpoch.diagnostics.failureStage
                      == PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                          None
                  && momentumPressureEpoch.acceptedState.accepted
                  && momentumPressureEpoch.acceptedState
                         .sourcePressureOperatorFingerprint
                      == predictionPressureOperator.fingerprint
                  && momentumPressureEpoch.acceptedState
                         .sourceOpeningFluxFingerprint
                      == momentumPressureEpoch
                             .sourcePredictedOpeningFluxFingerprint,
              "transported uniform predictor reaches a deterministic accepted pressure endpoint on every axis");
        checkNear(
            momentumPressureEpoch.diagnostics.pressureStep
                .kineticEnergyChangeJoules,
            0.0, 8.0e-13,
            "uniform transported pressure epoch has no spurious kinetic-energy change");
        validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
            momentumPressureEpoch, prediction, predictionPressureOperator,
            predictionBasePressureOperator, geometry, predictionSweep,
            predictionFragments, predictionTopology,
            predictionVolumeRates, definitions, predictionOpenings,
            zeroResistance, predictionBaseMetric, predictionMetric,
            momentumPressureSettings);

        const auto acceptedPressureFlow =
            capturePlanarPressureRegionFragmentOpeningVelocityState(
                momentumPressureEpoch.acceptedState,
                predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, predictionBaseMetric, predictionMetric);
        const auto secondTransport =
            advancePlanarPressureRegionFragmentOpeningMomentum(
                transport, targetMetric, acceptedPressureFlow,
                predictionMetric, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates);
        const auto repeatedSecondTransport =
            advancePlanarPressureRegionFragmentOpeningMomentum(
                transport, targetMetric, acceptedPressureFlow,
                predictionMetric, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates);
        check(secondTransport == repeatedSecondTransport
                  && secondTransport.diagnostics.accepted
                  && secondTransport.sourceStateFingerprint == 0
                  && secondTransport.sourceTransportFingerprint
                      == transport.fingerprint
                  && secondTransport.sourceMetricFingerprint
                      == targetMetric.fingerprint
                  && secondTransport.diagnostics
                         .momentumBeforeKilogramMetersPerSecond
                      == transport.diagnostics
                             .momentumAfterKilogramMetersPerSecond
                  && secondTransport.diagnostics
                         .momentumResidualNormKilogramMetersPerSecond
                      < 2.0e-13,
              "accepted pressure flow drives a deterministic second momentum transport on every axis");
        for (const auto& control : secondTransport.controls) {
            checkNear(control.velocityMetersPerSecond.x,
                      uniformVelocity.x, 3.0e-14,
                      "re-entrant momentum transport preserves uniform X velocity");
            checkNear(control.velocityMetersPerSecond.y,
                      uniformVelocity.y, 3.0e-14,
                      "re-entrant momentum transport preserves uniform Y velocity");
            checkNear(control.velocityMetersPerSecond.z,
                      uniformVelocity.z, 3.0e-14,
                      "re-entrant momentum transport preserves uniform Z velocity");
        }
        checkNear(
            secondTransport.diagnostics.advectiveKineticEnergyLossJoules,
            0.0, 6.0e-13,
            "uniform re-entrant transport has zero mixing loss");
        validatePlanarPressureRegionFragmentOpeningMomentumTransport(
            secondTransport, transport, targetMetric,
            acceptedPressureFlow, predictionMetric, geometry,
            predictionSweep, predictionFragments, predictionTopology,
            predictionVolumeRates);

        auto warmCurrent = predictionCurrent;
        warmCurrent[0].physicalPlaneCoordinateMeters += 0.05;
        warmCurrent[1].physicalPlaneCoordinateMeters += 0.05;
        const auto warmSweep = makePlanarPressureRegionSweepLedger(
            geometry, predictionCurrent, warmCurrent, 0.5);
        const auto warmFragments = buildPlanarPressureRegionFragments(
            geometry, warmSweep);
        const auto warmTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, warmSweep, warmFragments);
        const auto warmVolumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, warmSweep, warmFragments, warmTopology);
        const auto warmOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, warmSweep, warmFragments, warmTopology,
                definitions);
        const auto warmBaseMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, warmSweep, warmFragments, warmTopology);
        const auto warmMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, warmSweep, warmFragments, warmTopology,
                warmBaseMetric, definitions, warmOpenings);
        const auto warmPrediction =
            predictPlanarPressureRegionFragmentOpeningMomentum(
                secondTransport, predictionMetric, geometry, warmSweep,
                warmFragments, warmTopology, warmVolumeRates,
                definitions, warmOpenings, warmBaseMetric, warmMetric);
        const auto warmBasePressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, warmSweep, warmFragments, warmTopology);
        const auto warmPressureOperator =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                warmBasePressureOperator, geometry, warmSweep,
                warmFragments, warmTopology, definitions, warmOpenings);
        const auto pressureWarmStart =
            buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
                momentumPressureEpoch.acceptedState,
                predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, warmPressureOperator,
                warmBasePressureOperator, warmSweep, warmFragments,
                warmTopology, warmVolumeRates, definitions, warmOpenings);
        const auto repeatedPressureWarmStart =
            buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
                momentumPressureEpoch.acceptedState,
                predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, warmPressureOperator,
                warmBasePressureOperator, warmSweep, warmFragments,
                warmTopology, warmVolumeRates, definitions, warmOpenings);
        check(pressureWarmStart == repeatedPressureWarmStart
                  && pressureWarmStart.sourceAcceptedStateFingerprint
                      == momentumPressureEpoch.acceptedState.fingerprint
                  && pressureWarmStart.currentPressureOperatorFingerprint
                      == warmPressureOperator.fingerprint
                  && pressureWarmStart.pressureCorrectionPascals.size()
                      == warmFragments.fragments.size(),
              "transported pressure history maps deterministically without stale velocity state");
        validatePlanarPressureRegionFragmentOpeningPressureWarmStart(
            pressureWarmStart, momentumPressureEpoch.acceptedState,
            predictionPressureOperator,
            predictionBasePressureOperator, geometry, predictionSweep,
            predictionFragments, predictionTopology,
            predictionVolumeRates, definitions, predictionOpenings,
            zeroResistance, warmPressureOperator,
            warmBasePressureOperator, warmSweep, warmFragments,
            warmTopology, warmVolumeRates, definitions, warmOpenings);
        const auto warmMomentumPressureEpoch =
            acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                warmPrediction, pressureWarmStart, warmPressureOperator,
                warmBasePressureOperator, geometry, warmSweep,
                warmFragments, warmTopology, warmVolumeRates, definitions,
                warmOpenings, zeroResistance, warmBaseMetric, warmMetric,
                momentumPressureSettings);
        const auto repeatedWarmMomentumPressureEpoch =
            acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                warmPrediction, pressureWarmStart, warmPressureOperator,
                warmBasePressureOperator, geometry, warmSweep,
                warmFragments, warmTopology, warmVolumeRates, definitions,
                warmOpenings, zeroResistance, warmBaseMetric, warmMetric,
                momentumPressureSettings);
        check(warmMomentumPressureEpoch
                      == repeatedWarmMomentumPressureEpoch
                  && warmMomentumPressureEpoch.diagnostics.accepted
                  && warmMomentumPressureEpoch.diagnostics
                         .usedTransportedPrediction
                  && warmMomentumPressureEpoch.diagnostics
                         .usedWarmPressureStart
                  && !warmMomentumPressureEpoch.diagnostics
                          .usedColdPressureStart
                  && warmMomentumPressureEpoch
                         .sourcePressureWarmStartFingerprint
                      == pressureWarmStart.fingerprint,
              "re-entrant uniform momentum reaches a deterministic warm pressure epoch on every axis");
        validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
            warmMomentumPressureEpoch, warmPrediction, pressureWarmStart,
            warmPressureOperator, warmBasePressureOperator, geometry,
            warmSweep, warmFragments, warmTopology, warmVolumeRates,
            definitions, warmOpenings, zeroResistance, warmBaseMetric,
            warmMetric, momentumPressureSettings);

        const auto momentumCycle =
            advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                transport, targetMetric,
                momentumPressureEpoch.acceptedState,
                predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, predictionBaseMetric, predictionMetric,
                warmPressureOperator, warmBasePressureOperator, warmSweep,
                warmFragments, warmTopology, warmVolumeRates, definitions,
                warmOpenings, zeroResistance, warmBaseMetric, warmMetric,
                {}, momentumPressureSettings);
        const auto repeatedMomentumCycle =
            advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                transport, targetMetric,
                momentumPressureEpoch.acceptedState,
                predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, predictionBaseMetric, predictionMetric,
                warmPressureOperator, warmBasePressureOperator, warmSweep,
                warmFragments, warmTopology, warmVolumeRates, definitions,
                warmOpenings, zeroResistance, warmBaseMetric, warmMetric,
                {}, momentumPressureSettings);
        check(momentumCycle == repeatedMomentumCycle
                  && momentumCycle.diagnostics.accepted
                  && momentumCycle.diagnostics.usedReentrantTransport
                  && momentumCycle.diagnostics.usedPressureWarmStart
                  && momentumCycle.diagnostics.failureStage
                      == PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                          None
                  && momentumCycle.transport == secondTransport
                  && momentumCycle.acceptedState
                      == warmMomentumPressureEpoch.acceptedState,
              "atomic transported cycle reproduces both explicit uniform endpoints on every axis");
        validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
            momentumCycle, transport, targetMetric,
            momentumPressureEpoch.acceptedState,
            predictionPressureOperator,
            predictionBasePressureOperator, geometry, predictionSweep,
            predictionFragments, predictionTopology,
            predictionVolumeRates, definitions, predictionOpenings,
            zeroResistance, predictionBaseMetric, predictionMetric,
            warmPressureOperator, warmBasePressureOperator, warmSweep,
            warmFragments, warmTopology, warmVolumeRates, definitions,
            warmOpenings, zeroResistance, warmBaseMetric, warmMetric);
        const auto bootstrapMomentumCycle =
            advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                targetFlowState, targetMetric,
                momentumPressureEpoch.acceptedState,
                predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, predictionBaseMetric, predictionMetric,
                warmPressureOperator, warmBasePressureOperator, warmSweep,
                warmFragments, warmTopology, warmVolumeRates, definitions,
                warmOpenings, zeroResistance, warmBaseMetric, warmMetric,
                {}, momentumPressureSettings);
        check(bootstrapMomentumCycle.diagnostics.accepted
                  && bootstrapMomentumCycle.diagnostics
                         .usedInitialStateTransport
                  && !bootstrapMomentumCycle.diagnostics
                          .usedReentrantTransport,
              "initial velocity state bootstrap accepts on every axis");
        check(bootstrapMomentumCycle.sourceStateFingerprint
                      == targetFlowState.fingerprint
                  && bootstrapMomentumCycle.sourceTransportFingerprint == 0,
              "initial velocity state bootstrap retains exclusive source lineage");
        for (const auto& control : bootstrapMomentumCycle.transport.controls) {
            checkNear(control.velocityMetersPerSecond.x,
                      uniformVelocity.x, 4.0e-14,
                      "initial bootstrap preserves uniform X velocity");
            checkNear(control.velocityMetersPerSecond.y,
                      uniformVelocity.y, 4.0e-14,
                      "initial bootstrap preserves uniform Y velocity");
            checkNear(control.velocityMetersPerSecond.z,
                      uniformVelocity.z, 4.0e-14,
                      "initial bootstrap preserves uniform Z velocity");
        }
        check(bootstrapMomentumCycle.transport.diagnostics
                      .momentumResidualNormKilogramMetersPerSecond
                  < 2.0e-13
                  && bootstrapMomentumCycle.acceptedState.accepted
                  && bootstrapMomentumCycle.acceptedState
                         .correctedContinuityResidualMaximumCubicMetersPerSecond
                      <= bootstrapMomentumCycle.acceptedState
                             .continuityToleranceCubicMetersPerSecond,
              "initial velocity state bootstrap closes momentum and pressure");
        validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
            bootstrapMomentumCycle, targetFlowState, targetMetric,
            momentumPressureEpoch.acceptedState,
            predictionPressureOperator,
            predictionBasePressureOperator, geometry, predictionSweep,
            predictionFragments, predictionTopology,
            predictionVolumeRates, definitions, predictionOpenings,
            zeroResistance, predictionBaseMetric, predictionMetric,
            warmPressureOperator, warmBasePressureOperator, warmSweep,
            warmFragments, warmTopology, warmVolumeRates, definitions,
            warmOpenings, zeroResistance, warmBaseMetric, warmMetric);

        if (axis != GridFaceAxis::X) continue;
        const auto bootstrapCycleState =
            capturePlanarPressureRegionFragmentOpeningMomentumCycleState(
                bootstrapMomentumCycle, predictionMetric, warmMetric);
        validatePlanarPressureRegionFragmentOpeningMomentumCycleState(
            bootstrapCycleState, predictionVolumeRates,
            predictionMetric, warmPressureOperator,
            warmBasePressureOperator, geometry, warmSweep, warmFragments,
            warmTopology, warmVolumeRates, definitions, warmOpenings,
            zeroResistance, warmMetric);
        check(bootstrapCycleState.transport
                      == bootstrapMomentumCycle.transport
                  && bootstrapCycleState.acceptedState
                      == bootstrapMomentumCycle.acceptedState
                  && bootstrapCycleState.transportMetricFingerprint
                      == predictionMetric.fingerprint
                  && bootstrapCycleState.acceptedMetricFingerprint
                      == warmMetric.fingerprint
                  && bootstrapCycleState.predictedOpeningFluxFingerprint
                      == bootstrapCycleState.acceptedState
                             .sourceOpeningFluxFingerprint,
              "accepted cycle capture retains the exact staggered restart pair");

        PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError
            cyclePersistenceError;
        std::vector<std::uint8_t> cycleStateBytes;
        std::vector<std::uint8_t> repeatedCycleStateBytes;
        check(serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                  bootstrapCycleState, predictionVolumeRates,
                  predictionMetric, warmPressureOperator,
                  warmBasePressureOperator, geometry, warmSweep,
                  warmFragments, warmTopology, warmVolumeRates,
                  definitions, warmOpenings, zeroResistance, warmMetric,
                  cycleStateBytes, &cyclePersistenceError)
                  && cyclePersistenceError.code
                      == PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                          None
                  && serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                      bootstrapCycleState, predictionVolumeRates,
                      predictionMetric, warmPressureOperator,
                      warmBasePressureOperator, geometry, warmSweep,
                      warmFragments, warmTopology, warmVolumeRates,
                      definitions, warmOpenings, zeroResistance, warmMetric,
                      repeatedCycleStateBytes, &cyclePersistenceError)
                  && cycleStateBytes == repeatedCycleStateBytes,
              "opening momentum-cycle state persistence encodes deterministically");
        PlanarPressureRegionFragmentOpeningMomentumCycleState
            restoredBootstrapCycleState;
        const bool restoredBootstrapState =
            deserializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                cycleStateBytes, predictionVolumeRates, predictionMetric,
                warmPressureOperator, warmBasePressureOperator, geometry,
                warmSweep, warmFragments, warmTopology, warmVolumeRates,
                definitions, warmOpenings, zeroResistance, warmMetric,
                restoredBootstrapCycleState, &cyclePersistenceError);
        std::vector<std::uint8_t> restoredCycleStateBytes;
        check(restoredBootstrapState
                  && cyclePersistenceError.code
                      == PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                          None
                  && restoredBootstrapCycleState == bootstrapCycleState
                  && serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                      restoredBootstrapCycleState, predictionVolumeRates,
                      predictionMetric, warmPressureOperator,
                      warmBasePressureOperator, geometry, warmSweep,
                      warmFragments, warmTopology, warmVolumeRates,
                      definitions, warmOpenings, zeroResistance, warmMetric,
                      restoredCycleStateBytes, &cyclePersistenceError)
                  && restoredCycleStateBytes == cycleStateBytes,
              "opening momentum-cycle restart round trips bit-exactly");

        const auto expectCycleStatePersistenceRejected =
            [&](const std::vector<std::uint8_t>& candidateBytes,
                const PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode
                    expectedCode,
                const PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceLimits&
                    persistenceLimits = {}) {
                auto retainedState = bootstrapCycleState;
                const auto before = retainedState;
                PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceError
                    rejectedError;
                const bool decoded =
                    deserializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                        candidateBytes, predictionVolumeRates,
                        predictionMetric, warmPressureOperator,
                        warmBasePressureOperator, geometry, warmSweep,
                        warmFragments, warmTopology, warmVolumeRates,
                        definitions, warmOpenings, zeroResistance,
                        warmMetric, retainedState, &rejectedError,
                        persistenceLimits);
                return !decoded && rejectedError.code == expectedCode
                    && retainedState == before;
            };
        auto corruptCycleStateBytes = cycleStateBytes;
        corruptCycleStateBytes[0] ^= 0x01U;
        check(expectCycleStatePersistenceRejected(
                  corruptCycleStateBytes,
                  PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                      InvalidMagic),
              "opening momentum-cycle persistence rejects foreign magic transactionally");
        corruptCycleStateBytes = cycleStateBytes;
        corruptCycleStateBytes[4] = 2;
        corruptCycleStateBytes[5] = 0;
        check(expectCycleStatePersistenceRejected(
                  corruptCycleStateBytes,
                  PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                      UnsupportedVersion),
              "opening momentum-cycle persistence rejects unsupported protocol transactionally");
        corruptCycleStateBytes = cycleStateBytes;
        corruptCycleStateBytes.back() ^= 0x01U;
        check(expectCycleStatePersistenceRejected(
                  corruptCycleStateBytes,
                  PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                      ChecksumMismatch),
              "opening momentum-cycle persistence detects payload corruption");
        corruptCycleStateBytes = cycleStateBytes;
        constexpr std::size_t encodedCycleStateFingerprintOffset = 40;
        corruptCycleStateBytes[encodedCycleStateFingerprintOffset] ^= 0x01U;
        refreshOpeningMomentumCycleStatePersistenceChecksum(
            corruptCycleStateBytes);
        check(expectCycleStatePersistenceRejected(
                  corruptCycleStateBytes,
                  PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                      InvalidData),
              "opening momentum-cycle persistence rejects recomputed-checksum state corruption");
        corruptCycleStateBytes = cycleStateBytes;
        corruptCycleStateBytes.pop_back();
        check(expectCycleStatePersistenceRejected(
                  corruptCycleStateBytes,
                  PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                      Truncated),
              "opening momentum-cycle persistence rejects truncation transactionally");
        corruptCycleStateBytes = cycleStateBytes;
        corruptCycleStateBytes.push_back(0);
        check(expectCycleStatePersistenceRejected(
                  corruptCycleStateBytes,
                  PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                      TrailingData),
              "opening momentum-cycle persistence rejects trailing bytes transactionally");
        auto cyclePersistenceLimits =
            PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceLimits{};
        cyclePersistenceLimits.maximumTransportControls =
            bootstrapCycleState.transport.controls.size() - 1;
        check(expectCycleStatePersistenceRejected(
                  cycleStateBytes,
                  PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                      LimitExceeded,
                  cyclePersistenceLimits),
              "opening momentum-cycle persistence enforces transport record limits");
        auto foreignCycleState = bootstrapCycleState;
        const auto retainedForeignCycleState = foreignCycleState;
        check(!deserializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                  cycleStateBytes, predictionVolumeRates,
                  predictionMetric, predictionPressureOperator,
                  predictionBasePressureOperator, geometry,
                  predictionSweep, predictionFragments,
                  predictionTopology, predictionVolumeRates, definitions,
                  predictionOpenings, zeroResistance, warmMetric,
                  foreignCycleState, &cyclePersistenceError)
                  && cyclePersistenceError.code
                      == PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                          SourceMismatch
                  && foreignCycleState == retainedForeignCycleState,
              "opening momentum-cycle persistence rejects foreign pressure sources transactionally");
        cyclePersistenceLimits = {};
        cyclePersistenceLimits.maximumEncodedBytes =
            cycleStateBytes.size() - 1;
        std::vector<std::uint8_t> retainedCycleEncoding{1, 2, 3};
        const auto retainedCycleEncodingBefore = retainedCycleEncoding;
        check(!serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                  bootstrapCycleState, predictionVolumeRates,
                  predictionMetric, warmPressureOperator,
                  warmBasePressureOperator, geometry, warmSweep,
                  warmFragments, warmTopology, warmVolumeRates,
                  definitions, warmOpenings, zeroResistance, warmMetric,
                  retainedCycleEncoding, &cyclePersistenceError,
                  cyclePersistenceLimits)
                  && cyclePersistenceError.code
                      == PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                          LimitExceeded
                  && retainedCycleEncoding == retainedCycleEncodingBefore,
              "opening momentum-cycle persistence enforces encode limits transactionally");
        auto corruptCycleState = bootstrapCycleState;
        corruptCycleState.transport.controls[0]
            .velocityMetersPerSecond.x += 0.1;
        retainedCycleEncoding = retainedCycleEncodingBefore;
        check(!serializePlanarPressureRegionFragmentOpeningMomentumCycleState(
                  corruptCycleState, predictionVolumeRates,
                  predictionMetric, warmPressureOperator,
                  warmBasePressureOperator, geometry, warmSweep,
                  warmFragments, warmTopology, warmVolumeRates,
                  definitions, warmOpenings, zeroResistance, warmMetric,
                  retainedCycleEncoding, &cyclePersistenceError)
                  && cyclePersistenceError.code
                      == PlanarPressureRegionFragmentOpeningMomentumCycleStatePersistenceErrorCode::
                          InvalidData
                  && retainedCycleEncoding == retainedCycleEncodingBefore,
              "opening momentum-cycle persistence rejects corrupt in-memory transport transactionally");

        auto fourthCurrent = warmCurrent;
        fourthCurrent[0].physicalPlaneCoordinateMeters += 0.02;
        fourthCurrent[1].physicalPlaneCoordinateMeters += 0.02;
        const auto fourthSweep = makePlanarPressureRegionSweepLedger(
            geometry, warmCurrent, fourthCurrent, 0.5);
        const auto fourthFragments = buildPlanarPressureRegionFragments(
            geometry, fourthSweep);
        const auto fourthTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, fourthSweep, fourthFragments);
        const auto fourthVolumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, fourthSweep, fourthFragments, fourthTopology);
        const auto fourthOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, fourthSweep, fourthFragments, fourthTopology,
                definitions);
        const auto fourthBaseMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, fourthSweep, fourthFragments, fourthTopology);
        const auto fourthMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, fourthSweep, fourthFragments, fourthTopology,
                fourthBaseMetric, definitions, fourthOpenings);
        const auto fourthBasePressureOperator =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, fourthSweep, fourthFragments, fourthTopology);
        const auto fourthPressureOperator =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                fourthBasePressureOperator, geometry, fourthSweep,
                fourthFragments, fourthTopology, definitions,
                fourthOpenings);
        const auto postBootstrapCycle =
            advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                restoredBootstrapCycleState.transport, predictionMetric,
                restoredBootstrapCycleState.acceptedState,
                warmPressureOperator, warmBasePressureOperator, geometry,
                warmSweep, warmFragments, warmTopology, warmVolumeRates,
                definitions, warmOpenings, zeroResistance, warmBaseMetric,
                warmMetric, fourthPressureOperator,
                fourthBasePressureOperator, fourthSweep, fourthFragments,
                fourthTopology, fourthVolumeRates, definitions,
                fourthOpenings, zeroResistance, fourthBaseMetric,
                fourthMetric, {}, momentumPressureSettings);
        check(postBootstrapCycle.diagnostics.accepted
                  && postBootstrapCycle.diagnostics
                         .usedReentrantTransport
                  && !postBootstrapCycle.diagnostics
                          .usedInitialStateTransport
                  && postBootstrapCycle.sourceTransportFingerprint
                      == restoredBootstrapCycleState.transport.fingerprint
                  && postBootstrapCycle.acceptedState.accepted,
              "bootstrap endpoint pair feeds the next re-entrant atomic cycle");
        for (const auto& control : postBootstrapCycle.transport.controls) {
            checkNear(control.velocityMetersPerSecond.x,
                      uniformVelocity.x, 5.0e-14,
                      "post-bootstrap cycle preserves uniform X velocity");
            checkNear(control.velocityMetersPerSecond.y,
                      uniformVelocity.y, 5.0e-14,
                      "post-bootstrap cycle preserves uniform Y velocity");
            checkNear(control.velocityMetersPerSecond.z,
                      uniformVelocity.z, 5.0e-14,
                      "post-bootstrap cycle preserves uniform Z velocity");
        }
        validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
            postBootstrapCycle, restoredBootstrapCycleState.transport,
            predictionMetric, restoredBootstrapCycleState.acceptedState,
            warmPressureOperator, warmBasePressureOperator, geometry,
            warmSweep, warmFragments, warmTopology, warmVolumeRates,
            definitions, warmOpenings, zeroResistance, warmBaseMetric,
            warmMetric, fourthPressureOperator,
            fourthBasePressureOperator, fourthSweep, fourthFragments,
            fourthTopology, fourthVolumeRates, definitions,
            fourthOpenings, zeroResistance, fourthBaseMetric,
            fourthMetric);
        std::vector<double> nonuniformNormal(
            sourceMetric.dofs.size(), 0.0);
        std::vector<double> nonuniformMaterial(
            sourceMetric.dofs.size(), 0.0);
        std::vector<double> nonuniformRelative(
            sourceMetric.dofs.size(), 0.0);
        for (const auto& dof : sourceMetric.dofs) {
            const double base = dof.axis == GridFaceAxis::X
                ? uniformVelocity.x
                : dof.axis == GridFaceAxis::Y
                ? uniformVelocity.y : uniformVelocity.z;
            const double value = base
                + 0.01 * static_cast<double>(
                    static_cast<int>(dof.dofIndex % 5) - 2);
            nonuniformNormal[dof.dofIndex] = value;
            if (dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SharedRegionGrid) {
                nonuniformRelative[dof.dofIndex] = value;
            } else {
                nonuniformMaterial[dof.dofIndex] = value;
            }
        }
        const auto nonuniformSource =
            buildPlanarPressureRegionFragmentOpeningVelocityState(
                sourceMetric, nonuniformNormal, nonuniformMaterial,
                nonuniformRelative, 1.2);
        const auto mixed =
            advancePlanarPressureRegionFragmentOpeningMomentum(
                nonuniformSource, sourceMetric, targetFlowState,
                targetMetric, geometry, targetSweep, targetFragments,
                targetTopology, targetVolumeRates);
        check(mixed.diagnostics.accepted
                  && mixed.diagnostics.finite
                  && mixed.diagnostics.advectiveKineticEnergyLossJoules
                      > 0.0
                  && mixed.diagnostics
                         .momentumResidualNormKilogramMetersPerSecond
                      <= mixed.settings
                             .absoluteMomentumToleranceKilogramMetersPerSecond,
              "opening momentum transport conserves nonuniform momentum with dissipative donor mixing");
        validatePlanarPressureRegionFragmentOpeningMomentumTransport(
            mixed, nonuniformSource, sourceMetric, targetFlowState,
            targetMetric, geometry, targetSweep, targetFragments,
            targetTopology, targetVolumeRates);
        const auto mixedPrediction =
            predictPlanarPressureRegionFragmentOpeningMomentum(
                mixed, targetMetric, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                predictionBaseMetric, predictionMetric);
        check(mixedPrediction.diagnostics.finite
                  && mixedPrediction.diagnostics
                         .maximumEndpointNormalVelocityJumpMetersPerSecond
                      > 0.0
                  && mixedPrediction.predictedVelocityState
                         .maximumCollocatedSpeedMetersPerSecond
                      > 0.0,
              "nonuniform transported momentum predicts a finite current face state");
        validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
            mixedPrediction, mixed, targetMetric, geometry,
            predictionSweep, predictionFragments, predictionTopology,
            predictionVolumeRates, definitions, predictionOpenings,
            predictionBaseMetric, predictionMetric);

        auto breathingCurrent = previous;
        breathingCurrent[0].physicalPlaneCoordinateMeters -= 0.1;
        breathingCurrent[1].physicalPlaneCoordinateMeters += 0.1;
        const auto breathingSweep = makePlanarPressureRegionSweepLedger(
            geometry, previous, breathingCurrent, 0.5);
        const auto breathingFragments =
            buildPlanarPressureRegionFragments(geometry, breathingSweep);
        const auto breathingTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, breathingSweep, breathingFragments);
        const auto breathingVolumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, breathingSweep, breathingFragments,
                breathingTopology);
        const auto breathingWall = std::ranges::find_if(
            breathingTopology.links,
            [](const auto& link) {
                return link.kind
                        == PlanarPressureRegionFragmentFaceKind::
                            PressureLayerWall
                    && link.surfaceStableId == 10;
            });
        check(breathingWall != breathingTopology.links.end(),
              "opening momentum transport finds its breathing aperture wall");
        if (breathingWall != breathingTopology.links.end()) {
            const auto breathingOpenings =
                buildPlanarPressureRegionFragmentOpenings(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology, definitions);
            const auto breathingBaseMetric =
                buildPlanarPressureRegionFragmentVelocityMetric(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology);
            const auto breathingMetric =
                buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology, breathingBaseMetric, definitions,
                    breathingOpenings);
            const auto breathingPressureOperator =
                buildPlanarPressureRegionFragmentPressureOperator(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology);
            const double requiredOpeningFlow =
                -breathingVolumeRates
                     .components[breathingWall->minusComponentIndex]
                     .geometryVolumeChangeRateCubicMetersPerSecond;
            const std::vector<
                PlanarPressureRegionFragmentOpeningVelocitySample>
                breathingSamples{{
                    100,
                    requiredOpeningFlow / definitions[0].areaSquareMeters,
                }};
            const auto breathingFlux =
                buildPlanarPressureRegionFragmentOpeningFluxState(
                    geometry, breathingSweep, breathingFragments,
                    breathingTopology, definitions, breathingOpenings,
                    breathingSamples);
            PlanarPressureRegionFragmentPressureProjectionSettings
                projectionSettings;
            projectionSettings.densityKgPerCubicMeter = 1.2;
            projectionSettings.timeStepSeconds = 0.5;
            projectionSettings.pressureSolve
                .absoluteResidualTolerancePascalsMeters = 1.0e-13;
            projectionSettings.pressureSolve.relativeResidualTolerance =
                1.0e-12;
            projectionSettings.pressureSolve.maximumIterations = 300;
            std::vector<double> breathingLinkVelocity(
                breathingTopology.links.size(), 0.0);
            std::vector<double> breathingPressure(
                breathingFragments.fragments.size(), 0.0);
            const auto breathingProjection =
                projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                    breathingPressureOperator, geometry, breathingSweep,
                    breathingFragments, breathingTopology,
                    breathingVolumeRates, definitions, breathingOpenings,
                    breathingFlux, breathingSamples, breathingLinkVelocity,
                    breathingPressure, projectionSettings);
            check(breathingProjection.accepted,
                  "opening momentum transport receives a corrected breathing flow");
            std::vector<double> breathingNormal(
                breathingMetric.dofs.size(), 0.0);
            std::vector<double> breathingMaterial(
                breathingMetric.dofs.size(), 0.0);
            std::vector<double> breathingRelative(
                breathingMetric.dofs.size(), 0.0);
            for (const auto& dof : breathingMetric.dofs) {
                const auto& link = breathingTopology.links.at(
                    dof.sourceFaceLinkIndex);
                if (dof.kind
                    == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                        SharedRegionGrid) {
                    breathingRelative[dof.dofIndex] =
                        breathingLinkVelocity[link.linkIndex];
                } else {
                    breathingMaterial[dof.dofIndex] = breathingVolumeRates
                        .fragments[link.minusFragmentIndex]
                        .upperBoundaryVelocityMetersPerSecond;
                    if (dof.kind
                        == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                            OpeningPatch) {
                        breathingRelative[dof.dofIndex] =
                            breathingSamples[dof.sourceOpeningPatchIndex]
                                .relativeNormalVelocityMetersPerSecond;
                    }
                }
                breathingNormal[dof.dofIndex] =
                    breathingMaterial[dof.dofIndex]
                    + breathingRelative[dof.dofIndex];
            }
            const auto breathingFlowState =
                buildPlanarPressureRegionFragmentOpeningVelocityState(
                    breathingMetric, breathingNormal, breathingMaterial,
                    breathingRelative, 1.2);
            const auto breathingTransport =
                advancePlanarPressureRegionFragmentOpeningMomentum(
                    nonuniformSource, sourceMetric, breathingFlowState,
                    breathingMetric, geometry, breathingSweep,
                    breathingFragments, breathingTopology,
                    breathingVolumeRates);
            check(breathingTransport.diagnostics.accepted
                      && breathingTransport.diagnostics.openingDofCount == 1
                      && breathingTransport.diagnostics
                             .maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond
                          > 0.0
                      && breathingTransport.diagnostics
                             .advectiveKineticEnergyLossJoules
                          > 0.0
                      && breathingTransport.diagnostics
                             .momentumResidualNormKilogramMetersPerSecond
                          <= breathingTransport.settings
                                 .absoluteMomentumToleranceKilogramMetersPerSecond,
                  "opening momentum transport carries full vector momentum through a live aperture");
            validatePlanarPressureRegionFragmentOpeningMomentumTransport(
                breathingTransport, nonuniformSource, sourceMetric,
                breathingFlowState, breathingMetric, geometry,
                breathingSweep, breathingFragments, breathingTopology,
                breathingVolumeRates);

            auto breathingPredictionCurrent = breathingCurrent;
            breathingPredictionCurrent[0]
                .physicalPlaneCoordinateMeters -= 0.05;
            breathingPredictionCurrent[1]
                .physicalPlaneCoordinateMeters += 0.05;
            const auto breathingPredictionSweep =
                makePlanarPressureRegionSweepLedger(
                    geometry, breathingCurrent,
                    breathingPredictionCurrent, 0.5);
            const auto breathingPredictionFragments =
                buildPlanarPressureRegionFragments(
                    geometry, breathingPredictionSweep);
            const auto breathingPredictionTopology =
                buildPlanarPressureRegionFragmentTopology(
                    geometry, breathingPredictionSweep,
                    breathingPredictionFragments);
            const auto breathingPredictionVolumeRates =
                buildPlanarPressureRegionFragmentVolumeRates(
                    geometry, breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology);
            const auto breathingPredictionOpenings =
                buildPlanarPressureRegionFragmentOpenings(
                    geometry, breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology, definitions);
            const auto breathingPredictionBaseMetric =
                buildPlanarPressureRegionFragmentVelocityMetric(
                    geometry, breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology);
            const auto breathingPredictionMetric =
                buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                    geometry, breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionBaseMetric, definitions,
                    breathingPredictionOpenings);
            const auto breathingPrediction =
                predictPlanarPressureRegionFragmentOpeningMomentum(
                    breathingTransport, breathingMetric, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings,
                    breathingPredictionBaseMetric,
                    breathingPredictionMetric);
            const auto predictedAperture = std::ranges::find(
                breathingPrediction.predictedVelocityState.samples,
                PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    OpeningPatch,
                &PlanarPressureRegionFragmentOpeningVelocityStateSample::
                    kind);
            check(breathingPrediction.diagnostics.finite
                      && breathingPrediction.diagnostics
                             .maximumAbsolutePredictedOpeningRelativeVelocityMetersPerSecond
                          > 0.0
                      && predictedAperture
                          != breathingPrediction
                                 .predictedVelocityState.samples.end()
                      && predictedAperture
                             ->relativeNormalVelocityMetersPerSecond
                          != 0.0,
                  "transported breathing momentum predicts a live material-relative aperture flow");
            validatePlanarPressureRegionFragmentOpeningMomentumPrediction(
                breathingPrediction, breathingTransport, breathingMetric,
                geometry, breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates, definitions,
                breathingPredictionOpenings,
                breathingPredictionBaseMetric,
                breathingPredictionMetric);
            const auto breathingPredictionBasePressureOperator =
                buildPlanarPressureRegionFragmentPressureOperator(
                    geometry, breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology);
            const auto breathingPredictionPressureOperator =
                buildPlanarPressureRegionFragmentOpeningPressureOperator(
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology, definitions,
                    breathingPredictionOpenings);
            const auto breathingMomentumPressureEpoch =
                acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                    breathingPrediction,
                    breathingPredictionPressureOperator,
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings, zeroResistance,
                    breathingPredictionBaseMetric,
                    breathingPredictionMetric, momentumPressureSettings);
            check(breathingMomentumPressureEpoch.diagnostics.accepted
                      && breathingMomentumPressureEpoch.diagnostics
                             .pressureStep.projection.pressureSolve
                             .iterationCount
                          > 0
                      && breathingMomentumPressureEpoch.acceptedState.accepted
                      && breathingMomentumPressureEpoch.acceptedState
                             .correctedContinuityResidualMaximumCubicMetersPerSecond
                          <= breathingMomentumPressureEpoch
                                 .acceptedState
                                 .continuityToleranceCubicMetersPerSecond,
                  "transported breathing predictor reaches a constrained pressure endpoint");
            validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
                breathingMomentumPressureEpoch, breathingPrediction,
                breathingPredictionPressureOperator,
                breathingPredictionBasePressureOperator, geometry,
                breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates, definitions,
                breathingPredictionOpenings, zeroResistance,
                breathingPredictionBaseMetric,
                breathingPredictionMetric, momentumPressureSettings);

            const auto breathingAcceptedPressureFlow =
                capturePlanarPressureRegionFragmentOpeningVelocityState(
                    breathingMomentumPressureEpoch.acceptedState,
                    breathingPredictionPressureOperator,
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings, zeroResistance,
                    breathingPredictionBaseMetric,
                    breathingPredictionMetric);
            const auto breathingSecondTransport =
                advancePlanarPressureRegionFragmentOpeningMomentum(
                    breathingTransport, breathingMetric,
                    breathingAcceptedPressureFlow,
                    breathingPredictionMetric, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates);
            check(breathingSecondTransport.diagnostics.accepted
                      && breathingSecondTransport
                             .sourceTransportFingerprint
                          == breathingTransport.fingerprint
                      && breathingSecondTransport.diagnostics
                             .maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond
                          > 0.0
                      && breathingSecondTransport.diagnostics
                             .kineticEnergyAfterJoules
                          <= breathingSecondTransport.diagnostics
                                 .kineticEnergyBeforeJoules
                              + breathingSecondTransport.settings
                                    .absoluteEnergyToleranceJoules
                      && breathingSecondTransport.diagnostics
                             .momentumResidualNormKilogramMetersPerSecond
                          <= breathingSecondTransport.settings
                                 .absoluteMomentumToleranceKilogramMetersPerSecond,
                  "accepted breathing pressure flow advances exact prior transported momentum");
            validatePlanarPressureRegionFragmentOpeningMomentumTransport(
                breathingSecondTransport, breathingTransport,
                breathingMetric, breathingAcceptedPressureFlow,
                breathingPredictionMetric, geometry,
                breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates);

            auto breathingWarmCurrent = breathingPredictionCurrent;
            breathingWarmCurrent[0]
                .physicalPlaneCoordinateMeters -= 0.02;
            breathingWarmCurrent[1]
                .physicalPlaneCoordinateMeters += 0.02;
            const auto breathingWarmSweep =
                makePlanarPressureRegionSweepLedger(
                    geometry, breathingPredictionCurrent,
                    breathingWarmCurrent, 0.5);
            const auto breathingWarmFragments =
                buildPlanarPressureRegionFragments(
                    geometry, breathingWarmSweep);
            const auto breathingWarmTopology =
                buildPlanarPressureRegionFragmentTopology(
                    geometry, breathingWarmSweep,
                    breathingWarmFragments);
            const auto breathingWarmVolumeRates =
                buildPlanarPressureRegionFragmentVolumeRates(
                    geometry, breathingWarmSweep,
                    breathingWarmFragments,
                    breathingWarmTopology);
            const auto breathingWarmOpenings =
                buildPlanarPressureRegionFragmentOpenings(
                    geometry, breathingWarmSweep,
                    breathingWarmFragments,
                    breathingWarmTopology, definitions);
            const auto breathingWarmBaseMetric =
                buildPlanarPressureRegionFragmentVelocityMetric(
                    geometry, breathingWarmSweep,
                    breathingWarmFragments,
                    breathingWarmTopology);
            const auto breathingWarmMetric =
                buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                    geometry, breathingWarmSweep,
                    breathingWarmFragments,
                    breathingWarmTopology, breathingWarmBaseMetric,
                    definitions, breathingWarmOpenings);
            const auto breathingWarmPrediction =
                predictPlanarPressureRegionFragmentOpeningMomentum(
                    breathingSecondTransport,
                    breathingPredictionMetric, geometry,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, breathingWarmVolumeRates,
                    definitions, breathingWarmOpenings,
                    breathingWarmBaseMetric, breathingWarmMetric);
            const auto breathingWarmBasePressureOperator =
                buildPlanarPressureRegionFragmentPressureOperator(
                    geometry, breathingWarmSweep,
                    breathingWarmFragments,
                    breathingWarmTopology);
            const auto breathingWarmPressureOperator =
                buildPlanarPressureRegionFragmentOpeningPressureOperator(
                    breathingWarmBasePressureOperator, geometry,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, definitions,
                    breathingWarmOpenings);
            const auto breathingPressureWarmStart =
                buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
                    breathingMomentumPressureEpoch.acceptedState,
                    breathingPredictionPressureOperator,
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings, zeroResistance,
                    breathingWarmPressureOperator,
                    breathingWarmBasePressureOperator,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, breathingWarmVolumeRates,
                    definitions, breathingWarmOpenings);
            check(breathingPressureWarmStart
                          .maximumAbsolutePressureCorrectionPascals
                      > 0.0
                      && breathingPressureWarmStart
                             .sourceAcceptedStateFingerprint
                          == breathingMomentumPressureEpoch
                                 .acceptedState.fingerprint,
                  "breathing pressure history retains a nonzero transported warm gauge");
            validatePlanarPressureRegionFragmentOpeningPressureWarmStart(
                breathingPressureWarmStart,
                breathingMomentumPressureEpoch.acceptedState,
                breathingPredictionPressureOperator,
                breathingPredictionBasePressureOperator, geometry,
                breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates, definitions,
                breathingPredictionOpenings, zeroResistance,
                breathingWarmPressureOperator,
                breathingWarmBasePressureOperator,
                breathingWarmSweep, breathingWarmFragments,
                breathingWarmTopology, breathingWarmVolumeRates,
                definitions, breathingWarmOpenings);
            const auto breathingWarmMomentumPressureEpoch =
                acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                    breathingWarmPrediction,
                    breathingPressureWarmStart,
                    breathingWarmPressureOperator,
                    breathingWarmBasePressureOperator, geometry,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, breathingWarmVolumeRates,
                    definitions, breathingWarmOpenings, zeroResistance,
                    breathingWarmBaseMetric, breathingWarmMetric,
                    momentumPressureSettings);
            check(breathingWarmMomentumPressureEpoch.diagnostics.accepted
                      && breathingWarmMomentumPressureEpoch.diagnostics
                             .usedWarmPressureStart
                      && !breathingWarmMomentumPressureEpoch.diagnostics
                              .usedColdPressureStart
                      && breathingWarmMomentumPressureEpoch.acceptedState
                             .accepted
                      && breathingWarmMomentumPressureEpoch.acceptedState
                             .correctedContinuityResidualMaximumCubicMetersPerSecond
                          <= breathingWarmMomentumPressureEpoch
                                 .acceptedState
                                 .continuityToleranceCubicMetersPerSecond,
                  "re-entrant breathing momentum reaches a constrained warm pressure epoch");
            validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
                breathingWarmMomentumPressureEpoch,
                breathingWarmPrediction, breathingPressureWarmStart,
                breathingWarmPressureOperator,
                breathingWarmBasePressureOperator, geometry,
                breathingWarmSweep, breathingWarmFragments,
                breathingWarmTopology, breathingWarmVolumeRates,
                definitions, breathingWarmOpenings, zeroResistance,
                breathingWarmBaseMetric, breathingWarmMetric,
                momentumPressureSettings);
            const auto breathingMomentumCycle =
                advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                    breathingTransport, breathingMetric,
                    breathingMomentumPressureEpoch.acceptedState,
                    breathingPredictionPressureOperator,
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings, zeroResistance,
                    breathingPredictionBaseMetric,
                    breathingPredictionMetric,
                    breathingWarmPressureOperator,
                    breathingWarmBasePressureOperator,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, breathingWarmVolumeRates,
                    definitions, breathingWarmOpenings, zeroResistance,
                    breathingWarmBaseMetric, breathingWarmMetric, {},
                    momentumPressureSettings);
            check(breathingMomentumCycle.diagnostics.accepted
                      && breathingMomentumCycle.diagnostics
                             .usedReentrantTransport
                      && breathingMomentumCycle.diagnostics
                             .usedPressureWarmStart
                      && breathingMomentumCycle.transport
                          == breathingSecondTransport
                      && breathingMomentumCycle.acceptedState
                          == breathingWarmMomentumPressureEpoch.acceptedState,
                  "atomic transported cycle reproduces both explicit breathing endpoints");
            validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
                breathingMomentumCycle, breathingTransport,
                breathingMetric,
                breathingMomentumPressureEpoch.acceptedState,
                breathingPredictionPressureOperator,
                breathingPredictionBasePressureOperator, geometry,
                breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates, definitions,
                breathingPredictionOpenings, zeroResistance,
                breathingPredictionBaseMetric,
                breathingPredictionMetric,
                breathingWarmPressureOperator,
                breathingWarmBasePressureOperator,
                breathingWarmSweep, breathingWarmFragments,
                breathingWarmTopology, breathingWarmVolumeRates,
                definitions, breathingWarmOpenings, zeroResistance,
                breathingWarmBaseMetric, breathingWarmMetric);
            const auto breathingBootstrapCycle =
                advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                    breathingFlowState, breathingMetric,
                    breathingMomentumPressureEpoch.acceptedState,
                    breathingPredictionPressureOperator,
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings, zeroResistance,
                    breathingPredictionBaseMetric,
                    breathingPredictionMetric,
                    breathingWarmPressureOperator,
                    breathingWarmBasePressureOperator,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, breathingWarmVolumeRates,
                    definitions, breathingWarmOpenings, zeroResistance,
                    breathingWarmBaseMetric, breathingWarmMetric, {},
                    momentumPressureSettings);
            check(breathingBootstrapCycle.diagnostics.accepted
                      && breathingBootstrapCycle.diagnostics
                             .usedInitialStateTransport
                      && !breathingBootstrapCycle.diagnostics
                              .usedReentrantTransport
                      && breathingBootstrapCycle.diagnostics.transport
                             .maximumAbsoluteOpeningRelativeVolumeFlowRateCubicMetersPerSecond
                          > 0.0
                      && breathingBootstrapCycle.acceptedState.accepted,
                  "live-aperture velocity state bootstraps a breathing atomic cycle");
            validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
                breathingBootstrapCycle, breathingFlowState,
                breathingMetric,
                breathingMomentumPressureEpoch.acceptedState,
                breathingPredictionPressureOperator,
                breathingPredictionBasePressureOperator, geometry,
                breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates, definitions,
                breathingPredictionOpenings, zeroResistance,
                breathingPredictionBaseMetric,
                breathingPredictionMetric,
                breathingWarmPressureOperator,
                breathingWarmBasePressureOperator,
                breathingWarmSweep, breathingWarmFragments,
                breathingWarmTopology, breathingWarmVolumeRates,
                definitions, breathingWarmOpenings, zeroResistance,
                breathingWarmBaseMetric, breathingWarmMetric);

            auto truncatedMomentumPressureSettings =
                momentumPressureSettings;
            truncatedMomentumPressureSettings.projection.pressureSolve
                .absoluteResidualTolerancePascalsMeters = 1.0e-16;
            truncatedMomentumPressureSettings.projection.pressureSolve
                .relativeResidualTolerance = 0.0;
            truncatedMomentumPressureSettings.projection.pressureSolve
                .maximumIterations = 1;
            const auto retainedBreathingPrediction = breathingPrediction;
            const auto rejectedMomentumPressureEpoch =
                acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                    breathingPrediction,
                    breathingPredictionPressureOperator,
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings, zeroResistance,
                    breathingPredictionBaseMetric,
                    breathingPredictionMetric,
                    truncatedMomentumPressureSettings);
            check(!rejectedMomentumPressureEpoch.diagnostics.accepted
                      && rejectedMomentumPressureEpoch.diagnostics.failureStage
                          == PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                              PressureProjection
                      && rejectedMomentumPressureEpoch.diagnostics
                             .pressureStep.resistance.accepted
                      && !rejectedMomentumPressureEpoch.acceptedState.accepted
                      && rejectedMomentumPressureEpoch.acceptedState
                             .fingerprint
                          == 0
                      && breathingPrediction
                          == retainedBreathingPrediction,
                  "truncated transported pressure epoch publishes no partial endpoint");
            validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
                rejectedMomentumPressureEpoch, breathingPrediction,
                breathingPredictionPressureOperator,
                breathingPredictionBasePressureOperator, geometry,
                breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates, definitions,
                breathingPredictionOpenings, zeroResistance,
                breathingPredictionBaseMetric,
                breathingPredictionMetric,
                truncatedMomentumPressureSettings);

            const auto retainedBreathingWarmPrediction =
                breathingWarmPrediction;
            const auto retainedBreathingPressureWarmStart =
                breathingPressureWarmStart;
            const auto rejectedWarmMomentumPressureEpoch =
                acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                    breathingWarmPrediction,
                    breathingPressureWarmStart,
                    breathingWarmPressureOperator,
                    breathingWarmBasePressureOperator, geometry,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, breathingWarmVolumeRates,
                    definitions, breathingWarmOpenings, zeroResistance,
                    breathingWarmBaseMetric, breathingWarmMetric,
                    truncatedMomentumPressureSettings);
            check(!rejectedWarmMomentumPressureEpoch.diagnostics.accepted
                      && rejectedWarmMomentumPressureEpoch.diagnostics
                             .usedWarmPressureStart
                      && rejectedWarmMomentumPressureEpoch.diagnostics
                             .failureStage
                          == PlanarPressureRegionFragmentOpeningMomentumPressureEpochFailureStage::
                              PressureProjection
                      && !rejectedWarmMomentumPressureEpoch.acceptedState
                              .accepted
                      && rejectedWarmMomentumPressureEpoch.acceptedState
                             .fingerprint
                          == 0
                      && breathingWarmPrediction
                          == retainedBreathingWarmPrediction
                      && breathingPressureWarmStart
                          == retainedBreathingPressureWarmStart,
                  "truncated transported warm pressure epoch rolls back without mutating either source");
            validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResult(
                rejectedWarmMomentumPressureEpoch,
                breathingWarmPrediction, breathingPressureWarmStart,
                breathingWarmPressureOperator,
                breathingWarmBasePressureOperator, geometry,
                breathingWarmSweep, breathingWarmFragments,
                breathingWarmTopology, breathingWarmVolumeRates,
                definitions, breathingWarmOpenings, zeroResistance,
                breathingWarmBaseMetric, breathingWarmMetric,
                truncatedMomentumPressureSettings);
            const auto rejectedPressureMomentumCycle =
                advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                    breathingTransport, breathingMetric,
                    breathingMomentumPressureEpoch.acceptedState,
                    breathingPredictionPressureOperator,
                    breathingPredictionBasePressureOperator, geometry,
                    breathingPredictionSweep,
                    breathingPredictionFragments,
                    breathingPredictionTopology,
                    breathingPredictionVolumeRates, definitions,
                    breathingPredictionOpenings, zeroResistance,
                    breathingPredictionBaseMetric,
                    breathingPredictionMetric,
                    breathingWarmPressureOperator,
                    breathingWarmBasePressureOperator,
                    breathingWarmSweep, breathingWarmFragments,
                    breathingWarmTopology, breathingWarmVolumeRates,
                    definitions, breathingWarmOpenings, zeroResistance,
                    breathingWarmBaseMetric, breathingWarmMetric, {},
                    truncatedMomentumPressureSettings);
            check(!rejectedPressureMomentumCycle.diagnostics.accepted
                      && rejectedPressureMomentumCycle.diagnostics
                             .failureStage
                          == PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                              PressureProjection
                      && rejectedPressureMomentumCycle.diagnostics
                             .transport.accepted
                      && rejectedPressureMomentumCycle.diagnostics
                             .prediction.finite
                      && rejectedPressureMomentumCycle.transport.fingerprint
                          == 0
                      && rejectedPressureMomentumCycle.transport.controls.empty()
                      && rejectedPressureMomentumCycle.acceptedState.fingerprint
                          == 0,
                  "atomic transported cycle publishes neither endpoint after pressure rejection");
            validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
                rejectedPressureMomentumCycle, breathingTransport,
                breathingMetric,
                breathingMomentumPressureEpoch.acceptedState,
                breathingPredictionPressureOperator,
                breathingPredictionBasePressureOperator, geometry,
                breathingPredictionSweep,
                breathingPredictionFragments,
                breathingPredictionTopology,
                breathingPredictionVolumeRates, definitions,
                breathingPredictionOpenings, zeroResistance,
                breathingPredictionBaseMetric,
                breathingPredictionMetric,
                breathingWarmPressureOperator,
                breathingWarmBasePressureOperator,
                breathingWarmSweep, breathingWarmFragments,
                breathingWarmTopology, breathingWarmVolumeRates,
                definitions, breathingWarmOpenings, zeroResistance,
                breathingWarmBaseMetric, breathingWarmMetric);
        }

        auto substepSettings =
            PlanarPressureRegionFragmentOpeningMomentumTransportSettings{};
        substepSettings.maximumOutgoingCourantNumber = 1.0e-6;
        substepSettings.maximumSubsteps = 1;
        const auto substepRejected =
            advancePlanarPressureRegionFragmentOpeningMomentum(
                sourceState, sourceMetric, targetFlowState, targetMetric,
                geometry, targetSweep, targetFragments, targetTopology,
                targetVolumeRates, substepSettings);
        check(!substepRejected.diagnostics.accepted
                  && substepRejected.controls.empty()
                  && substepRejected.diagnostics.substepCount > 1
                  && substepRejected.diagnostics.failureStage
                      == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                          SubstepLimit,
              "opening momentum transport rejects a step beyond its subcycling limit");
        const auto rejectedTransportMomentumCycle =
            advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                transport, targetMetric,
                momentumPressureEpoch.acceptedState,
                predictionPressureOperator,
                predictionBasePressureOperator, geometry, predictionSweep,
                predictionFragments, predictionTopology,
                predictionVolumeRates, definitions, predictionOpenings,
                zeroResistance, predictionBaseMetric, predictionMetric,
                warmPressureOperator, warmBasePressureOperator, warmSweep,
                warmFragments, warmTopology, warmVolumeRates, definitions,
                warmOpenings, zeroResistance, warmBaseMetric, warmMetric,
                substepSettings, momentumPressureSettings);
        check(!rejectedTransportMomentumCycle.diagnostics.accepted
                  && rejectedTransportMomentumCycle.diagnostics.failureStage
                      == PlanarPressureRegionFragmentOpeningMomentumCycleFailureStage::
                          MomentumTransport
                  && !rejectedTransportMomentumCycle.diagnostics
                          .transport.accepted
                  && rejectedTransportMomentumCycle.predictionFingerprint == 0
                  && rejectedTransportMomentumCycle
                         .pressureWarmStartFingerprint
                      == 0
                  && rejectedTransportMomentumCycle.transport.fingerprint == 0
                  && rejectedTransportMomentumCycle.acceptedState.fingerprint
                      == 0,
              "atomic transported cycle publishes neither endpoint after transport rejection");
        validatePlanarPressureRegionFragmentOpeningMomentumCycleResult(
            rejectedTransportMomentumCycle, transport, targetMetric,
            momentumPressureEpoch.acceptedState,
            predictionPressureOperator,
            predictionBasePressureOperator, geometry, predictionSweep,
            predictionFragments, predictionTopology,
            predictionVolumeRates, definitions, predictionOpenings,
            zeroResistance, predictionBaseMetric, predictionMetric,
            warmPressureOperator, warmBasePressureOperator, warmSweep,
            warmFragments, warmTopology, warmVolumeRates, definitions,
            warmOpenings, zeroResistance, warmBaseMetric, warmMetric);
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentum(
                        substepRejected, targetMetric,
                        acceptedPressureFlow, predictionMetric, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates));
            },
            "re-entrant momentum transport rejects an unaccepted source transport");
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentum(
                        transport, sourceMetric, acceptedPressureFlow,
                        predictionMetric, geometry, predictionSweep,
                        predictionFragments, predictionTopology,
                        predictionVolumeRates));
            },
            "re-entrant momentum transport rejects a foreign source metric");
        auto corruptSourceTransport = transport;
        corruptSourceTransport.controls[0]
            .momentumKilogramMetersPerSecond.x += 0.1;
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentum(
                        corruptSourceTransport, targetMetric,
                        acceptedPressureFlow, predictionMetric, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates));
            },
            "re-entrant momentum transport rejects corrupted source controls");

        std::vector<double> incompatibleNormal(
            targetMetric.dofs.size(), 0.0);
        std::vector<double> incompatibleMaterial(
            targetMetric.dofs.size(), 0.0);
        std::vector<double> incompatibleRelative(
            targetMetric.dofs.size(), 0.0);
        for (const auto& dof : targetMetric.dofs) {
            if (dof.kind
                != PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    SharedRegionGrid) {
                const double material = dof.axis == GridFaceAxis::X
                    ? uniformVelocity.x
                    : dof.axis == GridFaceAxis::Y
                    ? uniformVelocity.y : uniformVelocity.z;
                incompatibleMaterial[dof.dofIndex] = material;
                incompatibleNormal[dof.dofIndex] = material;
            }
        }
        const auto incompatibleFlow =
            buildPlanarPressureRegionFragmentOpeningVelocityState(
                targetMetric, incompatibleNormal, incompatibleMaterial,
                incompatibleRelative, 1.2);
        const auto rejected =
            advancePlanarPressureRegionFragmentOpeningMomentum(
                sourceState, sourceMetric, incompatibleFlow, targetMetric,
                geometry, targetSweep, targetFragments, targetTopology,
                targetVolumeRates);
        check(!rejected.diagnostics.accepted
                  && rejected.controls.empty()
                  && rejected.diagnostics.failureStage
                      == PlanarPressureRegionFragmentOpeningMomentumTransportFailureStage::
                          FlowContinuity,
              "opening momentum transport rejects relative flow that violates moving continuity");
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentum(
                        sourceState, sourceMetric, sourceState, sourceMetric,
                        geometry, targetSweep, targetFragments,
                        targetTopology, targetVolumeRates));
            },
            "opening momentum transport rejects a target metric foreign to the current geometry");
        auto corruptTransport = transport;
        corruptTransport.controls[0].volumeCubicMeters += 0.1;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningMomentumTransportIntegrity(
                    corruptTransport);
            },
            "opening momentum transport rejects control corruption");
        auto transportLimits =
            PlanarPressureRegionFragmentOpeningMomentumTransportLimits{};
        transportLimits.maximumOwnedBytes = transport.ownedStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentum(
                        sourceState, sourceMetric, targetFlowState,
                        targetMetric, geometry, targetSweep, targetFragments,
                        targetTopology, targetVolumeRates, {},
                        transportLimits));
            },
            "opening momentum transport enforces its owned-storage limit before allocation");
        transportLimits = {};
        transportLimits.maximumWorkingBytes =
            transport.workingStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentum(
                        sourceState, sourceMetric, targetFlowState,
                        targetMetric, geometry, targetSweep, targetFragments,
                        targetTopology, targetVolumeRates, {},
                        transportLimits));
            },
            "opening momentum transport enforces its working-storage limit before allocation");

        auto corruptPrediction = prediction;
        corruptPrediction.predictedVelocityState.samples[0]
            .normalVelocityMetersPerSecond += 0.1;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningMomentumPredictionIntegrity(
                    corruptPrediction);
            },
            "opening momentum prediction rejects nested state corruption");
        expectRejected(
            [&] {
                static_cast<void>(
                    predictPlanarPressureRegionFragmentOpeningMomentum(
                        transport, sourceMetric, geometry, predictionSweep,
                        predictionFragments, predictionTopology,
                        predictionVolumeRates, definitions,
                        predictionOpenings, predictionBaseMetric,
                        predictionMetric));
            },
            "opening momentum prediction rejects a foreign transport endpoint metric");
        auto predictionLimits =
            PlanarPressureRegionFragmentOpeningMomentumPredictionLimits{};
        predictionLimits.maximumOwnedBytes =
            prediction.ownedStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    predictPlanarPressureRegionFragmentOpeningMomentum(
                        transport, targetMetric, geometry, predictionSweep,
                        predictionFragments, predictionTopology,
                        predictionVolumeRates, definitions,
                        predictionOpenings, predictionBaseMetric,
                        predictionMetric, predictionLimits));
            },
            "opening momentum prediction enforces its owned-storage limit before allocation");
        predictionLimits = {};
        predictionLimits.maximumWorkingBytes =
            prediction.workingStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    predictPlanarPressureRegionFragmentOpeningMomentum(
                        transport, targetMetric, geometry, predictionSweep,
                        predictionFragments, predictionTopology,
                        predictionVolumeRates, definitions,
                        predictionOpenings, predictionBaseMetric,
                        predictionMetric, predictionLimits));
            },
            "opening momentum prediction enforces its working-storage limit before allocation");

        auto corruptMomentumPressureEpoch = momentumPressureEpoch;
        corruptMomentumPressureEpoch.acceptedState
            .pressureCorrectionPascals[0] += 0.1;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningMomentumPressureEpochResultIntegrity(
                    corruptMomentumPressureEpoch);
            },
            "opening momentum pressure epoch rejects accepted-state corruption");
        auto momentumPressureLimits =
            PlanarPressureRegionFragmentOpeningMomentumPressureEpochLimits{};
        momentumPressureLimits.maximumOwnedBytes =
            momentumPressureEpoch.ownedStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                        prediction, predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        predictionBaseMetric, predictionMetric,
                        momentumPressureSettings, momentumPressureLimits));
            },
            "opening momentum pressure epoch enforces its owned-storage limit");
        momentumPressureLimits = {};
        momentumPressureLimits.maximumWorkingBytes =
            momentumPressureEpoch.workingStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                        prediction, predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        predictionBaseMetric, predictionMetric,
                        momentumPressureSettings, momentumPressureLimits));
            },
            "opening momentum pressure epoch enforces its working-storage limit");

        auto corruptPressureWarmStart = pressureWarmStart;
        corruptPressureWarmStart.pressureCorrectionPascals[0] += 0.1;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningPressureWarmStartIntegrity(
                    corruptPressureWarmStart);
            },
            "opening pressure warm start rejects correction corruption");
        expectRejected(
            [&] {
                static_cast<void>(
                    acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                        prediction, pressureWarmStart,
                        predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        predictionBaseMetric, predictionMetric,
                        momentumPressureSettings));
            },
            "transported pressure epoch rejects a warm start bound to a foreign current epoch");
        auto pressureWarmStartLimits =
            PlanarPressureRegionFragmentOpeningPressureWarmStartLimits{};
        pressureWarmStartLimits.maximumOwnedBytes =
            pressureWarmStart.ownedStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
                        momentumPressureEpoch.acceptedState,
                        predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        warmPressureOperator, warmBasePressureOperator,
                        warmSweep, warmFragments, warmTopology,
                        warmVolumeRates, definitions, warmOpenings,
                        pressureWarmStartLimits));
            },
            "opening pressure warm start enforces its owned-storage limit");
        pressureWarmStartLimits = {};
        pressureWarmStartLimits.maximumWorkingBytes =
            pressureWarmStart.workingStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    buildPlanarPressureRegionFragmentOpeningPressureWarmStart(
                        momentumPressureEpoch.acceptedState,
                        predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        warmPressureOperator, warmBasePressureOperator,
                        warmSweep, warmFragments, warmTopology,
                        warmVolumeRates, definitions, warmOpenings,
                        pressureWarmStartLimits));
            },
            "opening pressure warm start enforces its working-storage limit before allocation");
        momentumPressureLimits = {};
        momentumPressureLimits.warmStartLimits.maximumOwnedBytes =
            pressureWarmStart.ownedStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    acceptPlanarPressureRegionFragmentOpeningMomentumPressureEpoch(
                        warmPrediction, pressureWarmStart,
                        warmPressureOperator, warmBasePressureOperator,
                        geometry, warmSweep, warmFragments, warmTopology,
                        warmVolumeRates, definitions, warmOpenings,
                        zeroResistance, warmBaseMetric, warmMetric,
                        momentumPressureSettings, momentumPressureLimits));
            },
            "transported warm pressure epoch enforces its source-artifact limit");

        auto corruptMomentumCycle = momentumCycle;
        corruptMomentumCycle.transport.controls[0]
            .momentumKilogramMetersPerSecond.x += 0.1;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
                    corruptMomentumCycle);
            },
            "atomic transported cycle rejects nested endpoint corruption");
        auto ambiguousBootstrapMomentumCycle = bootstrapMomentumCycle;
        ambiguousBootstrapMomentumCycle.sourceTransportFingerprint =
            transport.fingerprint;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningMomentumCycleResultIntegrity(
                    ambiguousBootstrapMomentumCycle);
            },
            "atomic momentum bootstrap rejects ambiguous source lineage");
        auto corruptCapturedCycleState = bootstrapCycleState;
        corruptCapturedCycleState.acceptedState
            .pressureCorrectionPascals[0] += 0.1;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningMomentumCycleStateIntegrity(
                    corruptCapturedCycleState);
            },
            "opening momentum-cycle restart state rejects nested endpoint corruption");
        expectRejected(
            [&] {
                static_cast<void>(
                    capturePlanarPressureRegionFragmentOpeningMomentumCycleState(
                        rejectedTransportMomentumCycle,
                        predictionMetric, warmMetric));
            },
            "opening momentum-cycle restart capture rejects an unaccepted cycle");
        auto cycleStateLimits =
            PlanarPressureRegionFragmentOpeningMomentumCycleStateLimits{};
        cycleStateLimits.maximumOwnedBytes =
            bootstrapCycleState.ownedStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    capturePlanarPressureRegionFragmentOpeningMomentumCycleState(
                        bootstrapMomentumCycle, predictionMetric,
                        warmMetric, cycleStateLimits));
            },
            "opening momentum-cycle restart capture enforces endpoint storage limits");
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                        transport, sourceMetric,
                        momentumPressureEpoch.acceptedState,
                        predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        predictionBaseMetric, predictionMetric,
                        warmPressureOperator, warmBasePressureOperator,
                        warmSweep, warmFragments, warmTopology,
                        warmVolumeRates, definitions, warmOpenings,
                        zeroResistance, warmBaseMetric, warmMetric, {},
                        momentumPressureSettings));
            },
            "atomic transported cycle rejects a foreign source metric");
        auto momentumCycleLimits =
            PlanarPressureRegionFragmentOpeningMomentumCycleLimits{};
        momentumCycleLimits.maximumOwnedBytes =
            momentumCycle.ownedStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                        transport, targetMetric,
                        momentumPressureEpoch.acceptedState,
                        predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        predictionBaseMetric, predictionMetric,
                        warmPressureOperator, warmBasePressureOperator,
                        warmSweep, warmFragments, warmTopology,
                        warmVolumeRates, definitions, warmOpenings,
                        zeroResistance, warmBaseMetric, warmMetric, {},
                        momentumPressureSettings, momentumCycleLimits));
            },
            "atomic transported cycle enforces its endpoint storage limit");
        momentumCycleLimits = {};
        momentumCycleLimits.maximumWorkingBytes =
            momentumCycle.workingStorageBytes - 1;
        expectRejected(
            [&] {
                static_cast<void>(
                    advancePlanarPressureRegionFragmentOpeningMomentumCycle(
                        transport, targetMetric,
                        momentumPressureEpoch.acceptedState,
                        predictionPressureOperator,
                        predictionBasePressureOperator, geometry,
                        predictionSweep, predictionFragments,
                        predictionTopology, predictionVolumeRates,
                        definitions, predictionOpenings, zeroResistance,
                        predictionBaseMetric, predictionMetric,
                        warmPressureOperator, warmBasePressureOperator,
                        warmSweep, warmFragments, warmTopology,
                        warmVolumeRates, definitions, warmOpenings,
                        zeroResistance, warmBaseMetric, warmMetric, {},
                        momentumPressureSettings, momentumCycleLimits));
            },
            "atomic transported cycle enforces aggregate working storage");
    }
}

void testPlanarRegionalResistedOpeningPressureStep() {
    const auto geometry = grid();
    const auto previous = pocketLayers();
    auto current = previous;
    current[0].physicalPlaneCoordinateMeters -= 0.1;
    current[1].physicalPlaneCoordinateMeters += 0.1;
    const auto sweep = makePlanarPressureRegionSweepLedger(
        geometry, previous, current, 0.5);
    const auto fragments = buildPlanarPressureRegionFragments(
        geometry, sweep);
    const auto topology = buildPlanarPressureRegionFragmentTopology(
        geometry, sweep, fragments);
    const auto base = buildPlanarPressureRegionFragmentPressureOperator(
        geometry, sweep, fragments, topology);
    const auto volumeRates = buildPlanarPressureRegionFragmentVolumeRates(
        geometry, sweep, fragments, topology);
    const auto wall = std::ranges::find_if(
        topology.links,
        [](const auto& link) {
            return link.kind
                    == PlanarPressureRegionFragmentFaceKind::PressureLayerWall
                && link.surfaceStableId == 10;
        });
    check(wall != topology.links.end(),
          "resisted opening pressure step finds its intake wall");
    if (wall == topology.links.end()) return;
    const double patchArea = 0.5 * wall->areaSquareMeters;
    const std::vector<PlanarPressureRegionFragmentOpeningPatchDefinition>
        definitions{{100, 1000, wall->surfaceStableId, wall->axis,
            wall->i, wall->j, wall->k, wall->minusRegionStableId,
            wall->plusRegionStableId, patchArea}};
    const auto openings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, definitions);
    const auto pressureOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            base, geometry, sweep, fragments, topology, definitions,
            openings);
    const auto velocityMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, sweep, fragments, topology);
    const auto openingVelocityMetric =
        buildPlanarPressureRegionFragmentOpeningVelocityMetric(
            geometry, sweep, fragments, topology, velocityMetric,
            definitions, openings);
    const auto repeatedOpeningVelocityMetric =
        buildPlanarPressureRegionFragmentOpeningVelocityMetric(
            geometry, sweep, fragments, topology, velocityMetric,
            definitions, openings);
    check(openingVelocityMetric == repeatedOpeningVelocityMetric
              && openingVelocityMetric.version
                  == planarPressureRegionFragmentOpeningVelocityMetricVersion
              && openingVelocityMetric.fingerprint != 0
              && openingVelocityMetric.sourceBaseMetricFingerprint
                  == velocityMetric.fingerprint
              && openingVelocityMetric.sourceOpeningFingerprint
                  == openings.fingerprint
              && openingVelocityMetric.sharedRegionGridDofCount
                  == topology.sameRegionGridLinkCount
              && openingVelocityMetric.solidWallTraceDofCount
                  == 2 * topology.pressureLayerWallLinkCount
              && openingVelocityMetric.openingPatchDofCount
                  == openings.patches.size()
              && openingVelocityMetric.dofs.size()
                  == velocityMetric.dofs.size() + openings.patches.size()
              && openingVelocityMetric.fragments.size()
                  == fragments.fragments.size()
              && openingVelocityMetric.components.size()
                  == openings.connectedComponents.size()
              && openingVelocityMetric.ownedStorageBytes > 0
              && openingVelocityMetric.workingStorageBytes > 0,
          "opening velocity metric deterministically partitions the regional inertia basis");
    checkNear(
        openingVelocityMetric.totalOpeningAreaSquareMeters,
        patchArea, 2.0e-15,
        "opening velocity metric retains exact aperture area");
    checkNear(
        openingVelocityMetric.totalSolidWallAreaSquareMeters
            + openingVelocityMetric.totalOpeningAreaSquareMeters,
        openingVelocityMetric.totalPressureWallAreaSquareMeters,
        2.0e-15,
        "opening velocity metric closes the physical wall-area partition");
    const double domainVolume = geometry.cellVolumeCubicMeters()
        * static_cast<double>(geometry.cellCount());
    checkNear(openingVelocityMetric.dualVolumeByAxisCubicMeters.x,
              domainVolume, 2.0e-13,
              "opening velocity metric closes X domain volume");
    checkNear(openingVelocityMetric.dualVolumeByAxisCubicMeters.y,
              domainVolume, 2.0e-13,
              "opening velocity metric closes Y domain volume");
    checkNear(openingVelocityMetric.dualVolumeByAxisCubicMeters.z,
              domainVolume, 2.0e-13,
              "opening velocity metric closes Z domain volume");
    checkNear(
        std::max({
            std::abs(openingVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.x),
            std::abs(openingVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.y),
            std::abs(openingVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.z),
            openingVelocityMetric
                .maximumAbsoluteFragmentVolumeClosureResidualCubicMeters,
            openingVelocityMetric
                .maximumAbsoluteComponentVolumeClosureResidualCubicMeters}),
        0.0, 2.0e-13,
        "opening velocity metric closes every fragment and connected component");
    const auto openingDof = std::ranges::find_if(
        openingVelocityMetric.dofs,
        [](const auto& dof) {
            return dof.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    OpeningPatch;
        });
    check(openingDof != openingVelocityMetric.dofs.end()
              && openingDof->sourceOpeningPatchStableId == 100
              && openingDof->sourceFaceLinkStableId == wall->stableId
              && openingDof->ownerFragmentIndex
                  == wall->minusFragmentIndex
              && openingDof->oppositeFragmentIndex
                  == wall->plusFragmentIndex,
          "opening velocity metric binds the aperture to both fluid half-volumes");
    if (openingDof != openingVelocityMetric.dofs.end()) {
        checkNear(
            openingDof->dualVolumeCubicMeters,
            patchArea * wall->centerDistanceMeters, 2.0e-15,
            "opening velocity metric recovers the exact aperture plug volume");
    }
    validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
        openingVelocityMetric);
    validatePlanarPressureRegionFragmentOpeningVelocityMetric(
        openingVelocityMetric, geometry, sweep, fragments, topology,
        velocityMetric, definitions, openings);
    auto corruptOpeningVelocityMetric = openingVelocityMetric;
    corruptOpeningVelocityMetric.dofs[0].dualVolumeCubicMeters += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningVelocityMetricIntegrity(
                corruptOpeningVelocityMetric);
        },
        "opening velocity metric rejects DOF corruption");
    auto openingVelocityMetricLimits =
        PlanarPressureRegionFragmentOpeningVelocityMetricLimits{};
    openingVelocityMetricLimits.maximumOwnedBytes =
        openingVelocityMetric.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                    geometry, sweep, fragments, topology, velocityMetric,
                    definitions, openings, openingVelocityMetricLimits));
        },
        "opening velocity metric enforces its owned-storage limit");
    constexpr double manufacturedDensity = 1.2;
    const Vector3 manufacturedVelocity{2.0, -0.5, 0.25};
    std::vector<double> manufacturedNormal(
        openingVelocityMetric.dofs.size(), 0.0);
    std::vector<double> manufacturedMaterial(
        openingVelocityMetric.dofs.size(), 0.0);
    std::vector<double> manufacturedRelative(
        openingVelocityMetric.dofs.size(), 0.0);
    for (const auto& dof : openingVelocityMetric.dofs) {
        const double component = dof.axis == GridFaceAxis::X
            ? manufacturedVelocity.x
            : dof.axis == GridFaceAxis::Y
            ? manufacturedVelocity.y : manufacturedVelocity.z;
        manufacturedNormal[dof.dofIndex] = component;
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            manufacturedRelative[dof.dofIndex] = component;
        } else {
            manufacturedMaterial[dof.dofIndex] = component;
        }
    }
    const auto manufacturedOpeningVelocityState =
        buildPlanarPressureRegionFragmentOpeningVelocityState(
            openingVelocityMetric, manufacturedNormal,
            manufacturedMaterial, manufacturedRelative,
            manufacturedDensity);
    const auto repeatedManufacturedOpeningVelocityState =
        buildPlanarPressureRegionFragmentOpeningVelocityState(
            openingVelocityMetric, manufacturedNormal,
            manufacturedMaterial, manufacturedRelative,
            manufacturedDensity);
    const double manufacturedMass = manufacturedDensity * domainVolume;
    check(manufacturedOpeningVelocityState
                  == repeatedManufacturedOpeningVelocityState
              && manufacturedOpeningVelocityState.version
                  == planarPressureRegionFragmentOpeningVelocityStateVersion
              && manufacturedOpeningVelocityState.fingerprint != 0
              && !manufacturedOpeningVelocityState
                      .mappedFromAcceptedEndpoint
              && manufacturedOpeningVelocityState
                     .sourceAcceptedStateFingerprint == 0
              && manufacturedOpeningVelocityState
                     .sourceVolumeRateFingerprint == 0
              && manufacturedOpeningVelocityState.samples.size()
                  == openingVelocityMetric.dofs.size()
              && manufacturedOpeningVelocityState.fragments.size()
                  == openingVelocityMetric.fragments.size()
              && manufacturedOpeningVelocityState.components.size()
                  == openingVelocityMetric.components.size()
              && manufacturedOpeningVelocityState
                     .sharedRegionGridSampleCount
                  == openingVelocityMetric.sharedRegionGridDofCount
              && manufacturedOpeningVelocityState
                     .solidWallTraceSampleCount
                  == openingVelocityMetric.solidWallTraceDofCount
              && manufacturedOpeningVelocityState.openingPatchSampleCount
                  == openingVelocityMetric.openingPatchDofCount,
          "opening velocity state deterministically owns the complete partitioned inertia basis");
    checkNear(
        manufacturedOpeningVelocityState.physicalMassKilograms,
        manufacturedMass, 2.0e-13,
        "opening velocity state closes physical mass once");
    checkNear(
        manufacturedOpeningVelocityState.diagonalMassByAxisKilograms.x,
        manufacturedMass, 2.0e-13,
        "opening velocity state closes X diagonal mass");
    checkNear(
        manufacturedOpeningVelocityState.diagonalMassByAxisKilograms.y,
        manufacturedMass, 2.0e-13,
        "opening velocity state closes Y diagonal mass");
    checkNear(
        manufacturedOpeningVelocityState.diagonalMassByAxisKilograms.z,
        manufacturedMass, 2.0e-13,
        "opening velocity state closes Z diagonal mass");
    checkNear(
        manufacturedOpeningVelocityState
            .momentumKilogramMetersPerSecond.x,
        manufacturedMass * manufacturedVelocity.x, 4.0e-13,
        "opening velocity state preserves uniform X momentum");
    checkNear(
        manufacturedOpeningVelocityState
            .momentumKilogramMetersPerSecond.y,
        manufacturedMass * manufacturedVelocity.y, 2.0e-13,
        "opening velocity state preserves uniform Y momentum");
    checkNear(
        manufacturedOpeningVelocityState
            .momentumKilogramMetersPerSecond.z,
        manufacturedMass * manufacturedVelocity.z, 1.0e-13,
        "opening velocity state preserves uniform Z momentum");
    const double manufacturedEnergy = 0.5 * manufacturedMass
        * (manufacturedVelocity.x * manufacturedVelocity.x
           + manufacturedVelocity.y * manufacturedVelocity.y
           + manufacturedVelocity.z * manufacturedVelocity.z);
    checkNear(
        manufacturedOpeningVelocityState.diagonalKineticEnergyJoules,
        manufacturedEnergy, 6.0e-13,
        "opening velocity state preserves uniform diagonal energy");
    checkNear(
        manufacturedOpeningVelocityState.collocatedKineticEnergyJoules,
        manufacturedEnergy, 6.0e-13,
        "opening velocity state reconstructs uniform collocated energy");
    checkNear(
        manufacturedOpeningVelocityState.staggeringKineticEnergyJoules,
        0.0, 8.0e-13,
        "uniform opening flow has zero staggering energy");
    for (const auto& fragmentState
         : manufacturedOpeningVelocityState.fragments) {
        checkNear(fragmentState.collocatedVelocityMetersPerSecond.x,
                  manufacturedVelocity.x, 2.0e-14,
                  "uniform opening flow reconstructs fragment X velocity");
        checkNear(fragmentState.collocatedVelocityMetersPerSecond.y,
                  manufacturedVelocity.y, 2.0e-14,
                  "uniform opening flow reconstructs fragment Y velocity");
        checkNear(fragmentState.collocatedVelocityMetersPerSecond.z,
                  manufacturedVelocity.z, 2.0e-14,
                  "uniform opening flow reconstructs fragment Z velocity");
    }
    validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(
        manufacturedOpeningVelocityState);
    validatePlanarPressureRegionFragmentOpeningVelocityState(
        manufacturedOpeningVelocityState, openingVelocityMetric);
    auto corruptOpeningVelocityState = manufacturedOpeningVelocityState;
    corruptOpeningVelocityState.samples[0]
        .normalMomentumKilogramMetersPerSecond += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningVelocityStateIntegrity(
                corruptOpeningVelocityState);
        },
        "opening velocity state rejects sample corruption");
    auto openingVelocityStateLimits =
        PlanarPressureRegionFragmentOpeningVelocityStateLimits{};
    openingVelocityStateLimits.maximumOwnedBytes =
        manufacturedOpeningVelocityState.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                buildPlanarPressureRegionFragmentOpeningVelocityState(
                    openingVelocityMetric, manufacturedNormal,
                    manufacturedMaterial, manufacturedRelative,
                    manufacturedDensity, openingVelocityStateLimits));
        },
        "opening velocity state enforces its owned-storage limit before allocation");
    auto invalidNormal = manufacturedNormal;
    auto invalidRelative = manufacturedRelative;
    const auto solidDof = std::ranges::find_if(
        openingVelocityMetric.dofs,
        [](const auto& dof) {
            return dof.kind
                    == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                        SolidWallMinusTrace
                || dof.kind
                    == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                        SolidWallPlusTrace;
        });
    check(solidDof != openingVelocityMetric.dofs.end(),
          "opening velocity-state composition test finds a solid trace");
    if (solidDof != openingVelocityMetric.dofs.end()) {
        invalidNormal[solidDof->dofIndex] += 0.1;
        invalidRelative[solidDof->dofIndex] += 0.1;
        expectRejected(
            [&] {
                static_cast<void>(
                    buildPlanarPressureRegionFragmentOpeningVelocityState(
                        openingVelocityMetric, invalidNormal,
                        manufacturedMaterial, invalidRelative,
                        manufacturedDensity));
            },
            "opening velocity state rejects relative flow through retained solid area");
    }
    auto fullOpeningDefinitions = definitions;
    fullOpeningDefinitions[0].areaSquareMeters = wall->areaSquareMeters;
    const auto fullOpenings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, fullOpeningDefinitions);
    const auto fullOpeningVelocityMetric =
        buildPlanarPressureRegionFragmentOpeningVelocityMetric(
            geometry, sweep, fragments, topology, velocityMetric,
            fullOpeningDefinitions, fullOpenings);
    check(fullOpeningVelocityMetric.solidWallTraceDofCount
                  + 2
              == openingVelocityMetric.solidWallTraceDofCount
              && fullOpeningVelocityMetric.openingPatchDofCount == 1
              && fullOpeningVelocityMetric.dofs.size() + 1
                  == velocityMetric.dofs.size(),
          "fully open wall replaces both solid traces with one shared aperture degree");
    checkNear(
        fullOpeningVelocityMetric.totalSolidWallAreaSquareMeters
            + fullOpeningVelocityMetric.totalOpeningAreaSquareMeters,
        fullOpeningVelocityMetric.totalPressureWallAreaSquareMeters,
        2.0e-15,
        "fully open velocity metric retains exact wall-area closure");
    checkNear(
        std::max({
            std::abs(fullOpeningVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.x),
            std::abs(fullOpeningVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.y),
            std::abs(fullOpeningVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.z),
            fullOpeningVelocityMetric
                .maximumAbsoluteFragmentVolumeClosureResidualCubicMeters,
            fullOpeningVelocityMetric
                .maximumAbsoluteComponentVolumeClosureResidualCubicMeters}),
        0.0, 2.0e-13,
        "fully open velocity metric loses no fragment or component inertia");
    std::vector<double> fullOpeningNormal(
        fullOpeningVelocityMetric.dofs.size(), 0.0);
    std::vector<double> fullOpeningMaterial(
        fullOpeningVelocityMetric.dofs.size(), 0.0);
    std::vector<double> fullOpeningRelative(
        fullOpeningVelocityMetric.dofs.size(), 0.0);
    for (const auto& dof : fullOpeningVelocityMetric.dofs) {
        const double component = dof.axis == GridFaceAxis::X
            ? manufacturedVelocity.x
            : dof.axis == GridFaceAxis::Y
            ? manufacturedVelocity.y : manufacturedVelocity.z;
        fullOpeningNormal[dof.dofIndex] = component;
        if (dof.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            fullOpeningRelative[dof.dofIndex] = component;
        } else {
            fullOpeningMaterial[dof.dofIndex] = component;
        }
    }
    const auto fullOpeningUniformState =
        buildPlanarPressureRegionFragmentOpeningVelocityState(
            fullOpeningVelocityMetric, fullOpeningNormal,
            fullOpeningMaterial, fullOpeningRelative,
            manufacturedDensity);
    checkNear(fullOpeningUniformState.physicalMassKilograms,
              manufacturedOpeningVelocityState.physicalMassKilograms,
              2.0e-13,
              "fully open velocity state preserves physical mass");
    checkNear(
        fullOpeningUniformState.momentumKilogramMetersPerSecond.x,
        manufacturedOpeningVelocityState
            .momentumKilogramMetersPerSecond.x,
        4.0e-13,
        "fully open velocity state preserves uniform momentum");
    checkNear(fullOpeningUniformState.staggeringKineticEnergyJoules,
              0.0, 8.0e-13,
              "fully open velocity state preserves the uniform free stream");
    const double requiredFlow =
        -volumeRates.components[wall->minusComponentIndex]
             .geometryVolumeChangeRateCubicMetersPerSecond;
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        initialSamples{{100, requiredFlow / patchArea}};
    auto initialFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, definitions, openings,
        initialSamples);
    std::vector<double> initialVelocity(topology.links.size(), 0.0);
    std::vector<double> initialPressure(base.rows.size(), 0.0);
    PlanarPressureRegionFragmentPressureProjectionSettings seedSettings;
    seedSettings.densityKgPerCubicMeter = 1.2;
    seedSettings.timeStepSeconds = 0.5;
    seedSettings.absoluteContinuityToleranceCubicMetersPerSecond = 2.0e-11;
    seedSettings.relativeContinuityTolerance = 1.0e-10;
    seedSettings.pressureSolve.absoluteResidualTolerancePascalsMeters =
        1.0e-13;
    seedSettings.pressureSolve.relativeResidualTolerance = 0.0;
    seedSettings.pressureSolve.maximumIterations = 300;
    const auto seeded =
        projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
            base, geometry, sweep, fragments, topology, volumeRates,
            definitions, openings, initialFlux, initialSamples,
            initialVelocity, initialPressure, seedSettings);
    check(seeded.accepted,
          "resisted opening pressure step starts from a compatible state");
    if (!seeded.accepted) return;

    PlanarPressureRegionFragmentOpeningPressureStepSettings settings;
    settings.projection.densityKgPerCubicMeter = 1.2;
    settings.projection.timeStepSeconds = 0.5;
    settings.projection.absoluteContinuityToleranceCubicMetersPerSecond =
        2.0e-11;
    settings.projection.relativeContinuityTolerance = 1.0e-10;
    settings.projection.absoluteMomentumResidualToleranceKilogramMetersPerSecond =
        2.0e-12;
    settings.projection.relativeMomentumResidualTolerance = 1.0e-10;
    settings.projection.absoluteEnergyResidualToleranceJoules = 2.0e-11;
    settings.projection.relativeEnergyResidualTolerance = 1.0e-10;
    settings.projection.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-13;
    settings.projection.pressureSolve.relativeResidualTolerance = 0.0;
    settings.projection.pressureSolve.maximumIterations = 300;
    const std::vector<
        PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistance{{100, {8.0, 3.0}}};
    auto velocity = initialVelocity;
    auto samples = initialSamples;
    auto flux = initialFlux;
    auto pressure = initialPressure;
    auto repeatedVelocity = initialVelocity;
    auto repeatedSamples = initialSamples;
    auto repeatedFlux = initialFlux;
    auto repeatedPressure = initialPressure;
    const auto diagnostics =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance, velocity, samples,
            flux, pressure, settings);
    const auto repeated =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance,
            repeatedVelocity, repeatedSamples, repeatedFlux,
            repeatedPressure, settings);
    check(diagnostics == repeated && velocity == repeatedVelocity
              && samples == repeatedSamples && flux == repeatedFlux
              && pressure == repeatedPressure
              && diagnostics.accepted && diagnostics.finite
              && diagnostics.energyAccepted
              && diagnostics.resistance.accepted
              && diagnostics.projection.accepted
              && diagnostics.sourceOpeningFluxFingerprint
                  == initialFlux.fingerprint
              && diagnostics.resultOpeningFluxFingerprint == flux.fingerprint,
          "resisted opening pressure step is deterministic and transactional");
    check(diagnostics.resistance.dissipatedEnergyJoules > 0.0
              && diagnostics.resistance.patches[0]
                     .plugFlow.velocityAfterMetersPerSecond
                  < initialSamples[0]
                        .relativeNormalVelocityMetersPerSecond
              && diagnostics.projection
                     .predictedContinuityResidualMaximumCubicMetersPerSecond
                  > 0.0
              && diagnostics.projection
                     .correctedContinuityResidualMaximumCubicMetersPerSecond
                  <= diagnostics.projection
                         .continuityToleranceCubicMetersPerSecond,
          "resistance creates a deficit that pressure restores to continuity");
    checkNear(samples[0].relativeNormalVelocityMetersPerSecond,
              initialSamples[0].relativeNormalVelocityMetersPerSecond,
              4.0e-11,
              "resisted pressure step restores the required intake flow");
    checkNear(flux.patches[0].relativeVolumeFlowRateCubicMetersPerSecond,
              requiredFlow, 4.0e-11,
              "resisted pressure step publishes the required aperture flux");
    check(diagnostics.dissipatedEnergyJoules > 0.0
              && diagnostics.correctionKineticEnergyJoules > 0.0
              && diagnostics.authoredPressureWorkJoules == 0.0
              && std::abs(diagnostics.energyResidualJoules)
                  <= diagnostics.energyToleranceJoules,
          "resisted pressure step closes loss plus projection energy");
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            check(velocity[link.linkIndex] == 0.0,
                  "resisted pressure step retains solid-wall zero flow");
        }
    }

    const auto acceptedState =
        capturePlanarPressureRegionFragmentOpeningAcceptedState(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance, diagnostics,
            velocity, samples, flux, pressure, settings);
    const auto repeatedAcceptedState =
        capturePlanarPressureRegionFragmentOpeningAcceptedState(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance, diagnostics,
            velocity, samples, flux, pressure, settings);
    check(acceptedState == repeatedAcceptedState
              && acceptedState.version
                  == planarPressureRegionFragmentOpeningAcceptedStateVersion
              && acceptedState.fingerprint != 0
              && acceptedState.accepted
              && acceptedState.sourcePressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && acceptedState.sourceBasePressureOperatorFingerprint
                  == base.fingerprint
              && acceptedState.sourceOpeningFingerprint
                  == openings.fingerprint
              && acceptedState.sourceFragmentFingerprint
                  == fragments.fingerprint
              && acceptedState.sourceTopologyFingerprint
                  == topology.fingerprint
              && acceptedState.sourceVolumeRateFingerprint
                  == volumeRates.fingerprint
              && acceptedState.sourceOpeningFluxFingerprint
                  == initialFlux.fingerprint
              && acceptedState.resultOpeningFluxFingerprint
                  == flux.fingerprint
              && acceptedState.settings == settings
              && acceptedState
                     .orientedTopologyLinkVelocityMetersPerSecond
                  == velocity
              && acceptedState.openingVelocitySamples == samples
              && acceptedState.openingFlux == flux
              && acceptedState.pressureCorrectionPascals == pressure
              && acceptedState.kineticEnergyAfterJoules
                  == diagnostics.kineticEnergyAfterJoules
              && acceptedState.correctedContinuityResidualMaximumCubicMetersPerSecond
                  == diagnostics.projection
                         .correctedContinuityResidualMaximumCubicMetersPerSecond
              && acceptedState.maximumAbsoluteCorrectionVolumeMeanPascals
                  == diagnostics.projection.pressureSolve
                         .maximumAbsoluteCorrectionVolumeMeanPascals
              && acceptedState.pressureSolveIterationCount
                  == diagnostics.projection.pressureSolve.iterationCount,
          "opening accepted state atomically retains the complete pressure-step endpoint");
    validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
        acceptedState);
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        acceptedState, pressureOperator, base, geometry, sweep, fragments,
        topology, volumeRates, definitions, openings, resistance);
    const auto mappedOpeningVelocityState =
        capturePlanarPressureRegionFragmentOpeningVelocityState(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, velocityMetric, openingVelocityMetric);
    const auto repeatedMappedOpeningVelocityState =
        capturePlanarPressureRegionFragmentOpeningVelocityState(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, velocityMetric, openingVelocityMetric);
    check(mappedOpeningVelocityState == repeatedMappedOpeningVelocityState
              && mappedOpeningVelocityState.mappedFromAcceptedEndpoint
              && mappedOpeningVelocityState.sourceMetricFingerprint
                  == openingVelocityMetric.fingerprint
              && mappedOpeningVelocityState
                     .sourceAcceptedStateFingerprint
                  == acceptedState.fingerprint
              && mappedOpeningVelocityState.sourceVolumeRateFingerprint
                  == volumeRates.fingerprint
              && mappedOpeningVelocityState.densityKgPerCubicMeter
                  == settings.projection.densityKgPerCubicMeter
              && mappedOpeningVelocityState
                     .maximumAbsoluteVelocityCompositionResidualMetersPerSecond
                  == 0.0
              && mappedOpeningVelocityState.staggeringKineticEnergyJoules
                  >= -2.0e-13,
          "opening velocity state maps the accepted pressure endpoint into absolute fragment momentum");
    const auto mappedApertureSample = std::ranges::find_if(
        mappedOpeningVelocityState.samples,
        [](const auto& sample) {
            return sample.kind
                == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                    OpeningPatch;
        });
    check(mappedApertureSample
                  != mappedOpeningVelocityState.samples.end(),
          "mapped opening velocity state retains its aperture sample");
    if (mappedApertureSample
        != mappedOpeningVelocityState.samples.end()) {
        const double expectedMaterialVelocity = volumeRates
            .fragments[wall->minusFragmentIndex]
            .upperBoundaryVelocityMetersPerSecond;
        check(mappedApertureSample
                      ->materialNormalVelocityMetersPerSecond
                  == expectedMaterialVelocity
                  && mappedApertureSample
                         ->relativeNormalVelocityMetersPerSecond
                      == acceptedState.openingFlux.patches[0]
                             .relativeNormalVelocityMetersPerSecond
                  && mappedApertureSample->normalVelocityMetersPerSecond
                      == mappedApertureSample
                                 ->materialNormalVelocityMetersPerSecond
                          + mappedApertureSample
                                 ->relativeNormalVelocityMetersPerSecond,
              "mapped aperture velocity is material motion plus accepted relative flow");
    }
    for (const auto& sample : mappedOpeningVelocityState.samples) {
        if (sample.kind
            == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                SharedRegionGrid) {
            check(sample.materialNormalVelocityMetersPerSecond == 0.0,
                  "fixed-grid regional flow has zero material velocity");
        } else if (sample.kind
                   != PlanarPressureRegionFragmentOpeningVelocityDofKind::
                       OpeningPatch) {
            check(sample.relativeNormalVelocityMetersPerSecond == 0.0,
                  "retained solid traces have zero relative flow");
        }
    }
    validatePlanarPressureRegionFragmentOpeningAcceptedVelocityState(
        mappedOpeningVelocityState, acceptedState, pressureOperator, base,
        geometry, sweep, fragments, topology, volumeRates, definitions,
        openings, resistance, velocityMetric, openingVelocityMetric);
    openingVelocityStateLimits = {};
    openingVelocityStateLimits.maximumWorkingBytes =
        mappedOpeningVelocityState.workingStorageBytes - 1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningAcceptedVelocityState(
                mappedOpeningVelocityState, acceptedState,
                pressureOperator, base, geometry, sweep, fragments,
                topology, volumeRates, definitions, openings, resistance,
                velocityMetric, openingVelocityMetric,
                openingVelocityStateLimits);
        },
        "mapped opening velocity state enforces validation working storage");
    auto corruptAcceptedState = acceptedState;
    corruptAcceptedState.openingFlux.patches[0]
        .relativeNormalVelocityMetersPerSecond += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
                corruptAcceptedState);
        },
        "opening accepted state rejects nested flux corruption");
    corruptAcceptedState = acceptedState;
    const auto sharedLink = std::ranges::find_if(
        topology.links,
        [](const auto& link) {
            return link.kind
                == PlanarPressureRegionFragmentFaceKind::SameRegionGrid;
        });
    check(sharedLink != topology.links.end(),
          "opening accepted-state corruption test finds a shared link");
    if (sharedLink != topology.links.end()) {
        corruptAcceptedState
            .orientedTopologyLinkVelocityMetersPerSecond[
                sharedLink->linkIndex] += 0.1;
        expectRejected(
            [&] {
                validatePlanarPressureRegionFragmentOpeningAcceptedStateIntegrity(
                    corruptAcceptedState);
            },
            "opening accepted state rejects velocity corruption");
    }
    auto foreignResistance = resistance;
    foreignResistance[0].resistance.linearPascalSecondsPerMeter += 1.0;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningAcceptedState(
                acceptedState, pressureOperator, base, geometry, sweep,
                fragments, topology, volumeRates, definitions, openings,
                foreignResistance);
        },
        "opening accepted state rejects foreign resistance provenance");
    auto acceptedLimits =
        PlanarPressureRegionFragmentOpeningAcceptedStateLimits{};
    acceptedLimits.maximumOwnedBytes =
        acceptedState.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningAcceptedState(
                acceptedState, pressureOperator, base, geometry, sweep,
                fragments, topology, volumeRates, definitions, openings,
                resistance, acceptedLimits);
        },
        "opening accepted state enforces its aggregate owned-byte limit");

    auto nextLayers = current;
    nextLayers[0].physicalPlaneCoordinateMeters -= 0.05;
    nextLayers[1].physicalPlaneCoordinateMeters += 0.05;
    const auto nextSweep = makePlanarPressureRegionSweepLedger(
        geometry, current, nextLayers, 0.5);
    const auto nextFragments = buildPlanarPressureRegionFragments(
        geometry, nextSweep);
    const auto nextTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, nextSweep, nextFragments);
    const auto nextBase =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, nextSweep, nextFragments, nextTopology);
    const auto nextVolumeRates =
        buildPlanarPressureRegionFragmentVolumeRates(
            geometry, nextSweep, nextFragments, nextTopology);
    const auto nextOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, nextSweep, nextFragments, nextTopology,
            definitions);
    const auto nextPressureOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            nextBase, geometry, nextSweep, nextFragments, nextTopology,
            definitions, nextOpenings);
    const auto nextVelocityMetric =
        buildPlanarPressureRegionFragmentVelocityMetric(
            geometry, nextSweep, nextFragments, nextTopology);
    const auto nextOpeningVelocityMetric =
        buildPlanarPressureRegionFragmentOpeningVelocityMetric(
            geometry, nextSweep, nextFragments, nextTopology,
            nextVelocityMetric, definitions, nextOpenings);
    check(nextOpeningVelocityMetric.dofs.size()
                  == openingVelocityMetric.dofs.size()
              && nextOpeningVelocityMetric.fragments.size()
                  == openingVelocityMetric.fragments.size()
              && nextOpeningVelocityMetric.components.size()
                  == openingVelocityMetric.components.size(),
          "opening velocity metric retains topology-stable entity coverage");
    bool openingMetricChanged = false;
    for (std::size_t index = 0;
         index < openingVelocityMetric.dofs.size(); ++index) {
        const auto& previousDof = openingVelocityMetric.dofs[index];
        const auto& nextDof = nextOpeningVelocityMetric.dofs[index];
        check(previousDof.stableId == nextDof.stableId
                  && previousDof.kind == nextDof.kind
                  && previousDof.sourceFaceLinkStableId
                      == nextDof.sourceFaceLinkStableId
                  && previousDof.sourceOpeningPatchStableId
                      == nextDof.sourceOpeningPatchStableId,
              "opening velocity metric preserves stable DOF identity across motion");
        openingMetricChanged = openingMetricChanged
            || previousDof.dualVolumeCubicMeters
                != nextDof.dualVolumeCubicMeters;
    }
    check(openingMetricChanged,
          "opening velocity metric updates moving half-volumes while retaining exact global closure");
    checkNear(
        std::max({
            std::abs(nextOpeningVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.x),
            std::abs(nextOpeningVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.y),
            std::abs(nextOpeningVelocityMetric
                         .domainVolumeClosureResidualByAxisCubicMeters.z)}),
        0.0, 2.0e-13,
        "moving opening velocity metric retains global volume closure");
    validatePlanarPressureRegionFragmentOpeningVelocityMetric(
        nextOpeningVelocityMetric, geometry, nextSweep, nextFragments,
        nextTopology, nextVelocityMetric, definitions, nextOpenings);
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningVelocityMetric(
                openingVelocityMetric, geometry, nextSweep, nextFragments,
                nextTopology, nextVelocityMetric, definitions,
                nextOpenings);
        },
        "opening velocity metric rejects a foreign moving endpoint");
    openingVelocityMetricLimits = {};
    openingVelocityMetricLimits.maximumWorkingBytes =
        openingVelocityMetric.workingStorageBytes - 1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningVelocityMetric(
                openingVelocityMetric, geometry, sweep, fragments,
                topology, velocityMetric, definitions, openings,
                openingVelocityMetricLimits);
        },
        "opening velocity metric enforces validation working storage");
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningAcceptedVelocityState(
                mappedOpeningVelocityState, acceptedState,
                nextPressureOperator, nextBase, geometry, nextSweep,
                nextFragments, nextTopology, nextVolumeRates, definitions,
                nextOpenings, resistance, nextVelocityMetric,
                nextOpeningVelocityMetric);
        },
        "mapped opening velocity state rejects a foreign accepted endpoint geometry");
    const auto continuation =
        buildPlanarPressureRegionFragmentOpeningContinuation(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, nextPressureOperator, nextBase, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings);
    const auto repeatedContinuation =
        buildPlanarPressureRegionFragmentOpeningContinuation(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, nextPressureOperator, nextBase, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings);
    check(continuation == repeatedContinuation
              && continuation.version
                  == planarPressureRegionFragmentOpeningContinuationVersion
              && continuation.fingerprint != 0
              && continuation.sourceAcceptedStateFingerprint
                  == acceptedState.fingerprint
              && continuation.previousPressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && continuation.currentPressureOperatorFingerprint
                  == nextPressureOperator.fingerprint
              && continuation.currentTopologyFingerprint
                  == nextTopology.fingerprint
              && continuation.currentVolumeRateFingerprint
                  == nextVolumeRates.fingerprint
              && continuation.currentOpeningFluxFingerprint
                  == continuation.openingFlux.fingerprint
              && continuation.topologyLinkCount
                  == nextTopology.links.size()
              && continuation.openingPatchCount
                  == nextOpenings.patches.size()
              && continuation.pressureCorrectionCount
                  == nextFragments.fragments.size()
              && continuation.connectedComponentCount
                  == nextPressureOperator.components.size(),
          "opening continuation maps one consecutive stable epoch deterministically");
    for (const auto& nextLink : nextTopology.links) {
        const auto previousLink = std::ranges::find(
            topology.links, nextLink.stableId,
            &PlanarPressureRegionFragmentFaceLink::stableId);
        check(previousLink != topology.links.end()
                  && continuation
                         .orientedTopologyLinkVelocityMetersPerSecond[
                             nextLink.linkIndex]
                      == acceptedState
                             .orientedTopologyLinkVelocityMetersPerSecond[
                                 previousLink->linkIndex],
              "opening continuation preserves every stable link velocity");
    }
    check(continuation.openingVelocitySamples
                  == acceptedState.openingVelocitySamples
              && continuation.openingFlux.sourceOpeningFingerprint
                  == nextOpenings.fingerprint,
          "opening continuation preserves aperture velocity and rebuilds current flux");
    for (const auto& component : nextPressureOperator.components) {
        checkNear(
            fragmentCorrectionVolumeMean(
                nextPressureOperator, nextFragments, component,
                continuation.pressureCorrectionPascals),
            0.0, 2.0e-12,
            "opening continuation applies the current connected-volume gauge");
    }
    validatePlanarPressureRegionFragmentOpeningContinuationIntegrity(
        continuation);
    validatePlanarPressureRegionFragmentOpeningContinuation(
        continuation, acceptedState, pressureOperator, base, geometry,
        sweep, fragments, topology, volumeRates, definitions, openings,
        resistance, nextPressureOperator, nextBase, nextSweep,
        nextFragments, nextTopology, nextVolumeRates, definitions,
        nextOpenings);

    const auto continuationMomentumAudit =
        auditPlanarPressureRegionFragmentOpeningContinuationMomentum(
            continuation, acceptedState, pressureOperator, base, geometry,
            sweep, fragments, topology, volumeRates, definitions, openings,
            resistance, nextPressureOperator, nextBase, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings);
    const auto repeatedContinuationMomentumAudit =
        auditPlanarPressureRegionFragmentOpeningContinuationMomentum(
            continuation, acceptedState, pressureOperator, base, geometry,
            sweep, fragments, topology, volumeRates, definitions, openings,
            resistance, nextPressureOperator, nextBase, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings);
    check(continuationMomentumAudit == repeatedContinuationMomentumAudit
              && continuationMomentumAudit.version
                  == planarPressureRegionFragmentOpeningContinuationMomentumAuditVersion
              && continuationMomentumAudit.fingerprint != 0
              && continuationMomentumAudit.audited
              && continuationMomentumAudit.sourceContinuationFingerprint
                  == continuation.fingerprint
              && continuationMomentumAudit.sourceAcceptedStateFingerprint
                  == acceptedState.fingerprint
              && continuationMomentumAudit.previousTopologyFingerprint
                  == topology.fingerprint
              && continuationMomentumAudit.currentTopologyFingerprint
                  == nextTopology.fingerprint
              && continuationMomentumAudit.samples.size()
                  == topology.sameRegionGridLinkCount
                      + openings.patches.size()
              && continuationMomentumAudit.sameRegionGridSampleCount
                  == topology.sameRegionGridLinkCount
              && continuationMomentumAudit.openingPatchSampleCount
                  == openings.patches.size()
              && continuationMomentumAudit.metricChanged
              && continuationMomentumAudit.metricChangedSampleCount == 33
              && continuationMomentumAudit.warmCarryChangesMomentum
              && continuationMomentumAudit.momentumChangedSampleCount == 33
              && continuationMomentumAudit
                     .maximumAbsoluteDualVolumeChangeCubicMeters > 0.0
              && continuationMomentumAudit
                     .maximumAbsoluteMomentumChangeKilogramMetersPerSecond
                  > 0.0
              && continuationMomentumAudit.ownedStorageBytes > 0
              && continuationMomentumAudit.workingStorageBytes > 0,
          "opening continuation momentum audit exposes metric-induced warm-carry drift");
    checkNear(
        continuationMomentumAudit
            .maximumAbsoluteDualVolumeChangeCubicMeters,
        0.1, 2.0e-15,
        "opening continuation momentum audit retains maximum metric change");
    checkNear(
        continuationMomentumAudit.maximumRelativeDualVolumeChange,
        0.5, 3.0e-14,
        "opening continuation momentum audit retains relative metric change");
    checkNear(
        continuationMomentumAudit
            .maximumAbsoluteMomentumChangeKilogramMetersPerSecond,
        0.048676001341721321, 2.0e-14,
        "opening continuation momentum audit retains maximum local momentum drift");
    checkNear(
        continuationMomentumAudit.kineticEnergyChangeJoules,
        0.013419632753952637, 2.0e-14,
        "opening continuation momentum audit retains aggregate energy drift");
    for (const auto& sample : continuationMomentumAudit.samples) {
        check(sample.previousMomentumKilogramMetersPerSecond
                      == sample.previousMassKilograms
                          * sample.carriedVelocityMetersPerSecond
                  && sample.carriedMomentumKilogramMetersPerSecond
                      == sample.currentMassKilograms
                          * sample.carriedVelocityMetersPerSecond
                  && sample.momentumChangeKilogramMetersPerSecond
                      == sample.carriedMomentumKilogramMetersPerSecond
                          - sample.previousMomentumKilogramMetersPerSecond
                  && sample.kineticEnergyChangeJoules
                      == sample.carriedKineticEnergyJoules
                          - sample.previousKineticEnergyJoules,
              "opening continuation momentum audit retains exact per-degree arithmetic");
    }
    validatePlanarPressureRegionFragmentOpeningContinuationMomentumAuditIntegrity(
        continuationMomentumAudit);
    validatePlanarPressureRegionFragmentOpeningContinuationMomentumAudit(
        continuationMomentumAudit, continuation, acceptedState,
        pressureOperator, base, geometry, sweep, fragments, topology,
        volumeRates, definitions, openings, resistance,
        nextPressureOperator, nextBase, nextSweep, nextFragments,
        nextTopology, nextVolumeRates, definitions, nextOpenings);

    auto corruptMomentumAudit = continuationMomentumAudit;
    corruptMomentumAudit.samples[0].currentMassKilograms += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningContinuationMomentumAuditIntegrity(
                corruptMomentumAudit);
        },
        "opening continuation momentum audit rejects sample corruption");
    auto momentumAuditLimits =
        PlanarPressureRegionFragmentOpeningContinuationMomentumAuditLimits{};
    momentumAuditLimits.maximumOwnedBytes =
        continuationMomentumAudit.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                auditPlanarPressureRegionFragmentOpeningContinuationMomentum(
                    continuation, acceptedState, pressureOperator, base,
                    geometry, sweep, fragments, topology, volumeRates,
                    definitions, openings, resistance,
                    nextPressureOperator, nextBase, nextSweep,
                    nextFragments, nextTopology, nextVolumeRates,
                    definitions, nextOpenings, momentumAuditLimits));
        },
        "opening continuation momentum audit enforces its owned-storage limit");

    const auto stationarySweep = makePlanarPressureRegionSweepLedger(
        geometry, current, current, 0.5);
    const auto stationaryFragments = buildPlanarPressureRegionFragments(
        geometry, stationarySweep);
    const auto stationaryTopology =
        buildPlanarPressureRegionFragmentTopology(
            geometry, stationarySweep, stationaryFragments);
    const auto stationaryBase =
        buildPlanarPressureRegionFragmentPressureOperator(
            geometry, stationarySweep, stationaryFragments,
            stationaryTopology);
    const auto stationaryVolumeRates =
        buildPlanarPressureRegionFragmentVolumeRates(
            geometry, stationarySweep, stationaryFragments,
            stationaryTopology);
    const auto stationaryOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, stationarySweep, stationaryFragments,
            stationaryTopology, definitions);
    const auto stationaryPressureOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            stationaryBase, geometry, stationarySweep,
            stationaryFragments, stationaryTopology, definitions,
            stationaryOpenings);
    const auto stationaryContinuation =
        buildPlanarPressureRegionFragmentOpeningContinuation(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, stationaryPressureOperator, stationaryBase,
            stationarySweep, stationaryFragments, stationaryTopology,
            stationaryVolumeRates, definitions, stationaryOpenings);
    const auto stationaryMomentumAudit =
        auditPlanarPressureRegionFragmentOpeningContinuationMomentum(
            stationaryContinuation, acceptedState, pressureOperator, base,
            geometry, sweep, fragments, topology, volumeRates, definitions,
            openings, resistance, stationaryPressureOperator,
            stationaryBase, stationarySweep, stationaryFragments,
            stationaryTopology, stationaryVolumeRates, definitions,
            stationaryOpenings);
    check(!stationaryMomentumAudit.metricChanged
              && !stationaryMomentumAudit.warmCarryChangesMomentum
              && stationaryMomentumAudit.metricChangedSampleCount == 0
              && stationaryMomentumAudit.momentumChangedSampleCount == 0
              && stationaryMomentumAudit.previousMassByAxisKilograms
                  == stationaryMomentumAudit.currentMassByAxisKilograms
              && stationaryMomentumAudit
                     .previousMomentumKilogramMetersPerSecond
                  == stationaryMomentumAudit
                         .carriedMomentumKilogramMetersPerSecond
              && stationaryMomentumAudit
                     .momentumChangeKilogramMetersPerSecond
                  == Vector3{}
              && stationaryMomentumAudit.previousKineticEnergyJoules
                  == stationaryMomentumAudit.carriedKineticEnergyJoules
              && stationaryMomentumAudit.kineticEnergyChangeJoules == 0.0
              && stationaryMomentumAudit
                     .maximumAbsoluteDualVolumeChangeCubicMeters == 0.0,
          "opening continuation momentum audit has an exact stationary-geometry zero oracle");
    validatePlanarPressureRegionFragmentOpeningContinuationMomentumAudit(
        stationaryMomentumAudit, stationaryContinuation, acceptedState,
        pressureOperator, base, geometry, sweep, fragments, topology,
        volumeRates, definitions, openings, resistance,
        stationaryPressureOperator, stationaryBase, stationarySweep,
        stationaryFragments, stationaryTopology, stationaryVolumeRates,
        definitions, stationaryOpenings);
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningContinuationMomentumAudit(
                continuationMomentumAudit, stationaryContinuation,
                acceptedState, pressureOperator, base, geometry, sweep,
                fragments, topology, volumeRates, definitions, openings,
                resistance, stationaryPressureOperator, stationaryBase,
                stationarySweep, stationaryFragments, stationaryTopology,
                stationaryVolumeRates, definitions, stationaryOpenings);
        },
        "opening continuation momentum audit rejects a foreign current epoch");
    momentumAuditLimits = {};
    momentumAuditLimits.maximumWorkingBytes =
        continuationMomentumAudit.workingStorageBytes - 1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningContinuationMomentumAudit(
                continuationMomentumAudit, continuation, acceptedState,
                pressureOperator, base, geometry, sweep, fragments,
                topology, volumeRates, definitions, openings, resistance,
                nextPressureOperator, nextBase, nextSweep, nextFragments,
                nextTopology, nextVolumeRates, definitions, nextOpenings,
                momentumAuditLimits);
        },
        "opening continuation momentum audit enforces validation working storage");

    auto corruptContinuation = continuation;
    corruptContinuation.pressureCorrectionPascals[0] += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningContinuationIntegrity(
                corruptContinuation);
        },
        "opening continuation rejects pressure corruption");
    auto continuationLimits =
        PlanarPressureRegionFragmentOpeningContinuationLimits{};
    continuationLimits.maximumWorkingBytes =
        continuation.workingStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                buildPlanarPressureRegionFragmentOpeningContinuation(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, nextPressureOperator, nextBase, nextSweep,
                    nextFragments, nextTopology, nextVolumeRates,
                    definitions, nextOpenings, continuationLimits));
        },
        "opening continuation enforces its working-storage limit");
    continuationLimits = {};
    continuationLimits.maximumOwnedBytes =
        continuation.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                buildPlanarPressureRegionFragmentOpeningContinuation(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, nextPressureOperator, nextBase, nextSweep,
                    nextFragments, nextTopology, nextVolumeRates,
                    definitions, nextOpenings, continuationLimits));
        },
        "opening continuation enforces its owned-storage limit");
    expectRejected(
        [&] {
            static_cast<void>(
                buildPlanarPressureRegionFragmentOpeningContinuation(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, nextPressureOperator, nextBase, sweep,
                    fragments, topology, volumeRates, definitions, openings));
        },
        "opening continuation rejects nonconsecutive endpoint profiles");
    auto renamedDefinitions = definitions;
    renamedDefinitions[0].patchStableId += 1;
    const auto renamedOpenings =
        buildPlanarPressureRegionFragmentOpenings(
            geometry, nextSweep, nextFragments, nextTopology,
            renamedDefinitions);
    const auto renamedPressureOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            nextBase, geometry, nextSweep, nextFragments, nextTopology,
            renamedDefinitions, renamedOpenings);
    expectRejected(
        [&] {
            static_cast<void>(
                buildPlanarPressureRegionFragmentOpeningContinuation(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, renamedPressureOperator, nextBase,
                    nextSweep, nextFragments, nextTopology, nextVolumeRates,
                    renamedDefinitions, renamedOpenings));
        },
        "opening continuation rejects aperture appearance and retirement");

    const auto pressureEpoch =
        acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, nextPressureOperator, nextBase, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings, resistance, settings);
    const auto repeatedPressureEpoch =
        acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, nextPressureOperator, nextBase, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings, resistance, settings);
    check(pressureEpoch == repeatedPressureEpoch
              && pressureEpoch.version
                  == planarPressureRegionFragmentOpeningPressureEpochVersion
              && pressureEpoch.diagnostics.accepted
              && pressureEpoch.diagnostics.usedConsecutiveContinuation
              && pressureEpoch.diagnostics.failureStage
                  == PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
                      None
              && pressureEpoch.diagnostics.pressureStep.accepted
              && pressureEpoch.sourceAcceptedStateFingerprint
                  == acceptedState.fingerprint
              && pressureEpoch.continuationFingerprint
                  == continuation.fingerprint
              && pressureEpoch.continuationOpeningFluxFingerprint
                  == continuation.openingFlux.fingerprint
              && pressureEpoch.currentPressureOperatorFingerprint
                  == nextPressureOperator.fingerprint
              && pressureEpoch.currentVolumeRateFingerprint
                  == nextVolumeRates.fingerprint
              && pressureEpoch.acceptedState.accepted
              && pressureEpoch.acceptedState.sourcePressureOperatorFingerprint
                  == nextPressureOperator.fingerprint
              && pressureEpoch.acceptedState.sourceOpeningFluxFingerprint
                  == continuation.openingFlux.fingerprint
              && pressureEpoch.acceptedState.resultOpeningFluxFingerprint
                  == pressureEpoch.diagnostics.pressureStep
                         .resultOpeningFluxFingerprint,
          "opening pressure epoch atomically accepts one continued moving state");
    validatePlanarPressureRegionFragmentOpeningPressureEpochResultIntegrity(
        pressureEpoch);
    validatePlanarPressureRegionFragmentOpeningPressureEpochResult(
        pressureEpoch, acceptedState, pressureOperator, base, geometry,
        sweep, fragments, topology, volumeRates, definitions, openings,
        resistance, nextPressureOperator, nextBase, nextSweep,
        nextFragments, nextTopology, nextVolumeRates, definitions,
        nextOpenings, resistance, settings);

    PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError
        persistenceError;
    std::vector<std::uint8_t> persistentStateBytes;
    std::vector<std::uint8_t> repeatedPersistentStateBytes;
    check(serializePlanarPressureRegionFragmentOpeningAcceptedState(
              pressureEpoch.acceptedState, nextPressureOperator, nextBase,
              geometry, nextSweep, nextFragments, nextTopology,
              nextVolumeRates, definitions, nextOpenings, resistance,
              persistentStateBytes, &persistenceError)
              && persistenceError.code
                  == PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                      None
              && serializePlanarPressureRegionFragmentOpeningAcceptedState(
                  pressureEpoch.acceptedState, nextPressureOperator,
                  nextBase, geometry, nextSweep, nextFragments,
                  nextTopology, nextVolumeRates, definitions, nextOpenings,
                  resistance, repeatedPersistentStateBytes,
                  &persistenceError)
              && persistentStateBytes == repeatedPersistentStateBytes,
          "opening accepted-state persistence encodes deterministically");
    PlanarPressureRegionFragmentOpeningAcceptedState restoredPressureState;
    const bool restoredPressureStateAccepted =
        deserializePlanarPressureRegionFragmentOpeningAcceptedState(
            persistentStateBytes, nextPressureOperator, nextBase, geometry,
            nextSweep, nextFragments, nextTopology, nextVolumeRates,
            definitions, nextOpenings, resistance, restoredPressureState,
            &persistenceError);
    std::vector<std::uint8_t> restoredPersistentStateBytes;
    check(restoredPressureStateAccepted
              && persistenceError.code
                  == PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                      None
              && restoredPressureState == pressureEpoch.acceptedState
              && serializePlanarPressureRegionFragmentOpeningAcceptedState(
                  restoredPressureState, nextPressureOperator, nextBase,
                  geometry, nextSweep, nextFragments, nextTopology,
                  nextVolumeRates, definitions, nextOpenings, resistance,
                  restoredPersistentStateBytes, &persistenceError)
              && restoredPersistentStateBytes == persistentStateBytes,
          "opening accepted-state persistence rebuilds flux and round trips bit-exactly");

    const auto expectPersistenceRejected =
        [&](const std::vector<std::uint8_t>& candidateBytes,
            const PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode
                expectedCode,
            const PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits&
                persistenceLimits = {}) {
            auto retainedState = acceptedState;
            const auto before = retainedState;
            PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceError
                rejectedError;
            const bool decoded =
                deserializePlanarPressureRegionFragmentOpeningAcceptedState(
                    candidateBytes, nextPressureOperator, nextBase, geometry,
                    nextSweep, nextFragments, nextTopology, nextVolumeRates,
                    definitions, nextOpenings, resistance, retainedState,
                    &rejectedError, persistenceLimits);
            return !decoded && rejectedError.code == expectedCode
                && retainedState == before;
        };
    auto corruptPersistentStateBytes = persistentStateBytes;
    corruptPersistentStateBytes[0] ^= 0x01U;
    check(expectPersistenceRejected(
              corruptPersistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  InvalidMagic),
          "opening accepted-state persistence rejects foreign magic transactionally");
    corruptPersistentStateBytes = persistentStateBytes;
    corruptPersistentStateBytes[4] = 2;
    corruptPersistentStateBytes[5] = 0;
    check(expectPersistenceRejected(
              corruptPersistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  UnsupportedVersion),
          "opening accepted-state persistence rejects unsupported protocol transactionally");
    corruptPersistentStateBytes = persistentStateBytes;
    corruptPersistentStateBytes[6] = 1;
    check(expectPersistenceRejected(
              corruptPersistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  InvalidData),
          "opening accepted-state persistence rejects reserved envelope bits");
    corruptPersistentStateBytes = persistentStateBytes;
    corruptPersistentStateBytes.back() ^= 0x01U;
    check(expectPersistenceRejected(
              corruptPersistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  ChecksumMismatch),
          "opening accepted-state persistence detects payload corruption");
    corruptPersistentStateBytes = persistentStateBytes;
    constexpr std::size_t encodedStateFingerprintOffset = 40;
    corruptPersistentStateBytes[encodedStateFingerprintOffset] ^= 0x01U;
    refreshOpeningAcceptedStatePersistenceChecksum(
        corruptPersistentStateBytes);
    check(expectPersistenceRejected(
              corruptPersistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  InvalidData),
          "opening accepted-state persistence rejects recomputed-checksum state corruption");
    corruptPersistentStateBytes = persistentStateBytes;
    corruptPersistentStateBytes.pop_back();
    check(expectPersistenceRejected(
              corruptPersistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  Truncated),
          "opening accepted-state persistence rejects truncation transactionally");
    corruptPersistentStateBytes = persistentStateBytes;
    corruptPersistentStateBytes.push_back(0);
    check(expectPersistenceRejected(
              corruptPersistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  TrailingData),
          "opening accepted-state persistence rejects trailing bytes transactionally");
    auto persistenceLimits =
        PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceLimits{};
    persistenceLimits.maximumTopologyLinkVelocities =
        nextTopology.links.size() - 1;
    check(expectPersistenceRejected(
              persistentStateBytes,
              PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                  LimitExceeded,
              persistenceLimits),
          "opening accepted-state persistence enforces decode record limits");
    auto foreignDecodedState = acceptedState;
    const auto retainedForeignDecodedState = foreignDecodedState;
    check(!deserializePlanarPressureRegionFragmentOpeningAcceptedState(
              persistentStateBytes, pressureOperator, base, geometry, sweep,
              fragments, topology, volumeRates, definitions, openings,
              resistance, foreignDecodedState, &persistenceError)
              && persistenceError.code
                  == PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                      SourceMismatch
              && foreignDecodedState == retainedForeignDecodedState,
          "opening accepted-state persistence rejects a foreign source epoch transactionally");
    persistenceLimits = {};
    persistenceLimits.maximumEncodedBytes = persistentStateBytes.size() - 1;
    std::vector<std::uint8_t> retainedEncoding{1, 2, 3};
    const auto retainedEncodingBefore = retainedEncoding;
    check(!serializePlanarPressureRegionFragmentOpeningAcceptedState(
              pressureEpoch.acceptedState, nextPressureOperator, nextBase,
              geometry, nextSweep, nextFragments, nextTopology,
              nextVolumeRates, definitions, nextOpenings, resistance,
              retainedEncoding, &persistenceError, persistenceLimits)
              && persistenceError.code
                  == PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                      LimitExceeded
              && retainedEncoding == retainedEncodingBefore,
          "opening accepted-state persistence enforces encode limits transactionally");
    auto corruptPersistentState = pressureEpoch.acceptedState;
    corruptPersistentState.pressureCorrectionPascals[0] += 0.1;
    retainedEncoding = retainedEncodingBefore;
    check(!serializePlanarPressureRegionFragmentOpeningAcceptedState(
              corruptPersistentState, nextPressureOperator, nextBase,
              geometry, nextSweep, nextFragments, nextTopology,
              nextVolumeRates, definitions, nextOpenings, resistance,
              retainedEncoding, &persistenceError)
              && persistenceError.code
                  == PlanarPressureRegionFragmentOpeningAcceptedStatePersistenceErrorCode::
                      InvalidData
              && retainedEncoding == retainedEncodingBefore,
          "opening accepted-state persistence rejects corrupt source state transactionally");
    auto corruptPressureEpoch = pressureEpoch;
    corruptPressureEpoch.acceptedState.pressureCorrectionPascals[0] += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningPressureEpochResultIntegrity(
                corruptPressureEpoch);
        },
        "opening pressure epoch rejects accepted-endpoint corruption");
    auto foreignEpochResistance = resistance;
    foreignEpochResistance[0].resistance.linearPascalSecondsPerMeter += 1.0;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningPressureEpochResult(
                pressureEpoch, acceptedState, pressureOperator, base,
                geometry, sweep, fragments, topology, volumeRates,
                definitions, openings, resistance, nextPressureOperator,
                nextBase, nextSweep, nextFragments, nextTopology,
                nextVolumeRates, definitions, nextOpenings,
                foreignEpochResistance, settings);
        },
        "opening pressure epoch rejects foreign current resistance provenance");

    auto rejectedEpochSettings = settings;
    rejectedEpochSettings.projection.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-16;
    rejectedEpochSettings.projection.pressureSolve
        .relativeResidualTolerance = 0.0;
    rejectedEpochSettings.projection.pressureSolve.maximumIterations = 1;
    const auto retainedPreviousState = acceptedState;
    const auto rejectedPressureEpoch =
        acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance, nextPressureOperator, nextBase, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings, resistance, rejectedEpochSettings);
    check(!rejectedPressureEpoch.diagnostics.accepted
              && rejectedPressureEpoch.diagnostics.usedConsecutiveContinuation
              && rejectedPressureEpoch.diagnostics.failureStage
                  == PlanarPressureRegionFragmentOpeningPressureEpochFailureStage::
                      PressureProjection
              && rejectedPressureEpoch.diagnostics.pressureStep
                     .resistance.accepted
              && !rejectedPressureEpoch.diagnostics.pressureStep
                      .projection.accepted
              && !rejectedPressureEpoch.acceptedState.accepted
              && rejectedPressureEpoch.acceptedState.fingerprint == 0
              && acceptedState == retainedPreviousState,
          "rejected opening pressure epoch publishes no partial endpoint");
    validatePlanarPressureRegionFragmentOpeningPressureEpochResultIntegrity(
        rejectedPressureEpoch);
    validatePlanarPressureRegionFragmentOpeningPressureEpochResult(
        rejectedPressureEpoch, acceptedState, pressureOperator, base,
        geometry, sweep, fragments, topology, volumeRates, definitions,
        openings, resistance, nextPressureOperator, nextBase, nextSweep,
        nextFragments, nextTopology, nextVolumeRates, definitions,
        nextOpenings, resistance, rejectedEpochSettings);
    auto pressureEpochLimits =
        PlanarPressureRegionFragmentOpeningPressureEpochLimits{};
    pressureEpochLimits.maximumOwnedBytes =
        pressureEpoch.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, nextPressureOperator, nextBase, nextSweep,
                    nextFragments, nextTopology, nextVolumeRates,
                    definitions, nextOpenings, resistance, settings,
                    pressureEpochLimits));
        },
        "opening pressure epoch enforces its aggregate owned-storage limit");

    if (pressureEpoch.diagnostics.accepted
        && restoredPressureStateAccepted) {
        auto thirdLayers = nextLayers;
        thirdLayers[0].physicalPlaneCoordinateMeters -= 0.02;
        thirdLayers[1].physicalPlaneCoordinateMeters += 0.02;
        const auto thirdSweep = makePlanarPressureRegionSweepLedger(
            geometry, nextLayers, thirdLayers, 0.5);
        const auto thirdFragments = buildPlanarPressureRegionFragments(
            geometry, thirdSweep);
        const auto thirdTopology =
            buildPlanarPressureRegionFragmentTopology(
                geometry, thirdSweep, thirdFragments);
        const auto thirdBase =
            buildPlanarPressureRegionFragmentPressureOperator(
                geometry, thirdSweep, thirdFragments, thirdTopology);
        const auto thirdVolumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, thirdSweep, thirdFragments, thirdTopology);
        const auto thirdOpenings =
            buildPlanarPressureRegionFragmentOpenings(
                geometry, thirdSweep, thirdFragments, thirdTopology,
                definitions);
        const auto thirdPressureOperator =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                thirdBase, geometry, thirdSweep, thirdFragments,
                thirdTopology, definitions, thirdOpenings);
        const auto thirdPressureEpoch =
            acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
                pressureEpoch.acceptedState, nextPressureOperator,
                nextBase, geometry, nextSweep, nextFragments, nextTopology,
                nextVolumeRates, definitions, nextOpenings, resistance,
                thirdPressureOperator, thirdBase, thirdSweep,
                thirdFragments, thirdTopology, thirdVolumeRates,
                definitions, thirdOpenings, resistance, settings);
        const auto restoredThirdPressureEpoch =
            acceptPlanarPressureRegionFragmentOpeningPressureEpoch(
                restoredPressureState, nextPressureOperator, nextBase,
                geometry, nextSweep, nextFragments, nextTopology,
                nextVolumeRates, definitions, nextOpenings, resistance,
                thirdPressureOperator, thirdBase, thirdSweep,
                thirdFragments, thirdTopology, thirdVolumeRates,
                definitions, thirdOpenings, resistance, settings);
        check(thirdPressureEpoch == restoredThirdPressureEpoch
                  && thirdPressureEpoch.diagnostics.accepted
                  && thirdPressureEpoch.sourceAcceptedStateFingerprint
                      == pressureEpoch.acceptedState.fingerprint
                  && thirdPressureEpoch.acceptedState
                         .sourcePressureOperatorFingerprint
                      == thirdPressureOperator.fingerprint
                  && thirdPressureEpoch.acceptedState.fingerprint
                      != pressureEpoch.acceptedState.fingerprint,
              "restored opening state continues with bit-exact consecutive acceptance");
        validatePlanarPressureRegionFragmentOpeningPressureEpochResult(
            thirdPressureEpoch, pressureEpoch.acceptedState,
            nextPressureOperator, nextBase, geometry, nextSweep,
            nextFragments, nextTopology, nextVolumeRates, definitions,
            nextOpenings, resistance, thirdPressureOperator, thirdBase,
            thirdSweep, thirdFragments, thirdTopology, thirdVolumeRates,
            definitions, thirdOpenings, resistance, settings);
    }

    const auto openingPressureState =
        composePlanarPressureRegionFragmentOpeningPressureState(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance);
    const auto repeatedOpeningPressureState =
        composePlanarPressureRegionFragmentOpeningPressureState(
            acceptedState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance);
    check(openingPressureState == repeatedOpeningPressureState
              && openingPressureState.version
                  == planarPressureRegionFragmentOpeningPressureStateVersion
              && openingPressureState.fingerprint != 0
              && openingPressureState.accepted
              && openingPressureState.sourceAcceptedStateFingerprint
                  == acceptedState.fingerprint
              && openingPressureState.sourcePressureOperatorFingerprint
                  == pressureOperator.fingerprint
              && openingPressureState.sourceBasePressureOperatorFingerprint
                  == base.fingerprint
              && openingPressureState.sourceOpeningFingerprint
                  == openings.fingerprint
              && openingPressureState.sourceVolumeRateFingerprint
                  == volumeRates.fingerprint
              && !openingPressureState.staticGeometry
              && openingPressureState.usesMovingVolumeRates
              && openingPressureState.controls.size()
                  == fragments.fragments.size()
              && openingPressureState.walls.size()
                  == topology.pressureLayerWallLinkCount
              && openingPressureState.components.size()
                  == openings.connectedComponents.size()
              && openingPressureState.components.size() == 1
              && openingPressureState.correctionGeometryPressureWorkJoules
                  == acceptedState.geometryPressureWorkJoules,
          "opening pressure state composes the accepted connected-gauge correction deterministically");
    for (std::size_t index = 0;
         index < openingPressureState.controls.size(); ++index) {
        const auto& control = openingPressureState.controls[index];
        check(control.fragmentIndex == index
                  && control.authoredPressurePascals
                      == fragments.fragments[index].pressurePascals
                  && control.correctionPressurePascals
                      == acceptedState.pressureCorrectionPascals[index]
                  && control.totalPressurePascals
                      == control.authoredPressurePascals
                         + control.correctionPressurePascals,
              "opening pressure state retains authored and correction control pressure");
    }
    check(openingPressureState.maximumAbsoluteCorrectionGaugePascals
              <= 2.0e-12
              && std::abs(
                     openingPressureState.wallGeometryWorkResidualJoules)
                  <= openingPressureState
                         .settings.absoluteWorkResidualToleranceJoules
              && std::abs(
                     openingPressureState.pressureWorkSplitResidualJoules)
                  <= openingPressureState
                         .settings.absoluteWorkResidualToleranceJoules
              && openingPressureState.totalPressureImpulseOnSheetNewtonSeconds
                  == Vector3{
                      openingPressureState
                              .totalPressureForceOnSheetNewtons.x
                          * openingPressureState.timeStepSeconds,
                      openingPressureState
                              .totalPressureForceOnSheetNewtons.y
                          * openingPressureState.timeStepSeconds,
                      openingPressureState
                              .totalPressureForceOnSheetNewtons.z
                          * openingPressureState.timeStepSeconds},
          "opening pressure state closes connected gauges, wall work, and impulse");
    validatePlanarPressureRegionFragmentOpeningPressureStateIntegrity(
        openingPressureState);
    validatePlanarPressureRegionFragmentOpeningPressureState(
        openingPressureState, acceptedState, pressureOperator, base,
        geometry, sweep, fragments, topology, volumeRates, definitions,
        openings, resistance);
    const auto openingSurfaceLoads =
        capturePlanarPressureRegionFragmentSurfaceLoads(
            openingPressureState);
    const auto openingAdjustedSurfaceLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            openingSurfaceLoads, openingPressureState, geometry, sweep,
            fragments, topology, definitions, openings);
    check(openingSurfaceLoads.accepted
              && openingSurfaceLoads.sourcePressureStateFingerprint
                  == openingPressureState.fingerprint
              && openingSurfaceLoads.sourceTopologyFingerprint
                  == topology.fingerprint
              && !openingSurfaceLoads.staticGeometry
              && openingSurfaceLoads.usesMovingVolumeRates
              && openingSurfaceLoads.totalPressureImpulseOnSheetNewtonSeconds
                  == openingPressureState
                         .totalPressureImpulseOnSheetNewtonSeconds,
          "opening pressure state captures the full-wall load ledger with exact provenance");
    check(openingAdjustedSurfaceLoads.accepted
              && openingAdjustedSurfaceLoads.sourceSurfaceLoadFingerprint
                  == openingSurfaceLoads.fingerprint
              && openingAdjustedSurfaceLoads.sourcePressureStateFingerprint
                  == openingPressureState.fingerprint
              && openingAdjustedSurfaceLoads.sourceOpeningFingerprint
                  == openings.fingerprint
              && openingAdjustedSurfaceLoads
                         .maximumAbsoluteAreaPartitionResidualSquareMeters
                  <= 1.0e-12
              && openingAdjustedSurfaceLoads
                         .maximumAbsoluteForcePartitionResidualNewtons
                  <= 1.0e-12
              && openingAdjustedSurfaceLoads
                         .maximumAbsoluteImpulsePartitionResidualNewtonSeconds
                  <= 1.0e-12
              && openingAdjustedSurfaceLoads
                         .maximumAbsoluteMomentPartitionResidualNewtonMeters
                  <= 1.0e-12
              && std::abs(
                     openingAdjustedSurfaceLoads.workPartitionResidualJoules)
                  <= 1.0e-12,
          "opening pressure loads close removed-aperture and retained-solid ownership");
    validatePlanarPressureRegionFragmentSurfaceLoads(
        openingSurfaceLoads, openingPressureState);
    validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
        openingAdjustedSurfaceLoads, openingSurfaceLoads,
        openingPressureState, geometry, sweep, fragments, topology,
        definitions, openings);
    const auto openingLoadState =
        capturePlanarPressureRegionFragmentOpeningLoadState(
            acceptedState, openingPressureState, openingSurfaceLoads,
            openingAdjustedSurfaceLoads, pressureOperator, base, geometry,
            sweep, fragments, topology, volumeRates, definitions, openings,
            resistance);
    const auto repeatedOpeningLoadState =
        capturePlanarPressureRegionFragmentOpeningLoadState(
            acceptedState, openingPressureState, openingSurfaceLoads,
            openingAdjustedSurfaceLoads, pressureOperator, base, geometry,
            sweep, fragments, topology, volumeRates, definitions, openings,
            resistance);
    check(openingLoadState == repeatedOpeningLoadState
              && openingLoadState.version
                  == planarPressureRegionFragmentOpeningLoadStateVersion
              && openingLoadState.fingerprint != 0
              && openingLoadState.accepted
              && openingLoadState.sourceAcceptedStateFingerprint
                  == acceptedState.fingerprint
              && openingLoadState.sourcePressureStateFingerprint
                  == openingPressureState.fingerprint
              && openingLoadState.sourceSurfaceLoadFingerprint
                  == openingSurfaceLoads.fingerprint
              && openingLoadState.sourceOpeningSurfaceLoadFingerprint
                  == openingAdjustedSurfaceLoads.fingerprint
              && openingLoadState.acceptedFlow == acceptedState
              && openingLoadState.pressure == openingPressureState
              && openingLoadState.surfaceLoads == openingSurfaceLoads
              && openingLoadState.openingSurfaceLoads
                  == openingAdjustedSurfaceLoads
              && openingLoadState.ownedStorageBytes
                  == acceptedState.ownedStorageBytes
                      + openingPressureState.ownedStorageBytes
                      + openingSurfaceLoads.ownedStorageBytes
                      + openingAdjustedSurfaceLoads.ownedStorageBytes,
          "opening load state atomically owns the complete aperture endpoint");
    check(openingLoadState.fullWallAreaSquareMeters
              == openingAdjustedSurfaceLoads.totalWallAreaSquareMeters
              && openingLoadState.openingAreaSquareMeters
                  == openingAdjustedSurfaceLoads.totalOpeningAreaSquareMeters
              && openingLoadState.solidAreaSquareMeters
                  == openingAdjustedSurfaceLoads.totalSolidAreaSquareMeters
              && openingLoadState.fullWallPressureForceOnSheetNewtons
                  == openingAdjustedSurfaceLoads
                         .sourceTotalPressureForceOnSheetNewtons
              && openingLoadState.openingRemovedPressureForceOnSheetNewtons
                  == openingAdjustedSurfaceLoads
                         .openingRemovedTotalPressureForceOnSheetNewtons
              && openingLoadState.solidPressureForceOnSheetNewtons
                  == openingAdjustedSurfaceLoads
                         .solidTotalPressureForceOnSheetNewtons
              && openingLoadState.solidPressureImpulseOnSheetNewtonSeconds
                  == openingAdjustedSurfaceLoads
                         .solidTotalPressureImpulseOnSheetNewtonSeconds
              && openingLoadState.solidPressureMomentOnSheetNewtonMeters
                  == openingAdjustedSurfaceLoads
                         .solidTotalPressureMomentOnSheetNewtonMeters
              && openingLoadState.solidPressureWorkToSheetJoules
                  == openingAdjustedSurfaceLoads
                         .solidTotalPressureWorkToSheetJoules,
          "opening load state exposes exact full, removed, and retained ledgers");
    validatePlanarPressureRegionFragmentOpeningLoadStateIntegrity(
        openingLoadState);
    validatePlanarPressureRegionFragmentOpeningLoadState(
        openingLoadState, pressureOperator, base, geometry, sweep, fragments,
        topology, volumeRates, definitions, openings, resistance);
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningLoadState(
                openingLoadState, pressureOperator, base, geometry, sweep,
                fragments, topology, volumeRates, definitions, openings,
                foreignResistance);
        },
        "opening load state rejects foreign resistance provenance");
    auto corruptOpeningLoadState = openingLoadState;
    corruptOpeningLoadState.openingSurfaceLoads.tiles[0]
        .solidAreaSquareMeters += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningLoadStateIntegrity(
                corruptOpeningLoadState);
        },
        "opening load-state integrity rejects nested load corruption");
    corruptOpeningLoadState = openingLoadState;
    corruptOpeningLoadState.solidPressureForceOnSheetNewtons.x += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningLoadStateIntegrity(
                corruptOpeningLoadState);
        },
        "opening load-state integrity rejects aggregate corruption");
    auto openingLoadStateLimits =
        PlanarPressureRegionFragmentOpeningLoadStateLimits{};
    openingLoadStateLimits.maximumOwnedBytes =
        openingLoadState.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningLoadState(
                    acceptedState, openingPressureState,
                    openingSurfaceLoads, openingAdjustedSurfaceLoads,
                    pressureOperator, base, geometry, sweep, fragments,
                    topology, volumeRates, definitions, openings,
                    resistance, openingLoadStateLimits));
        },
        "opening load state enforces its aggregate owned-byte limit");
    openingLoadStateLimits = {};
    openingLoadStateLimits.surfaceLoadLimits.maximumTiles =
        openingAdjustedSurfaceLoads.tiles.size() - 1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningLoadState(
                openingLoadState, pressureOperator, base, geometry, sweep,
                fragments, topology, volumeRates, definitions, openings,
                resistance, openingLoadStateLimits);
        },
        "opening load state enforces nested load limits");
    auto corruptOpeningSurfaceLoads = openingSurfaceLoads;
    corruptOpeningSurfaceLoads.tiles[0]
        .totalPressureForceOnSheetNewtons.x += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentSurfaceLoads(
                corruptOpeningSurfaceLoads, openingPressureState);
        },
        "opening pressure surface-load validation rejects corruption");
    auto openingSurfaceLoadLimits =
        PlanarPressureRegionFragmentSurfaceLoadLimits{};
    openingSurfaceLoadLimits.maximumTiles =
        openingSurfaceLoads.tiles.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentSurfaceLoads(
                    openingPressureState, openingSurfaceLoadLimits));
        },
        "opening pressure surface loads enforce the tile limit");
    auto corruptOpeningPressureState = openingPressureState;
    corruptOpeningPressureState.controls[0].totalPressurePascals += 0.1;
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningPressureStateIntegrity(
                corruptOpeningPressureState);
        },
        "opening pressure-state integrity rejects control corruption");
    auto openingPressureLimits =
        PlanarPressureRegionFragmentOpeningPressureStateLimits{};
    openingPressureLimits.maximumControls =
        openingPressureState.controls.size() - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composePlanarPressureRegionFragmentOpeningPressureState(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, {}, openingPressureLimits));
        },
        "opening pressure state enforces its control limit");
    openingPressureLimits = {};
    openingPressureLimits.maximumOwnedBytes =
        openingPressureState.ownedStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composePlanarPressureRegionFragmentOpeningPressureState(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, {}, openingPressureLimits));
        },
        "opening pressure state enforces its owned-byte limit");
    openingPressureLimits = {};
    openingPressureLimits.maximumWorkingBytes =
        openingPressureState.workingStorageBytes - 1;
    expectRejected(
        [&] {
            static_cast<void>(
                composePlanarPressureRegionFragmentOpeningPressureState(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, {}, openingPressureLimits));
        },
        "opening pressure state enforces its working-byte limit");
    auto invalidOpeningPressureSettings =
        PlanarPressureRegionFragmentOpeningPressureStateSettings{};
    invalidOpeningPressureSettings.absolutePressureResidualTolerancePascals =
        0.0;
    invalidOpeningPressureSettings.relativePressureResidualTolerance = 0.0;
    expectRejected(
        [&] {
            static_cast<void>(
                composePlanarPressureRegionFragmentOpeningPressureState(
                    acceptedState, pressureOperator, base, geometry, sweep,
                    fragments, topology, volumeRates, definitions, openings,
                    resistance, invalidOpeningPressureSettings));
        },
        "opening pressure state rejects disabled pressure tolerances");

    const std::vector<
        PlanarPressureRegionFragmentOpeningPatchDefinition>
        splitDefinitions{
            {102, 1000, wall->surfaceStableId, wall->axis,
             wall->i, wall->j, wall->k, wall->minusRegionStableId,
             wall->plusRegionStableId, 0.4 * patchArea},
            {101, 1000, wall->surfaceStableId, wall->axis,
             wall->i, wall->j, wall->k, wall->minusRegionStableId,
             wall->plusRegionStableId, 0.6 * patchArea},
        };
    const auto splitOpenings = buildPlanarPressureRegionFragmentOpenings(
        geometry, sweep, fragments, topology, splitDefinitions);
    const auto splitPressureOperator =
        buildPlanarPressureRegionFragmentOpeningPressureOperator(
            base, geometry, sweep, fragments, topology,
            splitDefinitions, splitOpenings);
    std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
        splitSamples{
            {102, requiredFlow / patchArea},
            {101, requiredFlow / patchArea},
        };
    auto splitFlux = buildPlanarPressureRegionFragmentOpeningFluxState(
        geometry, sweep, fragments, topology, splitDefinitions,
        splitOpenings, splitSamples);
    std::vector<double> splitVelocity(topology.links.size(), 0.0);
    std::vector<double> splitPressure(base.rows.size(), 0.0);
    const auto splitSeeded =
        projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
            base, geometry, sweep, fragments, topology, volumeRates,
            splitDefinitions, splitOpenings, splitFlux, splitSamples,
            splitVelocity, splitPressure, seedSettings);
    const std::vector<
        PlanarPressureRegionFragmentOpeningResistanceDefinition>
        splitResistance{{102, {8.0, 3.0}}, {101, {8.0, 3.0}}};
    const auto splitDiagnostics =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            splitPressureOperator, base, geometry, sweep, fragments,
            topology, volumeRates, splitDefinitions, splitOpenings,
            splitResistance, splitVelocity, splitSamples, splitFlux,
            splitPressure, settings);
    check(splitSeeded.accepted && splitDiagnostics.accepted,
          "opening accepted-state canonical-order fixture advances");
    if (splitSeeded.accepted && splitDiagnostics.accepted) {
        auto reversedSplitSamples = splitSamples;
        std::ranges::reverse(reversedSplitSamples);
        const auto splitState =
            capturePlanarPressureRegionFragmentOpeningAcceptedState(
                splitPressureOperator, base, geometry, sweep, fragments,
                topology, volumeRates, splitDefinitions, splitOpenings,
                splitResistance, splitDiagnostics, splitVelocity,
                splitSamples, splitFlux, splitPressure, settings);
        const auto reorderedSplitState =
            capturePlanarPressureRegionFragmentOpeningAcceptedState(
                splitPressureOperator, base, geometry, sweep, fragments,
                topology, volumeRates, splitDefinitions, splitOpenings,
                splitResistance, splitDiagnostics, splitVelocity,
                reversedSplitSamples, splitFlux, splitPressure, settings);
        check(splitState == reorderedSplitState
                  && splitState.openingVelocitySamples.size() == 2
                  && splitState.openingVelocitySamples[0].patchStableId
                      == 101
                  && splitState.openingVelocitySamples[1].patchStableId
                      == 102,
              "opening accepted state canonicalizes caller sample order");
    }

    auto authoredDriveSettings = settings;
    authoredDriveSettings.useAuthoredPressureDrive = true;
    auto authoredDriveVelocity = initialVelocity;
    auto authoredDriveSamples = initialSamples;
    auto authoredDriveFlux = initialFlux;
    auto authoredDrivePressure = initialPressure;
    const auto authoredDrive =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance,
            authoredDriveVelocity, authoredDriveSamples,
            authoredDriveFlux, authoredDrivePressure,
            authoredDriveSettings);
    check(authoredDrive.accepted && authoredDrive.finite
              && authoredDrive.energyAccepted
              && authoredDrive.resistance.accepted
              && authoredDrive.resistance.usesAuthoredPressureDrive
              && authoredDrive.resistance.patches[0]
                     .drivingPressureRisePascals
                  == -wall->pressureJumpPascals
              && authoredDrive.authoredPressureWorkJoules
                  == authoredDrive.resistance.authoredPressureWorkJoules
              && authoredDrive.dissipatedEnergyJoules > 0.0
              && authoredDrive.correctionKineticEnergyJoules > 0.0
              && std::abs(authoredDrive.energyResidualJoules)
                  <= authoredDrive.energyToleranceJoules,
          "composed opening step closes authored drive, resistance, and projection energy");
    checkNear(
        authoredDriveSamples[0].relativeNormalVelocityMetersPerSecond,
        initialSamples[0].relativeNormalVelocityMetersPerSecond,
        4.0e-11,
        "projection restores required intake after authored pressure drive");
    check(authoredDrivePressure != pressure,
          "authored aperture drive changes the compensating correction pressure");
    for (const auto& link : topology.links) {
        if (link.kind
            == PlanarPressureRegionFragmentFaceKind::PressureLayerWall) {
            check(authoredDriveVelocity[link.linkIndex] == 0.0,
                  "authored-pressure composed step retains solid-wall zero flow");
        }
    }
    const auto authoredDriveState =
        capturePlanarPressureRegionFragmentOpeningAcceptedState(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance, authoredDrive,
            authoredDriveVelocity, authoredDriveSamples, authoredDriveFlux,
            authoredDrivePressure, authoredDriveSettings);
    check(authoredDriveState.accepted
              && authoredDriveState.settings.useAuthoredPressureDrive
              && authoredDriveState.authoredPressureWorkJoules
                  == authoredDrive.authoredPressureWorkJoules
              && authoredDriveState.pressureCorrectionPascals
                  == authoredDrivePressure
              && authoredDriveState.pressureCorrectionPascals
                  != acceptedState.pressureCorrectionPascals
              && authoredDriveState.resultOpeningFluxFingerprint
                  == authoredDriveFlux.fingerprint,
          "opening accepted state distinguishes authored-pressure continuation");
    validatePlanarPressureRegionFragmentOpeningAcceptedState(
        authoredDriveState, pressureOperator, base, geometry, sweep,
        fragments, topology, volumeRates, definitions, openings,
        resistance);
    const auto authoredDrivePressureState =
        composePlanarPressureRegionFragmentOpeningPressureState(
            authoredDriveState, pressureOperator, base, geometry, sweep,
            fragments, topology, volumeRates, definitions, openings,
            resistance);
    check(authoredDrivePressureState.accepted
              && authoredDrivePressureState.sourceAcceptedStateFingerprint
                  == authoredDriveState.fingerprint
              && authoredDrivePressureState.correctionGeometryPressureWorkJoules
                  == authoredDriveState.geometryPressureWorkJoules
              && authoredDrivePressureState.controls
                  != openingPressureState.controls
              && authoredDrivePressureState.correctionPressureForceOnSheetNewtons
                  != openingPressureState
                         .correctionPressureForceOnSheetNewtons
              && std::abs(
                     authoredDrivePressureState
                         .wallGeometryWorkResidualJoules)
                  <= authoredDrivePressureState
                         .settings.absoluteWorkResidualToleranceJoules,
          "opening pressure state follows the authored-drive correction without mixing sealed gauges");
    const auto authoredDriveSurfaceLoads =
        capturePlanarPressureRegionFragmentSurfaceLoads(
            authoredDrivePressureState);
    const auto authoredDriveOpeningAdjustedLoads =
        capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
            authoredDriveSurfaceLoads, authoredDrivePressureState, geometry,
            sweep, fragments, topology, definitions, openings);
    check(authoredDriveSurfaceLoads.sourcePressureStateFingerprint
              == authoredDrivePressureState.fingerprint
              && authoredDriveOpeningAdjustedLoads
                         .sourcePressureStateFingerprint
                  == authoredDrivePressureState.fingerprint
              && authoredDriveOpeningAdjustedLoads
                         .solidCorrectionPressureForceOnSheetNewtons
                  != openingAdjustedSurfaceLoads
                         .solidCorrectionPressureForceOnSheetNewtons,
          "authored aperture drive reaches the retained-solid load ledger");
    const auto authoredDriveLoadState =
        capturePlanarPressureRegionFragmentOpeningLoadState(
            authoredDriveState, authoredDrivePressureState,
            authoredDriveSurfaceLoads, authoredDriveOpeningAdjustedLoads,
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance);
    check(authoredDriveLoadState.accepted
              && authoredDriveLoadState.fingerprint
                  != openingLoadState.fingerprint
              && authoredDriveLoadState.solidPressureForceOnSheetNewtons
                  != openingLoadState.solidPressureForceOnSheetNewtons,
          "opening load state distinguishes authored-drive retained loads");
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningLoadState(
                    acceptedState, authoredDrivePressureState,
                    authoredDriveSurfaceLoads,
                    authoredDriveOpeningAdjustedLoads, pressureOperator,
                    base, geometry, sweep, fragments, topology, volumeRates,
                    definitions, openings, resistance));
        },
        "opening load state rejects mixed accepted-flow lineage");
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningSurfaceLoads(
                openingAdjustedSurfaceLoads, openingSurfaceLoads,
                authoredDrivePressureState, geometry, sweep, fragments,
                topology, definitions, openings);
        },
        "opening surface loads reject a foreign pressure endpoint");
    expectRejected(
        [&] {
            validatePlanarPressureRegionFragmentOpeningPressureState(
                openingPressureState, authoredDriveState, pressureOperator,
                base, geometry, sweep, fragments, topology, volumeRates,
                definitions, openings, resistance);
        },
        "opening pressure state rejects a foreign accepted endpoint");

    const std::vector<
        PlanarPressureRegionFragmentOpeningResistanceDefinition>
        zeroResistance{{100, {0.0, 0.0}}};
    auto zeroVelocity = initialVelocity;
    auto zeroSamples = initialSamples;
    auto zeroFlux = initialFlux;
    auto zeroPressure = initialPressure;
    const auto zero =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, zeroResistance, zeroVelocity,
            zeroSamples, zeroFlux, zeroPressure, settings);
    auto directVelocity = initialVelocity;
    auto directSamples = initialSamples;
    auto directFlux = initialFlux;
    auto directPressure = initialPressure;
    const auto direct =
        projectMovingPlanarPressureRegionFragmentVelocitiesWithPressureOpenings(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, directVelocity, directSamples,
            directFlux, directPressure, settings.projection);
    check(zero.accepted && direct.accepted
              && zero.resistance.zeroResistancePatchCount == 1
              && zero.resistance.dissipatedEnergyJoules == 0.0
              && zero.projection == direct
              && zeroVelocity == directVelocity
              && zeroSamples == directSamples
              && zeroFlux == directFlux
              && zeroPressure == directPressure,
          "zero-loss composed step is bit-exact with direct projection");

    auto truncatedSettings = settings;
    truncatedSettings.projection.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-16;
    truncatedSettings.projection.pressureSolve.relativeResidualTolerance = 0.0;
    truncatedSettings.projection.pressureSolve.maximumIterations = 1;
    auto rejectedVelocity = initialVelocity;
    auto rejectedSamples = initialSamples;
    auto rejectedFlux = initialFlux;
    auto rejectedPressure = initialPressure;
    const auto rejected =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance, rejectedVelocity,
            rejectedSamples, rejectedFlux, rejectedPressure,
            truncatedSettings);
    check(!rejected.accepted && rejected.finite
              && rejected.resistance.accepted
              && !rejected.projection.accepted
              && rejected.projection.pressureSolve.compatible
              && !rejected.projection.pressureSolve.converged
              && rejectedVelocity == initialVelocity
              && rejectedSamples == initialSamples
              && rejectedFlux == initialFlux
              && rejectedPressure == initialPressure,
          "failed resisted projection rolls back all composed fields");
    expectRejected(
        [&] {
            static_cast<void>(
                capturePlanarPressureRegionFragmentOpeningAcceptedState(
                    pressureOperator, base, geometry, sweep, fragments,
                    topology, volumeRates, definitions, openings,
                    resistance, rejected, rejectedVelocity,
                    rejectedSamples, rejectedFlux, rejectedPressure,
                    truncatedSettings));
        },
        "opening accepted state rejects an unaccepted pressure step");

    auto rejectedAuthoredVelocity = initialVelocity;
    auto rejectedAuthoredSamples = initialSamples;
    auto rejectedAuthoredFlux = initialFlux;
    auto rejectedAuthoredPressure = initialPressure;
    auto truncatedAuthoredSettings = truncatedSettings;
    truncatedAuthoredSettings.useAuthoredPressureDrive = true;
    const auto rejectedAuthored =
        advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
            pressureOperator, base, geometry, sweep, fragments, topology,
            volumeRates, definitions, openings, resistance,
            rejectedAuthoredVelocity, rejectedAuthoredSamples,
            rejectedAuthoredFlux, rejectedAuthoredPressure,
            truncatedAuthoredSettings);
    check(!rejectedAuthored.accepted
              && rejectedAuthored.resistance.accepted
              && rejectedAuthored.resistance.usesAuthoredPressureDrive
              && !rejectedAuthored.projection.accepted
              && rejectedAuthoredVelocity == initialVelocity
              && rejectedAuthoredSamples == initialSamples
              && rejectedAuthoredFlux == initialFlux
              && rejectedAuthoredPressure == initialPressure,
          "failed projection rolls back the completed authored-pressure stage");

    auto limits = PlanarPressureRegionFragmentOpeningPressureStepLimits{};
    limits.maximumWorkingBytes = 1;
    expectRejected([&] {
        auto candidateVelocity = initialVelocity;
        auto candidateSamples = initialSamples;
        auto candidateFlux = initialFlux;
        auto candidatePressure = initialPressure;
        static_cast<void>(
            advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
                pressureOperator, base, geometry, sweep, fragments, topology,
                volumeRates, definitions, openings, resistance,
                candidateVelocity, candidateSamples, candidateFlux,
                candidatePressure, settings, limits));
    }, "resisted opening pressure step enforces its working limit");
}

void testPlanarRegionalOpeningPressureStateAxes() {
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
        const auto base = buildPlanarPressureRegionFragmentPressureOperator(
            geometry, sweep, fragments, topology);
        const auto volumeRates =
            buildPlanarPressureRegionFragmentVolumeRates(
                geometry, sweep, fragments, topology);
        const auto wall = std::ranges::find_if(
            topology.links,
            [](const auto& link) {
                return link.kind
                        == PlanarPressureRegionFragmentFaceKind::
                            PressureLayerWall
                    && link.surfaceStableId == 10;
            });
        check(wall != topology.links.end(),
              "opening pressure state finds a wall on every axis");
        if (wall == topology.links.end()) continue;
        const double patchArea = 0.5 * wall->areaSquareMeters;
        const std::vector<
            PlanarPressureRegionFragmentOpeningPatchDefinition>
            definitions{{
                100, 1000, wall->surfaceStableId, wall->axis,
                wall->i, wall->j, wall->k, wall->minusRegionStableId,
                wall->plusRegionStableId, patchArea,
            }};
        const auto openings = buildPlanarPressureRegionFragmentOpenings(
            geometry, sweep, fragments, topology, definitions);
        const auto baseVelocityMetric =
            buildPlanarPressureRegionFragmentVelocityMetric(
                geometry, sweep, fragments, topology);
        const auto openingVelocityMetric =
            buildPlanarPressureRegionFragmentOpeningVelocityMetric(
                geometry, sweep, fragments, topology,
                baseVelocityMetric, definitions, openings);
        const auto pressureOperator =
            buildPlanarPressureRegionFragmentOpeningPressureOperator(
                base, geometry, sweep, fragments, topology, definitions,
                openings);
        const double requiredFlow =
            -volumeRates.components[wall->minusComponentIndex]
                 .geometryVolumeChangeRateCubicMetersPerSecond;
        std::vector<PlanarPressureRegionFragmentOpeningVelocitySample>
            samples{{100, requiredFlow / patchArea}};
        auto flux = buildPlanarPressureRegionFragmentOpeningFluxState(
            geometry, sweep, fragments, topology, definitions, openings,
            samples);
        std::vector<double> velocity(topology.links.size(), 0.0);
        std::vector<double> pressure(base.rows.size(), 0.0);
        PlanarPressureRegionFragmentPressureProjectionSettings seedSettings;
        seedSettings.densityKgPerCubicMeter = 1.2;
        seedSettings.timeStepSeconds = 0.5;
        seedSettings.absoluteContinuityToleranceCubicMetersPerSecond =
            2.0e-11;
        seedSettings.relativeContinuityTolerance = 1.0e-10;
        seedSettings.pressureSolve
            .absoluteResidualTolerancePascalsMeters = 1.0e-13;
        seedSettings.pressureSolve.relativeResidualTolerance = 0.0;
        seedSettings.pressureSolve.maximumIterations = 300;
        const auto seeded =
            projectMovingPlanarPressureRegionFragmentFaceVelocitiesWithOpenings(
                base, geometry, sweep, fragments, topology, volumeRates,
                definitions, openings, flux, samples, velocity, pressure,
                seedSettings);
        PlanarPressureRegionFragmentOpeningPressureStepSettings settings;
        settings.projection.densityKgPerCubicMeter = 1.2;
        settings.projection.timeStepSeconds = 0.5;
        settings.projection
            .absoluteContinuityToleranceCubicMetersPerSecond = 2.0e-11;
        settings.projection.relativeContinuityTolerance = 1.0e-10;
        settings.projection
            .absoluteMomentumResidualToleranceKilogramMetersPerSecond =
            2.0e-12;
        settings.projection.relativeMomentumResidualTolerance = 1.0e-10;
        settings.projection.absoluteEnergyResidualToleranceJoules = 2.0e-11;
        settings.projection.relativeEnergyResidualTolerance = 1.0e-10;
        settings.projection.pressureSolve
            .absoluteResidualTolerancePascalsMeters = 1.0e-13;
        settings.projection.pressureSolve.relativeResidualTolerance = 0.0;
        settings.projection.pressureSolve.maximumIterations = 300;
        const std::vector<
            PlanarPressureRegionFragmentOpeningResistanceDefinition>
            resistance{{100, {0.0, 0.0}}};
        const auto diagnostics =
            advanceMovingPlanarPressureRegionFragmentOpeningPressureStep(
                pressureOperator, base, geometry, sweep, fragments,
                topology, volumeRates, definitions, openings, resistance,
                velocity, samples, flux, pressure, settings);
        check(seeded.accepted && diagnostics.accepted,
              "opening pressure-state all-axis fixture advances");
        if (!seeded.accepted || !diagnostics.accepted) continue;
        const auto accepted =
            capturePlanarPressureRegionFragmentOpeningAcceptedState(
                pressureOperator, base, geometry, sweep, fragments,
                topology, volumeRates, definitions, openings, resistance,
                diagnostics, velocity, samples, flux, pressure, settings);
        const auto acceptedVelocityState =
            capturePlanarPressureRegionFragmentOpeningVelocityState(
                accepted, pressureOperator, base, geometry, sweep,
                fragments, topology, volumeRates, definitions, openings,
                resistance, baseVelocityMetric, openingVelocityMetric);
        const auto acceptedAperture = std::ranges::find_if(
            acceptedVelocityState.samples,
            [](const auto& sample) {
                return sample.kind
                    == PlanarPressureRegionFragmentOpeningVelocityDofKind::
                        OpeningPatch;
            });
        check(acceptedVelocityState.mappedFromAcceptedEndpoint
                  && acceptedVelocityState.sourceAcceptedStateFingerprint
                      == accepted.fingerprint
                  && acceptedAperture
                      != acceptedVelocityState.samples.end()
                  && acceptedAperture->axis == axis
                  && acceptedVelocityState
                         .maximumAbsoluteVelocityCompositionResidualMetersPerSecond
                      == 0.0,
              "accepted opening velocity state maps material and relative flow on every axis");
        validatePlanarPressureRegionFragmentOpeningAcceptedVelocityState(
            acceptedVelocityState, accepted, pressureOperator, base,
            geometry, sweep, fragments, topology, volumeRates,
            definitions, openings, resistance, baseVelocityMetric,
            openingVelocityMetric);
        const auto state =
            composePlanarPressureRegionFragmentOpeningPressureState(
                accepted, pressureOperator, base, geometry, sweep,
                fragments, topology, volumeRates, definitions, openings,
                resistance);
        const auto surfaceLoads =
            capturePlanarPressureRegionFragmentSurfaceLoads(state);
        const auto solidLoads =
            capturePlanarPressureRegionFragmentOpeningSurfaceLoads(
                surfaceLoads, state, geometry, sweep, fragments, topology,
                definitions, openings);
        const auto pressureWall = std::ranges::find_if(
            state.walls,
            [&](const auto& candidate) {
                return candidate.sourceFaceLinkStableId == wall->stableId;
            });
        check(pressureWall != state.walls.end(),
              "opening pressure state retains its source wall on every axis");
        if (pressureWall == state.walls.end()) continue;
        const auto solidTile = std::ranges::find_if(
            solidLoads.tiles,
            [&](const auto& candidate) {
                return candidate.sourceFaceLinkStableId == wall->stableId;
            });
        check(solidTile != solidLoads.tiles.end(),
              "opening pressure loads retain the aperture wall on every axis");
        if (solidTile == solidLoads.tiles.end()) continue;
        const Vector3 expectedAuthoredForce = axis == GridFaceAxis::X
            ? Vector3{-70.0 * wall->areaSquareMeters, 0.0, 0.0}
            : axis == GridFaceAxis::Y
            ? Vector3{0.0, -70.0 * wall->areaSquareMeters, 0.0}
            : Vector3{0.0, 0.0, -70.0 * wall->areaSquareMeters};
        check(state.accepted && state.components.size() == 1
                  && pressureWall->axis == axis
                  && pressureWall->authoredPressureForceOnSheetNewtons
                      == expectedAuthoredForce
                  && pressureWall->totalPressureImpulseOnSheetNewtonSeconds
                      == Vector3{
                          pressureWall
                                  ->totalPressureForceOnSheetNewtons.x
                              * state.timeStepSeconds,
                          pressureWall
                                  ->totalPressureForceOnSheetNewtons.y
                              * state.timeStepSeconds,
                          pressureWall
                                  ->totalPressureForceOnSheetNewtons.z
                              * state.timeStepSeconds}
                  && std::abs(state.wallGeometryWorkResidualJoules)
                      <= state.settings
                             .absoluteWorkResidualToleranceJoules,
              "opening pressure state closes force, impulse, and work on every axis");
        check(surfaceLoads.accepted && solidLoads.accepted
                  && !surfaceLoads.staticGeometry
                  && surfaceLoads.usesMovingVolumeRates
                  && solidTile->axis == axis
                  && solidTile->openingAreaSquareMeters == patchArea
                  && solidTile->solidAreaSquareMeters == patchArea
                  && solidLoads.sourcePressureStateFingerprint
                      == state.fingerprint,
              "opening-connected pressure partitions full-wall load on every axis");
        checkNear(
            solidTile->solidTotalPressureForceOnSheetNewtons.x,
            0.5 * pressureWall->totalPressureForceOnSheetNewtons.x,
            1.0e-12,
            "opening-connected solid X force uses retained area");
        checkNear(
            solidTile->solidTotalPressureForceOnSheetNewtons.y,
            0.5 * pressureWall->totalPressureForceOnSheetNewtons.y,
            1.0e-12,
            "opening-connected solid Y force uses retained area");
        checkNear(
            solidTile->solidTotalPressureForceOnSheetNewtons.z,
            0.5 * pressureWall->totalPressureForceOnSheetNewtons.z,
            1.0e-12,
            "opening-connected solid Z force uses retained area");
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
    testPlanarRegionalFragmentOpeningTopology();
    testPlanarRegionalFragmentOpeningFluxState();
    testPlanarRegionalFragmentVolumeRates();
    testPlanarRegionalFragmentVolumeRatesAxesAndRejection();
    testPlanarRegionalFragmentVelocityMetric();
    testPlanarRegionalFragmentVelocityMetricAxesAndRejection();
    testPlanarRegionalFragmentVelocityState();
    testPlanarRegionalFragmentVelocityStateMotionAndRejection();
    testPlanarRegionalFragmentPressureJumpEnergyAudit();
    testPlanarRegionalMovingPressureJumpEnergyAudit();
    testPlanarRegionalFragmentPressureOperator();
    testPlanarRegionalFragmentPressureOperatorAxesAndRejection();
    testPlanarRegionalFragmentOpeningPressureOperator();
    testPlanarRegionalFragmentPressureCorrectionSolve();
    testPlanarRegionalFragmentOpeningPressureCorrectionSolve();
    testPlanarRegionalFragmentPressureCorrectionRollback();
    testPlanarRegionalFragmentPressureProjection();
    testPlanarRegionalFragmentProjectionEnergyAudit();
    testPlanarRegionalFragmentProjectionEnergyAxesAndRejection();
    testPlanarRegionalFragmentPressureProjectionAxesAndRollback();
    testPlanarRegionalMovingFragmentPressureProjection();
    testPlanarRegionalMovingFragmentPressureProjectionAllAxes();
    testPlanarRegionalOpeningFluxPressureProjection();
    testPlanarRegionalPressureDrivenOpeningProjection();
    testPlanarRegionalOpeningResistance();
    testPlanarRegionalOpeningMomentumTransport();
    testPlanarRegionalResistedOpeningPressureStep();
    testPlanarRegionalOpeningPressureStateAxes();
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
