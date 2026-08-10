#include "fluid/planar_region_fragment_surface_load.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
            value >>= 8U;
        }
    }

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1U : 0U));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

void validateLimits(
    const PlanarPressureRegionFragmentSurfaceLoadLimits& limits) {
    if (limits.maximumTiles == 0 || limits.maximumSurfaces == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "planar regional surface-load limits are invalid");
    }
}

std::size_t checkedProduct(const std::size_t first,
                           const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "planar regional surface-load storage overflows");
    }
    return first * second;
}

std::size_t checkedSum(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "planar regional surface-load storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(const std::size_t tileCount,
                              const std::size_t surfaceCount) {
    return checkedSum(
        checkedProduct(
            tileCount,
            sizeof(PlanarPressureRegionFragmentSurfaceLoadTile)),
        checkedProduct(
            surfaceCount,
            sizeof(PlanarPressureRegionFragmentSurfaceLoadSummary)));
}

Vector3 scaledVector(const Vector3& value, const double scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Vector3 vectorSum(const Vector3& first, const Vector3& second) {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z,
    };
}

Vector3 vectorDifference(const Vector3& first, const Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

void addVector(Vector3& target, const Vector3& value) {
    target.x += value.x;
    target.y += value.y;
    target.z += value.z;
}

double maximumAbsoluteComponent(const Vector3& value) {
    return std::max({
        std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double closureTolerance(const double scale) {
    return std::max(
        1.0e-12,
        64.0 * std::numeric_limits<double>::epsilon()
            * std::max(1.0, std::abs(scale)));
}

void fingerprintVector(Fingerprint& fingerprint, const Vector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t ledgerFingerprint(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& ledger) {
    Fingerprint fingerprint;
    fingerprint.integer(ledger.version);
    fingerprint.integer(ledger.sourcePressureStateFingerprint);
    fingerprint.integer(ledger.sourceTopologyFingerprint);
    fingerprint.boolean(ledger.staticGeometry);
    fingerprint.boolean(ledger.usesMovingVolumeRates);
    fingerprint.real(ledger.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(ledger.tiles.size()));
    for (const auto& tile : ledger.tiles) {
        for (const std::size_t value : {
                 tile.tileIndex,
                 tile.sourcePressureWallIndex,
                 tile.sourceFaceLinkIndex}) {
            fingerprint.integer(static_cast<std::uint64_t>(value));
        }
        fingerprint.integer(tile.sourceFaceLinkStableId);
        fingerprint.integer(tile.surfaceStableId);
        fingerprint.enumeration(tile.axis);
        fingerprint.integer(tile.minusRegionStableId);
        fingerprint.integer(tile.plusRegionStableId);
        fingerprint.real(tile.areaSquareMeters);
        fingerprintVector(fingerprint, tile.wrappedCentroidMeters);
        fingerprintVector(fingerprint, tile.unitNormalMinusToPlus);
        fingerprintVector(
            fingerprint, tile.authoredPressureTractionOnSheetPascals);
        fingerprintVector(
            fingerprint, tile.correctionPressureTractionOnSheetPascals);
        fingerprintVector(
            fingerprint, tile.totalPressureTractionOnSheetPascals);
        fingerprintVector(
            fingerprint, tile.authoredPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, tile.correctionPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, tile.totalPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, tile.totalPressureImpulseOnSheetNewtonSeconds);
        fingerprint.real(tile.totalPressureWorkToSheetJoules);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        ledger.surfaces.size()));
    for (const auto& surface : ledger.surfaces) {
        fingerprint.integer(static_cast<std::uint64_t>(
            surface.surfaceIndex));
        fingerprint.integer(surface.surfaceStableId);
        fingerprint.enumeration(surface.axis);
        fingerprint.integer(surface.minusRegionStableId);
        fingerprint.integer(surface.plusRegionStableId);
        fingerprint.integer(static_cast<std::uint64_t>(surface.tileCount));
        fingerprint.real(surface.areaSquareMeters);
        fingerprintVector(
            fingerprint, surface.areaWeightedCentroidMeters);
        fingerprintVector(
            fingerprint, surface.authoredPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, surface.correctionPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint, surface.totalPressureForceOnSheetNewtons);
        fingerprintVector(
            fingerprint,
            surface.totalPressureImpulseOnSheetNewtonSeconds);
        fingerprint.real(surface.totalPressureWorkToSheetJoules);
    }
    fingerprint.real(ledger.totalAreaSquareMeters);
    fingerprintVector(
        fingerprint, ledger.authoredPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.correctionPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.totalPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, ledger.pressureForceSplitResidualNewtons);
    fingerprintVector(
        fingerprint, ledger.totalPressureImpulseOnSheetNewtonSeconds);
    fingerprintVector(
        fingerprint, ledger.impulseTimeIntegrationResidualNewtonSeconds);
    fingerprint.real(ledger.totalPressureWorkToSheetJoules);
    fingerprint.real(ledger.sourceWorkResidualJoules);
    fingerprint.real(
        ledger.maximumAbsoluteForceReconstructionResidualNewtons);
    fingerprint.real(
        ledger.maximumAbsoluteSurfaceAggregationResidualNewtons);
    fingerprint.boolean(ledger.accepted);
    fingerprint.integer(static_cast<std::uint64_t>(
        ledger.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        ledger.workingStorageBytes));
    return fingerprint.value();
}

PlanarPressureRegionFragmentSurfaceLoadLedger buildLedger(
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLimits& limits) {
    validateLimits(limits);
    validatePlanarPressureRegionFragmentPressureStateIntegrity(pressureState);
    if (pressureState.walls.size() > limits.maximumTiles) {
        throw std::length_error(
            "planar regional surface-load tile limit exceeded");
    }

    std::vector<std::uint64_t> surfaceIds;
    surfaceIds.reserve(pressureState.walls.size());
    for (const auto& wall : pressureState.walls) {
        surfaceIds.push_back(wall.surfaceStableId);
    }
    std::ranges::sort(surfaceIds);
    const auto uniqueEnd = std::ranges::unique(surfaceIds).begin();
    surfaceIds.erase(uniqueEnd, surfaceIds.end());
    if (surfaceIds.size() > limits.maximumSurfaces) {
        throw std::length_error(
            "planar regional surface-load surface limit exceeded");
    }

    PlanarPressureRegionFragmentSurfaceLoadLedger result;
    result.sourcePressureStateFingerprint = pressureState.fingerprint;
    result.sourceTopologyFingerprint =
        pressureState.sourceTopologyFingerprint;
    result.staticGeometry = pressureState.staticGeometry;
    result.usesMovingVolumeRates = pressureState.usesMovingVolumeRates;
    result.timeStepSeconds = pressureState.timeStepSeconds;
    result.ownedStorageBytes = ownedStorageBytes(
        pressureState.walls.size(), surfaceIds.size());
    result.workingStorageBytes = checkedProduct(
        surfaceIds.size(), sizeof(std::uint64_t));
    if (result.ownedStorageBytes > limits.maximumOwnedBytes
        || result.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "planar regional surface-load storage limit exceeded");
    }
    result.tiles.reserve(pressureState.walls.size());
    result.surfaces.resize(surfaceIds.size());
    for (std::size_t index = 0; index < surfaceIds.size(); ++index) {
        result.surfaces[index].surfaceIndex = index;
        result.surfaces[index].surfaceStableId = surfaceIds[index];
    }

    double maximumForceNewtons = 0.0;
    for (const auto& wall : pressureState.walls) {
        if (wall.surfaceStableId == 0
            || !std::isfinite(wall.areaSquareMeters)
            || !(wall.areaSquareMeters > 0.0)
            || !finiteVector(wall.wrappedCentroidMeters)
            || !finiteVector(wall.unitNormalMinusToPlus)
            || !finiteVector(wall.authoredPressureForceOnSheetNewtons)
            || !finiteVector(wall.correctionPressureForceOnSheetNewtons)
            || !finiteVector(wall.totalPressureForceOnSheetNewtons)
            || !std::isfinite(wall.totalPressureWorkToSheetJoules)) {
            throw std::invalid_argument(
                "planar regional surface-load source wall is invalid");
        }
        const Vector3 authoredTraction = scaledVector(
            wall.authoredPressureForceOnSheetNewtons,
            1.0 / wall.areaSquareMeters);
        const Vector3 correctionTraction = scaledVector(
            wall.correctionPressureForceOnSheetNewtons,
            1.0 / wall.areaSquareMeters);
        const Vector3 totalTraction = scaledVector(
            wall.totalPressureForceOnSheetNewtons,
            1.0 / wall.areaSquareMeters);
        const Vector3 reconstructedForce = scaledVector(
            totalTraction, wall.areaSquareMeters);
        const Vector3 forceReconstructionResidual = vectorDifference(
            reconstructedForce, wall.totalPressureForceOnSheetNewtons);
        const Vector3 impulse = scaledVector(
            wall.totalPressureForceOnSheetNewtons,
            pressureState.timeStepSeconds);
        maximumForceNewtons = std::max({
            maximumForceNewtons,
            maximumAbsoluteComponent(
                wall.authoredPressureForceOnSheetNewtons),
            maximumAbsoluteComponent(
                wall.correctionPressureForceOnSheetNewtons),
            maximumAbsoluteComponent(wall.totalPressureForceOnSheetNewtons)});
        result.maximumAbsoluteForceReconstructionResidualNewtons = std::max(
            result.maximumAbsoluteForceReconstructionResidualNewtons,
            maximumAbsoluteComponent(forceReconstructionResidual));
        if (maximumAbsoluteComponent(forceReconstructionResidual)
            > closureTolerance(maximumForceNewtons)) {
            throw std::invalid_argument(
                "planar regional surface-load traction reconstruction "
                "failed");
        }
        result.tiles.push_back({
            result.tiles.size(),
            wall.wallIndex,
            wall.sourceFaceLinkIndex,
            wall.sourceFaceLinkStableId,
            wall.surfaceStableId,
            wall.axis,
            wall.minusRegionStableId,
            wall.plusRegionStableId,
            wall.areaSquareMeters,
            wall.wrappedCentroidMeters,
            wall.unitNormalMinusToPlus,
            authoredTraction,
            correctionTraction,
            totalTraction,
            wall.authoredPressureForceOnSheetNewtons,
            wall.correctionPressureForceOnSheetNewtons,
            wall.totalPressureForceOnSheetNewtons,
            impulse,
            wall.totalPressureWorkToSheetJoules,
        });
        result.totalAreaSquareMeters += wall.areaSquareMeters;
        addVector(
            result.authoredPressureForceOnSheetNewtons,
            wall.authoredPressureForceOnSheetNewtons);
        addVector(
            result.correctionPressureForceOnSheetNewtons,
            wall.correctionPressureForceOnSheetNewtons);
        addVector(
            result.totalPressureForceOnSheetNewtons,
            wall.totalPressureForceOnSheetNewtons);
        addVector(
            result.totalPressureImpulseOnSheetNewtonSeconds, impulse);
        result.totalPressureWorkToSheetJoules +=
            wall.totalPressureWorkToSheetJoules;

        const auto surfacePosition = std::ranges::lower_bound(
            surfaceIds, wall.surfaceStableId);
        const std::size_t surfaceIndex = static_cast<std::size_t>(
            surfacePosition - surfaceIds.begin());
        auto& surface = result.surfaces[surfaceIndex];
        if (surface.tileCount == 0) {
            surface.axis = wall.axis;
            surface.minusRegionStableId = wall.minusRegionStableId;
            surface.plusRegionStableId = wall.plusRegionStableId;
        } else if (surface.axis != wall.axis
                   || surface.minusRegionStableId
                       != wall.minusRegionStableId
                   || surface.plusRegionStableId
                       != wall.plusRegionStableId) {
            throw std::invalid_argument(
                "planar regional surface-load authored surface is "
                "inconsistent");
        }
        ++surface.tileCount;
        surface.areaSquareMeters += wall.areaSquareMeters;
        addVector(
            surface.areaWeightedCentroidMeters,
            scaledVector(
                wall.wrappedCentroidMeters, wall.areaSquareMeters));
        addVector(
            surface.authoredPressureForceOnSheetNewtons,
            wall.authoredPressureForceOnSheetNewtons);
        addVector(
            surface.correctionPressureForceOnSheetNewtons,
            wall.correctionPressureForceOnSheetNewtons);
        addVector(
            surface.totalPressureForceOnSheetNewtons,
            wall.totalPressureForceOnSheetNewtons);
        addVector(
            surface.totalPressureImpulseOnSheetNewtonSeconds, impulse);
        surface.totalPressureWorkToSheetJoules +=
            wall.totalPressureWorkToSheetJoules;
    }

    Vector3 surfaceAuthoredForce;
    Vector3 surfaceCorrectionForce;
    Vector3 surfaceTotalForce;
    Vector3 surfaceTotalImpulse;
    double surfaceWork = 0.0;
    double surfaceArea = 0.0;
    for (auto& surface : result.surfaces) {
        if (surface.tileCount == 0
            || !(surface.areaSquareMeters > 0.0)) {
            throw std::invalid_argument(
                "planar regional surface-load surface is empty");
        }
        surface.areaWeightedCentroidMeters = scaledVector(
            surface.areaWeightedCentroidMeters,
            1.0 / surface.areaSquareMeters);
        surfaceArea += surface.areaSquareMeters;
        addVector(
            surfaceAuthoredForce,
            surface.authoredPressureForceOnSheetNewtons);
        addVector(
            surfaceCorrectionForce,
            surface.correctionPressureForceOnSheetNewtons);
        addVector(
            surfaceTotalForce, surface.totalPressureForceOnSheetNewtons);
        addVector(
            surfaceTotalImpulse,
            surface.totalPressureImpulseOnSheetNewtonSeconds);
        surfaceWork += surface.totalPressureWorkToSheetJoules;
    }

    result.pressureForceSplitResidualNewtons = vectorDifference(
        result.totalPressureForceOnSheetNewtons,
        vectorSum(
            result.authoredPressureForceOnSheetNewtons,
            result.correctionPressureForceOnSheetNewtons));
    result.impulseTimeIntegrationResidualNewtonSeconds = vectorDifference(
        result.totalPressureImpulseOnSheetNewtonSeconds,
        scaledVector(
            result.totalPressureForceOnSheetNewtons,
            result.timeStepSeconds));
    result.sourceWorkResidualJoules =
        result.totalPressureWorkToSheetJoules
        - pressureState.totalPressureWorkToSheetJoules;
    const Vector3 surfaceAuthoredResidual = vectorDifference(
        surfaceAuthoredForce,
        result.authoredPressureForceOnSheetNewtons);
    const Vector3 surfaceCorrectionResidual = vectorDifference(
        surfaceCorrectionForce,
        result.correctionPressureForceOnSheetNewtons);
    const Vector3 surfaceTotalResidual = vectorDifference(
        surfaceTotalForce, result.totalPressureForceOnSheetNewtons);
    const Vector3 surfaceImpulseResidual = vectorDifference(
        surfaceTotalImpulse,
        result.totalPressureImpulseOnSheetNewtonSeconds);
    result.maximumAbsoluteSurfaceAggregationResidualNewtons = std::max({
        maximumAbsoluteComponent(surfaceAuthoredResidual),
        maximumAbsoluteComponent(surfaceCorrectionResidual),
        maximumAbsoluteComponent(surfaceTotalResidual),
        maximumAbsoluteComponent(surfaceImpulseResidual)
            / result.timeStepSeconds});
    const Vector3 authoredSourceResidual = vectorDifference(
        result.authoredPressureForceOnSheetNewtons,
        pressureState.authoredPressureForceOnSheetNewtons);
    const Vector3 correctionSourceResidual = vectorDifference(
        result.correctionPressureForceOnSheetNewtons,
        pressureState.correctionPressureForceOnSheetNewtons);
    const Vector3 totalSourceResidual = vectorDifference(
        result.totalPressureForceOnSheetNewtons,
        pressureState.totalPressureForceOnSheetNewtons);
    const double forceTolerance = closureTolerance(maximumForceNewtons);
    const double impulseTolerance =
        forceTolerance * result.timeStepSeconds;
    const double workTolerance = closureTolerance(std::max(
        std::abs(result.totalPressureWorkToSheetJoules),
        std::abs(pressureState.totalPressureWorkToSheetJoules)));
    if (!finiteVector(result.pressureForceSplitResidualNewtons)
        || !finiteVector(result.impulseTimeIntegrationResidualNewtonSeconds)
        || maximumAbsoluteComponent(
               result.pressureForceSplitResidualNewtons)
            > forceTolerance
        || maximumAbsoluteComponent(
               result.impulseTimeIntegrationResidualNewtonSeconds)
            > impulseTolerance
        || maximumAbsoluteComponent(authoredSourceResidual) > forceTolerance
        || maximumAbsoluteComponent(correctionSourceResidual) > forceTolerance
        || maximumAbsoluteComponent(totalSourceResidual) > forceTolerance
        || maximumAbsoluteComponent(surfaceAuthoredResidual) > forceTolerance
        || maximumAbsoluteComponent(surfaceCorrectionResidual) > forceTolerance
        || maximumAbsoluteComponent(surfaceTotalResidual) > forceTolerance
        || maximumAbsoluteComponent(surfaceImpulseResidual)
            > impulseTolerance
        || std::abs(surfaceArea - result.totalAreaSquareMeters)
            > closureTolerance(result.totalAreaSquareMeters)
        || !std::isfinite(result.sourceWorkResidualJoules)
        || std::abs(result.sourceWorkResidualJoules) > workTolerance
        || std::abs(surfaceWork - result.totalPressureWorkToSheetJoules)
            > workTolerance) {
        throw std::invalid_argument(
            "planar regional surface-load aggregate closure failed");
    }

    result.accepted = true;
    result.fingerprint = ledgerFingerprint(result);
    return result;
}

} // namespace

PlanarPressureRegionFragmentSurfaceLoadLedger
capturePlanarPressureRegionFragmentSurfaceLoads(
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLimits& limits) {
    return buildLedger(pressureState, limits);
}

void validatePlanarPressureRegionFragmentSurfaceLoadLedgerIntegrity(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& ledger) {
    if (ledger.version != planarPressureRegionFragmentSurfaceLoadVersion
        || ledger.fingerprint == 0 || !ledger.accepted
        || ledger.sourcePressureStateFingerprint == 0
        || ledger.sourceTopologyFingerprint == 0
        || !std::isfinite(ledger.timeStepSeconds)
        || !(ledger.timeStepSeconds > 0.0)
        || ledger.fingerprint != ledgerFingerprint(ledger)) {
        throw std::invalid_argument(
            "planar regional surface-load ledger integrity is invalid");
    }
}

void validatePlanarPressureRegionFragmentSurfaceLoads(
    const PlanarPressureRegionFragmentSurfaceLoadLedger& ledger,
    const PlanarPressureRegionFragmentPressureState& pressureState,
    const PlanarPressureRegionFragmentSurfaceLoadLimits& limits) {
    validateLimits(limits);
    if (ledger.tiles.size() > limits.maximumTiles
        || ledger.surfaces.size() > limits.maximumSurfaces) {
        throw std::length_error(
            "planar regional surface-load validation limit exceeded");
    }
    if (ledger != buildLedger(pressureState, limits)) {
        throw std::invalid_argument(
            "planar regional surface-load ledger is corrupted");
    }
}

} // namespace simwing::fsi::fluid
