#include "scene_fluid_pressure_epoch.h"

#include <bit>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi {
namespace {

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    template<typename Unsigned>
    void integer(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        for (std::size_t index = 0; index < sizeof(value); ++index) {
            byte(static_cast<std::uint8_t>(value & 0xffU));
            value >>= 8U;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

bool checkedMultiply(const std::size_t first,
                     const std::size_t second,
                     std::size_t& result) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        return false;
    }
    result = first * second;
    return true;
}

template<typename Value>
void addVectorBytes(const std::vector<Value>& values,
                    std::size_t& total) {
    std::size_t bytes = 0;
    std::size_t next = 0;
    if (!checkedMultiply(values.size(), sizeof(Value), bytes)
        || !checkedAdd(total, bytes, next)) {
        throw std::length_error(
            "scene fluid pressure-epoch storage size overflows");
    }
    total = next;
}

void addOwnedBytes(const std::size_t bytes, std::size_t& total) {
    std::size_t next = 0;
    if (!checkedAdd(total, bytes, next)) {
        throw std::length_error(
            "scene fluid pressure-epoch storage size overflows");
    }
    total = next;
}

std::size_t ownedStorageBytes(const SceneFluidPressureEpoch& epoch) {
    std::size_t result = 0;
    addOwnedBytes(epoch.gridEpoch.ownedStorageBytes, result);
    addVectorBytes(epoch.openingCaps.caps, result);
    addVectorBytes(epoch.openingCaps.triangles, result);
    addOwnedBytes(epoch.openingQuadrature.ownedStorageBytes, result);
    addOwnedBytes(epoch.openingPatches.ownedStorageBytes, result);
    addVectorBytes(epoch.cellVolumes.cells, result);
    addVectorBytes(epoch.cellVolumes.cellRegionVolumes, result);
    addVectorBytes(epoch.cellVolumes.regionVolumes, result);
    addOwnedBytes(epoch.pressureControlVolumes.ownedStorageBytes, result);
    addOwnedBytes(epoch.pressureFaceLinks.ownedStorageBytes, result);
    addOwnedBytes(epoch.pressureOperator.ownedStorageBytes, result);
    return result;
}

std::uint64_t epochFingerprint(const SceneFluidPressureEpoch& epoch) {
    Fingerprint fingerprint;
    fingerprint.integer(epoch.version);
    fingerprint.integer(epoch.surfaceDefinitionFingerprint);
    fingerprint.integer(epoch.surfaceStateFingerprint);
    fingerprint.integer(epoch.structureDefinitionFingerprint);
    fingerprint.integer(epoch.regionConnectivityFingerprint);
    fingerprint.integer(epoch.acceptedStepCount);
    fingerprint.real(epoch.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(epoch.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(epoch.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(epoch.cellCounts.z));
    fingerprint.real(epoch.lowerMeters.x);
    fingerprint.real(epoch.lowerMeters.y);
    fingerprint.real(epoch.lowerMeters.z);
    fingerprint.real(epoch.upperMeters.x);
    fingerprint.real(epoch.upperMeters.y);
    fingerprint.real(epoch.upperMeters.z);
    fingerprint.integer(static_cast<std::uint64_t>(epoch.ownedStorageBytes));
    fingerprint.integer(epoch.gridEpoch.fingerprint);
    fingerprint.integer(epoch.openingCaps.fingerprint);
    fingerprint.integer(epoch.openingQuadrature.fingerprint);
    fingerprint.integer(epoch.openingPatches.fingerprint);
    fingerprint.integer(epoch.cellVolumes.fingerprint);
    fingerprint.integer(epoch.pressureControlVolumes.fingerprint);
    fingerprint.integer(epoch.pressureFaceLinks.fingerprint);
    fingerprint.integer(epoch.pressureOperator.fingerprint);
    return fingerprint.value();
}

void validateIdentity(const SceneFluidPressureEpoch& epoch,
                      const SceneFluidSurfaceDefinition& surface,
                      const SceneFluidSurfaceState& state,
                      const fluid::PeriodicCartesianGrid& grid,
                      const SceneFluidSurfaceTransfer& transfer,
                      const SceneFluidRegionConnectivity& connectivity) {
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidRegionConnectivity(connectivity, surface);
    const auto gridCellCounts = grid.cellCounts();
    const auto gridLowerMeters = grid.lowerMeters();
    const auto gridUpperMeters = grid.upperMeters();
    if (epoch.version != sceneFluidPressureEpochVersion
        || epoch.fingerprint == 0
        || epoch.surfaceDefinitionFingerprint != surface.fingerprint
        || epoch.surfaceStateFingerprint != state.fingerprint
        || epoch.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || epoch.structureDefinitionFingerprint
            != transfer.targetDefinitionFingerprint()
        || epoch.regionConnectivityFingerprint != connectivity.fingerprint
        || epoch.acceptedStepCount != state.acceptedStepCount
        || epoch.simulationTimeSeconds != state.simulationTimeSeconds
        || epoch.cellCounts != gridCellCounts
        || epoch.lowerMeters != gridLowerMeters
        || epoch.upperMeters != gridUpperMeters
        || epoch.ownedStorageBytes != ownedStorageBytes(epoch)
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "scene fluid pressure-epoch identity is invalid");
    }
}

} // namespace

SceneFluidPressureEpoch buildSceneFluidPressureEpoch(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureEpochSettings& settings,
    const SceneFluidPressureEpochLimits& limits) {
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidRegionConnectivity(connectivity, surface);

    SceneFluidPressureEpoch result;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = state.fingerprint;
    result.structureDefinitionFingerprint =
        state.structureDefinitionFingerprint;
    result.regionConnectivityFingerprint = connectivity.fingerprint;
    result.acceptedStepCount = state.acceptedStepCount;
    result.simulationTimeSeconds = state.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.gridEpoch = buildSceneFluidGridEpoch(
        surface, state, grid, transfer, settings.gridEpoch,
        limits.gridEpoch);
    result.openingCaps = buildSceneFluidOpeningCaps(
        surface, state, settings.cellVolumes.openingCaps,
        limits.cellVolumes.openingCaps);
    result.openingQuadrature = buildSceneFluidOpeningQuadrature(
        surface, state, result.openingCaps, limits.openingQuadrature);
    result.openingPatches = buildSceneFluidOpeningGridPatches(
        surface, state, result.openingCaps, result.openingQuadrature, grid,
        settings.openingPatches, limits.openingPatches);
    result.cellVolumes = buildSceneFluidCellVolumes(
        surface, state, grid, transfer, result.gridEpoch,
        settings.cellVolumes, limits.cellVolumes);
    result.pressureControlVolumes = buildSceneFluidPressureControlVolumes(
        surface, result.cellVolumes, connectivity,
        limits.pressureControlVolumes);
    result.pressureFaceLinks = buildSceneFluidPressureFaceLinks(
        surface, state, grid, transfer, result.gridEpoch,
        result.openingCaps, result.openingQuadrature,
        result.openingPatches, result.cellVolumes, connectivity,
        result.pressureControlVolumes, settings.faceLinks,
        limits.faceLinks);
    result.pressureOperator = buildSceneFluidPressureOperator(
        surface, state, grid, transfer, result.gridEpoch,
        result.openingCaps, result.openingQuadrature,
        result.openingPatches, result.cellVolumes, connectivity,
        result.pressureControlVolumes, result.pressureFaceLinks,
        limits.pressureOperator);
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumEpochBytes) {
        throw std::length_error(
            "scene fluid pressure epoch exceeds its aggregate byte limit");
    }
    result.fingerprint = epochFingerprint(result);
    validateIdentity(
        result, surface, state, grid, transfer, connectivity);
    return result;
}

void validateSceneFluidPressureEpoch(
    const SceneFluidPressureEpoch& epoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity) {
    validateIdentity(
        epoch, surface, state, grid, transfer, connectivity);
    validateSceneFluidGridEpoch(
        epoch.gridEpoch, surface, state, grid, transfer);
    validateSceneFluidOpeningCaps(epoch.openingCaps, surface, state);
    validateSceneFluidOpeningQuadrature(
        epoch.openingQuadrature, surface, state, epoch.openingCaps);
    validateSceneFluidOpeningGridPatches(
        epoch.openingPatches, surface, state, epoch.openingCaps,
        epoch.openingQuadrature, grid);
    validateSceneFluidCellVolumes(
        epoch.cellVolumes, surface, state, grid, transfer,
        epoch.gridEpoch);
    validateSceneFluidPressureControlVolumes(
        epoch.pressureControlVolumes, surface, epoch.cellVolumes,
        connectivity);
    validateSceneFluidPressureFaceLinks(
        epoch.pressureFaceLinks, surface, state, grid, transfer,
        epoch.gridEpoch, epoch.openingCaps, epoch.openingQuadrature,
        epoch.openingPatches, epoch.cellVolumes, connectivity,
        epoch.pressureControlVolumes);
    validateSceneFluidPressureOperator(
        epoch.pressureOperator, surface, state, grid, transfer,
        epoch.gridEpoch, epoch.openingCaps, epoch.openingQuadrature,
        epoch.openingPatches, epoch.cellVolumes, connectivity,
        epoch.pressureControlVolumes, epoch.pressureFaceLinks);
}

} // namespace simwing::fsi
