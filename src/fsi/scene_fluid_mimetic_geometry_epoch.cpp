#include "scene_fluid_mimetic_geometry_epoch.h"

#include <bit>
#include <cmath>
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
            value_ ^= static_cast<std::uint8_t>(value & 0xffU);
            value_ *= fnvPrime;
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
    std::uint64_t value_ = fnvOffsetBasis;
};

void addBytes(const std::size_t bytes, std::size_t& total) {
    if (bytes > std::numeric_limits<std::size_t>::max() - total) {
        throw std::length_error(
            "scene fluid mimetic geometry-epoch storage overflows");
    }
    total += bytes;
}

template<typename Value>
void addVectorBytes(const std::vector<Value>& values,
                    std::size_t& total) {
    if (values.size()
        > std::numeric_limits<std::size_t>::max() / sizeof(Value)) {
        throw std::length_error(
            "scene fluid mimetic geometry-epoch storage overflows");
    }
    addBytes(values.size() * sizeof(Value), total);
}

std::size_t ownedStorageBytes(
    const SceneFluidMimeticGeometryEpoch& epoch) {
    std::size_t result = 0;
    addBytes(epoch.gridEpoch.ownedStorageBytes, result);
    addVectorBytes(epoch.openingCaps.caps, result);
    addVectorBytes(epoch.openingCaps.triangles, result);
    addBytes(epoch.openingQuadrature.ownedStorageBytes, result);
    addBytes(epoch.openingPatches.ownedStorageBytes, result);
    addBytes(epoch.openingFaceCrossings.ownedStorageBytes, result);
    addBytes(epoch.cappedFacePartitions.ownedStorageBytes, result);
    addVectorBytes(epoch.cellVolumes.cells, result);
    addVectorBytes(epoch.cellVolumes.cellRegionVolumes, result);
    addVectorBytes(epoch.cellVolumes.regionVolumes, result);
    addBytes(epoch.pressureControlVolumes.ownedStorageBytes, result);
    addBytes(epoch.pressureFaceLinks.ownedStorageBytes, result);
    return result;
}

std::uint64_t epochFingerprint(
    const SceneFluidMimeticGeometryEpoch& epoch) {
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
    for (const std::uint64_t value : {
             epoch.gridEpoch.fingerprint,
             epoch.openingCaps.fingerprint,
             epoch.openingQuadrature.fingerprint,
             epoch.openingPatches.fingerprint,
             epoch.openingFaceCrossings.fingerprint,
             epoch.cappedFacePartitions.fingerprint,
             epoch.cellVolumes.fingerprint,
             epoch.pressureControlVolumes.fingerprint,
             epoch.pressureFaceLinks.fingerprint}) {
        fingerprint.integer(value);
    }
    return fingerprint.value();
}

void validateIdentity(
    const SceneFluidMimeticGeometryEpoch& epoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity) {
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidRegionConnectivity(connectivity, surface);
    if (epoch.version != sceneFluidMimeticGeometryEpochVersion
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
        || epoch.cellCounts != grid.cellCounts()
        || epoch.lowerMeters != grid.lowerMeters()
        || epoch.upperMeters != grid.upperMeters()
        || epoch.ownedStorageBytes != ownedStorageBytes(epoch)
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "scene fluid mimetic geometry-epoch identity is invalid");
    }
}

} // namespace

SceneFluidMimeticGeometryEpoch buildSceneFluidMimeticGeometryEpoch(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidMimeticGeometryEpochSettings& settings,
    const SceneFluidMimeticGeometryEpochLimits& limits) {
    if (limits.maximumEpochBytes == 0) {
        throw std::invalid_argument(
            "scene fluid mimetic geometry-epoch limit is invalid");
    }
    validateSceneFluidSurfaceState(surface, state);
    validateSceneFluidRegionConnectivity(connectivity, surface);

    SceneFluidMimeticGeometryEpoch result;
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
    result.openingFaceCrossings = buildSceneFluidOpeningFaceCrossings(
        surface, state, result.openingCaps, result.openingQuadrature,
        result.openingPatches, grid, limits.openingFaceCrossings);
    result.cappedFacePartitions = buildSceneFluidCappedFacePartitions(
        surface, state, grid, transfer, result.gridEpoch,
        result.openingCaps, result.openingQuadrature,
        result.openingPatches, result.openingFaceCrossings,
        settings.cappedFacePartitions, limits.cappedFacePartitions);
    result.cellVolumes = buildSceneFluidCellVolumes(
        surface, state, grid, transfer, result.gridEpoch,
        settings.cellVolumes, limits.cellVolumes);
    result.pressureControlVolumes = buildSceneFluidPressureControlVolumes(
        surface, result.cellVolumes, connectivity,
        limits.pressureControlVolumes);
    result.pressureFaceLinks = buildSceneFluidPressureFaceLinks(
        surface, state, grid, transfer, result.gridEpoch,
        result.openingCaps, result.openingQuadrature,
        result.openingPatches, result.openingFaceCrossings,
        result.cappedFacePartitions, result.cellVolumes, connectivity,
        result.pressureControlVolumes, settings.faceLinks,
        limits.faceLinks);
    result.ownedStorageBytes = ownedStorageBytes(result);
    if (result.ownedStorageBytes > limits.maximumEpochBytes) {
        throw std::length_error(
            "scene fluid mimetic geometry epoch exceeds its aggregate byte limit");
    }
    result.fingerprint = epochFingerprint(result);
    validateSceneFluidMimeticGeometryEpoch(
        result, surface, state, grid, transfer, connectivity);
    return result;
}

void validateSceneFluidMimeticGeometryEpochIntegrity(
    const SceneFluidMimeticGeometryEpoch& epoch) {
    validateSceneFluidCellVolumeIntegrity(epoch.cellVolumes);
    validateSceneFluidPressureControlVolumeIntegrity(
        epoch.pressureControlVolumes);
    validateSceneFluidPressureFaceLinkIntegrity(epoch.pressureFaceLinks);
    if (epoch.version != sceneFluidMimeticGeometryEpochVersion
        || epoch.fingerprint == 0
        || epoch.surfaceDefinitionFingerprint == 0
        || epoch.surfaceStateFingerprint == 0
        || epoch.structureDefinitionFingerprint == 0
        || epoch.regionConnectivityFingerprint == 0
        || !std::isfinite(epoch.simulationTimeSeconds)
        || epoch.simulationTimeSeconds < 0.0
        || epoch.cellCounts.x == 0 || epoch.cellCounts.y == 0
        || epoch.cellCounts.z == 0
        || !(epoch.upperMeters.x > epoch.lowerMeters.x)
        || !(epoch.upperMeters.y > epoch.lowerMeters.y)
        || !(epoch.upperMeters.z > epoch.lowerMeters.z)
        || epoch.gridEpoch.fingerprint == 0
        || epoch.gridEpoch.quadrature.surfaceDefinitionFingerprint
            != epoch.surfaceDefinitionFingerprint
        || epoch.gridEpoch.quadrature.surfaceStateFingerprint
            != epoch.surfaceStateFingerprint
        || epoch.gridEpoch.quadrature.structureDefinitionFingerprint
            != epoch.structureDefinitionFingerprint
        || epoch.gridEpoch.quadrature.acceptedStepCount
            != epoch.acceptedStepCount
        || epoch.gridEpoch.quadrature.simulationTimeSeconds
            != epoch.simulationTimeSeconds
        || epoch.openingCaps.surfaceStateFingerprint
            != epoch.surfaceStateFingerprint
        || epoch.openingQuadrature.openingCapFingerprint
            != epoch.openingCaps.fingerprint
        || epoch.openingPatches.openingQuadratureFingerprint
            != epoch.openingQuadrature.fingerprint
        || epoch.cellVolumes.gridEpochFingerprint
            != epoch.gridEpoch.fingerprint
        || epoch.pressureControlVolumes.cellVolumeFingerprint
            != epoch.cellVolumes.fingerprint
        || epoch.pressureControlVolumes.regionConnectivityFingerprint
            != epoch.regionConnectivityFingerprint
        || epoch.pressureFaceLinks.pressureControlVolumeFingerprint
            != epoch.pressureControlVolumes.fingerprint
        || epoch.pressureFaceLinks.gridEpochFingerprint
            != epoch.gridEpoch.fingerprint
        || epoch.ownedStorageBytes != ownedStorageBytes(epoch)
        || epoch.fingerprint != epochFingerprint(epoch)) {
        throw std::invalid_argument(
            "scene fluid mimetic geometry-epoch integrity is invalid");
    }
}

void validateSceneFluidMimeticGeometryEpoch(
    const SceneFluidMimeticGeometryEpoch& epoch,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidRegionConnectivity& connectivity) {
    validateIdentity(epoch, surface, state, grid, transfer, connectivity);
    validateSceneFluidGridEpoch(
        epoch.gridEpoch, surface, state, grid, transfer);
    validateSceneFluidOpeningCaps(epoch.openingCaps, surface, state);
    validateSceneFluidOpeningQuadrature(
        epoch.openingQuadrature, surface, state, epoch.openingCaps);
    validateSceneFluidOpeningGridPatches(
        epoch.openingPatches, surface, state, epoch.openingCaps,
        epoch.openingQuadrature, grid);
    validateSceneFluidOpeningFaceCrossings(
        epoch.openingFaceCrossings, surface, state, epoch.openingCaps,
        epoch.openingQuadrature, epoch.openingPatches, grid);
    validateSceneFluidCappedFacePartitions(
        epoch.cappedFacePartitions, surface, state, grid, transfer,
        epoch.gridEpoch, epoch.openingCaps, epoch.openingQuadrature,
        epoch.openingPatches, epoch.openingFaceCrossings);
    validateSceneFluidCellVolumes(
        epoch.cellVolumes, surface, state, grid, transfer,
        epoch.gridEpoch);
    validateSceneFluidPressureControlVolumes(
        epoch.pressureControlVolumes, surface, epoch.cellVolumes,
        connectivity);
    validateSceneFluidPressureFaceLinks(
        epoch.pressureFaceLinks, surface, state, grid, transfer,
        epoch.gridEpoch, epoch.openingCaps, epoch.openingQuadrature,
        epoch.openingPatches, epoch.openingFaceCrossings,
        epoch.cappedFacePartitions, epoch.cellVolumes, connectivity,
        epoch.pressureControlVolumes);
}

} // namespace simwing::fsi
