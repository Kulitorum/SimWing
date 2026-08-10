#include "scene_fluid_regional_pressure_sampling.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <compare>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
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

struct TileKey {
    StableId surfaceStableId = invalidStableId;
    std::size_t firstTransverseCoordinate = 0;
    std::size_t secondTransverseCoordinate = 0;

    auto operator<=>(const TileKey&) const = default;
};

struct CellCoordinates {
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;
};

bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

StructureVector3 converted(const fluid::Vector3& value) {
    return {value.x, value.y, value.z};
}

StructureVector3 subtract(const StructureVector3& first,
                          const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

StructureVector3 add(const StructureVector3& first,
                     const StructureVector3& second) {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

StructureVector3 scaled(const StructureVector3& value,
                        const double scale) {
    return {scale * value.x, scale * value.y, scale * value.z};
}

double dot(const StructureVector3& first,
           const StructureVector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

StructureVector3 cross(const StructureVector3& first,
                       const StructureVector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x,
    };
}

double maximumAbsolute(const StructureVector3& value) {
    return std::max({std::abs(value.x),
                     std::abs(value.y),
                     std::abs(value.z)});
}

double tolerance(const double scale) {
    return 16384.0 * std::numeric_limits<double>::epsilon()
        * std::max(std::abs(scale), 1.0);
}

std::size_t checkedMultiply(const std::size_t first,
                            const std::size_t second) {
    if (first != 0
        && second > std::numeric_limits<std::size_t>::max() / first) {
        throw std::length_error(
            "scene fluid regional pressure-sampling storage overflows");
    }
    return first * second;
}

std::size_t checkedAdd(const std::size_t first,
                       const std::size_t second) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        throw std::length_error(
            "scene fluid regional pressure-sampling storage overflows");
    }
    return first + second;
}

std::size_t ownedStorageBytes(
    const SceneFluidRegionalPressureSampleSet& samples) {
    return checkedAdd(
        checkedAdd(
            checkedMultiply(
                samples.bindings.size(),
                sizeof(SceneFluidRegionalPressureSampleBinding)),
            checkedMultiply(
                samples.pressures.size(),
                sizeof(SceneFluidQuadraturePressure))),
        checkedMultiply(
            samples.tiles.size(),
            sizeof(SceneFluidRegionalPressureTileCoverage)));
}

std::size_t workingStorageBytes(const std::size_t sampleCount,
                                const std::size_t tileCount) {
    const std::size_t perTile = checkedAdd(
        sizeof(std::pair<TileKey, std::size_t>),
        checkedAdd(sizeof(double), sizeof(std::size_t)));
    return checkedAdd(
        checkedMultiply(
            sampleCount, sizeof(SceneFluidQuadratureKinematics)),
        checkedMultiply(tileCount, perTile));
}

void fingerprintVector(Fingerprint& fingerprint,
                       const StructureVector3& value) {
    fingerprint.real(value.x);
    fingerprint.real(value.y);
    fingerprint.real(value.z);
}

std::uint64_t samplingFingerprint(
    const SceneFluidRegionalPressureSampleSet& samples) {
    Fingerprint fingerprint;
    fingerprint.integer(samples.version);
    fingerprint.integer(samples.regionalAcceptedStateFingerprint);
    fingerprint.integer(samples.regionalOpeningLoadStateFingerprint);
    fingerprint.integer(samples.regionalPressureStateFingerprint);
    fingerprint.integer(samples.regionalSurfaceLoadFingerprint);
    fingerprint.integer(samples.regionalTopologyFingerprint);
    fingerprint.integer(samples.quadratureFingerprint);
    fingerprint.integer(samples.surfaceDefinitionFingerprint);
    fingerprint.integer(samples.surfaceStateFingerprint);
    fingerprint.integer(samples.structureDefinitionFingerprint);
    fingerprint.integer(samples.acceptedStepCount);
    fingerprint.real(samples.simulationTimeSeconds);
    fingerprint.boolean(samples.openingAware);
    fingerprint.boolean(samples.staticGeometry);
    fingerprint.boolean(samples.usesMovingVolumeRates);
    fingerprint.real(samples.timeStepSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(samples.bindings.size()));
    for (const auto& binding : samples.bindings) {
        fingerprint.integer(static_cast<std::uint64_t>(binding.sampleIndex));
        fingerprint.integer(binding.stableId);
        fingerprint.integer(binding.triangleId);
        fingerprint.integer(static_cast<std::uint64_t>(
            binding.surfaceLoadTileIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            binding.sourcePressureWallIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            binding.sourceFaceLinkIndex));
        fingerprint.integer(binding.sourceFaceLinkStableId);
        fingerprint.integer(binding.surfaceStableId);
        fingerprint.integer(binding.negativeSideRegionId);
        fingerprint.integer(binding.positiveSideRegionId);
        fingerprint.real(binding.areaSquareMeters);
        fingerprint.real(binding.negativeSidePressurePascals);
        fingerprint.real(binding.positiveSidePressurePascals);
        fingerprint.real(binding.pressureDifferencePascals);
    }
    fingerprint.integer(static_cast<std::uint64_t>(samples.pressures.size()));
    for (const auto& pressure : samples.pressures) {
        fingerprint.integer(pressure.stableId);
        fingerprint.real(pressure.negativeSidePressurePascals);
        fingerprint.real(pressure.positiveSidePressurePascals);
    }
    fingerprint.integer(static_cast<std::uint64_t>(samples.tiles.size()));
    for (const auto& tile : samples.tiles) {
        fingerprint.integer(static_cast<std::uint64_t>(tile.tileIndex));
        fingerprint.integer(tile.sourceFaceLinkStableId);
        fingerprint.integer(tile.surfaceStableId);
        fingerprint.integer(static_cast<std::uint64_t>(tile.sampleCount));
        fingerprint.real(tile.sourceAreaSquareMeters);
        fingerprint.real(tile.sampledAreaSquareMeters);
        fingerprint.real(tile.areaResidualSquareMeters);
    }
    fingerprint.real(samples.sampledAreaSquareMeters);
    fingerprint.real(samples.maximumAbsoluteTileAreaResidualSquareMeters);
    fingerprintVector(fingerprint, samples.sampledPressureForceOnSheetNewtons);
    fingerprintVector(
        fingerprint, samples.sampledPressureMomentOnSheetNewtonMeters);
    fingerprint.real(samples.sampledPressurePowerToSheetWatts);
    fingerprintVector(fingerprint, samples.sourceForceResidualNewtons);
    fingerprintVector(fingerprint, samples.sourceMomentResidualNewtonMeters);
    fingerprint.real(samples.sourcePowerResidualWatts);
    fingerprint.integer(static_cast<std::uint64_t>(samples.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(samples.workingStorageBytes));
    return fingerprint.value();
}

std::size_t applicationOwnedStorageBytes(
    const SceneFluidRegionalPressureLoadApplication& application) {
    return checkedMultiply(
        application.nodeLoads.size(),
        sizeof(SceneFluidRegionalPressureAppliedNodeLoad));
}

std::size_t applicationWorkingStorageBytes(
    const std::size_t structureNodeCount) {
    return checkedMultiply(
        structureNodeCount, sizeof(StructureVector3));
}

std::uint64_t applicationFingerprint(
    const SceneFluidRegionalPressureLoadApplication& application) {
    Fingerprint fingerprint;
    fingerprint.integer(application.version);
    fingerprint.integer(application.sourceSamplingFingerprint);
    fingerprint.integer(application.sourceSurfaceStateFingerprint);
    fingerprint.integer(application.couplingSurfaceFingerprint);
    fingerprint.integer(application.targetDefinitionFingerprint);
    fingerprint.integer(application.acceptedStepCount);
    fingerprint.real(application.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        application.structureNodeCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        application.nodeLoads.size()));
    for (const auto& load : application.nodeLoads) {
        fingerprint.integer(static_cast<std::uint64_t>(load.loadIndex));
        fingerprint.integer(load.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(load.structureNode));
        fingerprintVector(fingerprint, load.priorPendingForceNewtons);
        fingerprintVector(fingerprint, load.appliedPressureForceNewtons);
        fingerprintVector(fingerprint, load.resultingPendingForceNewtons);
        fingerprintVector(fingerprint, load.applicationResidualNewtons);
    }
    fingerprintVector(fingerprint, application.priorPendingForceNewtons);
    fingerprintVector(fingerprint, application.appliedPressureForceNewtons);
    fingerprintVector(fingerprint, application.resultingPendingForceNewtons);
    fingerprintVector(fingerprint, application.applicationResidualNewtons);
    fingerprint.boolean(application.applied);
    fingerprint.integer(static_cast<std::uint64_t>(
        application.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        application.workingStorageBytes));
    return fingerprint.value();
}

void validateLimits(
    const SceneFluidRegionalPressureSamplingLimits& limits) {
    if (limits.maximumSamples == 0 || limits.maximumTiles == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "scene fluid regional pressure-sampling limits are invalid");
    }
}

void validateApplicationLimits(
    const SceneFluidRegionalPressureLoadApplicationLimits& limits) {
    if (limits.maximumNodeLoads == 0
        || limits.maximumStructureNodes == 0
        || limits.maximumOwnedBytes == 0
        || limits.maximumWorkingBytes == 0) {
        throw std::invalid_argument(
            "scene fluid regional pressure load-application limits are invalid");
    }
}

CellCoordinates cellCoordinates(
    const fluid::PeriodicCartesianGrid& grid,
    const std::size_t cellIndex) {
    const auto counts = grid.cellCounts();
    if (cellIndex >= grid.cellCount()) {
        throw std::invalid_argument(
            "scene fluid regional pressure sample references a foreign cell");
    }
    CellCoordinates result;
    result.i = cellIndex % counts.x;
    const std::size_t remaining = cellIndex / counts.x;
    result.j = remaining % counts.y;
    result.k = remaining / counts.y;
    return result;
}

TileKey tileKey(const StableId surfaceStableId,
                const fluid::GridFaceAxis axis,
                const CellCoordinates& coordinates) {
    switch (axis) {
    case fluid::GridFaceAxis::X:
        return {surfaceStableId, coordinates.j, coordinates.k};
    case fluid::GridFaceAxis::Y:
        return {surfaceStableId, coordinates.i, coordinates.k};
    case fluid::GridFaceAxis::Z:
        return {surfaceStableId, coordinates.i, coordinates.j};
    }
    throw std::invalid_argument(
        "scene fluid regional pressure sample has an invalid axis");
}

TileKey tileKey(
    const fluid::PlanarPressureRegionFragmentSurfaceLoadTile& tile,
    const fluid::PlanarPressureRegionFragmentFaceLink& link) {
    return tileKey(
        tile.surfaceStableId, tile.axis, {link.i, link.j, link.k});
}

double axisCoordinate(const StructureVector3& value,
                      const fluid::GridFaceAxis axis) {
    switch (axis) {
    case fluid::GridFaceAxis::X: return value.x;
    case fluid::GridFaceAxis::Y: return value.y;
    case fluid::GridFaceAxis::Z: return value.z;
    }
    throw std::invalid_argument(
        "scene fluid regional pressure sample has an invalid axis");
}

StructureVector3 currentTriangleNormal(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const StableId triangleId) {
    const auto triangleIndex = surface.mappings.triangleIndex(triangleId);
    if (!triangleIndex) {
        throw std::invalid_argument(
            "scene fluid regional pressure sample references a foreign triangle");
    }
    const auto& triangle = surface.triangles[*triangleIndex];
    const auto& first = state.vertices[triangle.vertexIndices[0]].positionMeters;
    const auto& second = state.vertices[triangle.vertexIndices[1]].positionMeters;
    const auto& third = state.vertices[triangle.vertexIndices[2]].positionMeters;
    const StructureVector3 firstEdge{
        second.x - first.x, second.y - first.y, second.z - first.z};
    const StructureVector3 secondEdge{
        third.x - first.x, third.y - first.y, third.z - first.z};
    const StructureVector3 raw = cross(firstEdge, secondEdge);
    const double magnitude = std::hypot(raw.x, raw.y, raw.z);
    if (!std::isfinite(magnitude) || !(magnitude > 0.0)) {
        throw std::invalid_argument(
            "scene fluid regional pressure sample has a degenerate triangle");
    }
    return scaled(raw, 1.0 / magnitude);
}

void validateSceneSources(
    const std::uint64_t sourceTopologyFingerprint,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature) {
    validateSceneFluidSurfaceDefinition(surface);
    validateSceneFluidSurfaceState(surface, surfaceState);
    validateSceneFluidQuadratureDefinition(quadrature);
    if (sourceTopologyFingerprint != topology.fingerprint
        || quadrature.surfaceDefinitionFingerprint != surface.fingerprint
        || quadrature.surfaceStateFingerprint != surfaceState.fingerprint
        || quadrature.structureDefinitionFingerprint
            != surfaceState.structureDefinitionFingerprint
        || quadrature.acceptedStepCount != surfaceState.acceptedStepCount
        || quadrature.simulationTimeSeconds
            != surfaceState.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid regional pressure-sampling source identity is invalid");
    }
}

void validateSources(
    const fluid::PlanarPressureRegionFragmentAcceptedState& acceptedState,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature) {
    fluid::validatePlanarPressureRegionFragmentAcceptedState(
        acceptedState, grid, sweep, fragments, topology, metric);
    validateSceneSources(
        acceptedState.sourceTopologyFingerprint, topology, surface,
        surfaceState, quadrature);
}

void validateOpeningSources(
    const fluid::PlanarPressureRegionFragmentOpeningLoadState& loadState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const fluid::PlanarPressureRegionFragmentOpeningLoadStateLimits& limits) {
    fluid::validatePlanarPressureRegionFragmentOpeningLoadState(
        loadState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, limits);
    validateSceneSources(
        loadState.sourceTopologyFingerprint, topology, surface, surfaceState,
        quadrature);
}

StructureVector3 fullWallMoment(
    const fluid::PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads) {
    StructureVector3 result;
    for (const auto& tile : surfaceLoads.tiles) {
        result = add(
            result,
            cross(converted(tile.wrappedCentroidMeters),
                  converted(tile.totalPressureForceOnSheetNewtons)));
    }
    return result;
}

template<typename PressureState>
SceneFluidRegionalPressureSampleSet buildSamples(
    const std::uint64_t acceptedStateFingerprint,
    const std::uint64_t openingLoadStateFingerprint,
    const PressureState& pressureState,
    const fluid::PlanarPressureRegionFragmentSurfaceLoadLedger& surfaceLoads,
    const fluid::PlanarPressureRegionFragmentOpeningSurfaceLoadLedger*
        openingSurfaceLoads,
    const bool staticGeometry,
    const bool usesMovingVolumeRates,
    const double timeStepSeconds,
    const StructureVector3& sourceForce,
    const StructureVector3& sourceMoment,
    const double sourceWorkToSheetJoules,
    const double sourceAreaSquareMeters,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSamplingLimits& limits) {
    validateLimits(limits);
    const std::size_t sampleCount = quadrature.points.size();
    const std::size_t tileCount = surfaceLoads.tiles.size();
    if (sampleCount == 0 || tileCount == 0
        || sampleCount > limits.maximumSamples
        || tileCount > limits.maximumTiles) {
        throw std::length_error(
            "scene fluid regional pressure-sampling count limit exceeded");
    }
    const std::size_t expectedOwnedBytes = checkedAdd(
        checkedAdd(
            checkedMultiply(
                sampleCount,
                sizeof(SceneFluidRegionalPressureSampleBinding)),
            checkedMultiply(
                sampleCount, sizeof(SceneFluidQuadraturePressure))),
        checkedMultiply(
            tileCount,
            sizeof(SceneFluidRegionalPressureTileCoverage)));
    const std::size_t expectedWorkingBytes = workingStorageBytes(
        sampleCount, tileCount);
    if (expectedOwnedBytes > limits.maximumOwnedBytes
        || expectedWorkingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "scene fluid regional pressure-sampling byte limit exceeded");
    }

    SceneFluidRegionalPressureSampleSet result;
    result.regionalAcceptedStateFingerprint = acceptedStateFingerprint;
    result.regionalOpeningLoadStateFingerprint =
        openingLoadStateFingerprint;
    result.regionalPressureStateFingerprint = pressureState.fingerprint;
    result.regionalSurfaceLoadFingerprint = surfaceLoads.fingerprint;
    result.regionalTopologyFingerprint = topology.fingerprint;
    result.quadratureFingerprint = quadrature.fingerprint;
    result.surfaceDefinitionFingerprint = surface.fingerprint;
    result.surfaceStateFingerprint = surfaceState.fingerprint;
    result.structureDefinitionFingerprint =
        surfaceState.structureDefinitionFingerprint;
    result.acceptedStepCount = surfaceState.acceptedStepCount;
    result.simulationTimeSeconds = surfaceState.simulationTimeSeconds;
    result.openingAware = openingSurfaceLoads != nullptr;
    result.staticGeometry = staticGeometry;
    result.usesMovingVolumeRates = usesMovingVolumeRates;
    result.timeStepSeconds = timeStepSeconds;
    result.bindings.reserve(sampleCount);
    result.pressures.reserve(sampleCount);
    result.tiles.resize(tileCount);

    std::vector<std::pair<TileKey, std::size_t>> tileByKey;
    tileByKey.reserve(tileCount);
    for (std::size_t index = 0; index < tileCount; ++index) {
        const auto& tile = surfaceLoads.tiles[index];
        if (tile.tileIndex != index
            || tile.sourceFaceLinkIndex >= topology.links.size()
            || tile.sourcePressureWallIndex
                >= pressureState.walls.size()) {
            throw std::invalid_argument(
                "scene fluid regional pressure tile provenance is invalid");
        }
        const auto& link = topology.links[tile.sourceFaceLinkIndex];
        const auto& wall =
            pressureState.walls[tile.sourcePressureWallIndex];
        if (link.kind
                != fluid::PlanarPressureRegionFragmentFaceKind::
                    PressureLayerWall
            || link.stableId != tile.sourceFaceLinkStableId
            || link.surfaceStableId != tile.surfaceStableId
            || link.axis != tile.axis
            || wall.sourceFaceLinkIndex != tile.sourceFaceLinkIndex
            || wall.sourceFaceLinkStableId != tile.sourceFaceLinkStableId) {
            throw std::invalid_argument(
                "scene fluid regional pressure tile is foreign to topology");
        }
        tileByKey.emplace_back(tileKey(tile, link), index);
        double sourceTileArea = tile.areaSquareMeters;
        if (openingSurfaceLoads != nullptr) {
            const auto& openingTile = openingSurfaceLoads->tiles[index];
            if (openingTile.tileIndex != index
                || openingTile.sourceSurfaceLoadTileIndex != index
                || openingTile.sourceFaceLinkStableId
                    != tile.sourceFaceLinkStableId
                || openingTile.surfaceStableId != tile.surfaceStableId) {
                throw std::invalid_argument(
                    "scene fluid regional opening load tile is foreign");
            }
            sourceTileArea = openingTile.solidAreaSquareMeters;
        }
        result.tiles[index] = {
            index,
            tile.sourceFaceLinkStableId,
            tile.surfaceStableId,
            0,
            sourceTileArea,
            0.0,
            0.0,
        };
    }
    std::ranges::sort(
        tileByKey, {}, &std::pair<TileKey, std::size_t>::first);
    if (std::ranges::adjacent_find(
            tileByKey,
            [](const auto& first, const auto& second) {
                return first.first == second.first;
            }) != tileByKey.end()) {
        throw std::invalid_argument(
            "scene fluid regional pressure tile key is ambiguous");
    }

    const auto kinematics = sampleSceneFluidQuadratureKinematics(
        surface, surfaceState, quadrature);
    if (kinematics.size() != sampleCount) {
        throw std::logic_error(
            "scene fluid regional quadrature kinematics changed size");
    }
    for (std::size_t index = 0; index < sampleCount; ++index) {
        const auto& point = quadrature.points[index];
        const auto& motion = kinematics[index];
        if (point.negativeSideCellIndex != point.positiveSideCellIndex
            || motion.stableId != point.stableId) {
            throw std::invalid_argument(
                "scene fluid regional pressure requires a subcell sheet patch");
        }
        const auto coordinates = cellCoordinates(
            grid, point.negativeSideCellIndex);
        const TileKey key = tileKey(
            point.sheetId, topology.profileAxis, coordinates);
        const auto found = std::ranges::lower_bound(
            tileByKey, key, {},
            &std::pair<TileKey, std::size_t>::first);
        if (found == tileByKey.end()) {
            throw std::invalid_argument(
                "scene fluid regional quadrature patch has no pressure tile");
        }
        if (found->first != key) {
            throw std::invalid_argument(
                "scene fluid regional quadrature patch has no pressure tile");
        }
        const std::size_t tileIndex = found->second;
        const auto& tile = surfaceLoads.tiles[tileIndex];
        const auto& wall =
            pressureState.walls[tile.sourcePressureWallIndex];
        if (point.negativeSideRegionId != tile.minusRegionStableId
            || point.positiveSideRegionId != tile.plusRegionStableId) {
            throw std::invalid_argument(
                "scene fluid regional quadrature has reversed pressure sides");
        }
        const StructureVector3 tileNormal = converted(
            tile.unitNormalMinusToPlus);
        const StructureVector3 normal = currentTriangleNormal(
            surface, surfaceState, point.triangleId);
        const StructureVector3 normalResidual = subtract(normal, tileNormal);
        const double coordinate = axisCoordinate(
            motion.positionMeters, tile.axis);
        const double sourceCoordinate = axisCoordinate(
            converted(tile.wrappedCentroidMeters), tile.axis);
        const double normalVelocity = dot(
            motion.velocityMetersPerSecond, tileNormal);
        const double geometryScale = std::max({
            std::abs(coordinate), std::abs(sourceCoordinate), 1.0});
        const double velocityScale = std::max({
            std::abs(normalVelocity),
            std::abs(wall.materialWallVelocityMetersPerSecond), 1.0});
        if (maximumAbsolute(normalResidual) > tolerance(1.0)
            || std::abs(coordinate - sourceCoordinate)
                > tolerance(geometryScale)
            || std::abs(normalVelocity
                        - wall.materialWallVelocityMetersPerSecond)
                > tolerance(velocityScale)) {
            throw std::invalid_argument(
                "scene fluid regional pressure patch geometry is incompatible");
        }

        const double negativePressure = wall.minusTotalPressurePascals;
        const double positivePressure = wall.plusTotalPressurePascals;
        const double difference = negativePressure - positivePressure;
        const double expectedDifference = dot(
            converted(tile.totalPressureTractionOnSheetPascals), tileNormal);
        if (!std::isfinite(difference)
            || std::abs(difference - expectedDifference)
                > tolerance(std::max({
                    std::abs(difference),
                    std::abs(expectedDifference), 1.0}))) {
            throw std::invalid_argument(
                "scene fluid regional pressure tile traction is inconsistent");
        }

        result.bindings.push_back({
            index,
            point.stableId,
            point.triangleId,
            tileIndex,
            tile.sourcePressureWallIndex,
            tile.sourceFaceLinkIndex,
            tile.sourceFaceLinkStableId,
            tile.surfaceStableId,
            tile.minusRegionStableId,
            tile.plusRegionStableId,
            point.areaSquareMeters,
            negativePressure,
            positivePressure,
            difference,
        });
        result.pressures.push_back({
            point.stableId, negativePressure, positivePressure});
        auto& coverage = result.tiles[tileIndex];
        ++coverage.sampleCount;
        coverage.sampledAreaSquareMeters += point.areaSquareMeters;
        result.sampledAreaSquareMeters += point.areaSquareMeters;

        const StructureVector3 force = scaled(
            normal, point.areaSquareMeters * difference);
        result.sampledPressureForceOnSheetNewtons = add(
            result.sampledPressureForceOnSheetNewtons, force);
        result.sampledPressureMomentOnSheetNewtonMeters = add(
            result.sampledPressureMomentOnSheetNewtonMeters,
            cross(motion.positionMeters, force));
        result.sampledPressurePowerToSheetWatts +=
            dot(force, motion.velocityMetersPerSecond);
    }

    for (std::size_t index = 0; index < tileCount; ++index) {
        auto& coverage = result.tiles[index];
        coverage.areaResidualSquareMeters =
            coverage.sampledAreaSquareMeters
            - coverage.sourceAreaSquareMeters;
        result.maximumAbsoluteTileAreaResidualSquareMeters = std::max(
            result.maximumAbsoluteTileAreaResidualSquareMeters,
            std::abs(coverage.areaResidualSquareMeters));
        const double areaScale = std::max({
            coverage.sampledAreaSquareMeters,
            coverage.sourceAreaSquareMeters, 1.0});
        const bool zeroSolidOpeningTile = result.openingAware
            && coverage.sourceAreaSquareMeters == 0.0;
        if ((zeroSolidOpeningTile
             && (coverage.sampleCount != 0
                 || coverage.sampledAreaSquareMeters != 0.0))
            || (!zeroSolidOpeningTile
                && (coverage.sampleCount == 0
                    || !(coverage.sourceAreaSquareMeters > 0.0)))
            || std::abs(coverage.areaResidualSquareMeters)
                > tolerance(areaScale)) {
            throw std::invalid_argument(
                "scene fluid regional pressure tile coverage is incomplete");
        }
    }

    const double sourcePower =
        sourceWorkToSheetJoules / timeStepSeconds;
    result.sourceForceResidualNewtons = subtract(
        result.sampledPressureForceOnSheetNewtons, sourceForce);
    result.sourceMomentResidualNewtonMeters = subtract(
        result.sampledPressureMomentOnSheetNewtonMeters, sourceMoment);
    result.sourcePowerResidualWatts =
        result.sampledPressurePowerToSheetWatts - sourcePower;
    const double forceScale = std::max({
        maximumAbsolute(result.sampledPressureForceOnSheetNewtons),
        maximumAbsolute(sourceForce), 1.0});
    const double momentScale = std::max({
        maximumAbsolute(result.sampledPressureMomentOnSheetNewtonMeters),
        maximumAbsolute(sourceMoment), 1.0});
    const double powerScale = std::max({
        std::abs(result.sampledPressurePowerToSheetWatts),
        std::abs(sourcePower), 1.0});
    if (std::abs(result.sampledAreaSquareMeters
                 - sourceAreaSquareMeters)
            > tolerance(std::max({
                result.sampledAreaSquareMeters,
                sourceAreaSquareMeters, 1.0}))
        || maximumAbsolute(result.sourceForceResidualNewtons)
            > tolerance(forceScale)
        || maximumAbsolute(result.sourceMomentResidualNewtonMeters)
            > tolerance(momentScale)
        || std::abs(result.sourcePowerResidualWatts)
            > tolerance(powerScale)) {
        throw std::invalid_argument(
            "scene fluid regional pressure sampling does not close its source ledger");
    }

    result.ownedStorageBytes = ownedStorageBytes(result);
    result.workingStorageBytes = expectedWorkingBytes;
    if (result.ownedStorageBytes != expectedOwnedBytes) {
        throw std::logic_error(
            "scene fluid regional pressure-sampling storage count changed");
    }
    result.fingerprint = samplingFingerprint(result);
    return result;
}

} // namespace

SceneFluidRegionalPressureSampleSet
sampleSceneFluidRegionalAcceptedPressure(
    const fluid::PlanarPressureRegionFragmentAcceptedState& acceptedState,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSamplingLimits& limits) {
    validateSources(
        acceptedState, grid, sweep, fragments, topology, metric,
        surface, surfaceState, quadrature);
    auto result = buildSamples(
        acceptedState.fingerprint, 0, acceptedState.pressure,
        acceptedState.surfaceLoads, nullptr, acceptedState.staticGeometry,
        acceptedState.usesMovingVolumeRates, acceptedState.timeStepSeconds,
        converted(acceptedState.pressureForceOnSheetNewtons),
        fullWallMoment(acceptedState.surfaceLoads),
        acceptedState.pressureWorkToSheetJoules,
        acceptedState.surfaceLoads.totalAreaSquareMeters, grid, topology,
        surface, surfaceState, quadrature, limits);
    validateSceneFluidRegionalPressureSampleIntegrity(result);
    return result;
}

SceneFluidRegionalPressureSampleSet
sampleSceneFluidRegionalOpeningPressure(
    const fluid::PlanarPressureRegionFragmentOpeningLoadState& loadState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const fluid::PlanarPressureRegionFragmentOpeningLoadStateLimits&
        loadStateLimits,
    const SceneFluidRegionalPressureSamplingLimits& limits) {
    validateOpeningSources(
        loadState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, surface, surfaceState, quadrature,
        loadStateLimits);
    auto result = buildSamples(
        0, loadState.fingerprint, loadState.pressure,
        loadState.surfaceLoads, &loadState.openingSurfaceLoads,
        loadState.staticGeometry, loadState.usesMovingVolumeRates,
        loadState.timeStepSeconds,
        converted(loadState.solidPressureForceOnSheetNewtons),
        converted(loadState.solidPressureMomentOnSheetNewtonMeters),
        loadState.solidPressureWorkToSheetJoules,
        loadState.solidAreaSquareMeters, grid, topology, surface,
        surfaceState, quadrature, limits);
    validateSceneFluidRegionalPressureSampleIntegrity(result);
    return result;
}

void validateSceneFluidRegionalPressureSampleIntegrity(
    const SceneFluidRegionalPressureSampleSet& samples) {
    if (samples.version != sceneFluidRegionalPressureSamplingVersion
        || samples.fingerprint == 0
        || (samples.openingAware
                ? samples.regionalAcceptedStateFingerprint != 0
                    || samples.regionalOpeningLoadStateFingerprint == 0
                : samples.regionalAcceptedStateFingerprint == 0
                    || samples.regionalOpeningLoadStateFingerprint != 0)
        || samples.regionalPressureStateFingerprint == 0
        || samples.regionalSurfaceLoadFingerprint == 0
        || samples.regionalTopologyFingerprint == 0
        || samples.quadratureFingerprint == 0
        || samples.surfaceDefinitionFingerprint == 0
        || samples.surfaceStateFingerprint == 0
        || samples.structureDefinitionFingerprint == 0
        || !std::isfinite(samples.simulationTimeSeconds)
        || !std::isfinite(samples.timeStepSeconds)
        || !(samples.timeStepSeconds > 0.0)
        || (!samples.openingAware
            && samples.staticGeometry == samples.usesMovingVolumeRates)
        || (samples.openingAware && !samples.usesMovingVolumeRates)
        || samples.bindings.empty()
        || samples.bindings.size() != samples.pressures.size()
        || samples.tiles.empty()
        || !std::isfinite(samples.sampledAreaSquareMeters)
        || !(samples.sampledAreaSquareMeters > 0.0)
        || !std::isfinite(
            samples.maximumAbsoluteTileAreaResidualSquareMeters)
        || !finite(samples.sampledPressureForceOnSheetNewtons)
        || !finite(samples.sampledPressureMomentOnSheetNewtonMeters)
        || !std::isfinite(samples.sampledPressurePowerToSheetWatts)
        || !finite(samples.sourceForceResidualNewtons)
        || !finite(samples.sourceMomentResidualNewtonMeters)
        || !std::isfinite(samples.sourcePowerResidualWatts)
        || samples.ownedStorageBytes != ownedStorageBytes(samples)
        || samples.workingStorageBytes
            != workingStorageBytes(
                samples.bindings.size(), samples.tiles.size())
        || samples.fingerprint != samplingFingerprint(samples)) {
        throw std::invalid_argument(
            "scene fluid regional pressure-sampling integrity is invalid");
    }
    std::vector<std::size_t> sampleCounts(samples.tiles.size(), 0);
    std::vector<double> sampledAreas(samples.tiles.size(), 0.0);
    double sampledArea = 0.0;
    for (std::size_t index = 0; index < samples.bindings.size(); ++index) {
        const auto& binding = samples.bindings[index];
        const auto& pressure = samples.pressures[index];
        if (binding.sampleIndex != index || binding.stableId == 0
            || binding.triangleId == invalidStableId
            || binding.surfaceLoadTileIndex >= samples.tiles.size()
            || binding.sourceFaceLinkStableId == 0
            || binding.surfaceStableId == invalidStableId
            || binding.negativeSideRegionId == invalidStableId
            || binding.positiveSideRegionId == invalidStableId
            || !std::isfinite(binding.areaSquareMeters)
            || !(binding.areaSquareMeters > 0.0)
            || !std::isfinite(binding.negativeSidePressurePascals)
            || !std::isfinite(binding.positiveSidePressurePascals)
            || !std::isfinite(binding.pressureDifferencePascals)
            || binding.pressureDifferencePascals
                != binding.negativeSidePressurePascals
                    - binding.positiveSidePressurePascals
            || pressure.stableId != binding.stableId
            || pressure.negativeSidePressurePascals
                != binding.negativeSidePressurePascals
            || pressure.positiveSidePressurePascals
                != binding.positiveSidePressurePascals) {
            throw std::invalid_argument(
                "scene fluid regional pressure sample row is invalid");
        }
        const auto& tile = samples.tiles[binding.surfaceLoadTileIndex];
        if (tile.sourceFaceLinkStableId
                != binding.sourceFaceLinkStableId
            || tile.surfaceStableId != binding.surfaceStableId) {
            throw std::invalid_argument(
                "scene fluid regional pressure sample tile binding is invalid");
        }
        ++sampleCounts[binding.surfaceLoadTileIndex];
        sampledAreas[binding.surfaceLoadTileIndex] +=
            binding.areaSquareMeters;
        sampledArea += binding.areaSquareMeters;
    }
    double maximumAreaResidual = 0.0;
    for (std::size_t index = 0; index < samples.tiles.size(); ++index) {
        const auto& tile = samples.tiles[index];
        const bool zeroSolidOpeningTile = samples.openingAware
            && tile.sourceAreaSquareMeters == 0.0;
        if (tile.tileIndex != index || tile.sourceFaceLinkStableId == 0
            || tile.surfaceStableId == invalidStableId
            || !std::isfinite(tile.sourceAreaSquareMeters)
            || !std::isfinite(tile.sampledAreaSquareMeters)
            || (zeroSolidOpeningTile
                    ? tile.sampleCount != 0
                        || tile.sampledAreaSquareMeters != 0.0
                    : tile.sampleCount == 0
                        || !(tile.sourceAreaSquareMeters > 0.0)
                        || !(tile.sampledAreaSquareMeters > 0.0))
            || !std::isfinite(tile.areaResidualSquareMeters)
            || tile.areaResidualSquareMeters
                != tile.sampledAreaSquareMeters
                    - tile.sourceAreaSquareMeters) {
            throw std::invalid_argument(
                "scene fluid regional pressure tile coverage is invalid");
        }
        if (tile.sampleCount != sampleCounts[index]
            || tile.sampledAreaSquareMeters != sampledAreas[index]) {
            throw std::invalid_argument(
                "scene fluid regional pressure tile coverage does not match samples");
        }
        maximumAreaResidual = std::max(
            maximumAreaResidual, std::abs(tile.areaResidualSquareMeters));
    }
    if (samples.sampledAreaSquareMeters != sampledArea
        || samples.maximumAbsoluteTileAreaResidualSquareMeters
            != maximumAreaResidual) {
        throw std::invalid_argument(
            "scene fluid regional pressure aggregate coverage is invalid");
    }
}

void validateSceneFluidRegionalAcceptedPressureSamples(
    const SceneFluidRegionalPressureSampleSet& samples,
    const fluid::PlanarPressureRegionFragmentAcceptedState& acceptedState,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVelocityMetric& metric,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSamplingLimits& limits) {
    validateSources(
        acceptedState, grid, sweep, fragments, topology, metric,
        surface, surfaceState, quadrature);
    validateSceneFluidRegionalPressureSampleIntegrity(samples);
    if (samples != buildSamples(
                       acceptedState.fingerprint, 0, acceptedState.pressure,
                       acceptedState.surfaceLoads, nullptr,
                       acceptedState.staticGeometry,
                       acceptedState.usesMovingVolumeRates,
                       acceptedState.timeStepSeconds,
                       converted(acceptedState.pressureForceOnSheetNewtons),
                       fullWallMoment(acceptedState.surfaceLoads),
                       acceptedState.pressureWorkToSheetJoules,
                       acceptedState.surfaceLoads.totalAreaSquareMeters,
                       grid, topology, surface, surfaceState, quadrature,
                       limits)) {
        throw std::invalid_argument(
            "scene fluid regional pressure sample payload is invalid");
    }
}

void validateSceneFluidRegionalOpeningPressureSamples(
    const SceneFluidRegionalPressureSampleSet& samples,
    const fluid::PlanarPressureRegionFragmentOpeningLoadState& loadState,
    const fluid::PlanarPressureRegionFragmentOpeningPressureOperator&
        pressureOperator,
    const fluid::PlanarPressureRegionFragmentPressureOperator&
        basePressureOperator,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::PlanarPressureRegionSweepLedger& sweep,
    const fluid::PlanarPressureRegionFragmentSet& fragments,
    const fluid::PlanarPressureRegionFragmentTopology& topology,
    const fluid::PlanarPressureRegionFragmentVolumeRateSet& volumeRates,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningPatchDefinition>
        openingDefinitions,
    const fluid::PlanarPressureRegionFragmentOpeningSet& openings,
    const std::span<
        const fluid::PlanarPressureRegionFragmentOpeningResistanceDefinition>
        resistanceDefinitions,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& surfaceState,
    const SceneFluidQuadratureDefinition& quadrature,
    const fluid::PlanarPressureRegionFragmentOpeningLoadStateLimits&
        loadStateLimits,
    const SceneFluidRegionalPressureSamplingLimits& limits) {
    validateOpeningSources(
        loadState, pressureOperator, basePressureOperator, grid, sweep,
        fragments, topology, volumeRates, openingDefinitions, openings,
        resistanceDefinitions, surface, surfaceState, quadrature,
        loadStateLimits);
    validateSceneFluidRegionalPressureSampleIntegrity(samples);
    if (samples != buildSamples(
                       0, loadState.fingerprint, loadState.pressure,
                       loadState.surfaceLoads,
                       &loadState.openingSurfaceLoads,
                       loadState.staticGeometry,
                       loadState.usesMovingVolumeRates,
                       loadState.timeStepSeconds,
                       converted(loadState.solidPressureForceOnSheetNewtons),
                       converted(
                           loadState
                               .solidPressureMomentOnSheetNewtonMeters),
                       loadState.solidPressureWorkToSheetJoules,
                       loadState.solidAreaSquareMeters, grid, topology,
                       surface, surfaceState, quadrature, limits)) {
        throw std::invalid_argument(
            "scene fluid regional opening-pressure sample payload is invalid");
    }
}

ConservativeTransferResult
evaluateSceneFluidRegionalAcceptedPressureQuadrature(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSampleSet& samples,
    const ConservativeTransferSettings& settings) {
    validateSceneFluidRegionalPressureSampleIntegrity(samples);
    if (samples.quadratureFingerprint != quadrature.fingerprint
        || samples.surfaceDefinitionFingerprint != surface.fingerprint
        || samples.surfaceStateFingerprint != state.fingerprint
        || samples.structureDefinitionFingerprint
            != state.structureDefinitionFingerprint
        || samples.acceptedStepCount != state.acceptedStepCount
        || samples.simulationTimeSeconds != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid regional pressure traction binding is invalid");
    }
    auto result = evaluateSceneFluidPressureQuadrature(
        surface, state, transfer, quadrature, samples.pressures, settings);
    const auto& diagnostics = result.diagnostics();
    const StructureVector3 forceResidual = subtract(
        diagnostics.integratedSurfaceForceNewtons,
        samples.sampledPressureForceOnSheetNewtons);
    const StructureVector3 momentResidual = subtract(
        diagnostics.integratedSurfaceMomentNewtonMeters,
        samples.sampledPressureMomentOnSheetNewtonMeters);
    const double powerResidual =
        diagnostics.integratedSurfacePowerWatts
        - samples.sampledPressurePowerToSheetWatts;
    if (maximumAbsolute(forceResidual)
            > tolerance(std::max({
                maximumAbsolute(diagnostics.integratedSurfaceForceNewtons),
                maximumAbsolute(samples.sampledPressureForceOnSheetNewtons),
                1.0}))
        || maximumAbsolute(momentResidual)
            > tolerance(std::max({
                maximumAbsolute(diagnostics.integratedSurfaceMomentNewtonMeters),
                maximumAbsolute(samples.sampledPressureMomentOnSheetNewtonMeters),
                1.0}))
        || std::abs(powerResidual)
            > tolerance(std::max({
                std::abs(diagnostics.integratedSurfacePowerWatts),
                std::abs(samples.sampledPressurePowerToSheetWatts), 1.0}))) {
        throw std::invalid_argument(
            "scene fluid regional pressure transfer changed the sampled ledger");
    }
    return result;
}

SceneFluidRegionalPressureLoadApplication
applySceneFluidRegionalAcceptedPressureLoads(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSampleSet& samples,
    Structure& target,
    const ConservativeTransferSettings& transferSettings,
    const SceneFluidRegionalPressureLoadApplicationLimits& limits) {
    validateApplicationLimits(limits);
    validateSceneFluidRegionalPressureSampleIntegrity(samples);
    const std::size_t nodeLoadCount = transfer.nodes().size();
    const std::size_t structureNodeCount = target.definition().nodes.size();
    if (nodeLoadCount == 0 || structureNodeCount == 0
        || nodeLoadCount > limits.maximumNodeLoads
        || structureNodeCount > limits.maximumStructureNodes) {
        throw std::length_error(
            "scene fluid regional pressure load-application count limit exceeded");
    }
    const std::size_t expectedOwnedBytes = checkedMultiply(
        nodeLoadCount,
        sizeof(SceneFluidRegionalPressureAppliedNodeLoad));
    const std::size_t expectedWorkingBytes =
        applicationWorkingStorageBytes(structureNodeCount);
    if (expectedOwnedBytes > limits.maximumOwnedBytes
        || expectedWorkingBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "scene fluid regional pressure load-application byte limit exceeded");
    }
    if (target.definitionFingerprint()
            != transfer.targetDefinitionFingerprint()
        || target.acceptedStepCount() != state.acceptedStepCount
        || target.simulationTimeSeconds() != state.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid regional pressure target epoch is stale");
    }
    const auto sampledKinematics = transfer.kinematics(state);
    const auto targetKinematics =
        transfer.conservativeTransfer().captureKinematics(target);
    if (sampledKinematics != targetKinematics) {
        throw std::invalid_argument(
            "scene fluid regional pressure target kinematics are stale");
    }
    const auto transferred =
        evaluateSceneFluidRegionalAcceptedPressureQuadrature(
            surface, state, transfer, quadrature, samples,
            transferSettings);
    const auto transferredLoads = transferred.nodeLoads();
    if (transferredLoads.size() != nodeLoadCount) {
        throw std::logic_error(
            "scene fluid regional pressure transfer changed node count");
    }

    const StructureCheckpoint before = target.checkpoint();
    if (before.definitionFingerprint != target.definitionFingerprint()
        || before.acceptedStepCount != state.acceptedStepCount
        || before.simulationTimeSeconds != state.simulationTimeSeconds
        || before.nodes.size() != structureNodeCount
        || before.pendingExternalForcesNewtons.size()
            != structureNodeCount) {
        throw std::invalid_argument(
            "scene fluid regional pressure target checkpoint is incompatible");
    }
    std::vector<StructureVector3> expectedPending =
        before.pendingExternalForcesNewtons;

    SceneFluidRegionalPressureLoadApplication application;
    application.sourceSamplingFingerprint = samples.fingerprint;
    application.sourceSurfaceStateFingerprint = state.fingerprint;
    application.couplingSurfaceFingerprint =
        transfer.couplingSurfaceFingerprint();
    application.targetDefinitionFingerprint =
        target.definitionFingerprint();
    application.acceptedStepCount = state.acceptedStepCount;
    application.simulationTimeSeconds = state.simulationTimeSeconds;
    application.structureNodeCount = structureNodeCount;
    application.nodeLoads.reserve(nodeLoadCount);
    for (const auto& force : before.pendingExternalForcesNewtons) {
        application.priorPendingForceNewtons = add(
            application.priorPendingForceNewtons, force);
    }
    for (std::size_t index = 0; index < nodeLoadCount; ++index) {
        const auto& load = transferredLoads[index];
        if (load.structureNode >= structureNodeCount
            || !finite(load.forceNewtons)) {
            throw std::invalid_argument(
                "scene fluid regional pressure node load is invalid");
        }
        const StructureVector3 prior = expectedPending[load.structureNode];
        const StructureVector3 resulting = add(prior, load.forceNewtons);
        if (!finite(resulting)) {
            throw std::overflow_error(
                "scene fluid regional pressure pending load is not finite");
        }
        expectedPending[load.structureNode] = resulting;
        application.appliedPressureForceNewtons = add(
            application.appliedPressureForceNewtons, load.forceNewtons);
        application.nodeLoads.push_back({
            index,
            load.stableId,
            load.structureNode,
            prior,
            load.forceNewtons,
            resulting,
            subtract(resulting, add(prior, load.forceNewtons)),
        });
    }
    for (const auto& force : expectedPending) {
        application.resultingPendingForceNewtons = add(
            application.resultingPendingForceNewtons, force);
    }
    application.applicationResidualNewtons = subtract(
        application.resultingPendingForceNewtons,
        add(application.priorPendingForceNewtons,
            application.appliedPressureForceNewtons));
    if (!finite(application.priorPendingForceNewtons)
        || !finite(application.appliedPressureForceNewtons)
        || !finite(application.resultingPendingForceNewtons)
        || !finite(application.applicationResidualNewtons)) {
        throw std::overflow_error(
            "scene fluid regional pressure load-application ledger is not finite");
    }
    application.ownedStorageBytes = applicationOwnedStorageBytes(application);
    application.workingStorageBytes = expectedWorkingBytes;
    if (application.ownedStorageBytes != expectedOwnedBytes) {
        throw std::logic_error(
            "scene fluid regional pressure load-application storage changed");
    }

    try {
        transfer.addLoadsTo(target, transferred);
        const StructureCheckpoint after = target.checkpoint();
        if (after.definitionFingerprint != before.definitionFingerprint
            || after.acceptedStepCount != before.acceptedStepCount
            || after.simulationTimeSeconds != before.simulationTimeSeconds
            || after.nodes != before.nodes
            || after.pendingExternalForcesNewtons != expectedPending
            || after.lastAppliedExternalForceNewtons
                != before.lastAppliedExternalForceNewtons) {
            throw std::logic_error(
                "scene fluid regional pressure load application changed non-load state");
        }
        application.applied = true;
        application.fingerprint = applicationFingerprint(application);
        validateSceneFluidRegionalPressureLoadApplicationIntegrity(
            application);
        validateSceneFluidRegionalPressureLoadApplication(
            application, surface, state, transfer, quadrature, samples,
            transferSettings, limits);
        return application;
    } catch (...) {
        target.restore(before);
        throw;
    }
}

void validateSceneFluidRegionalPressureLoadApplicationIntegrity(
    const SceneFluidRegionalPressureLoadApplication& application) {
    if (application.version
            != sceneFluidRegionalPressureLoadApplicationVersion
        || application.fingerprint == 0
        || application.sourceSamplingFingerprint == 0
        || application.sourceSurfaceStateFingerprint == 0
        || application.couplingSurfaceFingerprint == 0
        || application.targetDefinitionFingerprint == 0
        || !std::isfinite(application.simulationTimeSeconds)
        || application.structureNodeCount == 0
        || application.nodeLoads.empty()
        || !finite(application.priorPendingForceNewtons)
        || !finite(application.appliedPressureForceNewtons)
        || !finite(application.resultingPendingForceNewtons)
        || !finite(application.applicationResidualNewtons)
        || !application.applied
        || application.ownedStorageBytes
            != applicationOwnedStorageBytes(application)
        || application.workingStorageBytes
            != applicationWorkingStorageBytes(
                application.structureNodeCount)
        || application.fingerprint != applicationFingerprint(application)) {
        throw std::invalid_argument(
            "scene fluid regional pressure load-application integrity is invalid");
    }
    StructureVector3 appliedForce;
    std::uint64_t previousStableId = 0;
    std::size_t previousStructureNode = 0;
    bool havePrevious = false;
    for (std::size_t index = 0;
         index < application.nodeLoads.size(); ++index) {
        const auto& load = application.nodeLoads[index];
        const StructureVector3 reconstructed = add(
            load.priorPendingForceNewtons,
            load.appliedPressureForceNewtons);
        const StructureVector3 residual = subtract(
            load.resultingPendingForceNewtons, reconstructed);
        if (load.loadIndex != index || load.stableId == 0
            || load.structureNode >= application.structureNodeCount
            || !finite(load.priorPendingForceNewtons)
            || !finite(load.appliedPressureForceNewtons)
            || !finite(load.resultingPendingForceNewtons)
            || !finite(load.applicationResidualNewtons)
            || load.resultingPendingForceNewtons != reconstructed
            || load.applicationResidualNewtons != residual
            || (havePrevious
                && (load.stableId <= previousStableId
                    || load.structureNode <= previousStructureNode))) {
            throw std::invalid_argument(
                "scene fluid regional pressure applied node load is invalid");
        }
        havePrevious = true;
        previousStableId = load.stableId;
        previousStructureNode = load.structureNode;
        appliedForce = add(appliedForce, load.appliedPressureForceNewtons);
    }
    const StructureVector3 aggregateResidual = subtract(
        application.resultingPendingForceNewtons,
        add(application.priorPendingForceNewtons,
            application.appliedPressureForceNewtons));
    const double forceScale = std::max({
        maximumAbsolute(application.priorPendingForceNewtons),
        maximumAbsolute(application.appliedPressureForceNewtons),
        maximumAbsolute(application.resultingPendingForceNewtons), 1.0});
    if (appliedForce != application.appliedPressureForceNewtons
        || aggregateResidual != application.applicationResidualNewtons
        || maximumAbsolute(application.applicationResidualNewtons)
            > tolerance(forceScale)) {
        throw std::invalid_argument(
            "scene fluid regional pressure load-application closure is invalid");
    }
}

void validateSceneFluidRegionalPressureLoadApplication(
    const SceneFluidRegionalPressureLoadApplication& application,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidRegionalPressureSampleSet& samples,
    const ConservativeTransferSettings& transferSettings,
    const SceneFluidRegionalPressureLoadApplicationLimits& limits) {
    validateApplicationLimits(limits);
    validateSceneFluidRegionalPressureSampleIntegrity(samples);
    validateSceneFluidRegionalPressureLoadApplicationIntegrity(application);
    if (application.nodeLoads.size() > limits.maximumNodeLoads
        || application.structureNodeCount > limits.maximumStructureNodes
        || application.ownedStorageBytes > limits.maximumOwnedBytes
        || application.workingStorageBytes > limits.maximumWorkingBytes) {
        throw std::length_error(
            "scene fluid regional pressure load-application validation limit exceeded");
    }
    if (application.sourceSamplingFingerprint != samples.fingerprint
        || application.sourceSurfaceStateFingerprint != state.fingerprint
        || application.couplingSurfaceFingerprint
            != transfer.couplingSurfaceFingerprint()
        || application.targetDefinitionFingerprint
            != transfer.targetDefinitionFingerprint()
        || application.acceptedStepCount != state.acceptedStepCount
        || application.simulationTimeSeconds != state.simulationTimeSeconds
        || application.nodeLoads.size() != transfer.nodes().size()) {
        throw std::invalid_argument(
            "scene fluid regional pressure load application is foreign to its source");
    }
    const auto transferred =
        evaluateSceneFluidRegionalAcceptedPressureQuadrature(
            surface, state, transfer, quadrature, samples,
            transferSettings);
    const auto transferredLoads = transferred.nodeLoads();
    if (transferredLoads.size() != application.nodeLoads.size()) {
        throw std::invalid_argument(
            "scene fluid regional pressure load application changed node count");
    }
    for (std::size_t index = 0;
         index < transferredLoads.size(); ++index) {
        const auto& expected = transferredLoads[index];
        const auto& actual = application.nodeLoads[index];
        if (actual.loadIndex != index
            || actual.stableId != expected.stableId
            || actual.structureNode != expected.structureNode
            || actual.appliedPressureForceNewtons
                != expected.forceNewtons) {
            throw std::invalid_argument(
                "scene fluid regional pressure applied node load is foreign to its source");
        }
    }
}

} // namespace simwing::fsi
