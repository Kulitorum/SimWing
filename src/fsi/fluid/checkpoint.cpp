#include "fluid/checkpoint.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi::fluid {
namespace {

constexpr std::uint64_t fnvOffset = 14695981039346656037ull;
constexpr std::uint64_t fnvPrime = 1099511628211ull;

void hashValue(std::uint64_t& hash, const std::uint64_t value) noexcept {
    for (unsigned int byte = 0; byte < 8; ++byte) {
        hash ^= (value >> (byte * 8)) & 0xffu;
        hash *= fnvPrime;
    }
}

void hashValue(std::uint64_t& hash, const double value) noexcept {
    hashValue(hash, std::bit_cast<std::uint64_t>(value));
}

void hashVector(std::uint64_t& hash, const Vector3& value) noexcept {
    hashValue(hash, value.x);
    hashValue(hash, value.y);
    hashValue(hash, value.z);
}

bool finite(const Vector3& value) noexcept {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Vector3 add(const Vector3& first, const Vector3& second) noexcept {
    return {first.x + second.x,
            first.y + second.y,
            first.z + second.z};
}

Vector3 scale(const Vector3& value, const double factor) noexcept {
    return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(const Vector3& first, const Vector3& second) noexcept {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

double length(const Vector3& value) noexcept {
    return std::hypot(value.x, value.y, value.z);
}

bool near(const double first, const double second) noexcept {
    return std::abs(first - second)
        <= 1.0e-10
            + 1.0e-11 * std::max(std::abs(first), std::abs(second));
}

bool near(const Vector3& first, const Vector3& second) noexcept {
    return near(length(add(first, scale(second, -1.0))), 0.0);
}

Vector3 axialVelocity(const GridFaceAxis axis,
                      const double value) noexcept {
    switch (axis) {
    case GridFaceAxis::X:
        return {value, 0.0, 0.0};
    case GridFaceAxis::Y:
        return {0.0, value, 0.0};
    case GridFaceAxis::Z:
        return {0.0, 0.0, value};
    }
    return {};
}

double constrainedVelocity(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocity,
    const GridFaceMovingInterface& face) {
    const std::size_t index = grid.cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case GridFaceAxis::X:
        return velocity.xFaces()[index];
    case GridFaceAxis::Y:
        return velocity.yFaces()[index];
    case GridFaceAxis::Z:
        return velocity.zFaces()[index];
    }
    return std::numeric_limits<double>::quiet_NaN();
}

bool gridMetadataMatches(const PeriodicCartesianGrid& grid,
                         const GridCellCounts cellCounts,
                         const Vector3& lowerMeters,
                         const Vector3& upperMeters) noexcept {
    return grid.cellCounts() == cellCounts
        && grid.lowerMeters() == lowerMeters
        && grid.upperMeters() == upperMeters;
}

std::uint64_t topologyFingerprint(
    const PeriodicCartesianGrid& grid,
    const FaceAlignedMovingInterface& interfaces) noexcept {
    std::uint64_t hash = fnvOffset;
    const GridCellCounts counts = grid.cellCounts();
    hashValue(hash, static_cast<std::uint64_t>(
        faceAlignedMovingInterfaceVersion));
    hashValue(hash, static_cast<std::uint64_t>(counts.x));
    hashValue(hash, static_cast<std::uint64_t>(counts.y));
    hashValue(hash, static_cast<std::uint64_t>(counts.z));
    hashVector(hash, grid.lowerMeters());
    hashVector(hash, grid.upperMeters());
    hashValue(hash, static_cast<std::uint64_t>(interfaces.faceCount()));
    for (const auto& face : interfaces.faces()) {
        hashValue(hash, face.surfaceStableId);
        hashValue(hash, face.minusRegionStableId);
        hashValue(hash, face.plusRegionStableId);
        hashValue(hash, static_cast<std::uint64_t>(face.axis));
        hashValue(hash, static_cast<std::uint64_t>(face.i));
        hashValue(hash, static_cast<std::uint64_t>(face.j));
        hashValue(hash, static_cast<std::uint64_t>(face.k));
    }
    for (const std::uint64_t stableId : interfaces.regionStableIds()) {
        hashValue(hash, stableId);
    }
    for (const std::uint64_t stableId : interfaces.cellRegionStableIds()) {
        hashValue(hash, stableId);
    }
    return hash;
}

bool finiteProjection(const ProjectionDiagnostics& diagnostics) noexcept {
    const std::array values{
        diagnostics.compatibilityDivergencePerSecond,
        diagnostics.initialResidualPascalsPerSquareMeter,
        diagnostics.finalResidualPascalsPerSquareMeter,
        diagnostics.divergenceL2BeforePerSecond,
        diagnostics.divergenceL2AfterPerSecond,
        diagnostics.divergenceMaximumBeforePerSecond,
        diagnostics.divergenceMaximumAfterPerSecond,
        diagnostics.kineticEnergyBeforeJoules,
        diagnostics.kineticEnergyAfterJoules,
        diagnostics.pressureMeanPascals,
        diagnostics.pressureJumpSourceCompatibilityPascalsPerSquareMeter,
    };
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

bool diagnosticsMatch(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocityMetersPerSecond,
    const FaceAlignedMovingInterface& interfaces,
    const MovingInterfaceProjectionDiagnostics& diagnostics) noexcept {
    if (!diagnostics.projection.converged
        || !diagnostics.finite
        || diagnostics.interfaceVersion != interfaces.version()
        || diagnostics.interfaceFaceCount != interfaces.faceCount()
        || diagnostics.fluidRegionCount != interfaces.regionCount()
        || diagnostics.faces.size() != interfaces.faceCount()
        || diagnostics.regions.size() != interfaces.regionCount()
        || !finiteProjection(diagnostics.projection)
        || !std::isfinite(
            diagnostics.maximumAbsoluteRegionVolumeRateCubicMetersPerSecond)
        || !std::isfinite(
            diagnostics.maximumNormalVelocityErrorMetersPerSecond)
        || !finite(diagnostics.totalPressureForceNewtons)
        || !finite(diagnostics.totalPressureImpulseNewtonSeconds)
        || !std::isfinite(diagnostics.totalPressurePowerWatts)
        || !std::isfinite(diagnostics.totalPressureWorkJoules)
        || !finite(diagnostics.totalConstraintReactionForceNewtons)
        || !finite(diagnostics.totalConstraintReactionImpulseNewtonSeconds)
        || !std::isfinite(diagnostics.totalConstraintReactionPowerWatts)
        || !std::isfinite(diagnostics.totalConstraintReactionWorkJoules)) {
        return false;
    }

    const auto regionIds = interfaces.regionStableIds();
    for (std::size_t index = 0; index < regionIds.size(); ++index) {
        const auto& region = diagnostics.regions[index];
        const std::size_t expectedCellCount = static_cast<std::size_t>(
            std::ranges::count(
                interfaces.cellRegionStableIds(), regionIds[index]));
        if (region.stableId != regionIds[index]
            || region.cellCount == 0
            || region.cellCount != expectedCellCount
            || !std::isfinite(
                region.compatibilityVolumeRateCubicMetersPerSecond)
            || !std::isfinite(region.pressureMeanBeforePascals)
            || !std::isfinite(region.pressureMeanAfterPascals)) {
            return false;
        }
    }

    struct SurfaceAggregate {
        std::size_t faceCount = 0;
        double areaSquareMeters = 0.0;
        Vector3 pressureForceNewtons;
        Vector3 directConstraintForceNewtons;
        Vector3 constraintReactionForceNewtons;
        double pressurePowerWatts = 0.0;
        double constraintReactionPowerWatts = 0.0;
    };
    std::map<std::uint64_t, SurfaceAggregate> aggregates;
    for (std::size_t index = 0; index < diagnostics.faces.size(); ++index) {
        const auto& face = diagnostics.faces[index];
        const auto& source = interfaces.faces()[index];
        if (face.surfaceStableId != source.surfaceStableId
            || face.minusRegionStableId != source.minusRegionStableId
            || face.plusRegionStableId != source.plusRegionStableId
            || face.axis != source.axis
            || face.i != source.i || face.j != source.j || face.k != source.k
            || face.normalVelocityMetersPerSecond
                != source.normalVelocityMetersPerSecond
            || constrainedVelocity(
                grid, velocityMetersPerSecond, source)
                != source.normalVelocityMetersPerSecond
            || !finite(face.lowerCornerMeters)
            || !finite(face.upperCornerMeters)
            || !std::isfinite(face.areaSquareMeters)
            || face.areaSquareMeters <= 0.0
            || !std::isfinite(face.normalVelocityMetersPerSecond)
            || !finite(face.pressureTractionPascals)
            || !finite(face.pressureForceNewtons)
            || !std::isfinite(face.pressurePowerWatts)
            || !finite(face.directConstraintForceNewtons)
            || !finite(face.constraintReactionTractionPascals)
            || !finite(face.constraintReactionForceNewtons)
            || !std::isfinite(face.constraintReactionPowerWatts)
            || !near(
                face.pressureForceNewtons,
                scale(face.pressureTractionPascals,
                      face.areaSquareMeters))
            || !near(
                face.constraintReactionForceNewtons,
                add(face.pressureForceNewtons,
                    face.directConstraintForceNewtons))
            || !near(
                face.constraintReactionForceNewtons,
                scale(face.constraintReactionTractionPascals,
                      face.areaSquareMeters))
            || !near(
                face.pressurePowerWatts,
                dot(face.pressureForceNewtons,
                    axialVelocity(face.axis,
                                  face.normalVelocityMetersPerSecond)))
            || !near(
                face.constraintReactionPowerWatts,
                dot(face.constraintReactionForceNewtons,
                    axialVelocity(face.axis,
                                  face.normalVelocityMetersPerSecond)))) {
            return false;
        }
        auto& aggregate = aggregates[face.surfaceStableId];
        ++aggregate.faceCount;
        aggregate.areaSquareMeters += face.areaSquareMeters;
        aggregate.pressureForceNewtons = add(
            aggregate.pressureForceNewtons, face.pressureForceNewtons);
        aggregate.directConstraintForceNewtons = add(
            aggregate.directConstraintForceNewtons,
            face.directConstraintForceNewtons);
        aggregate.constraintReactionForceNewtons = add(
            aggregate.constraintReactionForceNewtons,
            face.constraintReactionForceNewtons);
        aggregate.pressurePowerWatts += face.pressurePowerWatts;
        aggregate.constraintReactionPowerWatts +=
            face.constraintReactionPowerWatts;
    }

    std::vector<std::uint64_t> expectedSurfaceIds;
    expectedSurfaceIds.reserve(interfaces.faceCount());
    for (const auto& face : interfaces.faces()) {
        expectedSurfaceIds.push_back(face.surfaceStableId);
    }
    std::ranges::sort(expectedSurfaceIds);
    const auto uniqueEnd = std::ranges::unique(expectedSurfaceIds).begin();
    expectedSurfaceIds.erase(uniqueEnd, expectedSurfaceIds.end());
    if (diagnostics.surfaces.size() != expectedSurfaceIds.size()) {
        return false;
    }
    Vector3 totalPressureForce;
    Vector3 totalPressureImpulse;
    Vector3 totalConstraintReactionForce;
    Vector3 totalConstraintReactionImpulse;
    double totalPressurePower = 0.0;
    double totalPressureWork = 0.0;
    double totalConstraintReactionPower = 0.0;
    double totalConstraintReactionWork = 0.0;
    for (std::size_t index = 0; index < diagnostics.surfaces.size(); ++index) {
        const auto& surface = diagnostics.surfaces[index];
        if (surface.stableId != expectedSurfaceIds[index]) {
            return false;
        }
        const auto foundAggregate = aggregates.find(surface.stableId);
        if (foundAggregate == aggregates.end()) {
            return false;
        }
        const auto& aggregate = foundAggregate->second;
        if (surface.faceCount == 0
            || surface.faceCount != aggregate.faceCount
            || !std::isfinite(surface.areaSquareMeters)
            || surface.areaSquareMeters <= 0.0
            || !finite(surface.pressureForceNewtons)
            || !finite(surface.pressureImpulseNewtonSeconds)
            || !std::isfinite(
                surface.maximumPressureTractionDeviationPascals)
            || !std::isfinite(surface.pressurePowerWatts)
            || !std::isfinite(surface.pressureWorkJoules)
            || !finite(surface.directConstraintForceNewtons)
            || !finite(surface.constraintReactionForceNewtons)
            || !finite(surface.constraintReactionImpulseNewtonSeconds)
            || !std::isfinite(
                surface.maximumConstraintReactionTractionDeviationPascals)
            || !std::isfinite(surface.constraintReactionPowerWatts)
            || !std::isfinite(surface.constraintReactionWorkJoules)
            || !near(surface.areaSquareMeters,
                     aggregate.areaSquareMeters)
            || !near(surface.pressureForceNewtons,
                     aggregate.pressureForceNewtons)
            || !near(surface.directConstraintForceNewtons,
                     aggregate.directConstraintForceNewtons)
            || !near(surface.constraintReactionForceNewtons,
                     aggregate.constraintReactionForceNewtons)
            || !near(surface.pressurePowerWatts,
                     aggregate.pressurePowerWatts)
            || !near(surface.constraintReactionPowerWatts,
                     aggregate.constraintReactionPowerWatts)) {
            return false;
        }
        totalPressureForce = add(
            totalPressureForce, surface.pressureForceNewtons);
        totalPressureImpulse = add(
            totalPressureImpulse, surface.pressureImpulseNewtonSeconds);
        totalConstraintReactionForce = add(
            totalConstraintReactionForce,
            surface.constraintReactionForceNewtons);
        totalConstraintReactionImpulse = add(
            totalConstraintReactionImpulse,
            surface.constraintReactionImpulseNewtonSeconds);
        totalPressurePower += surface.pressurePowerWatts;
        totalPressureWork += surface.pressureWorkJoules;
        totalConstraintReactionPower +=
            surface.constraintReactionPowerWatts;
        totalConstraintReactionWork +=
            surface.constraintReactionWorkJoules;
    }
    return near(diagnostics.totalPressureForceNewtons, totalPressureForce)
        && near(diagnostics.totalPressureImpulseNewtonSeconds,
                totalPressureImpulse)
        && near(diagnostics.totalPressurePowerWatts, totalPressurePower)
        && near(diagnostics.totalPressureWorkJoules, totalPressureWork)
        && near(diagnostics.totalConstraintReactionForceNewtons,
                totalConstraintReactionForce)
        && near(diagnostics.totalConstraintReactionImpulseNewtonSeconds,
                totalConstraintReactionImpulse)
        && near(diagnostics.totalConstraintReactionPowerWatts,
                totalConstraintReactionPower)
        && near(diagnostics.totalConstraintReactionWorkJoules,
                totalConstraintReactionWork);
}

void validateAcceptedState(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocityMetersPerSecond,
    const CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const MovingInterfaceProjectionDiagnostics& diagnostics) {
    if (!velocityMetersPerSecond.matches(grid)
        || !pressurePascals.matches(grid)
        || !interfaces.matches(grid)) {
        throw std::invalid_argument(
            "fluid checkpoint state does not match its grid");
    }
    if (!isFinite(velocityMetersPerSecond)
        || !isFinite(pressurePascals)
        || !diagnosticsMatch(
            grid, velocityMetersPerSecond, interfaces, diagnostics)) {
        throw std::invalid_argument(
            "fluid checkpoint requires accepted finite projection state");
    }
}

} // namespace

struct MovingInterfaceFluidCheckpoint::Detail {
    GridCellCounts cellCounts;
    Vector3 lowerMeters;
    Vector3 upperMeters;
    MacVelocityField velocityMetersPerSecond;
    CellScalarField pressurePascals;
    FaceAlignedMovingInterface interfaces;
    MovingInterfaceProjectionDiagnostics diagnostics;
    std::uint64_t topologyFingerprint = 0;
};

MovingInterfaceFluidCheckpoint checkpointMovingInterfaceFluidState(
    const PeriodicCartesianGrid& grid,
    const MacVelocityField& velocityMetersPerSecond,
    const CellScalarField& pressurePascals,
    const FaceAlignedMovingInterface& interfaces,
    const MovingInterfaceProjectionDiagnostics& diagnostics) {
    validateAcceptedState(
        grid, velocityMetersPerSecond, pressurePascals,
        interfaces, diagnostics);

    MovingInterfaceFluidCheckpoint result;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.scalarSampleCount = grid.cellCount();
    result.topologyFingerprint = topologyFingerprint(grid, interfaces);
    auto detail = std::make_shared<MovingInterfaceFluidCheckpoint::Detail>(
        MovingInterfaceFluidCheckpoint::Detail{
            result.cellCounts,
            result.lowerMeters,
            result.upperMeters,
            velocityMetersPerSecond,
            pressurePascals,
            interfaces,
            diagnostics,
            result.topologyFingerprint,
        });
    result.detail = std::move(detail);
    return result;
}

MovingInterfaceFluidState restoreMovingInterfaceFluidState(
    const PeriodicCartesianGrid& grid,
    const MovingInterfaceFluidCheckpoint& checkpoint) {
    if (checkpoint.version != movingInterfaceFluidCheckpointVersion
        || !checkpoint.detail
        || checkpoint.scalarSampleCount != grid.cellCount()
        || !gridMetadataMatches(
            grid, checkpoint.cellCounts,
            checkpoint.lowerMeters, checkpoint.upperMeters)
        || checkpoint.cellCounts != checkpoint.detail->cellCounts
        || checkpoint.lowerMeters != checkpoint.detail->lowerMeters
        || checkpoint.upperMeters != checkpoint.detail->upperMeters
        || checkpoint.topologyFingerprint
            != checkpoint.detail->topologyFingerprint
        || checkpoint.topologyFingerprint
            != topologyFingerprint(grid, checkpoint.detail->interfaces)) {
        throw std::invalid_argument(
            "fluid checkpoint version, grid, or topology binding is invalid");
    }
    validateAcceptedState(
        grid,
        checkpoint.detail->velocityMetersPerSecond,
        checkpoint.detail->pressurePascals,
        checkpoint.detail->interfaces,
        checkpoint.detail->diagnostics);
    return {
        checkpoint.detail->velocityMetersPerSecond,
        checkpoint.detail->pressurePascals,
        checkpoint.detail->interfaces,
        checkpoint.detail->diagnostics,
    };
}

} // namespace simwing::fsi::fluid
