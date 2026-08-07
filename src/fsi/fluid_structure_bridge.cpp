#include "fluid_structure_bridge.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace simwing::fsi {
namespace {

[[nodiscard]] StructureVector3 toStructure(const fluid::Vector3& value) {
    return {value.x, value.y, value.z};
}

[[nodiscard]] StructureVector3 subtract(const StructureVector3& first,
                                        const StructureVector3& second) {
    return {first.x - second.x,
            first.y - second.y,
            first.z - second.z};
}

[[nodiscard]] StructureVector3 scale(const StructureVector3& value,
                                     const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

[[nodiscard]] double length(const StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

[[nodiscard]] bool finite(const StructureVector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

[[nodiscard]] double combinedTolerance(const double absoluteTolerance,
                                       const double relativeTolerance,
                                       const double firstMagnitude,
                                       const double secondMagnitude) {
    return absoluteTolerance
        + relativeTolerance * std::max(firstMagnitude, secondMagnitude);
}

void validateSettings(const UniformFluidStructureBridgeSettings& settings) {
    const double values[] = {
        settings.maximumPressureTractionDeviationPascals,
        settings.absoluteAreaToleranceSquareMeters,
        settings.relativeAreaTolerance,
        settings.absoluteForceToleranceNewtons,
        settings.relativeForceTolerance,
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance,
    };
    if (!std::ranges::all_of(values, [](const double value) {
            return std::isfinite(value) && value >= 0.0;
        })) {
        throw std::invalid_argument(
            "uniform fluid-structure bridge tolerances must be finite and nonnegative");
    }
}

[[nodiscard]] const fluid::MovingInterfaceSurfaceDiagnostics& findSurface(
    const fluid::MovingInterfaceProjectionDiagnostics& diagnostics,
    const std::uint64_t stableId) {
    const auto found = std::lower_bound(
        diagnostics.surfaces.begin(), diagnostics.surfaces.end(), stableId,
        [](const fluid::MovingInterfaceSurfaceDiagnostics& surface,
           const std::uint64_t id) { return surface.stableId < id; });
    if (found == diagnostics.surfaces.end() || found->stableId != stableId) {
        throw std::invalid_argument(
            "fluid diagnostics do not contain the bound surface stable ID");
    }
    if (std::next(found) != diagnostics.surfaces.end()
        && std::next(found)->stableId == stableId) {
        throw std::invalid_argument(
            "fluid diagnostics contain a duplicate surface stable ID");
    }
    return *found;
}

} // namespace

UniformFluidStructureTransferResult::UniformFluidStructureTransferResult(
    ConservativeTransferResult transferResult,
    UniformFluidStructureBridgeDiagnostics diagnostics)
    : transferResult_(std::move(transferResult)),
      diagnostics_(std::move(diagnostics)) {}

const ConservativeTransferResult&
UniformFluidStructureTransferResult::transferResult() const noexcept {
    return transferResult_;
}

const UniformFluidStructureBridgeDiagnostics&
UniformFluidStructureTransferResult::diagnostics() const noexcept {
    return diagnostics_;
}

UniformFluidStructureBridge::UniformFluidStructureBridge(
    const Structure& target,
    const std::uint64_t fluidSurfaceStableId,
    std::vector<CouplingSurfaceNodeDefinition> nodes,
    std::vector<CouplingSurfaceTriangleDefinition> triangles)
    : fluidSurfaceStableId_(fluidSurfaceStableId),
      transfer_(target, std::move(nodes), std::move(triangles)) {
    if (fluidSurfaceStableId_ == 0) {
        throw std::invalid_argument(
            "uniform fluid-structure bridge surface stable ID must be nonzero");
    }
}

std::uint64_t UniformFluidStructureBridge::fluidSurfaceStableId() const
    noexcept {
    return fluidSurfaceStableId_;
}

const ConservativeSurfaceTransfer&
UniformFluidStructureBridge::transfer() const noexcept {
    return transfer_;
}

UniformFluidStructureTransferResult UniformFluidStructureBridge::evaluate(
    const fluid::MovingInterfaceProjectionDiagnostics& fluidDiagnostics,
    const std::span<const CouplingNodeKinematics> nodeKinematics,
    const UniformFluidStructureBridgeSettings& settings) const {
    validateSettings(settings);
    if (fluidDiagnostics.interfaceVersion
            != fluid::faceAlignedMovingInterfaceVersion
        || !fluidDiagnostics.projection.converged
        || !fluidDiagnostics.finite) {
        throw std::invalid_argument(
            "uniform fluid-structure bridge requires accepted finite fluid diagnostics");
    }
    if (!std::ranges::is_sorted(
            fluidDiagnostics.surfaces, {},
            &fluid::MovingInterfaceSurfaceDiagnostics::stableId)) {
        throw std::invalid_argument(
            "fluid surface diagnostics must be in canonical stable-ID order");
    }

    const auto& source = findSurface(
        fluidDiagnostics, fluidSurfaceStableId_);
    if (source.faceCount == 0
        || !std::isfinite(source.areaSquareMeters)
        || !(source.areaSquareMeters > 0.0)
        || !finite(toStructure(source.pressureForceNewtons))
        || !std::isfinite(source.pressurePowerWatts)
        || !std::isfinite(
            source.maximumPressureTractionDeviationPascals)) {
        throw std::invalid_argument(
            "bound fluid surface diagnostics are invalid");
    }
    if (source.maximumPressureTractionDeviationPascals
            > settings.maximumPressureTractionDeviationPascals) {
        throw std::invalid_argument(
            "bound fluid surface pressure traction is not uniform enough");
    }

    const StructureVector3 sourceForce =
        toStructure(source.pressureForceNewtons);
    const StructureVector3 uniformTraction = scale(
        sourceForce, 1.0 / source.areaSquareMeters);
    std::vector<CouplingTriangleTraction> tractions;
    tractions.reserve(transfer_.triangles().size());
    for (const auto& triangle : transfer_.triangles()) {
        tractions.push_back({triangle.stableId, uniformTraction});
    }
    ConservativeTransferResult transferred = transfer_.evaluate(
        nodeKinematics, tractions, settings.transfer);
    const auto& target = transferred.diagnostics();

    UniformFluidStructureBridgeDiagnostics diagnostics;
    diagnostics.fluidSurfaceStableId = fluidSurfaceStableId_;
    diagnostics.fluidFaceCount = source.faceCount;
    diagnostics.structureTriangleCount = target.triangleCount;
    diagnostics.fluidAreaSquareMeters = source.areaSquareMeters;
    diagnostics.structureAreaSquareMeters = target.surfaceAreaSquareMeters;
    diagnostics.areaResidualSquareMeters =
        target.surfaceAreaSquareMeters - source.areaSquareMeters;
    diagnostics.maximumPressureTractionDeviationPascals =
        source.maximumPressureTractionDeviationPascals;
    diagnostics.uniformPressureTractionPascals = uniformTraction;
    diagnostics.fluidPressureForceNewtons = sourceForce;
    diagnostics.structureSurfaceForceNewtons =
        target.integratedSurfaceForceNewtons;
    diagnostics.forceResidualNewtons = subtract(
        target.integratedSurfaceForceNewtons, sourceForce);
    diagnostics.forceResidualNormNewtons = length(
        diagnostics.forceResidualNewtons);
    diagnostics.fluidPressurePowerWatts = source.pressurePowerWatts;
    diagnostics.structureSurfacePowerWatts =
        target.integratedSurfacePowerWatts;
    diagnostics.powerResidualWatts =
        target.integratedSurfacePowerWatts - source.pressurePowerWatts;
    diagnostics.finite = std::isfinite(diagnostics.areaResidualSquareMeters)
        && finite(uniformTraction)
        && finite(diagnostics.forceResidualNewtons)
        && std::isfinite(diagnostics.forceResidualNormNewtons)
        && std::isfinite(diagnostics.powerResidualWatts)
        && target.finite;
    if (!diagnostics.finite) {
        throw std::invalid_argument(
            "uniform fluid-structure bridge produced non-finite diagnostics");
    }

    const double areaTolerance = combinedTolerance(
        settings.absoluteAreaToleranceSquareMeters,
        settings.relativeAreaTolerance,
        source.areaSquareMeters, target.surfaceAreaSquareMeters);
    const double forceTolerance = combinedTolerance(
        settings.absoluteForceToleranceNewtons,
        settings.relativeForceTolerance,
        length(sourceForce), length(target.integratedSurfaceForceNewtons));
    const double powerTolerance = combinedTolerance(
        settings.absolutePowerToleranceWatts,
        settings.relativePowerTolerance,
        std::abs(source.pressurePowerWatts),
        std::abs(target.integratedSurfacePowerWatts));
    if (std::abs(diagnostics.areaResidualSquareMeters) > areaTolerance) {
        throw std::invalid_argument(
            "fluid and structural coupling surface areas do not match");
    }
    if (diagnostics.forceResidualNormNewtons > forceTolerance) {
        throw std::invalid_argument(
            "fluid and structural pressure-force ledgers do not match");
    }
    if (std::abs(diagnostics.powerResidualWatts) > powerTolerance) {
        throw std::invalid_argument(
            "fluid and structural interface-power ledgers do not match");
    }
    return {std::move(transferred), std::move(diagnostics)};
}

} // namespace simwing::fsi
