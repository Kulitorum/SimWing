#include "scene_fluid_pressure_coupling.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace simwing::fsi {
namespace {

double vectorNorm(const StructureVector3& value) {
    return std::sqrt(value.x * value.x
                     + value.y * value.y
                     + value.z * value.z);
}

std::size_t couplingInterfaceBytes(const std::size_t nodeCount) {
    constexpr std::size_t bytesPerNode =
        3 * sizeof(double) + sizeof(CouplingNodeLoad);
    if (nodeCount > std::numeric_limits<std::size_t>::max()
                        / bytesPerNode) {
        throw std::length_error(
            "scene pressure coupling interface storage size overflows");
    }
    return nodeCount * bytesPerNode;
}

std::vector<double> flattenLoads(
    const ConservativeTransferResult& transfer) {
    if (transfer.nodeLoads().size()
        > std::numeric_limits<std::size_t>::max() / 3) {
        throw std::length_error(
            "scene pressure coupling interface size overflows");
    }
    std::vector<double> result;
    result.reserve(3 * transfer.nodeLoads().size());
    for (const auto& load : transfer.nodeLoads()) {
        result.push_back(load.forceNewtons.x);
        result.push_back(load.forceNewtons.y);
        result.push_back(load.forceNewtons.z);
    }
    return result;
}

std::vector<CouplingNodeLoad> expandLoads(
    const SceneFluidSurfaceTransfer& transfer,
    const std::span<const double> values) {
    const auto nodes = transfer.nodes();
    if (nodes.size() > std::numeric_limits<std::size_t>::max() / 3
        || values.size() != 3 * nodes.size()) {
        throw std::invalid_argument(
            "scene pressure coupling interface size is invalid");
    }
    std::vector<CouplingNodeLoad> result;
    result.reserve(nodes.size());
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        result.push_back({
            nodes[index].stableId,
            nodes[index].structureNode,
            {values[3 * index], values[3 * index + 1],
             values[3 * index + 2]},
        });
    }
    return result;
}

void replaceTractionResidual(
    CouplingResidualNorms& residuals,
    const std::span<const CouplingNodeLoad> candidate,
    const std::span<const CouplingNodeLoad> appliedGuess) {
    if (candidate.size() != appliedGuess.size()) {
        throw std::invalid_argument(
            "scene pressure coupling load residual size is invalid");
    }
    residuals.tractionNewtons = 0.0;
    residuals.tractionReferenceNewtons = 0.0;
    for (std::size_t index = 0; index < candidate.size(); ++index) {
        if (candidate[index].stableId != appliedGuess[index].stableId
            || candidate[index].structureNode
                != appliedGuess[index].structureNode) {
            throw std::invalid_argument(
                "scene pressure coupling load residual binding is invalid");
        }
        const StructureVector3 difference{
            candidate[index].forceNewtons.x
                - appliedGuess[index].forceNewtons.x,
            candidate[index].forceNewtons.y
                - appliedGuess[index].forceNewtons.y,
            candidate[index].forceNewtons.z
                - appliedGuess[index].forceNewtons.z,
        };
        residuals.tractionNewtons = std::max(
            residuals.tractionNewtons, vectorNorm(difference));
        residuals.tractionReferenceNewtons = std::max({
            residuals.tractionReferenceNewtons,
            vectorNorm(candidate[index].forceNewtons),
            vectorNorm(appliedGuess[index].forceNewtons),
        });
    }
}

std::vector<SceneFluidQuadraturePressure> zeroPressures(
    const SceneFluidQuadratureDefinition& quadrature) {
    std::vector<SceneFluidQuadraturePressure> result;
    result.reserve(quadrature.points.size());
    for (const auto& point : quadrature.points) {
        result.push_back({point.stableId, 0.0, 0.0});
    }
    return result;
}

std::vector<SceneFluidQuadratureTraction> combinedTractions(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidQuadratureDefinition& quadrature,
    const std::span<const SceneFluidQuadraturePressure> pressures,
    const std::span<const SceneFluidQuadratureTraction> wallTractions) {
    auto result = buildSceneFluidPressureTractions(
        surface, state, quadrature, pressures);
    if (wallTractions.empty()) {
        return result;
    }
    if (wallTractions.size() != result.size()) {
        throw std::invalid_argument(
            "scene pressure coupling wall traction size is invalid");
    }
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (wallTractions[index].stableId != result[index].stableId) {
            throw std::invalid_argument(
                "scene pressure coupling wall traction binding is invalid");
        }
        result[index].tractionPascals.x +=
            wallTractions[index].tractionPascals.x;
        result[index].tractionPascals.y +=
            wallTractions[index].tractionPascals.y;
        result[index].tractionPascals.z +=
            wallTractions[index].tractionPascals.z;
    }
    return result;
}

bool finite(const SceneFluidPressureCouplingStepDiagnostics& diagnostics) {
    return diagnostics.iteration.convergence.finite
        && diagnostics.iteration.relaxation.finite
        && diagnostics.structure.finite
        && diagnostics.pressureProjection.finite
        && (!diagnostics.usesRegionWall
            || (diagnostics.regionWall.finite
                && diagnostics.regionWall.accepted))
        && (!diagnostics.usesRegionRebase
            || diagnostics.regionRebase.finite)
        && (!diagnostics.usesMimeticPressureAudit
            || (diagnostics.mimeticPressureAudit.accepted
                && diagnostics.mimeticPressureAudit.pressureSolve
                    .reducedTraceSolve.finite
                && diagnostics.mimeticPressureComparison.finite
                && diagnostics.mimeticPressureOwnerTransition.fingerprint
                    != 0
                && diagnostics.mimeticPressureOwnerTransition
                       .comparisonFingerprint
                    == diagnostics.mimeticPressureComparisonFingerprint))
        && diagnostics.pressureTransfer.finite
        && diagnostics.totalFluidTransfer.finite
        && std::isfinite(diagnostics.interfaceForceClosureNewtons)
        && std::isfinite(diagnostics.interfaceForceReferenceNewtons);
}

bool sameControlTopology(
    const SceneFluidPressureControlVolumeSet& first,
    const SceneFluidPressureControlVolumeSet& second) {
    if (first.controlVolumes.size() != second.controlVolumes.size()) {
        return false;
    }
    for (std::size_t index = 0;
         index < first.controlVolumes.size(); ++index) {
        if (first.controlVolumes[index].stableId
            != second.controlVolumes[index].stableId) {
            return false;
        }
    }
    return true;
}

std::vector<std::uint64_t> controlStableIds(
    const SceneFluidPressureControlVolumeSet& controls) {
    std::vector<std::uint64_t> result;
    result.reserve(controls.controlVolumes.size());
    for (const auto& control : controls.controlVolumes) {
        result.push_back(control.stableId);
    }
    return result;
}

bool sameControlTopology(
    const std::span<const std::uint64_t> stableIds,
    const SceneFluidPressureControlVolumeSet& controls) {
    if (stableIds.size() != controls.controlVolumes.size()) {
        return false;
    }
    for (std::size_t index = 0; index < stableIds.size(); ++index) {
        if (stableIds[index] != controls.controlVolumes[index].stableId) {
            return false;
        }
    }
    return true;
}

bool sameVector(const StructureVector3& first,
                const StructureVector3& second) {
    return first.x == second.x
        && first.y == second.y
        && first.z == second.z;
}

bool sameSettings(const SceneFluidPressureCouplingSettings& first,
                  const SceneFluidPressureCouplingSettings& second) {
    return first.structure.timeStepSeconds
               == second.structure.timeStepSeconds
        && first.structure.substeps == second.structure.substeps
        && first.structure.constraintIterations
               == second.structure.constraintIterations
        && first.structure.cableConstraintSweepPairs
               == second.structure.cableConstraintSweepPairs
        && sameVector(
            first.structure.gravityMetersPerSecondSquared,
            second.structure.gravityMetersPerSecondSquared)
        && first.structure.velocityDampingPerSecond
               == second.structure.velocityDampingPerSecond
        && sameVector(
            first.structure.dampingReferenceVelocityMetersPerSecond,
            second.structure.dampingReferenceVelocityMetersPerSecond)
        && first.structure.workerThreads == second.structure.workerThreads
        && first.relaxation == second.relaxation
        && first.convergence == second.convergence
        && first.pressureEpoch == second.pressureEpoch
        && first.pressureProjection == second.pressureProjection
        && first.regionWall == second.regionWall
        && sameVector(
            first.transfer.momentReferenceMeters,
            second.transfer.momentReferenceMeters)
        && first.transfer.minimumTriangleAreaSquareMeters
               == second.transfer.minimumTriangleAreaSquareMeters
        && first.transfer.minimumQuadratureAreaSquareMeters
               == second.transfer.minimumQuadratureAreaSquareMeters
        && first.transfer.barycentricTolerance
               == second.transfer.barycentricTolerance;
}

} // namespace

SceneFluidPressureCouplingSettings::SceneFluidPressureCouplingSettings() {
    pressureProjection.timeStepSeconds = structure.timeStepSeconds;
    regionWall.timeStepSeconds = structure.timeStepSeconds;
}

SceneFluidPressureCoupling::SceneFluidPressureCoupling(
    SceneFluidSurfaceDefinition surface,
    SceneStructureMappings structureMappings,
    const Structure& target,
    fluid::PeriodicCartesianGrid grid,
    const SceneFluidPressureCouplingSettings& settings,
    const SceneFluidPressureCouplingLimits& limits,
    const SceneFluidMimeticPressureAuditConfiguration&
        mimeticPressureAudit)
    : surface_(std::move(surface)),
      structureMappings_(std::move(structureMappings)),
      grid_(std::move(grid)),
      connectivity_(buildSceneFluidRegionConnectivity(
          surface_, limits.connectivity)),
      transfer_(surface_, structureMappings_, target),
      macroCoupling_(transfer_.conservativeTransfer()),
      settings_(settings),
      limits_(limits),
      mimeticPressureAuditConfiguration_(mimeticPressureAudit),
      acceptedSurfaceState_(captureSceneFluidSurfaceState(
          surface_, structureMappings_, target)),
      acceptedPressureEpoch_(buildSceneFluidPressureEpoch(
          surface_, acceptedSurfaceState_, grid_, transfer_, connectivity_,
          settings_.pressureEpoch, limits_.pressureEpoch)),
      acceptedPressurePascals_(
          acceptedPressureEpoch_.pressureOperator.rows.size(), 0.0) {
    if (settings_.structure.timeStepSeconds
            != settings_.pressureProjection.timeStepSeconds) {
        throw std::invalid_argument(
            "scene pressure coupling Structure and pressure time steps differ");
    }
    if (mimeticPressureAuditConfiguration_.enabled
        && mimeticPressureAuditConfiguration_.settings
                .densityKgPerCubicMeter
            != settings_.pressureProjection.densityKgPerCubicMeter) {
        throw std::invalid_argument(
            "scene pressure coupling graph and mimetic densities differ");
    }
    if (mimeticPressureAuditConfiguration_.enabled
        && mimeticPressureAuditConfiguration_.settings.timeStepSeconds
            != settings_.structure.timeStepSeconds) {
        throw std::invalid_argument(
            "scene pressure coupling graph and mimetic time steps differ");
    }
    if (mimeticPressureAuditConfiguration_.enabled) {
        validateSceneFluidPressureOwnerTransitionPolicy(
            mimeticPressureAuditConfiguration_.ownerTransitionPolicy);
    }
    if (transfer_.nodes().size() > limits_.maximumCouplingNodes
        || couplingInterfaceBytes(transfer_.nodes().size())
            > limits_.maximumInterfaceBytes) {
        throw std::length_error(
            "scene pressure coupling interface exceeds its limits");
    }
    const auto pressures = zeroPressures(
        acceptedPressureEpoch_.gridEpoch.quadrature);
    acceptedPressureTransfer_.emplace(
        evaluateSceneFluidPressureQuadrature(
            surface_, acceptedSurfaceState_, transfer_,
            acceptedPressureEpoch_.gridEpoch.quadrature, pressures,
            settings_.transfer));
    acceptedPressureOnlyTransfer_ = acceptedPressureTransfer_;
}

const fluid::PeriodicCartesianGrid&
SceneFluidPressureCoupling::grid() const noexcept {
    return grid_;
}

const SceneFluidSurfaceState&
SceneFluidPressureCoupling::acceptedSurfaceState() const noexcept {
    return acceptedSurfaceState_;
}

const SceneFluidPressureEpoch&
SceneFluidPressureCoupling::acceptedPressureEpoch() const noexcept {
    return acceptedPressureEpoch_;
}

const ConservativeTransferResult&
SceneFluidPressureCoupling::acceptedPressureTransfer() const noexcept {
    return *acceptedPressureTransfer_;
}

const ConservativeTransferResult&
SceneFluidPressureCoupling::acceptedPressureOnlyTransfer() const noexcept {
    return *acceptedPressureOnlyTransfer_;
}

const SceneFluidPressureProjection*
SceneFluidPressureCoupling::acceptedPressureProjection() const noexcept {
    return acceptedPressureProjection_
        ? &*acceptedPressureProjection_ : nullptr;
}

const SceneFluidPressureSampleSet*
SceneFluidPressureCoupling::acceptedPressureSamples() const noexcept {
    return acceptedPressureSamples_
        ? &*acceptedPressureSamples_ : nullptr;
}

const SceneFluidAcceptedWallTractionSet*
SceneFluidPressureCoupling::acceptedWallTractions() const noexcept {
    return acceptedWallTractions_ ? &*acceptedWallTractions_ : nullptr;
}

const SceneFluidMimeticPressureAuditEndpoint*
SceneFluidPressureCoupling::acceptedMimeticPressureAudit() const noexcept {
    return acceptedMimeticPressureAudit_
        ? &*acceptedMimeticPressureAudit_ : nullptr;
}

const SceneFluidPressureShadowComparison*
SceneFluidPressureCoupling::acceptedMimeticPressureComparison()
    const noexcept {
    return acceptedMimeticPressureComparison_
        ? &*acceptedMimeticPressureComparison_ : nullptr;
}

const SceneFluidPressureOwnerTransitionDecision*
SceneFluidPressureCoupling::acceptedMimeticPressureOwnerTransition()
    const noexcept {
    return acceptedMimeticPressureOwnerTransition_
        ? &*acceptedMimeticPressureOwnerTransition_ : nullptr;
}

std::uint64_t
SceneFluidPressureCoupling::mimeticPressureAuditSettingsFingerprint()
    const {
    return mimeticPressureAuditConfiguration_.enabled
        ? sceneFluidMimeticPressureAuditSettingsFingerprint(
            mimeticPressureAuditConfiguration_.settings)
        : 0;
}

SceneFluidMimeticPressureAuditTopology
SceneFluidPressureCoupling::rebuildMimeticPressureAuditTopology(
    const Structure& target,
    const StructureCheckpoint& structureCheckpoint) const {
    if (!mimeticPressureAuditConfiguration_.enabled
        || target.definitionFingerprint()
            != transfer_.targetDefinitionFingerprint()) {
        throw std::invalid_argument(
            "scene pressure coupling mimetic restart topology is unavailable");
    }
    Structure restoredStructure(target.definition());
    restoredStructure.restore(structureCheckpoint);
    const auto restoredState = captureSceneFluidSurfaceState(
        surface_, structureMappings_, restoredStructure);
    const auto restoredEpoch = buildSceneFluidPressureEpoch(
        surface_, restoredState, grid_, transfer_, connectivity_,
        settings_.pressureEpoch, limits_.pressureEpoch);
    return buildSceneFluidMimeticPressureAuditTopology(
        surface_, restoredState, grid_, restoredEpoch,
        mimeticPressureAuditConfiguration_.settings,
        mimeticPressureAuditConfiguration_.limits);
}

SceneFluidPressureMacVelocityCollapse
SceneFluidPressureCoupling::acceptedPressureCorrectedMacVelocity() const {
    if (!acceptedPressureProjection_) {
        throw std::logic_error(
            "scene pressure coupling has no accepted pressure velocity");
    }
    const auto& projection = *acceptedPressureProjection_;
    validateSceneFluidPressureProjectionIntegrity(projection);
    const auto& faceLinks = acceptedPressureEpoch_.pressureFaceLinks;
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    const auto& patches = acceptedPressureEpoch_.openingPatches;
    if (!projection.diagnostics.accepted
        || projection.pressureFaceLinkFingerprint != faceLinks.fingerprint
        || projection.acceptedStepCount != faceLinks.acceptedStepCount
        || projection.simulationTimeSeconds
            != faceLinks.simulationTimeSeconds
        || projection.cellCounts != grid_.cellCounts()
        || projection.lowerMeters != grid_.lowerMeters()
        || projection.upperMeters != grid_.upperMeters()
        || projection.links.size() != faceLinks.links.size()) {
        throw std::invalid_argument(
            "scene pressure corrected MAC velocity identity is invalid");
    }
    std::map<std::uint64_t, const SceneFluidOpeningGridPatch*>
        openingPatches;
    for (const auto& patch : patches.patches) {
        if (patch.stableId == 0
            || !openingPatches.emplace(patch.stableId, &patch).second) {
            throw std::invalid_argument(
                "scene pressure corrected MAC opening identity is invalid");
        }
    }

    SceneFluidPressureMacVelocityCollapse result(grid_);
    auto& diagnostics = result.diagnostics;
    diagnostics.pressureProjectionFingerprint = projection.fingerprint;
    diagnostics.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    diagnostics.openingPatchFingerprint = patches.fingerprint;
    diagnostics.acceptedStepCount = projection.acceptedStepCount;
    diagnostics.simulationTimeSeconds = projection.simulationTimeSeconds;
    diagnostics.faceCount = faceLinks.faces.size();
    diagnostics.linkCount = faceLinks.links.size();

    const auto axisComponent = [](const Vec3& value,
                                  const fluid::GridFaceAxis axis) {
        switch (axis) {
        case fluid::GridFaceAxis::X:
            return value.x;
        case fluid::GridFaceAxis::Y:
            return value.y;
        case fluid::GridFaceAxis::Z:
            return value.z;
        }
        throw std::invalid_argument(
            "scene pressure corrected MAC face axis is invalid");
    };
    const auto absoluteFlow = [&](
        const SceneFluidPressureFaceLink& source,
        const SceneFluidPressureProjectedLink& projected,
        const fluid::GridFaceAxis axis) {
        double flow = projected
            .correctedRelativeVolumeFlowRateCubicMetersPerSecond;
        if (source.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            const auto found = openingPatches.find(
                source.openingPatchStableId);
            if (found == openingPatches.end()) {
                throw std::invalid_argument(
                    "scene pressure corrected MAC opening patch is missing");
            }
            const auto& patch = *found->second;
            const double normal = axisComponent(
                patch.unitNormalNegativeToPositive, axis);
            if (patch.openingId != source.openingId
                || patch.areaSquareMeters != source.areaSquareMeters
                || std::abs(std::abs(normal) - 1.0) > 1.0e-10) {
                throw std::invalid_argument(
                    "scene pressure corrected MAC opening patch is foreign");
            }
            flow += (normal > 0.0 ? 1.0 : -1.0)
                * patch.surfaceSweepRateCubicMetersPerSecond;
        }
        return flow;
    };
    const auto embeddedAbsoluteFlow = [&](
        const SceneFluidPressureFaceLink& source,
        const SceneFluidPressureProjectedLink& projected) {
        if (source.geometryKind
                != SceneFluidPressureLinkGeometryKind::EmbeddedOpening
            || source.kind
                != SceneFluidPressureFaceLinkKind::AuthoredOpening
            || source.faceIndex != invalidSceneFluidPressureFaceIndex) {
            throw std::invalid_argument(
                "scene pressure corrected MAC embedded link is invalid");
        }
        const auto found = openingPatches.find(
            source.openingPatchStableId);
        if (found == openingPatches.end()) {
            throw std::invalid_argument(
                "scene pressure corrected MAC embedded patch is missing");
        }
        const auto& patch = *found->second;
        const bool forwardRegions =
            patch.negativeSideRegionId == source.minusRegionId
            && patch.positiveSideRegionId == source.plusRegionId;
        const bool reverseRegions =
            patch.negativeSideRegionId == source.plusRegionId
            && patch.positiveSideRegionId == source.minusRegionId;
        const double alignment =
            patch.unitNormalNegativeToPositive.x
                * source.unitNormalMinusToPlus.x
            + patch.unitNormalNegativeToPositive.y
                * source.unitNormalMinusToPlus.y
            + patch.unitNormalNegativeToPositive.z
                * source.unitNormalMinusToPlus.z;
        const double orientation = forwardRegions ? 1.0 : -1.0;
        if (patch.openingId != source.openingId
            || patch.areaSquareMeters != source.areaSquareMeters
            || (!forwardRegions && !reverseRegions)
            || orientation * alignment < 1.0 - 1.0e-10) {
            throw std::invalid_argument(
                "scene pressure corrected MAC embedded patch is foreign");
        }
        return projected
                .correctedRelativeVolumeFlowRateCubicMetersPerSecond
            + orientation * patch.surfaceSweepRateCubicMetersPerSecond;
    };

    for (const auto& face : faceLinks.faces) {
        if (!(face.faceAreaSquareMeters > 0.0)
            || face.linkCount == 0
            || face.firstLink > faceLinks.links.size()
            || face.linkCount
                > faceLinks.links.size() - face.firstLink) {
            throw std::invalid_argument(
                "scene pressure corrected MAC face is unresolved");
        }
        if (face.linkCount > 1) {
            ++diagnostics.multiLinkFaceCount;
        }
        double totalFlow = 0.0;
        for (std::size_t offset = 0; offset < face.linkCount; ++offset) {
            const std::size_t index = face.firstLink + offset;
            const auto& source = faceLinks.links[index];
            const auto& projected = projection.links[index];
            if (source.linkIndex != index
                || source.faceIndex != face.faceIndex
                || projected.linkIndex != index
                || projected.stableId != source.stableId
                || projected.faceIndex != face.faceIndex
                || projected.kind != source.kind
                || projected.minusControlVolumeIndex
                    != source.minusControlVolumeIndex
                || projected.plusControlVolumeIndex
                    != source.plusControlVolumeIndex
                || projected.openingPatchStableId
                    != source.openingPatchStableId
                || !(source.areaSquareMeters > 0.0)) {
                throw std::invalid_argument(
                    "scene pressure corrected MAC link binding is invalid");
            }
            if (source.kind
                == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
                ++diagnostics.openingLinkCount;
            }
            totalFlow += absoluteFlow(source, projected, face.axis);
        }
        const double collapsedVelocity =
            totalFlow / face.faceAreaSquareMeters;
        if (!std::isfinite(collapsedVelocity)) {
            throw std::overflow_error(
                "scene pressure corrected MAC velocity is non-finite");
        }
        diagnostics.maximumAbsoluteVelocityMetersPerSecond = std::max(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond,
            std::abs(collapsedVelocity));
        diagnostics.maximumVolumeFlowClosureCubicMetersPerSecond = std::max(
            diagnostics.maximumVolumeFlowClosureCubicMetersPerSecond,
            std::abs(collapsedVelocity * face.faceAreaSquareMeters
                     - totalFlow));
        for (std::size_t offset = 0; offset < face.linkCount; ++offset) {
            const std::size_t index = face.firstLink + offset;
            const auto& source = faceLinks.links[index];
            const double subfaceVelocity = absoluteFlow(
                source, projection.links[index], face.axis)
                / source.areaSquareMeters;
            diagnostics.maximumSubfaceVelocityDeviationMetersPerSecond =
                std::max(
                    diagnostics
                        .maximumSubfaceVelocityDeviationMetersPerSecond,
                    std::abs(subfaceVelocity - collapsedVelocity));
        }

        const std::size_t cellIndex = grid_.cellIndex(
            face.i, face.j, face.k);
        switch (face.axis) {
        case fluid::GridFaceAxis::X:
            result.velocityMetersPerSecond.xFaces()[cellIndex] =
                collapsedVelocity;
            break;
        case fluid::GridFaceAxis::Y:
            result.velocityMetersPerSecond.yFaces()[cellIndex] =
                collapsedVelocity;
            break;
        case fluid::GridFaceAxis::Z:
            result.velocityMetersPerSecond.zFaces()[cellIndex] =
                collapsedVelocity;
            break;
        }
    }
    for (std::size_t index = 0; index < faceLinks.links.size(); ++index) {
        const auto& source = faceLinks.links[index];
        if (source.geometryKind
            != SceneFluidPressureLinkGeometryKind::EmbeddedOpening) {
            continue;
        }
        const auto& projected = projection.links[index];
        if (source.linkIndex != index
            || projected.linkIndex != index
            || projected.stableId != source.stableId
            || projected.faceIndex != invalidSceneFluidPressureFaceIndex
            || projected.kind != source.kind
            || projected.minusControlVolumeIndex
                != source.minusControlVolumeIndex
            || projected.plusControlVolumeIndex
                != source.plusControlVolumeIndex
            || projected.openingPatchStableId
                != source.openingPatchStableId
            || !(source.areaSquareMeters > 0.0)
            || !std::isfinite(embeddedAbsoluteFlow(source, projected))) {
            throw std::invalid_argument(
                "scene pressure corrected MAC embedded binding is invalid");
        }
        ++diagnostics.openingLinkCount;
        ++diagnostics.embeddedOpeningLinkCount;
    }
    if (diagnostics.embeddedOpeningLinkCount
        != faceLinks.embeddedOpeningLinkCount) {
        throw std::invalid_argument(
            "scene pressure corrected MAC embedded count is invalid");
    }
    diagnostics.finite = fluid::isFinite(result.velocityMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumAbsoluteVelocityMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumSubfaceVelocityDeviationMetersPerSecond)
        && std::isfinite(
            diagnostics.maximumVolumeFlowClosureCubicMetersPerSecond);
    if (!diagnostics.finite) {
        throw std::overflow_error(
            "scene pressure corrected MAC diagnostics are non-finite");
    }
    return result;
}

const SceneFluidPressureCouplingSettings&
SceneFluidPressureCoupling::settings() const noexcept {
    return settings_;
}

SceneFluidPressureCouplingCheckpoint
SceneFluidPressureCoupling::checkpoint(const Structure& target) const {
    const auto currentState = captureSceneFluidSurfaceState(
        surface_, structureMappings_, target);
    if (currentState != acceptedSurfaceState_) {
        throw std::invalid_argument(
            "scene pressure coupling checkpoint Structure state is foreign");
    }
    SceneFluidPressureCouplingCheckpoint result;
    result.surfaceDefinitionFingerprint = surface_.fingerprint;
    result.couplingSurfaceFingerprint =
        transfer_.couplingSurfaceFingerprint();
    result.structureDefinitionFingerprint =
        transfer_.targetDefinitionFingerprint();
    result.cellCounts = grid_.cellCounts();
    result.lowerMeters = grid_.lowerMeters();
    result.upperMeters = grid_.upperMeters();
    result.settings = settings_;
    result.structure = target.checkpoint();
    result.pressureProjection = acceptedPressureProjection_;
    result.wallTractions = acceptedWallTractions_;
    result.mimeticPressureAuditSettingsFingerprint =
        mimeticPressureAuditSettingsFingerprint();
    if (acceptedMimeticPressureAudit_) {
        result.mimeticPressureState = acceptedMimeticPressureAudit_
            ->pressureEpoch.acceptedPressureState;
    } else if (acceptedMimeticPressureAuditWarmState_) {
        result.mimeticPressureState =
            acceptedMimeticPressureAuditWarmState_
                ->acceptedPressureState;
    }
    if ((result.mimeticPressureAuditSettingsFingerprint != 0)
            != mimeticPressureAuditConfiguration_.enabled
        || (result.mimeticPressureState.has_value()
            != (mimeticPressureAuditConfiguration_.enabled
                && acceptedPressureProjection_.has_value()))) {
        throw std::logic_error(
            "scene pressure coupling mimetic checkpoint state is incomplete");
    }
    return result;
}

void SceneFluidPressureCoupling::restore(
    Structure& target,
    const SceneFluidPressureCouplingCheckpoint& checkpointValue) {
    if (checkpointValue.version
            != sceneFluidPressureCouplingCheckpointVersion
        || checkpointValue.surfaceDefinitionFingerprint
            != surface_.fingerprint
        || checkpointValue.couplingSurfaceFingerprint
            != transfer_.couplingSurfaceFingerprint()
        || checkpointValue.structureDefinitionFingerprint
            != transfer_.targetDefinitionFingerprint()
        || checkpointValue.cellCounts != grid_.cellCounts()
        || checkpointValue.lowerMeters != grid_.lowerMeters()
        || checkpointValue.upperMeters != grid_.upperMeters()
        || !sameSettings(checkpointValue.settings, settings_)
        || checkpointValue.mimeticPressureAuditSettingsFingerprint
            != mimeticPressureAuditSettingsFingerprint()
        || target.definitionFingerprint()
            != transfer_.targetDefinitionFingerprint()) {
        throw std::invalid_argument(
            "scene pressure coupling checkpoint identity is invalid");
    }

    Structure restoredStructure(target.definition());
    restoredStructure.restore(checkpointValue.structure);
    auto restoredState = captureSceneFluidSurfaceState(
        surface_, structureMappings_, restoredStructure);
    auto restoredEpoch = buildSceneFluidPressureEpoch(
        surface_, restoredState, grid_, transfer_, connectivity_,
        settings_.pressureEpoch, limits_.pressureEpoch);
    std::vector<double> restoredPressure(
        restoredEpoch.pressureOperator.rows.size(), 0.0);
    std::optional<ConservativeTransferResult> restoredTransfer;
    std::optional<ConservativeTransferResult> restoredPressureOnlyTransfer;
    std::optional<SceneFluidPressureProjection> restoredProjection;
    std::optional<SceneFluidPressureSampleSet> restoredSamples;
    std::optional<SceneFluidAcceptedWallTractionSet> restoredWallTractions;
    std::optional<SceneFluidMimeticPressureAuditWarmState>
        restoredMimeticPressureWarmState;
    if (checkpointValue.pressureProjection) {
        validateSceneFluidPressureProjectionIntegrity(
            *checkpointValue.pressureProjection);
        if (checkpointValue.pressureProjection->acceptedStepCount
                != restoredState.acceptedStepCount
            || checkpointValue.pressureProjection->simulationTimeSeconds
                != restoredState.simulationTimeSeconds
            || checkpointValue.pressureProjection
                   ->pressureControlVolumeFingerprint
                != restoredEpoch.pressureControlVolumes.fingerprint
            || checkpointValue.pressureProjection
                   ->pressureFaceLinkFingerprint
                != restoredEpoch.pressureFaceLinks.fingerprint
            || checkpointValue.pressureProjection
                   ->pressureOperatorFingerprint
                != restoredEpoch.pressureOperator.fingerprint) {
            throw std::invalid_argument(
                "scene pressure coupling checkpoint projection is foreign");
        }
        const auto wallFingerprint = checkpointValue.pressureProjection
            ->regionWallExchangeFingerprint;
        if ((wallFingerprint != 0)
                != checkpointValue.wallTractions.has_value()) {
            throw std::invalid_argument(
                "scene pressure coupling checkpoint wall state is incomplete");
        }
        if (checkpointValue.wallTractions) {
            validateSceneFluidAcceptedWallTractions(
                *checkpointValue.wallTractions,
                restoredEpoch.gridEpoch.quadrature, wallFingerprint);
            restoredWallTractions = checkpointValue.wallTractions;
        }
        auto samples = sampleSceneFluidProjectedPressure(
            restoredEpoch.gridEpoch.quadrature,
            restoredEpoch.pressureControlVolumes,
            *checkpointValue.pressureProjection,
            limits_.pressureSampling);
        const auto tractions = combinedTractions(
            surface_, restoredState, restoredEpoch.gridEpoch.quadrature,
            samples.pressures,
            restoredWallTractions
                ? std::span<const SceneFluidQuadratureTraction>(
                    restoredWallTractions->tractions)
                : std::span<const SceneFluidQuadratureTraction>{});
        restoredPressureOnlyTransfer.emplace(
            evaluateSceneFluidPressureQuadrature(
                surface_, restoredState, transfer_,
                restoredEpoch.gridEpoch.quadrature, samples.pressures,
                settings_.transfer));
        restoredTransfer.emplace(evaluateSceneFluidQuadrature(
            transfer_, restoredState, restoredEpoch.gridEpoch.quadrature,
            tractions, settings_.transfer));
        restoredSamples = std::move(samples);
        restoredPressure =
            checkpointValue.pressureProjection->pressurePascals;
        restoredProjection = checkpointValue.pressureProjection;
    } else {
        if (checkpointValue.wallTractions
            || restoredState.acceptedStepCount != 0
            || restoredState.simulationTimeSeconds != 0.0) {
            throw std::invalid_argument(
                "scene pressure coupling noninitial checkpoint lacks pressure");
        }
        const auto pressures = zeroPressures(
            restoredEpoch.gridEpoch.quadrature);
        restoredTransfer.emplace(evaluateSceneFluidPressureQuadrature(
            surface_, restoredState, transfer_,
            restoredEpoch.gridEpoch.quadrature, pressures,
            settings_.transfer));
        restoredPressureOnlyTransfer = restoredTransfer;
    }
    const bool expectsMimeticPressureState =
        mimeticPressureAuditConfiguration_.enabled
        && checkpointValue.pressureProjection.has_value();
    if (checkpointValue.mimeticPressureState.has_value()
            != expectsMimeticPressureState) {
        throw std::invalid_argument(
            "scene pressure coupling mimetic checkpoint state is incomplete");
    }
    if (checkpointValue.mimeticPressureState) {
        auto topology = buildSceneFluidMimeticPressureAuditTopology(
            surface_, restoredState, grid_, restoredEpoch,
            mimeticPressureAuditConfiguration_.settings,
            mimeticPressureAuditConfiguration_.limits);
        restoredMimeticPressureWarmState.emplace(
            bindSceneFluidMimeticPressureAuditWarmState(
                std::move(topology),
                *checkpointValue.mimeticPressureState,
                mimeticPressureAuditConfiguration_.limits));
    }

    target.restore(checkpointValue.structure);
    acceptedSurfaceState_ = std::move(restoredState);
    acceptedPressureEpoch_ = std::move(restoredEpoch);
    acceptedPressurePascals_ = std::move(restoredPressure);
    acceptedPressureTransfer_ = std::move(restoredTransfer);
    acceptedPressureOnlyTransfer_ =
        std::move(restoredPressureOnlyTransfer);
    acceptedPressureProjection_ = std::move(restoredProjection);
    acceptedPressureSamples_ = std::move(restoredSamples);
    acceptedWallTractions_ = std::move(restoredWallTractions);
    acceptedMimeticPressureAudit_.reset();
    acceptedMimeticPressureAuditWarmState_ =
        std::move(restoredMimeticPressureWarmState);
    acceptedMimeticPressureComparison_.reset();
    acceptedMimeticPressureOwnerTransition_.reset();
}

SceneFluidPressureCouplingStepDiagnostics
SceneFluidPressureCoupling::advance(
    Structure& target,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond) {
    return advanceImpl(target, predictedVelocityMetersPerSecond, nullptr);
}

SceneFluidPressureCouplingStepDiagnostics
SceneFluidPressureCoupling::advance(
    Structure& target,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidRegionTransport& transportedRegionMomentum) {
    return advanceImpl(
        target, predictedVelocityMetersPerSecond,
        &transportedRegionMomentum);
}

SceneFluidPressureCouplingStepDiagnostics
SceneFluidPressureCoupling::advanceImpl(
    Structure& target,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidRegionTransport* const transportedRegionMomentum) {
    const auto currentAcceptedState = captureSceneFluidSurfaceState(
        surface_, structureMappings_, target);
    if (currentAcceptedState != acceptedSurfaceState_
        || !predictedVelocityMetersPerSecond.matches(grid_)) {
        throw std::invalid_argument(
            "scene pressure coupling advance baseline is invalid");
    }
    if (transportedRegionMomentum != nullptr) {
        validateSceneFluidRegionTransportIntegrity(
            *transportedRegionMomentum);
        if (!acceptedPressureProjection_
            || !transportedRegionMomentum->diagnostics.accepted
            || transportedRegionMomentum->pressureProjectionFingerprint
                != acceptedPressureProjection_->fingerprint
            || transportedRegionMomentum->pressureFaceLinkFingerprint
                != acceptedPressureEpoch_.pressureFaceLinks.fingerprint
            || transportedRegionMomentum->acceptedStepCount
                != acceptedSurfaceState_.acceptedStepCount
            || transportedRegionMomentum->sourceSimulationTimeSeconds
                != acceptedSurfaceState_.simulationTimeSeconds
            || transportedRegionMomentum->targetSimulationTimeSeconds
                != acceptedSurfaceState_.simulationTimeSeconds
                    + settings_.structure.timeStepSeconds
            || transportedRegionMomentum->densityKgPerCubicMeter
                != settings_.pressureProjection.densityKgPerCubicMeter) {
            throw std::invalid_argument(
                "scene pressure coupling region transport baseline is foreign");
        }
    }
    const StructureCheckpoint structureBaseline = target.checkpoint();
    const auto baselineKinematics = transfer_.kinematics(
        acceptedSurfaceState_);
    auto previousKinematics = baselineKinematics;
    ConservativeTransferResult previousTraction =
        *acceptedPressureTransfer_;
    std::vector<double> warmPressure = acceptedPressurePascals_;
    std::vector<std::uint64_t> warmPressureStableIds = controlStableIds(
        acceptedPressureEpoch_.pressureControlVolumes);
    StrongCouplingIteration iteration(
        transfer_.couplingSurfaceFingerprint(),
        flattenLoads(*acceptedPressureTransfer_),
        settings_.relaxation, settings_.convergence);

    SceneFluidPressureCouplingStepDiagnostics diagnostics;
    diagnostics.previousPressureEpochFingerprint =
        acceptedPressureEpoch_.fingerprint;
    try {
        for (;;) {
            if (diagnostics.solverRunCount != 0) {
                target.restore(structureBaseline);
            }
            ++diagnostics.solverRunCount;
            const auto appliedEndLoads = expandLoads(
                transfer_, iteration.currentInterface());
            const auto structureDiagnostics =
                macroCoupling_.advanceStructureWithEndpointLoads(
                    target, *acceptedPressureTransfer_, appliedEndLoads,
                    settings_.structure);
            if (!structureDiagnostics.finite) {
                throw std::runtime_error(
                    "scene pressure coupling Structure step was not accepted");
            }

            auto currentState = captureSceneFluidSurfaceState(
                surface_, structureMappings_, target);
            auto currentEpoch = buildSceneFluidPressureEpoch(
                surface_, currentState, grid_, transfer_, connectivity_,
                settings_.pressureEpoch, limits_.pressureEpoch);
            const auto topologyTransition =
                buildSceneFluidPressureTopologyTransition(
                    acceptedPressureEpoch_.pressureControlVolumes,
                    acceptedPressureEpoch_.pressureFaceLinks,
                    currentEpoch.pressureControlVolumes,
                    currentEpoch.pressureFaceLinks,
                    limits_.topologyTransition);
            const auto rates = buildSceneFluidPressureVolumeRates(
                acceptedPressureEpoch_.cellVolumes,
                currentEpoch.cellVolumes,
                currentEpoch.pressureControlVolumes,
                topologyTransition,
                limits_.volumeRates);
            const auto openingFlux = evaluateSceneFluidOpeningFlux(
                surface_, currentState, currentEpoch.openingCaps,
                currentEpoch.openingQuadrature,
                currentEpoch.openingPatches, grid_,
                predictedVelocityMetersPerSecond, limits_.openingFlux);
            auto projectionSettings = settings_.pressureProjection;
            projectionSettings.timeStepSeconds = rates.durationSeconds;
            std::vector<double> rebasedWarmPressure;
            std::span<const double> projectionWarmPressure = warmPressure;
            if (!sameControlTopology(
                    warmPressureStableIds,
                    currentEpoch.pressureControlVolumes)) {
                rebasedWarmPressure = rebaseSceneFluidPressureWarmStart(
                    acceptedPressureEpoch_.pressureControlVolumes,
                    currentEpoch.pressureControlVolumes,
                    topologyTransition,
                    acceptedPressurePascals_);
                projectionWarmPressure = rebasedWarmPressure;
            }
            SceneFluidPressureProjection projection;
            std::optional<SceneFluidRegionWallExchange> wallExchange;
            std::optional<SceneFluidRegionRebase> regionRebase;
            if (transportedRegionMomentum != nullptr) {
                if (!sameControlTopology(
                        acceptedPressureEpoch_.pressureControlVolumes,
                        currentEpoch.pressureControlVolumes)) {
                    regionRebase.emplace(rebaseSceneFluidRegionTransport(
                        *transportedRegionMomentum,
                        acceptedPressureEpoch_.pressureControlVolumes,
                        currentEpoch.pressureControlVolumes,
                        topologyTransition,
                        limits_.regionRebase));
                }
                auto wallSettings = settings_.regionWall;
                wallSettings.timeStepSeconds = rates.durationSeconds;
                if (regionRebase) {
                    wallExchange.emplace(
                        exchangeSceneFluidRegionWallMomentum(
                            *regionRebase, grid_,
                            currentEpoch.pressureControlVolumes, surface_,
                            currentState,
                            currentEpoch.gridEpoch.quadrature,
                            wallSettings, limits_.regionWall));
                } else {
                    wallExchange.emplace(
                        exchangeSceneFluidRegionWallMomentum(
                            *transportedRegionMomentum, grid_,
                            currentEpoch.pressureControlVolumes, surface_,
                            currentState,
                            currentEpoch.gridEpoch.quadrature,
                            wallSettings, limits_.regionWall));
                }
                if (!wallExchange->diagnostics.accepted) {
                    throw std::runtime_error(
                        "scene pressure coupling wall exchange was not accepted");
                }
                const auto regionPrediction =
                    predictSceneFluidRegionLinkFlows(
                        *wallExchange, grid_,
                        currentEpoch.pressureControlVolumes,
                        currentEpoch.pressureFaceLinks, openingFlux,
                        limits_.regionLinkFlow);
                projection = projectSceneFluidPressureLinkFlows(
                    surface_, currentState, grid_, transfer_,
                    currentEpoch.gridEpoch, currentEpoch.openingCaps,
                    currentEpoch.openingQuadrature,
                    currentEpoch.openingPatches, openingFlux,
                    predictedVelocityMetersPerSecond, regionPrediction,
                    currentEpoch.cellVolumes, connectivity_,
                    currentEpoch.pressureControlVolumes,
                    currentEpoch.pressureFaceLinks,
                    currentEpoch.pressureOperator, rates,
                    projectionWarmPressure,
                    projectionSettings, limits_.pressureProjection);
            } else {
                projection = projectSceneFluidPressureLinkFlows(
                    surface_, currentState, grid_, transfer_,
                    currentEpoch.gridEpoch, currentEpoch.openingCaps,
                    currentEpoch.openingQuadrature,
                    currentEpoch.openingPatches, openingFlux,
                    predictedVelocityMetersPerSecond,
                    currentEpoch.cellVolumes, connectivity_,
                    currentEpoch.pressureControlVolumes,
                    currentEpoch.pressureFaceLinks,
                    currentEpoch.pressureOperator, rates,
                    projectionWarmPressure,
                    projectionSettings, limits_.pressureProjection);
            }
            if (!projection.diagnostics.accepted) {
                throw std::runtime_error(
                    "scene pressure coupling projection was not accepted");
            }
            auto samples = sampleSceneFluidProjectedPressure(
                currentEpoch.gridEpoch.quadrature,
                currentEpoch.pressureControlVolumes, projection,
                limits_.pressureSampling);
            std::optional<SceneFluidAcceptedWallTractionSet>
                acceptedWallCandidate;
            if (wallExchange) {
                acceptedWallCandidate.emplace(
                    captureSceneFluidAcceptedWallTractions(*wallExchange));
            }
            auto pressureOnlyTraction =
                evaluateSceneFluidPressureQuadrature(
                    surface_, currentState, transfer_,
                    currentEpoch.gridEpoch.quadrature, samples.pressures,
                    settings_.transfer);
            const auto tractions = combinedTractions(
                surface_, currentState, currentEpoch.gridEpoch.quadrature,
                samples.pressures,
                acceptedWallCandidate
                    ? std::span<const SceneFluidQuadratureTraction>(
                        acceptedWallCandidate->tractions)
                    : std::span<const SceneFluidQuadratureTraction>{});
            auto currentTraction = evaluateSceneFluidQuadrature(
                transfer_, currentState,
                currentEpoch.gridEpoch.quadrature, tractions,
                settings_.transfer);
            const auto currentKinematics = transfer_.kinematics(currentState);
            auto residuals = macroCoupling_.measureResiduals(
                baselineKinematics, previousKinematics,
                currentKinematics, previousTraction, currentTraction);
            replaceTractionResidual(
                residuals, currentTraction.nodeLoads(), appliedEndLoads);
            const auto candidate = flattenLoads(currentTraction);
            const auto iterationResult = iteration.advance(
                candidate, residuals);

            diagnostics.iteration = iterationResult;
            diagnostics.structure = structureDiagnostics;
            diagnostics.pressureProjection = projection.diagnostics;
            diagnostics.usesRegionWall = wallExchange.has_value();
            if (wallExchange) {
                diagnostics.regionWall = wallExchange->diagnostics;
            }
            diagnostics.usesRegionRebase = regionRebase.has_value();
            if (regionRebase) {
                diagnostics.regionRebase = regionRebase->diagnostics;
            }
            diagnostics.pressureTransfer =
                pressureOnlyTraction.diagnostics();
            diagnostics.totalFluidTransfer = currentTraction.diagnostics();
            diagnostics.interfaceForceClosureNewtons =
                residuals.tractionNewtons;
            diagnostics.interfaceForceReferenceNewtons =
                residuals.tractionReferenceNewtons;
            diagnostics.currentPressureEpochFingerprint =
                currentEpoch.fingerprint;
            std::optional<SceneFluidMimeticPressureAuditEndpoint>
                mimeticPressureAuditCandidate;
            std::optional<SceneFluidPressureShadowComparison>
                mimeticPressureComparisonCandidate;
            std::optional<SceneFluidPressureOwnerTransitionDecision>
                mimeticPressureOwnerTransitionCandidate;
            if (iterationResult.status
                    == StrongCouplingIterationStatus::Converged
                && mimeticPressureAuditConfiguration_.enabled) {
                const auto& auditSettings =
                    mimeticPressureAuditConfiguration_.settings;
                const auto& auditLimits =
                    mimeticPressureAuditConfiguration_.limits;
                if (acceptedMimeticPressureAudit_) {
                    if (wallExchange) {
                        mimeticPressureAuditCandidate.emplace(
                            buildSceneFluidMimeticPressureAuditEndpoint(
                                surface_, currentState, grid_, currentEpoch,
                                openingFlux, *wallExchange, rates,
                                topologyTransition,
                                *acceptedMimeticPressureAudit_,
                                auditSettings, auditLimits));
                    } else {
                        mimeticPressureAuditCandidate.emplace(
                            buildSceneFluidMimeticPressureAuditEndpoint(
                                surface_, currentState, grid_, currentEpoch,
                                openingFlux,
                                predictedVelocityMetersPerSecond, rates,
                                topologyTransition,
                                *acceptedMimeticPressureAudit_,
                                auditSettings, auditLimits));
                    }
                } else if (acceptedMimeticPressureAuditWarmState_) {
                    if (wallExchange) {
                        mimeticPressureAuditCandidate.emplace(
                            buildSceneFluidMimeticPressureAuditEndpoint(
                                surface_, currentState, grid_, currentEpoch,
                                openingFlux, *wallExchange, rates,
                                topologyTransition,
                                *acceptedMimeticPressureAuditWarmState_,
                                auditSettings, auditLimits));
                    } else {
                        mimeticPressureAuditCandidate.emplace(
                            buildSceneFluidMimeticPressureAuditEndpoint(
                                surface_, currentState, grid_, currentEpoch,
                                openingFlux,
                                predictedVelocityMetersPerSecond, rates,
                                topologyTransition,
                                *acceptedMimeticPressureAuditWarmState_,
                                auditSettings, auditLimits));
                    }
                } else if (wallExchange) {
                    mimeticPressureAuditCandidate.emplace(
                        buildSceneFluidMimeticPressureAuditEndpoint(
                            surface_, currentState, grid_, currentEpoch,
                            openingFlux, *wallExchange, rates,
                            topologyTransition, auditSettings,
                            auditLimits));
                } else {
                    mimeticPressureAuditCandidate.emplace(
                        buildSceneFluidMimeticPressureAuditEndpoint(
                            surface_, currentState, grid_, currentEpoch,
                            openingFlux,
                            predictedVelocityMetersPerSecond, rates,
                            topologyTransition, auditSettings,
                            auditLimits));
                }
                diagnostics.usesMimeticPressureAudit = true;
                diagnostics.mimeticPressureAuditFingerprint =
                    mimeticPressureAuditCandidate->fingerprint;
                diagnostics.mimeticPressureAudit =
                    mimeticPressureAuditCandidate->pressureEpoch.diagnostics;
                mimeticPressureComparisonCandidate.emplace(
                    compareSceneFluidPressureShadow(
                        surface_, currentState, transfer_,
                        currentEpoch.gridEpoch.quadrature, projection,
                        samples,
                        mimeticPressureAuditCandidate->controlCells,
                        mimeticPressureAuditCandidate->pressureSources,
                        mimeticPressureAuditCandidate->pressureEpoch
                            .acceptedPressureSamples,
                        settings_.transfer,
                        mimeticPressureAuditConfiguration_
                            .comparisonLimits));
                diagnostics.mimeticPressureComparisonFingerprint =
                    mimeticPressureComparisonCandidate->fingerprint;
                diagnostics.mimeticPressureComparison =
                    mimeticPressureComparisonCandidate->diagnostics;
                mimeticPressureOwnerTransitionCandidate.emplace(
                    decideSceneFluidPressureOwnerTransition(
                        *mimeticPressureComparisonCandidate,
                        mimeticPressureAuditConfiguration_
                            .ownerTransitionPolicy));
                diagnostics.mimeticPressureOwnerTransition =
                    *mimeticPressureOwnerTransitionCandidate;
            }
            diagnostics.finite = finite(diagnostics);
            if (!diagnostics.finite) {
                throw std::runtime_error(
                    "scene pressure coupling diagnostics are non-finite");
            }

            if (iterationResult.status
                == StrongCouplingIterationStatus::Converged) {
                diagnostics.accepted = true;
                acceptedSurfaceState_ = std::move(currentState);
                acceptedPressureEpoch_ = std::move(currentEpoch);
                acceptedPressurePascals_ = projection.pressurePascals;
                acceptedPressureTransfer_ = std::move(currentTraction);
                acceptedPressureOnlyTransfer_ =
                    std::move(pressureOnlyTraction);
                acceptedPressureProjection_ = std::move(projection);
                acceptedPressureSamples_ = std::move(samples);
                acceptedWallTractions_ = std::move(acceptedWallCandidate);
                if (mimeticPressureAuditCandidate) {
                    acceptedMimeticPressureAudit_ =
                        std::move(mimeticPressureAuditCandidate);
                    acceptedMimeticPressureAuditWarmState_.reset();
                    acceptedMimeticPressureComparison_ =
                        std::move(mimeticPressureComparisonCandidate);
                    acceptedMimeticPressureOwnerTransition_ =
                        std::move(mimeticPressureOwnerTransitionCandidate);
                }
                return diagnostics;
            }
            if (iterationResult.status
                == StrongCouplingIterationStatus::Exhausted) {
                target.restore(structureBaseline);
                return diagnostics;
            }

            previousKinematics = currentKinematics;
            previousTraction = std::move(currentTraction);
            warmPressure = std::move(projection.pressurePascals);
            warmPressureStableIds = controlStableIds(
                currentEpoch.pressureControlVolumes);
        }
    } catch (...) {
        target.restore(structureBaseline);
        throw;
    }
}

} // namespace simwing::fsi
