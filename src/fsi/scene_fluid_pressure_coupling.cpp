#include "scene_fluid_pressure_coupling.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

bool finite(const SceneFluidPressureCouplingStepDiagnostics& diagnostics) {
    return diagnostics.iteration.convergence.finite
        && diagnostics.iteration.relaxation.finite
        && diagnostics.structure.finite
        && diagnostics.pressureProjection.finite
        && diagnostics.pressureTransfer.finite
        && std::isfinite(diagnostics.interfaceForceClosureNewtons)
        && std::isfinite(diagnostics.interfaceForceReferenceNewtons);
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
}

SceneFluidPressureCoupling::SceneFluidPressureCoupling(
    SceneFluidSurfaceDefinition surface,
    SceneStructureMappings structureMappings,
    const Structure& target,
    fluid::PeriodicCartesianGrid grid,
    const SceneFluidPressureCouplingSettings& settings,
    const SceneFluidPressureCouplingLimits& limits)
    : surface_(std::move(surface)),
      structureMappings_(std::move(structureMappings)),
      grid_(std::move(grid)),
      connectivity_(buildSceneFluidRegionConnectivity(
          surface_, limits.connectivity)),
      transfer_(surface_, structureMappings_, target),
      macroCoupling_(transfer_.conservativeTransfer()),
      settings_(settings),
      limits_(limits),
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

const SceneFluidPressureProjection*
SceneFluidPressureCoupling::acceptedPressureProjection() const noexcept {
    return acceptedPressureProjection_
        ? &*acceptedPressureProjection_ : nullptr;
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
    std::optional<SceneFluidPressureProjection> restoredProjection;
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
        auto samples = sampleSceneFluidProjectedPressure(
            restoredEpoch.gridEpoch.quadrature,
            restoredEpoch.pressureControlVolumes,
            *checkpointValue.pressureProjection,
            limits_.pressureSampling);
        restoredTransfer.emplace(
            evaluateSceneFluidProjectedPressureQuadrature(
                surface_, restoredState, transfer_,
                restoredEpoch.gridEpoch.quadrature, samples,
                settings_.transfer));
        restoredPressure =
            checkpointValue.pressureProjection->pressurePascals;
        restoredProjection = checkpointValue.pressureProjection;
    } else {
        if (restoredState.acceptedStepCount != 0
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
    }

    target.restore(checkpointValue.structure);
    acceptedSurfaceState_ = std::move(restoredState);
    acceptedPressureEpoch_ = std::move(restoredEpoch);
    acceptedPressurePascals_ = std::move(restoredPressure);
    acceptedPressureTransfer_ = std::move(restoredTransfer);
    acceptedPressureProjection_ = std::move(restoredProjection);
}

SceneFluidPressureCouplingStepDiagnostics
SceneFluidPressureCoupling::advance(
    Structure& target,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond) {
    const auto currentAcceptedState = captureSceneFluidSurfaceState(
        surface_, structureMappings_, target);
    if (currentAcceptedState != acceptedSurfaceState_
        || !predictedVelocityMetersPerSecond.matches(grid_)) {
        throw std::invalid_argument(
            "scene pressure coupling advance baseline is invalid");
    }
    const StructureCheckpoint structureBaseline = target.checkpoint();
    const auto baselineKinematics = transfer_.kinematics(
        acceptedSurfaceState_);
    auto previousKinematics = baselineKinematics;
    ConservativeTransferResult previousTraction =
        *acceptedPressureTransfer_;
    std::vector<double> warmPressure = acceptedPressurePascals_;
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
            const auto rates = buildSceneFluidPressureVolumeRates(
                acceptedPressureEpoch_.cellVolumes,
                currentEpoch.cellVolumes,
                currentEpoch.pressureControlVolumes,
                limits_.volumeRates);
            const auto openingFlux = evaluateSceneFluidOpeningFlux(
                surface_, currentState, currentEpoch.openingCaps,
                currentEpoch.openingQuadrature,
                currentEpoch.openingPatches, grid_,
                predictedVelocityMetersPerSecond, limits_.openingFlux);
            auto projection = projectSceneFluidPressureLinkFlows(
                surface_, currentState, grid_, transfer_,
                currentEpoch.gridEpoch, currentEpoch.openingCaps,
                currentEpoch.openingQuadrature,
                currentEpoch.openingPatches, openingFlux,
                predictedVelocityMetersPerSecond, currentEpoch.cellVolumes,
                connectivity_, currentEpoch.pressureControlVolumes,
                currentEpoch.pressureFaceLinks,
                currentEpoch.pressureOperator, rates, warmPressure,
                settings_.pressureProjection,
                limits_.pressureProjection);
            if (!projection.diagnostics.accepted) {
                throw std::runtime_error(
                    "scene pressure coupling projection was not accepted");
            }
            auto samples = sampleSceneFluidProjectedPressure(
                currentEpoch.gridEpoch.quadrature,
                currentEpoch.pressureControlVolumes, projection,
                limits_.pressureSampling);
            auto currentTraction =
                evaluateSceneFluidProjectedPressureQuadrature(
                    surface_, currentState, transfer_,
                    currentEpoch.gridEpoch.quadrature, samples,
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
            diagnostics.pressureTransfer = currentTraction.diagnostics();
            diagnostics.interfaceForceClosureNewtons =
                residuals.tractionNewtons;
            diagnostics.interfaceForceReferenceNewtons =
                residuals.tractionReferenceNewtons;
            diagnostics.currentPressureEpochFingerprint =
                currentEpoch.fingerprint;
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
                acceptedPressureProjection_ = std::move(projection);
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
        }
    } catch (...) {
        target.restore(structureBaseline);
        throw;
    }
}

} // namespace simwing::fsi
