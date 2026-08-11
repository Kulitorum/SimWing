#include "frozen_scene_pressure_case.h"

#include "fluid/evolution.h"
#include "fluid_frame.h"
#include "scene_fluid_capped_face_partition.h"
#include "scene_fluid_cell_volume.h"
#include "scene_fluid_mimetic_geometry_epoch.h"
#include "scene_fluid_mimetic_geometry_epoch_transition.h"
#include "scene_fluid_mimetic_pressure_audit.h"
#include "scene_fluid_mimetic_pressure_flow.h"
#include "scene_fluid_mimetic_pressure_sampling.h"
#include "scene_fluid_opening_cap.h"
#include "scene_fluid_opening_face_crossing.h"
#include "scene_fluid_opening_flux.h"
#include "scene_fluid_opening_patch.h"
#include "scene_fluid_opening_quadrature.h"
#include "scene_fluid_pressure_control_volume.h"
#include "scene_fluid_pressure_face_link.h"
#include "scene_fluid_pressure_traction.h"
#include "scene_fluid_pressure_volume_rate.h"
#include "scene_fluid_region_connectivity.h"
#include "scene_fluid_region_momentum.h"
#include "scene_fluid_region_rebase.h"
#include "scene_fluid_region_transport.h"
#include "scene_fluid_region_wall.h"
#include "scene_fluid_surface.h"
#include "scene_fluid_surface_transfer.h"
#include "scene_structure.h"
#include "structure_frame.h"
#include "structure_rest_audit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace simwing::fsi {
namespace {

void validateSettings(const FrozenScenePressureCaseSettings& settings) {
    const auto counts = settings.cellCounts;
    if (counts.x < 2 || counts.y < 2 || counts.z < 2
        || !std::isfinite(settings.domainPaddingMeters)
        || !(settings.domainPaddingMeters > 0.0)
        || (settings.useExplicitDomain
            && (!std::isfinite(settings.lowerMeters.x)
                || !std::isfinite(settings.lowerMeters.y)
                || !std::isfinite(settings.lowerMeters.z)
                || !std::isfinite(settings.upperMeters.x)
                || !std::isfinite(settings.upperMeters.y)
                || !std::isfinite(settings.upperMeters.z)
                || !(settings.upperMeters.x > settings.lowerMeters.x)
                || !(settings.upperMeters.y > settings.lowerMeters.y)
                || !(settings.upperMeters.z > settings.lowerMeters.z)))
        || !std::isfinite(settings.backgroundWindMetersPerSecond.x)
        || !std::isfinite(settings.backgroundWindMetersPerSecond.y)
        || !std::isfinite(settings.backgroundWindMetersPerSecond.z)
        || !std::isfinite(settings.windRampSeconds)
        || settings.windRampSeconds < 0.0
        || settings.windRampSeconds > 60.0
        || (settings.useCorrectedTraceFlowContinuation
            && settings.useRegionalTransportFlowPrediction)
        || (settings.useMovingGeometryFsi
            && (settings.useCorrectedTraceFlowContinuation
                || settings.useRegionalTransportFlowPrediction))
        || !std::isfinite(
            settings.diagnosticPerturbationSpeedMetersPerSecond)
        || settings.diagnosticPerturbationSpeedMetersPerSecond < 0.0
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || !std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)) {
        throw std::invalid_argument(
            "frozen scene pressure settings are invalid");
    }
}

std::string assemblyError(const SceneStructureAssembly& assembly) {
    if (assembly.diagnostics.empty()) {
        return "scene has no supported structural assembly";
    }
    return assembly.diagnostics.front().message;
}

std::string surfaceError(const SceneFluidSurfaceAssembly& surface) {
    if (surface.diagnostics.empty()) {
        return "scene has no supported fluid surface";
    }
    return surface.diagnostics.front().message;
}

fluid::PeriodicCartesianGrid makeGrid(
    const SceneFluidSurfaceState& state,
    const FrozenScenePressureCaseSettings& settings) {
    if (state.vertices.empty()) {
        throw std::invalid_argument(
            "frozen scene pressure surface is empty");
    }
    if (settings.useExplicitDomain) {
        return {
            settings.cellCounts, settings.lowerMeters,
            settings.upperMeters,
        };
    }
    fluid::Vector3 lower{
        state.vertices.front().positionMeters.x,
        state.vertices.front().positionMeters.y,
        state.vertices.front().positionMeters.z,
    };
    fluid::Vector3 upper = lower;
    for (const auto& vertex : state.vertices) {
        lower.x = std::min(lower.x, vertex.positionMeters.x);
        lower.y = std::min(lower.y, vertex.positionMeters.y);
        lower.z = std::min(lower.z, vertex.positionMeters.z);
        upper.x = std::max(upper.x, vertex.positionMeters.x);
        upper.y = std::max(upper.y, vertex.positionMeters.y);
        upper.z = std::max(upper.z, vertex.positionMeters.z);
    }
    lower.x -= settings.domainPaddingMeters;
    lower.y -= settings.domainPaddingMeters;
    lower.z -= settings.domainPaddingMeters;
    upper.x += settings.domainPaddingMeters;
    upper.y += settings.domainPaddingMeters;
    upper.z += settings.domainPaddingMeters;
    return {settings.cellCounts, lower, upper};
}

fluid::MacVelocityField makePredictor(
    const fluid::PeriodicCartesianGrid& grid,
    const FrozenScenePressureCaseSettings& settings,
    const double rampFraction) {
    fluid::MacVelocityField result(grid);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const std::size_t index = grid.cellIndex(i, j, k);
                const double sample = static_cast<double>(index + 1);
                result.xFaces()[index] =
                    rampFraction
                        * settings.backgroundWindMetersPerSecond.x
                    + settings.diagnosticPerturbationSpeedMetersPerSecond
                        * 0.01 * sample;
                result.yFaces()[index] =
                    rampFraction
                        * settings.backgroundWindMetersPerSecond.y
                    + settings.diagnosticPerturbationSpeedMetersPerSecond
                        * 0.02 * sample;
                result.zFaces()[index] =
                    rampFraction
                        * settings.backgroundWindMetersPerSecond.z
                    - settings.diagnosticPerturbationSpeedMetersPerSecond
                        * 0.015 * sample;
            }
        }
    }
    return result;
}

double windRampFraction(
    const FrozenScenePressureCaseSettings& settings,
    const std::size_t sampleNumber) {
    if (settings.windRampSeconds == 0.0) {
        return 1.0;
    }
    return std::min(
        static_cast<double>(sampleNumber) * settings.timeStepSeconds
            / settings.windRampSeconds,
        1.0);
}

double meanXVelocity(const fluid::MacVelocityField& velocity) {
    double sum = 0.0;
    for (const double value : velocity.xFaces()) {
        sum += value;
    }
    return sum / static_cast<double>(velocity.xFaces().size());
}

struct BuiltCase {
    Scene scene;
    SceneFluidSurfaceAssembly surface;
    SceneStructureAssembly assembly;
    Structure structure;
    SceneFluidSurfaceState surfaceState;
    SceneFluidSurfaceTransfer transfer;
    fluid::PeriodicCartesianGrid grid;
    SceneFluidRegionConnectivity connectivity;
    SceneFluidMimeticGeometryEpoch geometryEpoch;
    SceneFluidMimeticPressureAuditSettings pressureSettings;
    SceneFluidMimeticPressureAuditEndpoint pressure;
    SceneFluidMimeticCorrectedTraceFlow correctedFlow;
    SceneFluidMimeticMacVelocityCollapse correctedMac;
    SceneFluidRegionMomentumState regionalMomentum;
    ConservativeTransferResult pressureTransfer;
    ConservativeTransferResult totalFluidTransfer;
    viewer::StructureFrameMapping frameMapping;
    std::optional<SceneFluidMimeticTraceFlowContinuation>
        traceFlowContinuation;
    std::optional<SceneFluidMimeticRegionTransportFlowPrediction>
        regionalTransportFlowPrediction;
    std::optional<SceneFluidRegionTransport> regionalTransport;
    std::optional<SceneFluidRegionWallExchange> regionWall;
    std::optional<ConservativeTransferResult> wallTransfer;
    std::optional<StructureDiagnostics> structureStep;
    std::optional<fluid::PeriodicFlowStrangSubcyclingDiagnostics> bulkFlow;
    std::optional<fluid::ProjectionDiagnostics> bulkProjection;
    double meanStreamwiseVelocityBeforePumpMetersPerSecond = 0.0;
    double streamwisePumpIncrementMetersPerSecond = 0.0;
    double windRampFraction = 0.0;
    double maximumGeometryDisplacementMeters = 0.0;
    bool pressureControlTopologyStable = false;
    std::size_t flowAdvanceCount = 0;
    std::size_t geometryAdvanceCount = 0;
};

struct FluidLoadTransfers {
    ConservativeTransferResult pressure;
    std::optional<ConservativeTransferResult> wall;
    ConservativeTransferResult total;
};

StructureStepSettings makeStructureStepSettings(
    const StructureDefinition& definition,
    const FrozenScenePressureCaseSettings& settings) {
    StructureStepSettings result;
    result.timeStepSeconds = settings.timeStepSeconds;
    result.gravityMetersPerSecondSquared = {};
    if (definition.suspension) {
        result.constraintIterations =
            definition.suspension->solverIterations;
    }
    if (settings.useMovingGeometryFsi) {
        // The authoritative tessellation contains much smaller membrane
        // charts than the analytic harnesses. Resolve their structural time
        // scale explicitly instead of altering authored mass or stiffness at
        // the scene boundary.
        result.substeps = 8;
    }
    return result;
}

FluidLoadTransfers evaluateFluidLoads(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidQuadratureDefinition& quadrature,
    const SceneFluidMimeticPressureSampleSet& pressureSamples,
    const SceneFluidAcceptedWallTractionSet* wallTractions = nullptr) {
    auto pressure = evaluateSceneFluidMimeticPressureQuadrature(
        surface, state, transfer, quadrature, pressureSamples);
    if (wallTractions == nullptr) {
        return {pressure, std::nullopt, std::move(pressure)};
    }
    validateSceneFluidAcceptedWallTractions(
        *wallTractions, quadrature,
        wallTractions->wallExchangeFingerprint);
    auto pressureTractions = buildSceneFluidPressureTractions(
        surface, state, quadrature, pressureSamples.pressures);
    if (pressureTractions.size() != wallTractions->tractions.size()) {
        throw std::invalid_argument(
            "moving scene wall traction size is invalid");
    }
    for (std::size_t index = 0; index < pressureTractions.size(); ++index) {
        const auto& wall = wallTractions->tractions[index];
        auto& combined = pressureTractions[index];
        if (combined.stableId != wall.stableId) {
            throw std::invalid_argument(
                "moving scene wall traction binding is invalid");
        }
        combined.tractionPascals.x += wall.tractionPascals.x;
        combined.tractionPascals.y += wall.tractionPascals.y;
        combined.tractionPascals.z += wall.tractionPascals.z;
    }
    auto wall = evaluateSceneFluidQuadrature(
        transfer, state, quadrature, wallTractions->tractions);
    auto total = evaluateSceneFluidQuadrature(
        transfer, state, quadrature, pressureTractions);
    return {std::move(pressure), std::move(wall), std::move(total)};
}

BuiltCase buildCase(
    Scene scene,
    const FrozenScenePressureCaseSettings& settings) {
    validateSettings(settings);
    const auto validation = validateScene(scene);
    if (!validation.ok()) {
        throw std::invalid_argument(
            "frozen scene pressure input is not a valid scene-v2 payload");
    }
    auto surface = assembleSceneFluidSurface(scene);
    if (!surface.ok()) {
        throw std::invalid_argument(
            "frozen scene pressure surface failed: "
            + surfaceError(surface));
    }
    SceneStructureSettings structureAssemblySettings;
    if (settings.useMovingGeometryFsi) {
        structureAssemblySettings.suspensionSolverIterations = 16;
    }
    auto assembly = assembleSceneStructure(
        scene, {}, structureAssemblySettings);
    if (!assembly.ok()) {
        throw std::invalid_argument(
            "frozen scene pressure structure failed: "
            + assemblyError(assembly));
    }
    Structure structure(assembly.definition);
    if (settings.useMovingGeometryFsi) {
        const auto restStepSettings = makeStructureStepSettings(
            assembly.definition, settings);
        const auto restAudit = auditStructureRestState(
            structure, restStepSettings);
        if (!restAudit.stationary) {
            std::ostringstream message;
            message.precision(17);
            message
                << "moving scene initial Structure rest audit rejected:"
                << " node=" << restAudit.maximumDisplacementNode;
            if (restAudit.maximumDisplacementNode
                < assembly.mappings.nodeVertexIds.size()) {
                message << " vertex-id="
                        << assembly.mappings.nodeVertexIds[
                               restAudit.maximumDisplacementNode];
            }
            message
                << " dx=" << restAudit.maximumNodeDisplacementMeters
                << " m rms-dx=" << restAudit.rmsNodeDisplacementMeters
                << " m node-speed="
                << restAudit.maximumNodeSpeedMetersPerSecond
                << " m/s payload-dx="
                << restAudit.payloadDisplacementMeters
                << " m payload-rotation="
                << restAudit.payloadRotationRadians
                << " rad harness-dx="
                << restAudit.maximumHarnessDisplacementMeters
                << " m initial-membrane-edge-mismatch="
                << restAudit.maximumInitialMembraneEdgeMismatchMeters
                << " m membrane-index="
                << restAudit.maximumInitialMembraneEdgeMismatchIndex;
            if (restAudit.maximumInitialMembraneEdgeMismatchIndex
                < assembly.mappings.membraneTriangleIds.size()) {
                const auto triangleId =
                    assembly.mappings.membraneTriangleIds[
                        restAudit.maximumInitialMembraneEdgeMismatchIndex];
                message << " triangle-id="
                        << triangleId;
                const auto triangle = std::ranges::find(
                    scene.triangles, triangleId, &Triangle::id);
                if (triangle != scene.triangles.end()) {
                    message
                        << " sheet-id=" << triangle->sheetId
                        << " role=" << static_cast<unsigned>(triangle->role)
                        << " vertices=" << triangle->vertexIds[0] << ','
                        << triangle->vertexIds[1] << ','
                        << triangle->vertexIds[2];
                }
            }
            message
                << " initial-line-extension="
                << restAudit.maximumInitialSuspensionExtensionMeters
                << " m line-id="
                << restAudit
                       .maximumInitialSuspensionExtensionSegmentStableId
                << " line-residual="
                << restAudit.maximumSuspensionResidualMeters << " m";
            throw std::runtime_error(message.str());
        }
    }
    auto state = captureSceneFluidSurfaceState(
        surface.definition, assembly.mappings, structure);
    SceneFluidSurfaceTransfer transfer(
        surface.definition, assembly.mappings, structure);
    auto grid = makeGrid(state, settings);
    auto connectivity = buildSceneFluidRegionConnectivity(
        surface.definition);
    auto geometryEpoch = buildSceneFluidMimeticGeometryEpoch(
        surface.definition, state, grid, transfer, connectivity);
    const double initialRampFraction = windRampFraction(settings, 1);
    const auto predictor = makePredictor(
        grid, settings, initialRampFraction);
    const auto openingFlux = evaluateSceneFluidOpeningFlux(
        surface.definition, state, geometryEpoch.openingCaps,
        geometryEpoch.openingQuadrature, geometryEpoch.openingPatches,
        grid, predictor);
    SceneFluidMimeticPressureAuditSettings pressureSettings;
    pressureSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    pressureSettings.timeStepSeconds = settings.timeStepSeconds;
    pressureSettings.pressureSolve
        .absoluteResidualTolerancePascalsMeters = 1.0e-10;
    pressureSettings.pressureSolve.relativeResidualTolerance = 1.0e-5;
    pressureSettings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters = 1.0e-8;
    pressureSettings.pressureSolve.maximumIterations = 4000;
    auto pressure = buildSceneFluidMimeticPressureAuditEndpoint(
        surface.definition, state, grid, geometryEpoch.gridEpoch,
        geometryEpoch.openingCaps, geometryEpoch.openingQuadrature,
        geometryEpoch.openingPatches,
        geometryEpoch.pressureControlVolumes,
        geometryEpoch.pressureFaceLinks, openingFlux, predictor,
        pressureSettings);
    if (!pressure.pressureEpoch.diagnostics.accepted) {
        throw std::runtime_error(
            "frozen scene mimetic pressure solve was not accepted");
    }
    auto correctedFlow = correctSceneFluidMimeticTraceFlows(pressure);
    if (!correctedFlow.accepted) {
        throw std::runtime_error(
            "frozen scene mimetic pressure correction was not accepted");
    }
    auto correctedMac = collapseSceneFluidMimeticCorrectedMacVelocity(
        correctedFlow, geometryEpoch.pressureFaceLinks,
        geometryEpoch.openingPatches, grid);
    auto regionalMomentum = reconstructSceneFluidRegionMomentumState(
        grid, geometryEpoch.pressureControlVolumes,
        geometryEpoch.pressureFaceLinks, geometryEpoch.openingPatches,
        pressure.controlCells, correctedFlow,
        correctedMac.velocityMetersPerSecond);
    auto fluidLoads = evaluateFluidLoads(
        surface.definition, state, transfer,
        geometryEpoch.gridEpoch.quadrature,
        pressure.pressureEpoch.acceptedPressureSamples);
    auto frameMapping = viewer::makeStructureFrameMapping(
        scene, assembly, structure);
    BuiltCase result{
        std::move(scene), std::move(surface), std::move(assembly),
        std::move(structure), std::move(state), std::move(transfer),
        std::move(grid), std::move(connectivity),
        std::move(geometryEpoch),
        std::move(pressureSettings), std::move(pressure),
        std::move(correctedFlow), std::move(correctedMac),
        std::move(regionalMomentum),
        std::move(fluidLoads.pressure),
        std::move(fluidLoads.total), std::move(frameMapping),
    };
    result.windRampFraction = initialRampFraction;
    return result;
}

void advanceFixedGeometryFlow(
    BuiltCase& built,
    const FrozenScenePressureCaseSettings& settings) {
    const auto& geometry = built.geometryEpoch;
    auto candidateVelocity = built.correctedMac.velocityMetersPerSecond;
    std::optional<SceneFluidMimeticTraceFlowPrediction> previousBulkBaseline;
    if (settings.useCorrectedTraceFlowContinuation) {
        const auto previousOpeningFlux = evaluateSceneFluidOpeningFlux(
            built.surface.definition, built.surfaceState,
            geometry.openingCaps, geometry.openingQuadrature,
            geometry.openingPatches, built.grid,
            built.correctedMac.velocityMetersPerSecond);
        previousBulkBaseline = sampleSceneFluidMimeticTraceFlows(
            built.pressure.controlCells,
            built.pressure.fullTraceSystem,
            geometry.pressureFaceLinks, previousOpeningFlux, built.grid,
            built.correctedMac.velocityMetersPerSecond);
    }
    const double meanBeforePump = meanXVelocity(candidateVelocity);
    const double candidateRampFraction = windRampFraction(
        settings, built.flowAdvanceCount + 2);
    const double pumpIncrement =
        candidateRampFraction
            * settings.backgroundWindMetersPerSecond.x
        - meanBeforePump;
    for (double& value : candidateVelocity.xFaces()) {
        value += pumpIncrement;
    }
    fluid::CellScalarField candidateBulkPressure(built.grid);
    fluid::ProjectionSettings projectionSettings;
    projectionSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    projectionSettings.timeStepSeconds = settings.timeStepSeconds;
    const auto candidateBulkProjection = fluid::projectVelocity(
        built.grid, candidateVelocity, candidateBulkPressure,
        projectionSettings);
    if (!candidateBulkProjection.converged) {
        throw std::runtime_error(
            "frozen scene bulk continuation projection did not converge");
    }
    fluid::PeriodicFlowStrangSubcyclingSettings bulkSettings;
    bulkSettings.flow.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    bulkSettings.flow.timeStepSeconds = settings.timeStepSeconds;
    bulkSettings.flow.advectionReconstruction =
        fluid::VariableMacReconstruction::DonorCell;
    bulkSettings.maximumSubsteps = 1024;
    auto candidateBulkFlow =
        fluid::advancePeriodicFlowStrangSspRk2Subcycled(
            built.grid, candidateVelocity, candidateBulkPressure,
            bulkSettings);
    if (!candidateBulkFlow.accepted) {
        throw std::runtime_error(
            "frozen scene bulk flow rejected its requested interval at stage "
            + std::to_string(static_cast<unsigned>(
                candidateBulkFlow.failureStage))
            + " (initial divergence "
            + std::to_string(candidateBulkFlow.initialDivergenceL2PerSecond)
            + ", final divergence "
            + std::to_string(candidateBulkFlow.finalDivergenceL2PerSecond)
            + ")");
    }
    SceneFluidRegionTransportSettings regionalTransportSettings;
    regionalTransportSettings.timeStepSeconds = settings.timeStepSeconds;
    regionalTransportSettings.maximumSubsteps = 1024;
    auto candidateRegionalTransport = advanceSceneFluidRegionMomentum(
        built.regionalMomentum, geometry.pressureFaceLinks,
        built.pressure.controlCells, built.correctedFlow, built.grid,
        built.correctedMac.velocityMetersPerSecond,
        candidateVelocity, regionalTransportSettings);
    if (!candidateRegionalTransport.diagnostics.accepted) {
        throw std::runtime_error(
            "frozen scene regional momentum transport was not accepted");
    }
    const auto candidateOpeningFlux = evaluateSceneFluidOpeningFlux(
        built.surface.definition, built.surfaceState,
        geometry.openingCaps, geometry.openingQuadrature,
        geometry.openingPatches, built.grid, candidateVelocity);
    std::optional<SceneFluidMimeticTraceFlowContinuation>
        candidateTraceFlowContinuation;
    std::optional<SceneFluidMimeticRegionTransportFlowPrediction>
        candidateRegionalTransportFlowPrediction;
    SceneFluidMimeticPressureAuditEndpoint candidatePressure;
    if (settings.useCorrectedTraceFlowContinuation) {
        const auto currentBulkBaseline = sampleSceneFluidMimeticTraceFlows(
            built.pressure.controlCells,
            built.pressure.fullTraceSystem,
            geometry.pressureFaceLinks, candidateOpeningFlux, built.grid,
            candidateVelocity);
        candidateTraceFlowContinuation =
            continueSceneFluidMimeticTraceFlowsFixedTopology(
                built.correctedFlow, *previousBulkBaseline,
                currentBulkBaseline);
        candidatePressure =
            advanceSceneFluidMimeticPressureAuditFixedTopology(
                built.pressure, geometry.gridEpoch.quadrature,
                geometry.pressureControlVolumes,
                geometry.pressureFaceLinks,
                *candidateTraceFlowContinuation,
                built.pressureSettings);
    } else if (settings.useRegionalTransportFlowPrediction) {
        const auto currentBulkBaseline = sampleSceneFluidMimeticTraceFlows(
            built.pressure.controlCells,
            built.pressure.fullTraceSystem,
            geometry.pressureFaceLinks, candidateOpeningFlux, built.grid,
            candidateVelocity);
        candidateRegionalTransportFlowPrediction =
            predictSceneFluidMimeticTraceFlowsFromRegionTransportFixedTopology(
                built.pressure.controlCells,
                built.pressure.fullTraceSystem,
                geometry.pressureFaceLinks, candidateOpeningFlux,
                built.correctedFlow, candidateRegionalTransport,
                currentBulkBaseline);
        candidatePressure =
            advanceSceneFluidMimeticPressureAuditFixedTopology(
                built.pressure, geometry.gridEpoch.quadrature,
                geometry.pressureControlVolumes,
                geometry.pressureFaceLinks,
                *candidateRegionalTransportFlowPrediction,
                built.pressureSettings);
    } else {
        candidatePressure =
            advanceSceneFluidMimeticPressureAuditFixedTopology(
                built.pressure, geometry.gridEpoch.quadrature,
                geometry.pressureControlVolumes,
                geometry.pressureFaceLinks,
                candidateOpeningFlux, built.grid, candidateVelocity,
                built.pressureSettings);
    }
    if (!candidatePressure.pressureEpoch.diagnostics.accepted) {
        throw std::runtime_error(
            "frozen scene continued pressure solve was not accepted");
    }
    auto candidateCorrectedFlow =
        correctSceneFluidMimeticTraceFlows(candidatePressure);
    if (!candidateCorrectedFlow.accepted) {
        throw std::runtime_error(
            "frozen scene continued pressure correction was not accepted");
    }
    auto candidateCorrectedMac =
        collapseSceneFluidMimeticCorrectedMacVelocity(
            candidateCorrectedFlow, geometry.pressureFaceLinks,
            geometry.openingPatches, built.grid);
    auto candidateRegionalMomentum =
        reconstructSceneFluidRegionMomentumState(
            built.grid, geometry.pressureControlVolumes,
            geometry.pressureFaceLinks, geometry.openingPatches,
            candidatePressure.controlCells,
            candidateCorrectedFlow,
            candidateCorrectedMac.velocityMetersPerSecond);
    auto candidatePressureTransfer =
        evaluateSceneFluidMimeticPressureQuadrature(
            built.surface.definition, built.surfaceState, built.transfer,
            geometry.gridEpoch.quadrature,
            candidatePressure.pressureEpoch.acceptedPressureSamples);
    if (!candidatePressureTransfer.diagnostics().finite) {
        throw std::runtime_error(
            "frozen scene continued pressure transfer was not accepted");
    }

    built.pressure = std::move(candidatePressure);
    built.correctedFlow = std::move(candidateCorrectedFlow);
    built.correctedMac = std::move(candidateCorrectedMac);
    built.regionalMomentum = std::move(candidateRegionalMomentum);
    built.pressureTransfer = std::move(candidatePressureTransfer);
    built.traceFlowContinuation =
        std::move(candidateTraceFlowContinuation);
    built.regionalTransportFlowPrediction =
        std::move(candidateRegionalTransportFlowPrediction);
    built.regionalTransport = std::move(candidateRegionalTransport);
    built.bulkFlow = std::move(candidateBulkFlow);
    built.bulkProjection = candidateBulkProjection;
    built.meanStreamwiseVelocityBeforePumpMetersPerSecond = meanBeforePump;
    built.streamwisePumpIncrementMetersPerSecond = pumpIncrement;
    built.windRampFraction = candidateRampFraction;
    ++built.flowAdvanceCount;
}

void advanceMovingGeometryFsi(
    BuiltCase& built,
    const FrozenScenePressureCaseSettings& settings,
    const StructureStepSettings& structureSettings) {
    auto candidateVelocity = built.correctedMac.velocityMetersPerSecond;
    const double meanBeforePump = meanXVelocity(candidateVelocity);
    const double candidateRampFraction = windRampFraction(
        settings, built.flowAdvanceCount + 2);
    const double pumpIncrement =
        candidateRampFraction
            * settings.backgroundWindMetersPerSecond.x
        - meanBeforePump;
    for (double& value : candidateVelocity.xFaces()) {
        value += pumpIncrement;
    }
    fluid::CellScalarField candidateBulkPressure(built.grid);
    fluid::ProjectionSettings projectionSettings;
    projectionSettings.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    projectionSettings.timeStepSeconds = settings.timeStepSeconds;
    const auto candidateBulkProjection = fluid::projectVelocity(
        built.grid, candidateVelocity, candidateBulkPressure,
        projectionSettings);
    if (!candidateBulkProjection.converged) {
        throw std::runtime_error(
            "moving scene bulk continuation projection did not converge");
    }
    fluid::PeriodicFlowStrangSubcyclingSettings bulkSettings;
    bulkSettings.flow.densityKgPerCubicMeter =
        settings.densityKgPerCubicMeter;
    bulkSettings.flow.timeStepSeconds = settings.timeStepSeconds;
    bulkSettings.flow.advectionReconstruction =
        fluid::VariableMacReconstruction::DonorCell;
    bulkSettings.maximumSubsteps = 1024;
    auto candidateBulkFlow =
        fluid::advancePeriodicFlowStrangSspRk2Subcycled(
            built.grid, candidateVelocity, candidateBulkPressure,
            bulkSettings);
    if (!candidateBulkFlow.accepted) {
        throw std::runtime_error(
            "moving scene bulk flow rejected its requested interval");
    }
    SceneFluidRegionTransportSettings regionalTransportSettings;
    regionalTransportSettings.timeStepSeconds = settings.timeStepSeconds;
    regionalTransportSettings.maximumSubsteps = 1024;
    auto candidateRegionalTransport = advanceSceneFluidRegionMomentum(
        built.regionalMomentum,
        built.geometryEpoch.pressureFaceLinks,
        built.pressure.controlCells, built.correctedFlow, built.grid,
        built.correctedMac.velocityMetersPerSecond, candidateVelocity,
        regionalTransportSettings);
    if (!candidateRegionalTransport.diagnostics.accepted) {
        throw std::runtime_error(
            "moving scene regional momentum transport was not accepted");
    }

    const auto preStepForce = built.totalFluidTransfer.diagnostics()
        .transferredNodalForceNewtons;
    const double preStepForceNormNewtons = std::hypot(
        preStepForce.x, preStepForce.y, preStepForce.z);
    const auto structureCheckpoint = built.structure.checkpoint();
    try {
        built.transfer.addLoadsTo(
            built.structure, built.totalFluidTransfer);
        auto structureStep = built.structure.step(structureSettings);
        if (!structureStep.finite) {
            throw std::runtime_error(
                "moving scene Structure step was not accepted");
        }
        auto currentSurfaceState = captureSceneFluidSurfaceState(
            built.surface.definition, built.assembly.mappings,
            built.structure);
        if (currentSurfaceState.vertices.size()
            != built.surfaceState.vertices.size()) {
            throw std::logic_error(
                "moving scene surface topology changed during Structure step");
        }
        double maximumDisplacement = 0.0;
        for (std::size_t index = 0;
             index < currentSurfaceState.vertices.size(); ++index) {
            const auto& previous = built.surfaceState.vertices[index];
            const auto& current = currentSurfaceState.vertices[index];
            if (previous.id != current.id) {
                throw std::logic_error(
                    "moving scene surface identity changed during Structure step");
            }
            const double dx = current.positionMeters.x
                - previous.positionMeters.x;
            const double dy = current.positionMeters.y
                - previous.positionMeters.y;
            const double dz = current.positionMeters.z
                - previous.positionMeters.z;
            maximumDisplacement = std::max(
                maximumDisplacement,
                std::sqrt(dx * dx + dy * dy + dz * dz));
        }
        const auto acceptedGridEpoch = [&] {
            try {
                return buildSceneFluidGridEpoch(
                    built.surface.definition, currentSurfaceState,
                    built.grid, built.transfer);
            } catch (const std::exception& exception) {
                std::ostringstream message;
                message.precision(17);
                message
                    << "moving scene geometry rebuild rejected after"
                    << " structure-dx=" << maximumDisplacement << " m"
                    << " line-residual="
                    << structureStep.maximumSuspensionResidualMeters
                    << " m fluid-resultant=" << preStepForceNormNewtons
                    << " N: " << exception.what();
                throw std::runtime_error(message.str());
            }
        }();
        auto geometryTransition = [&] {
            try {
                return buildSceneFluidMimeticGeometryEpochTransition(
                    built.geometryEpoch, built.surface.definition,
                    built.surfaceState, currentSurfaceState, built.grid,
                    built.transfer, built.connectivity, acceptedGridEpoch);
            } catch (const std::exception& exception) {
                std::ostringstream message;
                message.precision(17);
                message
                    << "moving scene geometry transition rejected after"
                    << " structure-dx=" << maximumDisplacement << " m"
                    << " line-residual="
                    << structureStep.maximumSuspensionResidualMeters
                    << " m fluid-resultant=" << preStepForceNormNewtons
                    << " N: " << exception.what();
                throw std::runtime_error(message.str());
            }
        }();
        auto volumeRates = buildSceneFluidPressureVolumeRates(
            built.geometryEpoch.cellVolumes,
            geometryTransition.currentGeometryEpoch.cellVolumes,
            geometryTransition.currentGeometryEpoch.pressureControlVolumes,
            geometryTransition.topologyTransition);
        auto regionalRebase = rebaseSceneFluidRegionTransport(
            candidateRegionalTransport,
            built.geometryEpoch.pressureControlVolumes,
            geometryTransition.currentGeometryEpoch.pressureControlVolumes,
            geometryTransition.topologyTransition);
        const auto openingFlux = evaluateSceneFluidOpeningFlux(
            built.surface.definition, currentSurfaceState,
            geometryTransition.currentGeometryEpoch.openingCaps,
            geometryTransition.currentGeometryEpoch.openingQuadrature,
            geometryTransition.currentGeometryEpoch.openingPatches,
            built.grid, candidateVelocity);
        SceneFluidRegionWallSettings wallSettings;
        wallSettings.timeStepSeconds = volumeRates.durationSeconds;
        auto regionWall = exchangeSceneFluidRegionWallMomentum(
            regionalRebase, built.grid,
            geometryTransition.currentGeometryEpoch.pressureControlVolumes,
            built.surface.definition, currentSurfaceState,
            geometryTransition.currentGeometryEpoch.gridEpoch.quadrature,
            wallSettings);
        if (!regionWall.diagnostics.accepted) {
            throw std::runtime_error(
                "moving scene material-wall exchange was not accepted");
        }
        auto candidatePressure =
            buildSceneFluidMimeticPressureAuditEndpoint(
                built.surface.definition, currentSurfaceState, built.grid,
                geometryTransition.currentGeometryEpoch, openingFlux,
                regionWall, volumeRates,
                geometryTransition.topologyTransition, built.pressure,
                built.pressureSettings);
        if (!candidatePressure.pressureEpoch.diagnostics.accepted) {
            throw std::runtime_error(
                "moving scene continued pressure solve was not accepted");
        }
        auto candidateCorrectedFlow =
            correctSceneFluidMimeticTraceFlows(candidatePressure);
        if (!candidateCorrectedFlow.accepted) {
            throw std::runtime_error(
                "moving scene pressure correction was not accepted");
        }
        auto candidateCorrectedMac =
            collapseSceneFluidMimeticCorrectedMacVelocity(
                candidateCorrectedFlow,
                geometryTransition.currentGeometryEpoch.pressureFaceLinks,
                geometryTransition.currentGeometryEpoch.openingPatches,
                built.grid);
        auto candidateRegionalMomentum =
            reconstructSceneFluidRegionMomentumState(
                built.grid,
                geometryTransition.currentGeometryEpoch
                    .pressureControlVolumes,
                geometryTransition.currentGeometryEpoch.pressureFaceLinks,
                geometryTransition.currentGeometryEpoch.openingPatches,
                candidatePressure.controlCells, candidateCorrectedFlow,
                candidateCorrectedMac.velocityMetersPerSecond);
        auto acceptedWall = captureSceneFluidAcceptedWallTractions(
            regionWall);
        auto candidateFluidLoads = evaluateFluidLoads(
            built.surface.definition, currentSurfaceState, built.transfer,
            geometryTransition.currentGeometryEpoch.gridEpoch.quadrature,
            candidatePressure.pressureEpoch.acceptedPressureSamples,
            &acceptedWall);
        if (!candidateFluidLoads.pressure.diagnostics().finite
            || !candidateFluidLoads.wall
            || !candidateFluidLoads.wall->diagnostics().finite
            || !candidateFluidLoads.total.diagnostics().finite) {
            throw std::runtime_error(
                "moving scene fluid load transfer was not accepted");
        }

        built.surfaceState = std::move(currentSurfaceState);
        built.pressureControlTopologyStable =
            geometryTransition.controlVolumeTopologyStable;
        built.geometryEpoch =
            std::move(geometryTransition.currentGeometryEpoch);
        built.pressure = std::move(candidatePressure);
        built.correctedFlow = std::move(candidateCorrectedFlow);
        built.correctedMac = std::move(candidateCorrectedMac);
        built.regionalMomentum = std::move(candidateRegionalMomentum);
        built.pressureTransfer =
            std::move(candidateFluidLoads.pressure);
        built.wallTransfer = std::move(candidateFluidLoads.wall);
        built.totalFluidTransfer = std::move(candidateFluidLoads.total);
        built.regionalTransport = std::move(candidateRegionalTransport);
        built.regionWall = std::move(regionWall);
        built.structureStep = structureStep;
        built.bulkFlow = std::move(candidateBulkFlow);
        built.bulkProjection = candidateBulkProjection;
        built.meanStreamwiseVelocityBeforePumpMetersPerSecond =
            meanBeforePump;
        built.streamwisePumpIncrementMetersPerSecond = pumpIncrement;
        built.windRampFraction = candidateRampFraction;
        built.maximumGeometryDisplacementMeters = maximumDisplacement;
        ++built.flowAdvanceCount;
        ++built.geometryAdvanceCount;
    } catch (...) {
        built.structure.restore(structureCheckpoint);
        throw;
    }
}

FrozenScenePressureCaseDiagnostics makeDiagnostics(
    const BuiltCase& built) {
    FrozenScenePressureCaseDiagnostics result;
    result.gridCellCounts = built.grid.cellCounts();
    result.gridLowerMeters = built.grid.lowerMeters();
    result.gridUpperMeters = built.grid.upperMeters();
    result.windRampFraction = built.windRampFraction;
    result.usesCorrectedTraceFlowContinuation =
        built.traceFlowContinuation.has_value();
    result.usesRegionalTransportFlowPrediction =
        built.regionalTransportFlowPrediction.has_value();
    result.usesMovingGeometryFsi = built.geometryAdvanceCount > 0;
    result.usesConsecutivePressureWarmStart =
        built.pressure.usesConsecutiveWarmStart;
    result.usesRegionWallPrediction =
        built.pressure.usesRegionWallPrediction;
    result.pressureControlTopologyStable =
        built.pressureControlTopologyStable;
    result.geometryAdvanceCount = built.geometryAdvanceCount;
    result.maximumGeometryDisplacementMeters =
        built.maximumGeometryDisplacementMeters;
    if (built.structureStep) {
        result.maximumSuspensionResidualMeters =
            built.structureStep->maximumSuspensionResidualMeters;
    }
    if (built.traceFlowContinuation) {
        result.maximumCarriedTraceCorrectionCubicMetersPerSecond =
            built.traceFlowContinuation
                ->maximumAbsoluteCarriedCorrectionCubicMetersPerSecond;
        result.maximumTraceBulkIncrementCubicMetersPerSecond =
            built.traceFlowContinuation
                ->maximumAbsoluteBulkIncrementCubicMetersPerSecond;
    }
    if (built.regionalTransportFlowPrediction) {
        result
            .maximumRegionalTransportFlowDifferenceFromBulkBaselineCubicMetersPerSecond =
            built.regionalTransportFlowPrediction
                ->maximumAbsoluteFlowDifferenceFromBulkBaselineCubicMetersPerSecond;
    }
    result.pressureControlCount =
        built.pressure.controlCells.controlCells.size();
    result.sharedTraceCount =
        built.pressure.fullTraceSystem.sharedTraceCount;
    result.pressureIterationCount = built.pressure.pressureEpoch
        .diagnostics.pressureSolve.reducedTraceSolve.iterationCount;
    result.extrapolatedZeroVolumePressureSideCount =
        built.pressure.pressureEpoch.acceptedPressureSamples
            .extrapolatedZeroVolumeSideCount;
    result.maximumPressureExtrapolationDistanceMeters =
        built.pressure.pressureEpoch.acceptedPressureSamples
            .maximumExtrapolationDistanceMeters;
    result.maximumAbsolutePressureDifferencePascals =
        built.pressure.pressureEpoch.acceptedPressureSamples
            .maximumAbsolutePressureDifferencePascals;
    result.maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond =
        built.pressure.pressureSources
            .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond;
    result.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond =
        built.correctedFlow
            .maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond;
    result.correctedContinuityToleranceCubicMetersPerSecond =
        built.correctedFlow
            .correctedContinuityToleranceCubicMetersPerSecond;
    result.maximumCollapsedMacVelocityMetersPerSecond =
        built.correctedMac.diagnostics
            .maximumAbsoluteVelocityMetersPerSecond;
    result.maximumCollapsedSubfaceVelocityDeviationMetersPerSecond =
        built.correctedMac.diagnostics
            .maximumSubfaceVelocityDeviationMetersPerSecond;
    result.embeddedOpeningTraceCount =
        built.correctedMac.diagnostics.embeddedOpeningTraceCount;
    result.maximumRegionalVelocityMetersPerSecond =
        built.regionalMomentum.diagnostics
            .maximumAbsoluteVelocityMetersPerSecond;
    result.maximumRegionalLinkVelocityResidualMetersPerSecond =
        built.regionalMomentum.diagnostics
            .maximumLinkNormalVelocityResidualMetersPerSecond;
    result.regionalKineticEnergyJoules =
        built.regionalMomentum.diagnostics.kineticEnergyJoules;
    result.flowAdvanceCount = built.flowAdvanceCount;
    if (built.regionalTransport) {
        result.regionalTransportSubstepCount =
            built.regionalTransport->diagnostics.substepCount;
        result.regionalTransportMaximumVelocityChangeMetersPerSecond =
            built.regionalTransport->diagnostics
                .maximumVelocityChangeMetersPerSecond;
        result.regionalTransportMomentumResidualKilogramMetersPerSecond =
            built.regionalTransport->diagnostics
                .momentumResidualNormKilogramMetersPerSecond;
        result.regionalTransportAdvectiveEnergyLossJoules =
            built.regionalTransport->diagnostics
                .advectiveEnergyLossJoules;
        result.regionalTransportViscousEnergyLossJoules =
            built.regionalTransport->diagnostics.viscousEnergyLossJoules;
    }
    result.meanStreamwiseVelocityBeforePumpMetersPerSecond =
        built.meanStreamwiseVelocityBeforePumpMetersPerSecond;
    result.streamwisePumpIncrementMetersPerSecond =
        built.streamwisePumpIncrementMetersPerSecond;
    if (built.bulkFlow) {
        result.bulkFlowSubstepCount =
            built.bulkFlow->completedSubstepCount;
        result.bulkFlowMaximumVelocityChangeMetersPerSecond =
            built.bulkFlow->maximumVelocityChangeMetersPerSecond;
    }
    if (built.bulkProjection) {
        result.bulkProjectionDivergenceBeforePerSecond =
            built.bulkProjection->divergenceL2BeforePerSecond;
        result.bulkProjectionDivergenceAfterPerSecond =
            built.bulkProjection->divergenceL2AfterPerSecond;
    }
    const auto& transfer = built.pressureTransfer.diagnostics();
    result.pressureForceNewtons = transfer.transferredNodalForceNewtons;
    result.pressureMomentNewtonMeters =
        transfer.transferredNodalMomentNewtonMeters;
    result.transferForceResidualNewtons = transfer.forceResidualNormNewtons;
    result.transferMomentResidualNewtonMeters =
        transfer.momentResidualNormNewtonMeters;
    const auto& totalTransfer = built.totalFluidTransfer.diagnostics();
    result.totalFluidForceNewtons =
        totalTransfer.transferredNodalForceNewtons;
    result.totalFluidTransferForceResidualNewtons =
        totalTransfer.forceResidualNormNewtons;
    result.totalFluidTransferMomentResidualNewtonMeters =
        totalTransfer.momentResidualNormNewtonMeters;
    if (built.wallTransfer) {
        result.wallForceNewtons = built.wallTransfer->diagnostics()
            .transferredNodalForceNewtons;
    }
    if (built.regionWall) {
        result.wallMomentumResidualKilogramMetersPerSecond =
            built.regionWall->diagnostics
                .momentumResidualNormKilogramMetersPerSecond;
    }
    result.finite = built.pressure.pressureEpoch.diagnostics.accepted
        && built.correctedFlow.accepted
        && built.correctedMac.diagnostics.finite
        && transfer.finite
        && std::isfinite(
            result.maximumPressureExtrapolationDistanceMeters)
        && std::isfinite(result.windRampFraction)
        && result.windRampFraction >= 0.0
        && result.windRampFraction <= 1.0
        && std::isfinite(result.maximumAbsolutePressureDifferencePascals)
        && std::isfinite(
            result.maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond)
        && std::isfinite(
            result.maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond)
        && std::isfinite(
            result.correctedContinuityToleranceCubicMetersPerSecond)
        && std::isfinite(result.maximumCollapsedMacVelocityMetersPerSecond)
        && std::isfinite(
            result.maximumCollapsedSubfaceVelocityDeviationMetersPerSecond)
        && built.regionalMomentum.diagnostics.finite
        && totalTransfer.finite
        && (!built.wallTransfer
            || built.wallTransfer->diagnostics().finite)
        && (!built.regionWall
            || built.regionWall->diagnostics.accepted)
        && (!built.structureStep || built.structureStep->finite)
        && std::isfinite(result.maximumRegionalVelocityMetersPerSecond)
        && std::isfinite(
            result.maximumRegionalLinkVelocityResidualMetersPerSecond)
        && std::isfinite(result.regionalKineticEnergyJoules)
        && (!built.regionalTransport
            || built.regionalTransport->diagnostics.accepted)
        && std::isfinite(
            result.regionalTransportMaximumVelocityChangeMetersPerSecond)
        && std::isfinite(
            result.regionalTransportMomentumResidualKilogramMetersPerSecond)
        && std::isfinite(
            result.regionalTransportAdvectiveEnergyLossJoules)
        && std::isfinite(
            result.regionalTransportViscousEnergyLossJoules)
        && std::isfinite(
            result.maximumCarriedTraceCorrectionCubicMetersPerSecond)
        && std::isfinite(
            result.maximumTraceBulkIncrementCubicMetersPerSecond)
        && std::isfinite(
            result
                .maximumRegionalTransportFlowDifferenceFromBulkBaselineCubicMetersPerSecond)
        && std::isfinite(
            result.bulkFlowMaximumVelocityChangeMetersPerSecond)
        && std::isfinite(
            result.meanStreamwiseVelocityBeforePumpMetersPerSecond)
        && std::isfinite(result.streamwisePumpIncrementMetersPerSecond)
        && std::isfinite(result.bulkProjectionDivergenceBeforePerSecond)
        && std::isfinite(result.bulkProjectionDivergenceAfterPerSecond)
        && std::isfinite(result.transferForceResidualNewtons)
        && std::isfinite(result.transferMomentResidualNewtonMeters)
        && std::isfinite(result.maximumGeometryDisplacementMeters)
        && std::isfinite(result.maximumSuspensionResidualMeters)
        && std::isfinite(
            result.wallMomentumResidualKilogramMetersPerSecond)
        && std::isfinite(
            result.totalFluidTransferForceResidualNewtons)
        && std::isfinite(
            result.totalFluidTransferMomentResidualNewtonMeters);
    return result;
}

void appendFluidDiagnostics(
    viewer::DiagnosticFrame& frame,
    const fluid::PeriodicCartesianGrid& grid,
    const fluid::MacVelocityField& velocity) {
    const auto fluidFields = viewer::buildPeriodicFluidCellFields(
        grid, velocity);
    const std::size_t structureVertexCount = frame.vertices.size();
    const std::size_t fluidVertexCount = grid.cellCount();
    if (fluidFields.velocityMetersPerSecond.size() != fluidVertexCount
        || fluidFields.speedMetersPerSecond.size() != fluidVertexCount
        || fluidFields.divergencePerSecond.size() != fluidVertexCount
        || fluidFields.vorticityPerSecond.size() != fluidVertexCount
        || fluidFields.vorticityMagnitudePerSecond.size()
            != fluidVertexCount) {
        throw std::logic_error(
            "frozen scene fluid viewer fields are incomplete");
    }
    for (auto& field : frame.scalarFields) {
        if (field.association == viewer::FieldAssociation::Vertex) {
            field.values.resize(
                structureVertexCount + fluidVertexCount, 0.0);
        }
    }
    for (auto& field : frame.vectorFields) {
        if (field.association == viewer::FieldAssociation::Vertex) {
            field.values.resize(
                structureVertexCount + fluidVertexCount, {});
        }
    }

    std::set<std::uint64_t> usedIds;
    std::uint64_t candidateId = 0;
    for (const auto& vertex : frame.vertices) {
        usedIds.insert(vertex.stableId);
        candidateId = std::max(candidateId, vertex.stableId);
    }
    const auto nextAvailableId = [&]() {
        do {
            candidateId = candidateId
                    == std::numeric_limits<std::uint64_t>::max()
                ? 1 : candidateId + 1;
        } while (candidateId == 0 || usedIds.contains(candidateId));
        usedIds.insert(candidateId);
        return candidateId;
    };
    frame.vertices.reserve(structureVertexCount + fluidVertexCount);
    const auto counts = grid.cellCounts();
    for (std::size_t k = 0; k < counts.z; ++k) {
        for (std::size_t j = 0; j < counts.y; ++j) {
            for (std::size_t i = 0; i < counts.x; ++i) {
                const auto position = grid.cellCenterMeters(i, j, k);
                frame.vertices.push_back({
                    nextAvailableId(),
                    {position.x, position.y, position.z},
                });
            }
        }
    }

    const auto scalarValues = [&](const std::vector<double>& values) {
        std::vector<double> result(
            structureVertexCount + fluidVertexCount, 0.0);
        std::copy(
            values.begin(), values.end(),
            result.begin() + static_cast<std::ptrdiff_t>(
                structureVertexCount));
        return result;
    };
    const auto vectorValues = [&](const std::vector<viewer::Vec3d>& values) {
        std::vector<viewer::Vec3d> result(
            structureVertexCount + fluidVertexCount);
        std::copy(
            values.begin(), values.end(),
            result.begin() + static_cast<std::ptrdiff_t>(
                structureVertexCount));
        return result;
    };
    frame.scalarFields.push_back({
        "frozen_scene.fluid_speed", "m/s",
        viewer::FieldAssociation::Vertex,
        scalarValues(fluidFields.speedMetersPerSecond),
    });
    frame.scalarFields.push_back({
        "frozen_scene.fluid_divergence", "1/s",
        viewer::FieldAssociation::Vertex,
        scalarValues(fluidFields.divergencePerSecond),
    });
    frame.scalarFields.push_back({
        "frozen_scene.fluid_vorticity_magnitude", "1/s",
        viewer::FieldAssociation::Vertex,
        scalarValues(fluidFields.vorticityMagnitudePerSecond),
    });
    frame.vectorFields.push_back({
        "frozen_scene.fluid_velocity", "m/s",
        viewer::FieldAssociation::Vertex,
        vectorValues(fluidFields.velocityMetersPerSecond),
    });
    frame.vectorFields.push_back({
        "frozen_scene.fluid_vorticity", "1/s",
        viewer::FieldAssociation::Vertex,
        vectorValues(fluidFields.vorticityPerSecond),
    });
}

} // namespace

struct FrozenScenePressureCase::Implementation {
    Implementation(
        BuiltCase builtValue,
        const FrozenScenePressureCaseSettings& settingsValue)
        : built(std::move(builtValue)), settings(settingsValue),
          diagnostics(makeDiagnostics(built)) {
        stepSettings = makeStructureStepSettings(
            built.assembly.definition, settings);
        diagnostics.backgroundWindMetersPerSecond =
            settings.backgroundWindMetersPerSecond;
        diagnostics.windRampSeconds = settings.windRampSeconds;
        diagnostics.diagnosticPerturbationSpeedMetersPerSecond =
            settings.diagnosticPerturbationSpeedMetersPerSecond;
        diagnostics.usesMovingGeometryFsi =
            settings.useMovingGeometryFsi;
        if (!diagnostics.finite) {
            throw std::runtime_error(
                "frozen scene pressure diagnostics are non-finite");
        }
    }

    BuiltCase built;
    FrozenScenePressureCaseSettings settings;
    FrozenScenePressureCaseDiagnostics diagnostics;
    StructureStepSettings stepSettings;
    std::uint64_t acceptedStepCount = 0;
    double simulationTimeSeconds = 0.0;
};

FrozenScenePressureCase::FrozenScenePressureCase(
    Scene scene,
    const FrozenScenePressureCaseSettings& settings)
    : implementation_(std::make_unique<Implementation>(
          buildCase(std::move(scene), settings), settings)) {}

FrozenScenePressureCase::~FrozenScenePressureCase() = default;
FrozenScenePressureCase::FrozenScenePressureCase(
    FrozenScenePressureCase&&) noexcept = default;
FrozenScenePressureCase& FrozenScenePressureCase::operator=(
    FrozenScenePressureCase&&) noexcept = default;

viewer::TraceHeader FrozenScenePressureCase::traceHeader() const {
    return {
        implementation_->built.scene.metadata.designChecksum,
        implementation_->settings.useMovingGeometryFsi
            ? movingScenePressureSolverId
            : frozenScenePressureSolverId,
    };
}

viewer::DiagnosticFrame FrozenScenePressureCase::advance() {
    auto& state = *implementation_;
    if (state.settings.useMovingGeometryFsi) {
        advanceMovingGeometryFsi(
            state.built, state.settings, state.stepSettings);
        state.diagnostics = makeDiagnostics(state.built);
        state.diagnostics.usesMovingGeometryFsi = true;
        state.diagnostics.backgroundWindMetersPerSecond =
            state.settings.backgroundWindMetersPerSecond;
        state.diagnostics.windRampSeconds =
            state.settings.windRampSeconds;
        state.diagnostics.diagnosticPerturbationSpeedMetersPerSecond =
            state.settings.diagnosticPerturbationSpeedMetersPerSecond;
    } else if (state.acceptedStepCount > 0) {
        advanceFixedGeometryFlow(state.built, state.settings);
        state.diagnostics = makeDiagnostics(state.built);
        state.diagnostics.backgroundWindMetersPerSecond =
            state.settings.backgroundWindMetersPerSecond;
        state.diagnostics.windRampSeconds =
            state.settings.windRampSeconds;
        state.diagnostics.diagnosticPerturbationSpeedMetersPerSecond =
            state.settings.diagnosticPerturbationSpeedMetersPerSecond;
    }
    const std::uint64_t nextStep = state.acceptedStepCount + 1;
    const double nextTime = state.simulationTimeSeconds
        + state.settings.timeStepSeconds;
    viewer::StructureFrameContext context;
    context.sceneChecksum =
        state.built.scene.metadata.designChecksum;
    context.solverCommit = state.settings.useMovingGeometryFsi
        ? movingScenePressureSolverId
        : frozenScenePressureSolverId;
    context.timeStepSeconds = state.settings.timeStepSeconds;
    context.couplingResiduals.fluid = state.diagnostics
        .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond;
    context.couplingResiduals.tractionNewtons =
        state.settings.useMovingGeometryFsi
        ? state.diagnostics.totalFluidTransferForceResidualNewtons
        : state.diagnostics.transferForceResidualNewtons;
    const auto& acceptedTransfer = state.settings.useMovingGeometryFsi
        ? state.built.totalFluidTransfer
        : state.built.pressureTransfer;
    context.conservation.interfaceForceResidualNewtons = {
        acceptedTransfer.diagnostics().forceResidualNewtons.x,
        acceptedTransfer.diagnostics().forceResidualNewtons.y,
        acceptedTransfer.diagnostics().forceResidualNewtons.z,
    };
    context.conservation.interfaceMomentResidualNewtonMetres = {
        acceptedTransfer.diagnostics()
            .momentResidualNewtonMeters.x,
        acceptedTransfer.diagnostics()
            .momentResidualNewtonMeters.y,
        acceptedTransfer.diagnostics()
            .momentResidualNewtonMeters.z,
    };
    auto frame = viewer::buildStructureFrame(
        state.built.structure, state.built.frameMapping, context);
    frame.step = nextStep;
    frame.simulationTimeSeconds = nextTime;
    appendFluidDiagnostics(
        frame, state.built.grid,
        state.built.correctedMac.velocityMetersPerSecond);

    std::map<StableId, std::size_t> triangleIndices;
    for (std::size_t index = 0;
         index < state.built.assembly.mappings.triangleIds.size();
         ++index) {
        triangleIndices.emplace(
            state.built.assembly.mappings.triangleIds[index], index);
    }
    std::vector<double> trianglePressure(
        state.built.structure.definition().triangles.size(), 0.0);
    std::vector<double> triangleArea(trianglePressure.size(), 0.0);
    const auto& quadrature =
        state.built.geometryEpoch.gridEpoch.quadrature;
    const auto& samples = state.built.pressure.pressureEpoch
        .acceptedPressureSamples;
    if (samples.bindings.size() != quadrature.points.size()) {
        throw std::logic_error(
            "frozen scene pressure sample count changed");
    }
    for (std::size_t index = 0; index < quadrature.points.size(); ++index) {
        const auto found = triangleIndices.find(
            quadrature.points[index].triangleId);
        if (found == triangleIndices.end()) {
            throw std::logic_error(
                "frozen scene pressure triangle is missing");
        }
        const double area = quadrature.points[index].areaSquareMeters;
        trianglePressure[found->second] +=
            samples.bindings[index].pressureDifferencePascals * area;
        triangleArea[found->second] += area;
    }
    for (std::size_t index = 0; index < trianglePressure.size(); ++index) {
        if (!(triangleArea[index] > 0.0)) {
            throw std::logic_error(
                "frozen scene triangle has no pressure quadrature");
        }
        trianglePressure[index] /= triangleArea[index];
    }
    frame.scalarFields.push_back({
        "frozen_scene.pressure_jump", "Pa",
        viewer::FieldAssociation::Triangle, std::move(trianglePressure),
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_pressure_jump", "Pa",
        viewer::FieldAssociation::Global,
        {state.diagnostics.maximumAbsolutePressureDifferencePascals},
    });
    frame.scalarFields.push_back({
        "frozen_scene.pressure_iterations", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.pressureIterationCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.pressure_controls", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.pressureControlCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.shared_traces", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.sharedTraceCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.component_continuity_residual", "m^3/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .maximumAbsoluteComponentContinuityResidualCubicMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.corrected_continuity_residual", "m^3/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .maximumAbsoluteCorrectedContinuityResidualCubicMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.corrected_continuity_tolerance", "m^3/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.correctedContinuityToleranceCubicMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_collapsed_mac_velocity", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.maximumCollapsedMacVelocityMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_subface_velocity_deviation", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .maximumCollapsedSubfaceVelocityDeviationMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_regional_velocity", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.maximumRegionalVelocityMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_regional_link_velocity_residual", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .maximumRegionalLinkVelocityResidualMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.regional_kinetic_energy", "J",
        viewer::FieldAssociation::Global,
        {state.diagnostics.regionalKineticEnergyJoules},
    });
    frame.scalarFields.push_back({
        "frozen_scene.regional_transport_substeps", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(
            state.diagnostics.regionalTransportSubstepCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.regional_transport_velocity_change", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .regionalTransportMaximumVelocityChangeMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.regional_transport_momentum_residual", "kg*m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .regionalTransportMomentumResidualKilogramMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.regional_transport_advective_loss", "J",
        viewer::FieldAssociation::Global,
        {state.diagnostics.regionalTransportAdvectiveEnergyLossJoules},
    });
    frame.scalarFields.push_back({
        "frozen_scene.regional_transport_viscous_loss", "J",
        viewer::FieldAssociation::Global,
        {state.diagnostics.regionalTransportViscousEnergyLossJoules},
    });
    frame.scalarFields.push_back({
        "frozen_scene.regional_transport_flow_prediction", "1",
        viewer::FieldAssociation::Global,
        {state.diagnostics.usesRegionalTransportFlowPrediction ? 1.0 : 0.0},
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_regional_transport_flow_delta", "m^3/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .maximumRegionalTransportFlowDifferenceFromBulkBaselineCubicMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.embedded_opening_traces", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.embeddedOpeningTraceCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.corrected_trace_continuation", "1",
        viewer::FieldAssociation::Global,
        {state.diagnostics.usesCorrectedTraceFlowContinuation ? 1.0 : 0.0},
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_carried_trace_correction", "m^3/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .maximumCarriedTraceCorrectionCubicMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.maximum_trace_bulk_increment", "m^3/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.maximumTraceBulkIncrementCubicMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.flow_advances", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.flowAdvanceCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.grid_cells", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(
             state.diagnostics.gridCellCounts.x
             * state.diagnostics.gridCellCounts.y
             * state.diagnostics.gridCellCounts.z)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.initial_perturbation", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics
             .diagnosticPerturbationSpeedMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.wind_ramp_fraction", "1",
        viewer::FieldAssociation::Global,
        {state.diagnostics.windRampFraction},
    });
    frame.scalarFields.push_back({
        "frozen_scene.wind_ramp_duration", "s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.windRampSeconds},
    });
    frame.scalarFields.push_back({
        "frozen_scene.bulk_substeps", "1",
        viewer::FieldAssociation::Global,
        {static_cast<double>(state.diagnostics.bulkFlowSubstepCount)},
    });
    frame.scalarFields.push_back({
        "frozen_scene.bulk_velocity_change", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.bulkFlowMaximumVelocityChangeMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.streamwise_pump_increment", "m/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.streamwisePumpIncrementMetersPerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.bulk_projection_divergence_before", "1/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.bulkProjectionDivergenceBeforePerSecond},
    });
    frame.scalarFields.push_back({
        "frozen_scene.bulk_projection_divergence_after", "1/s",
        viewer::FieldAssociation::Global,
        {state.diagnostics.bulkProjectionDivergenceAfterPerSecond},
    });

    std::vector<viewer::Vec3d> nodalPressureForces(
        frame.vertices.size());
    for (const auto& load : state.built.pressureTransfer.nodeLoads()) {
        nodalPressureForces[load.structureNode] = {
            load.forceNewtons.x, load.forceNewtons.y,
            load.forceNewtons.z,
        };
    }
    frame.vectorFields.push_back({
        "frozen_scene.nodal_pressure_force", "N",
        viewer::FieldAssociation::Vertex,
        std::move(nodalPressureForces),
    });
    frame.vectorFields.push_back({
        "frozen_scene.total_pressure_force", "N",
        viewer::FieldAssociation::Global,
        {{state.diagnostics.pressureForceNewtons.x,
          state.diagnostics.pressureForceNewtons.y,
          state.diagnostics.pressureForceNewtons.z}},
    });
    frame.vectorFields.push_back({
        "frozen_scene.total_pressure_moment", "N m",
        viewer::FieldAssociation::Global,
        {{state.diagnostics.pressureMomentNewtonMeters.x,
          state.diagnostics.pressureMomentNewtonMeters.y,
          state.diagnostics.pressureMomentNewtonMeters.z}},
    });
    frame.vectorFields.push_back({
        "frozen_scene.background_wind", "m/s",
        viewer::FieldAssociation::Global,
        {{state.settings.backgroundWindMetersPerSecond.x,
          state.settings.backgroundWindMetersPerSecond.y,
          state.settings.backgroundWindMetersPerSecond.z}},
    });
    frame.vectorFields.push_back({
        "frozen_scene.current_target_wind", "m/s",
        viewer::FieldAssociation::Global,
        {{state.diagnostics.windRampFraction
              * state.settings.backgroundWindMetersPerSecond.x,
          state.diagnostics.windRampFraction
              * state.settings.backgroundWindMetersPerSecond.y,
          state.diagnostics.windRampFraction
              * state.settings.backgroundWindMetersPerSecond.z}},
    });
    if (state.settings.useMovingGeometryFsi) {
        frame.scalarFields.push_back({
            "moving_scene.geometry_advances", "1",
            viewer::FieldAssociation::Global,
            {static_cast<double>(state.diagnostics.geometryAdvanceCount)},
        });
        frame.scalarFields.push_back({
            "moving_scene.maximum_geometry_displacement", "m",
            viewer::FieldAssociation::Global,
            {state.diagnostics.maximumGeometryDisplacementMeters},
        });
        frame.scalarFields.push_back({
            "moving_scene.maximum_suspension_residual", "m",
            viewer::FieldAssociation::Global,
            {state.diagnostics.maximumSuspensionResidualMeters},
        });
        frame.scalarFields.push_back({
            "moving_scene.consecutive_pressure_warm_start", "1",
            viewer::FieldAssociation::Global,
            {state.diagnostics.usesConsecutivePressureWarmStart
                 ? 1.0 : 0.0},
        });
        frame.scalarFields.push_back({
            "moving_scene.wall_momentum_residual", "kg*m/s",
            viewer::FieldAssociation::Global,
            {state.diagnostics
                 .wallMomentumResidualKilogramMetersPerSecond},
        });
        std::vector<viewer::Vec3d> nodalWallForces(
            frame.vertices.size());
        if (state.built.wallTransfer) {
            for (const auto& load :
                 state.built.wallTransfer->nodeLoads()) {
                nodalWallForces[load.structureNode] = {
                    load.forceNewtons.x, load.forceNewtons.y,
                    load.forceNewtons.z,
                };
            }
        }
        frame.vectorFields.push_back({
            "moving_scene.nodal_wall_force", "N",
            viewer::FieldAssociation::Vertex,
            std::move(nodalWallForces),
        });
        std::vector<viewer::Vec3d> nodalTotalFluidForces(
            frame.vertices.size());
        for (const auto& load :
             state.built.totalFluidTransfer.nodeLoads()) {
            nodalTotalFluidForces[load.structureNode] = {
                load.forceNewtons.x, load.forceNewtons.y,
                load.forceNewtons.z,
            };
        }
        frame.vectorFields.push_back({
            "moving_scene.nodal_total_fluid_force", "N",
            viewer::FieldAssociation::Vertex,
            std::move(nodalTotalFluidForces),
        });
        frame.vectorFields.push_back({
            "moving_scene.total_fluid_force", "N",
            viewer::FieldAssociation::Global,
            {{state.diagnostics.totalFluidForceNewtons.x,
              state.diagnostics.totalFluidForceNewtons.y,
              state.diagnostics.totalFluidForceNewtons.z}},
        });
    }

    viewer::ProtocolError error;
    if (!viewer::validateFrame(frame, &error)) {
        throw std::runtime_error(
            "frozen scene pressure frame is invalid: " + error.message);
    }
    state.acceptedStepCount = nextStep;
    state.simulationTimeSeconds = nextTime;
    return frame;
}

const StructureStepSettings&
FrozenScenePressureCase::stepSettings() const noexcept {
    return implementation_->stepSettings;
}

std::uint64_t FrozenScenePressureCase::acceptedStepCount() const noexcept {
    return implementation_->acceptedStepCount;
}

double FrozenScenePressureCase::simulationTimeSeconds() const noexcept {
    return implementation_->simulationTimeSeconds;
}

const FrozenScenePressureCaseDiagnostics&
FrozenScenePressureCase::diagnostics() const noexcept {
    return implementation_->diagnostics;
}

} // namespace simwing::fsi
