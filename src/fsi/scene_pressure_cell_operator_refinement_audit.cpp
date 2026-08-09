#include "scene_pressure_cell_operator_refinement_audit.h"

#include "scene_pressure_cell_operator_phase_audit.h"

#include "scene_fluid_mimetic_condensed_trace_system.h"
#include "scene_fluid_mimetic_control_cell.h"
#include "scene_fluid_mimetic_trace_system.h"
#include "scene_fluid_surface_transfer.h"
#include "scene_pressure_cell_geometry.h"

#include <algorithm>
#include <bit>
#include <cmath>
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

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

std::size_t gridCellCount(const fluid::GridCellCounts counts) {
    std::size_t xy = 0;
    std::size_t xyz = 0;
    if (!checkedMultiply(counts.x, counts.y, xy)
        || !checkedMultiply(xy, counts.z, xyz)) {
        throw std::length_error(
            "scene pressure-cell refinement grid count overflows");
    }
    return xyz;
}

std::size_t storageBytes(
    const ScenePressureCellOperatorRefinementAudit& audit) {
    if (audit.samples.size() > std::numeric_limits<std::size_t>::max()
            / sizeof(ScenePressureCellOperatorRefinementSample)) {
        throw std::length_error(
            "scene pressure-cell refinement storage overflows");
    }
    std::size_t result = audit.samples.size()
        * sizeof(ScenePressureCellOperatorRefinementSample);
    for (const auto& sample : audit.samples) {
        if (sample.response.ownedStorageBytes
            > std::numeric_limits<std::size_t>::max() - result) {
            throw std::length_error(
                "scene pressure-cell refinement storage overflows");
        }
        result += sample.response.ownedStorageBytes;
    }
    return result;
}

void fingerprintResponseSettings(
    Fingerprint& fingerprint,
    const SceneFluidPressureOperatorResponseAuditSettings& settings) {
    fingerprint.real(settings.manufacturedPressureL2Pascals);
    fingerprint.real(
        settings.graphSolve.absoluteResidualTolerancePascalsMeters);
    fingerprint.real(settings.graphSolve.relativeResidualTolerance);
    fingerprint.real(settings.graphSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        settings.graphSolve.maximumIterations));
    fingerprint.real(
        settings.shadowSolve.absoluteResidualTolerancePascalsMeters);
    fingerprint.real(settings.shadowSolve.relativeResidualTolerance);
    fingerprint.real(settings.shadowSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        settings.shadowSolve.maximumIterations));
}

bool validPhaseComponent(const double phase) {
    return std::isfinite(phase)
        && phase >= -0.5 && phase < 0.5
        && (phase != 0.0 || !std::signbit(phase));
}

void validateSettings(
    const ScenePressureCellOperatorRefinementAuditSettings& settings) {
    if (!validPhaseComponent(settings.gridPhaseFraction.x)
        || !validPhaseComponent(settings.gridPhaseFraction.y)
        || !validPhaseComponent(settings.gridPhaseFraction.z)) {
        throw std::invalid_argument(
            "scene pressure-cell refinement grid phase is invalid");
    }
}

std::uint64_t productFingerprint(
    const ScenePressureCellOperatorRefinementAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.structureDefinitionFingerprint);
    fingerprintResponseSettings(fingerprint, audit.settings.response);
    fingerprint.real(audit.settings.gridPhaseFraction.x);
    fingerprint.real(audit.settings.gridPhaseFraction.y);
    fingerprint.real(audit.settings.gridPhaseFraction.z);
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(audit.samples.size()));
    for (const auto& sample : audit.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.sampleIndex));
        fingerprint.integer(static_cast<std::uint64_t>(sample.cellCounts.x));
        fingerprint.integer(static_cast<std::uint64_t>(sample.cellCounts.y));
        fingerprint.integer(static_cast<std::uint64_t>(sample.cellCounts.z));
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.gridCellCount));
        fingerprint.real(sample.gridLowerMeters.x);
        fingerprint.real(sample.gridLowerMeters.y);
        fingerprint.real(sample.gridLowerMeters.z);
        fingerprint.real(sample.cellSpacingMeters.x);
        fingerprint.real(sample.cellSpacingMeters.y);
        fingerprint.real(sample.cellSpacingMeters.z);
        fingerprint.real(sample.intakeAreaSquareMeters);
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.controlVolumeCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.fullTraceCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            sample.reducedTraceCount));
        fingerprint.real(sample.graphConductanceMeters);
        fingerprint.real(sample.shadowConductanceMeters);
        fingerprint.real(sample.graphToShadowConductanceRatio);
        fingerprint.real(sample.normalizedGraphConductance);
        fingerprint.real(sample.normalizedShadowConductance);
        fingerprint.integer(sample.response.fingerprint);
    }
    return fingerprint.value();
}

ScenePressureCellOperatorRefinementSample buildSample(
    const std::size_t sampleIndex,
    const fluid::GridCellCounts counts,
    const ScenePressureCellOperatorRefinementAuditSettings& settings,
    const ScenePressureCellOperatorRefinementAuditLimits& limits) {
    const Scene scene = makeScenePressureCellRefinementGeometry();
    const auto surface = assembleSceneFluidSurface(scene);
    const auto assembly = makeScenePressureCellAssembly(scene);
    Structure structure(assembly.definition);
    const SceneFluidSurfaceTransfer transfer(
        surface.definition, assembly.mappings, structure);
    const auto state = captureSceneFluidSurfaceState(
        surface.definition, assembly.mappings, structure);
    const double cellSpacingMeters = 4.0
        / static_cast<double>(counts.x);
    const fluid::Vector3 gridLowerMeters{
        settings.gridPhaseFraction.x * cellSpacingMeters,
        settings.gridPhaseFraction.y * cellSpacingMeters,
        settings.gridPhaseFraction.z * cellSpacingMeters,
    };
    const auto grid = makeScenePressureCellGrid(counts, gridLowerMeters);
    const auto epoch = buildSceneFluidGridEpoch(
        surface.definition, state, grid, transfer);
    const auto caps = buildSceneFluidOpeningCaps(
        surface.definition, state);
    const auto openingQuadrature = buildSceneFluidOpeningQuadrature(
        surface.definition, state, caps);
    const auto openingPatches = buildSceneFluidOpeningGridPatches(
        surface.definition, state, caps, openingQuadrature, grid);
    const auto openingFaceCrossings =
        buildSceneFluidOpeningFaceCrossings(
            surface.definition, state, caps, openingQuadrature,
            openingPatches, grid);
    const auto cappedFacePartitions =
        buildSceneFluidCappedFacePartitions(
            surface.definition, state, grid, transfer, epoch, caps,
            openingQuadrature, openingPatches, openingFaceCrossings);
    const auto volumes = buildSceneFluidCellVolumes(
        surface.definition, state, grid, transfer, epoch);
    const auto connectivity = buildSceneFluidRegionConnectivity(
        surface.definition);
    const auto pressureVolumes = buildSceneFluidPressureControlVolumes(
        surface.definition, volumes, connectivity);
    const auto faceLinks = buildSceneFluidPressureFaceLinks(
        surface.definition, state, grid, transfer, epoch, caps,
        openingQuadrature, openingPatches, openingFaceCrossings,
        cappedFacePartitions, volumes, connectivity, pressureVolumes);
    const auto graphOperator = buildSceneFluidPressureOperator(
        surface.definition, state, grid, transfer, epoch, caps,
        openingQuadrature, openingPatches, volumes, connectivity,
        pressureVolumes, faceLinks);
    const auto controlCells = buildSceneFluidMimeticControlCells(
        surface.definition, state, grid, epoch, caps, openingQuadrature,
        openingPatches, pressureVolumes, faceLinks);
    const auto fullTraceSystem = buildSceneFluidMimeticTraceSystem(
        controlCells);
    const auto condensedTraceSystem =
        buildSceneFluidMimeticCondensedTraceSystem(fullTraceSystem);
    auto response = auditSceneFluidPressureOperatorResponses(
        graphOperator, controlCells, fullTraceSystem,
        condensedTraceSystem, {}, settings.response, limits.response);
    const auto& regionMode = response.modes.back();
    if (regionMode.kind
            != SceneFluidPressureOperatorResponseModeKind::RegionContrast
        || !regionMode.hasTwoTerminalConductance) {
        throw std::runtime_error(
            "scene pressure-cell refinement has no two-terminal response");
    }

    ScenePressureCellOperatorRefinementSample result;
    result.sampleIndex = sampleIndex;
    result.cellCounts = counts;
    result.gridCellCount = gridCellCount(counts);
    result.gridLowerMeters = grid.lowerMeters();
    result.cellSpacingMeters = grid.cellSpacingMeters();
    result.intakeAreaSquareMeters = openingQuadrature.totalAreaSquareMeters;
    result.controlVolumeCount = controlCells.controlCells.size();
    result.fullTraceCount = fullTraceSystem.traces.size();
    result.reducedTraceCount = condensedTraceSystem.traces.size();
    result.graphConductanceMeters =
        regionMode.graphTwoTerminalConductanceMeters;
    result.shadowConductanceMeters =
        regionMode.shadowTwoTerminalConductanceMeters;
    result.graphToShadowConductanceRatio =
        regionMode.graphToShadowTwoTerminalConductanceRatio;
    result.normalizedGraphConductance = result.graphConductanceMeters
        * result.cellSpacingMeters.x / result.intakeAreaSquareMeters;
    result.normalizedShadowConductance = result.shadowConductanceMeters
        * result.cellSpacingMeters.x / result.intakeAreaSquareMeters;
    result.response = std::move(response);
    return result;
}

void validatePhaseAuditSettings(
    const ScenePressureCellOperatorPhaseAuditSettings& settings) {
    const auto counts = settings.cellCounts;
    if (counts.x < 2 || counts.x != counts.y || counts.x != counts.z) {
        throw std::invalid_argument(
            "scene pressure-cell phase audit requires an isotropic grid");
    }
    static_cast<void>(gridCellCount(counts));
}

bool validIncompleteFaceOwnership(
    const SceneFluidPressureIncompleteFaceOwnership& rejection) {
    std::size_t resolved = 0;
    std::size_t unresolved = 0;
    return checkedAdd(
               rejection.resolvedFullFaceCount,
               rejection.resolvedPartitionFaceCount, resolved)
        && checkedAdd(
               resolved, rejection.resolvedOpeningFaceCount, resolved)
        && checkedAdd(
               rejection.unresolvedActiveFaceCount,
               rejection.unresolvedCappedFaceCount, unresolved)
        && checkedAdd(
               unresolved, rejection.unresolvedAmbiguousFaceCount,
               unresolved)
        && checkedAdd(
               unresolved, rejection.unresolvedOpeningFaceCount,
               unresolved)
        && checkedAdd(
               unresolved,
               rejection.unresolvedEmbeddedOpeningPatchCount,
               unresolved)
        && rejection.faceCount != 0
        && resolved <= rejection.faceCount
        && rejection.unresolvedActiveFaceCount <= rejection.faceCount
        && rejection.unresolvedCappedFaceCount <= rejection.faceCount
        && rejection.unresolvedAmbiguousFaceCount <= rejection.faceCount
        && rejection.unresolvedOpeningFaceCount <= rejection.faceCount
        && (unresolved != 0 || resolved != rejection.faceCount);
}

std::size_t phaseAuditStorageBytes(
    const ScenePressureCellOperatorPhaseAudit& audit) {
    if (audit.samples.size() > std::numeric_limits<std::size_t>::max()
            / sizeof(ScenePressureCellOperatorPhaseSample)) {
        throw std::length_error(
            "scene pressure-cell phase audit storage overflows");
    }
    std::size_t result = audit.samples.size()
        * sizeof(ScenePressureCellOperatorPhaseSample);
    for (const auto& sample : audit.samples) {
        if (!sample.acceptedAudit.has_value()) {
            continue;
        }
        if (sample.acceptedAudit->ownedStorageBytes
            > std::numeric_limits<std::size_t>::max() - result) {
            throw std::length_error(
                "scene pressure-cell phase audit storage overflows");
        }
        result += sample.acceptedAudit->ownedStorageBytes;
    }
    return result;
}

ScenePressureCellOperatorPhaseStatistics phaseStatistics(
    const ScenePressureCellOperatorPhaseAudit& audit) {
    ScenePressureCellOperatorPhaseStatistics result;
    bool firstAccepted = true;
    double graphSum = 0.0;
    double shadowSum = 0.0;
    for (const auto& phaseSample : audit.samples) {
        if (phaseSample.status
                == ScenePressureCellOperatorPhaseSampleStatus::
                    RejectedIncompleteFaceOwnership) {
            ++result.rejectedIncompleteFaceOwnershipSampleCount;
            continue;
        }
        if (!phaseSample.acceptedAudit.has_value()) {
            throw std::invalid_argument(
                "scene pressure-cell phase sample has no accepted audit");
        }
        ++result.acceptedSampleCount;
        const auto& sample = phaseSample.acceptedAudit->samples.front();
        if (firstAccepted) {
            result.minimumNormalizedGraphConductance =
                sample.normalizedGraphConductance;
            result.maximumNormalizedGraphConductance =
                sample.normalizedGraphConductance;
            result.minimumNormalizedShadowConductance =
                sample.normalizedShadowConductance;
            result.maximumNormalizedShadowConductance =
                sample.normalizedShadowConductance;
            firstAccepted = false;
        } else {
            result.minimumNormalizedGraphConductance = std::min(
                result.minimumNormalizedGraphConductance,
                sample.normalizedGraphConductance);
            result.maximumNormalizedGraphConductance = std::max(
                result.maximumNormalizedGraphConductance,
                sample.normalizedGraphConductance);
            result.minimumNormalizedShadowConductance = std::min(
                result.minimumNormalizedShadowConductance,
                sample.normalizedShadowConductance);
            result.maximumNormalizedShadowConductance = std::max(
                result.maximumNormalizedShadowConductance,
                sample.normalizedShadowConductance);
        }
        graphSum += sample.normalizedGraphConductance;
        shadowSum += sample.normalizedShadowConductance;
    }
    if (result.acceptedSampleCount == 0) {
        return result;
    }
    const double acceptedCount =
        static_cast<double>(result.acceptedSampleCount);
    result.meanNormalizedGraphConductance = graphSum / acceptedCount;
    result.meanNormalizedShadowConductance = shadowSum / acceptedCount;
    double graphSquaredDeviation = 0.0;
    double shadowSquaredDeviation = 0.0;
    for (const auto& phaseSample : audit.samples) {
        if (!phaseSample.acceptedAudit.has_value()) {
            continue;
        }
        const auto& sample = phaseSample.acceptedAudit->samples.front();
        const double graphDeviation = sample.normalizedGraphConductance
            - result.meanNormalizedGraphConductance;
        const double shadowDeviation = sample.normalizedShadowConductance
            - result.meanNormalizedShadowConductance;
        graphSquaredDeviation += graphDeviation * graphDeviation;
        shadowSquaredDeviation += shadowDeviation * shadowDeviation;
    }
    result.graphCoefficientOfVariation = std::sqrt(
        graphSquaredDeviation / acceptedCount)
        / result.meanNormalizedGraphConductance;
    result.shadowCoefficientOfVariation = std::sqrt(
        shadowSquaredDeviation / acceptedCount)
        / result.meanNormalizedShadowConductance;
    return result;
}

std::uint64_t phaseAuditFingerprint(
    const ScenePressureCellOperatorPhaseAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.structureDefinitionFingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.settings.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.settings.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.settings.cellCounts.z));
    fingerprintResponseSettings(fingerprint, audit.settings.response);
    fingerprint.integer(static_cast<std::uint64_t>(audit.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.statistics.acceptedSampleCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.statistics.rejectedIncompleteFaceOwnershipSampleCount));
    fingerprint.real(audit.statistics.minimumNormalizedGraphConductance);
    fingerprint.real(audit.statistics.maximumNormalizedGraphConductance);
    fingerprint.real(audit.statistics.meanNormalizedGraphConductance);
    fingerprint.real(audit.statistics.graphCoefficientOfVariation);
    fingerprint.real(audit.statistics.minimumNormalizedShadowConductance);
    fingerprint.real(audit.statistics.maximumNormalizedShadowConductance);
    fingerprint.real(audit.statistics.meanNormalizedShadowConductance);
    fingerprint.real(audit.statistics.shadowCoefficientOfVariation);
    fingerprint.integer(static_cast<std::uint64_t>(audit.samples.size()));
    for (const auto& sample : audit.samples) {
        fingerprint.integer(static_cast<std::uint64_t>(sample.sampleIndex));
        fingerprint.real(sample.gridPhaseFraction.x);
        fingerprint.real(sample.gridPhaseFraction.y);
        fingerprint.real(sample.gridPhaseFraction.z);
        fingerprint.integer(static_cast<std::uint8_t>(sample.status));
        fingerprint.integer(sample.faceOwnershipRejection.has_value()
            ? 1ULL : 0ULL);
        if (sample.faceOwnershipRejection.has_value()) {
            const auto& rejection = *sample.faceOwnershipRejection;
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.faceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.resolvedFullFaceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.resolvedPartitionFaceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.resolvedOpeningFaceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.unresolvedActiveFaceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.unresolvedCappedFaceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.unresolvedAmbiguousFaceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.unresolvedOpeningFaceCount));
            fingerprint.integer(static_cast<std::uint64_t>(
                rejection.unresolvedEmbeddedOpeningPatchCount));
        }
        fingerprint.integer(sample.acceptedAudit.has_value()
            ? sample.acceptedAudit->fingerprint : 0ULL);
    }
    return fingerprint.value();
}

} // namespace

ScenePressureCellOperatorRefinementAudit
auditScenePressureCellOperatorRefinement(
    const std::span<const fluid::GridCellCounts> resolutions,
    const ScenePressureCellOperatorRefinementAuditSettings& settings,
    const ScenePressureCellOperatorRefinementAuditLimits& limits) {
    if (resolutions.empty() || resolutions.size() > limits.maximumSamples) {
        throw std::length_error(
            "scene pressure-cell refinement sample limit exceeded");
    }
    validateSettings(settings);

    ScenePressureCellOperatorRefinementAudit result;
    result.settings = settings;
    result.samples.reserve(resolutions.size());
    std::size_t previousCount = 0;
    for (std::size_t index = 0; index < resolutions.size(); ++index) {
        const auto counts = resolutions[index];
        const std::size_t count = gridCellCount(counts);
        if (counts.x < 2 || counts.x != counts.y || counts.x != counts.z
            || count <= previousCount) {
            throw std::invalid_argument(
                "scene pressure-cell refinement resolutions must be "
                "strictly increasing isotropic grids");
        }
        if (count > limits.maximumGridCellsPerSample) {
            throw std::length_error(
                "scene pressure-cell refinement grid limit exceeded");
        }
        auto sample = buildSample(index, counts, settings, limits);
        if (index == 0) {
            result.structureDefinitionFingerprint =
                sample.response.structureDefinitionFingerprint;
        } else if (sample.response.structureDefinitionFingerprint
                   != result.structureDefinitionFingerprint) {
            throw std::logic_error(
                "scene pressure-cell refinement structure changed");
        }
        result.samples.push_back(std::move(sample));
        result.ownedStorageBytes = storageBytes(result);
        if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
            throw std::length_error(
                "scene pressure-cell refinement byte limit exceeded");
        }
        previousCount = count;
    }
    result.fingerprint = productFingerprint(result);
    validateScenePressureCellOperatorRefinementAuditIntegrity(result);
    return result;
}

void validateScenePressureCellOperatorRefinementAuditIntegrity(
    const ScenePressureCellOperatorRefinementAudit& audit) {
    if (audit.version
            != scenePressureCellOperatorRefinementAuditVersion
        || audit.fingerprint == 0
        || audit.structureDefinitionFingerprint == 0
        || audit.samples.empty()
        || audit.ownedStorageBytes != storageBytes(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell refinement integrity is invalid");
    }
    validateSettings(audit.settings);
    std::size_t previousCount = 0;
    double intakeAreaSquareMeters = 0.0;
    for (std::size_t index = 0; index < audit.samples.size(); ++index) {
        const auto& sample = audit.samples[index];
        validateSceneFluidPressureOperatorResponseAuditIntegrity(
            sample.response);
        const std::size_t count = gridCellCount(sample.cellCounts);
        const double expectedSpacingMeters = 4.0
            / static_cast<double>(sample.cellCounts.x);
        const fluid::Vector3 expectedLowerMeters{
            audit.settings.gridPhaseFraction.x * expectedSpacingMeters,
            audit.settings.gridPhaseFraction.y * expectedSpacingMeters,
            audit.settings.gridPhaseFraction.z * expectedSpacingMeters,
        };
        const auto& regionMode = sample.response.modes.back();
        if (sample.sampleIndex != index
            || sample.cellCounts.x < 2
            || sample.cellCounts.x != sample.cellCounts.y
            || sample.cellCounts.x != sample.cellCounts.z
            || count <= previousCount || sample.gridCellCount != count
            || sample.gridLowerMeters != expectedLowerMeters
            || !std::isfinite(sample.cellSpacingMeters.x)
            || !std::isfinite(sample.cellSpacingMeters.y)
            || !std::isfinite(sample.cellSpacingMeters.z)
            || !(sample.cellSpacingMeters.x > 0.0)
            || sample.cellSpacingMeters.x != expectedSpacingMeters
            || sample.cellSpacingMeters.x != sample.cellSpacingMeters.y
            || sample.cellSpacingMeters.x != sample.cellSpacingMeters.z
            || !(sample.intakeAreaSquareMeters > 0.0)
            || (index != 0
                && sample.intakeAreaSquareMeters
                    != intakeAreaSquareMeters)
            || sample.controlVolumeCount
                != sample.response.controlVolumeCount
            || sample.fullTraceCount == 0 || sample.reducedTraceCount == 0
            || sample.reducedTraceCount > sample.fullTraceCount
            || sample.response.structureDefinitionFingerprint
                != audit.structureDefinitionFingerprint
            || sample.response.settings != audit.settings.response
            || regionMode.kind
                != SceneFluidPressureOperatorResponseModeKind::RegionContrast
            || !regionMode.hasTwoTerminalConductance
            || sample.graphConductanceMeters
                != regionMode.graphTwoTerminalConductanceMeters
            || sample.shadowConductanceMeters
                != regionMode.shadowTwoTerminalConductanceMeters
            || sample.graphToShadowConductanceRatio
                != regionMode.graphToShadowTwoTerminalConductanceRatio
            || sample.normalizedGraphConductance
                != sample.graphConductanceMeters
                    * sample.cellSpacingMeters.x
                    / sample.intakeAreaSquareMeters
            || sample.normalizedShadowConductance
                != sample.shadowConductanceMeters
                    * sample.cellSpacingMeters.x
                    / sample.intakeAreaSquareMeters
            || !std::isfinite(sample.normalizedGraphConductance)
            || !std::isfinite(sample.normalizedShadowConductance)
            || !(sample.normalizedGraphConductance > 0.0)
            || !(sample.normalizedShadowConductance > 0.0)) {
            throw std::invalid_argument(
                "scene pressure-cell refinement sample is invalid");
        }
        if (index == 0) {
            intakeAreaSquareMeters = sample.intakeAreaSquareMeters;
        }
        previousCount = count;
    }
    if (audit.fingerprint != productFingerprint(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell refinement fingerprint is invalid");
    }
}

ScenePressureCellOperatorPhaseAudit
auditScenePressureCellOperatorGridPhases(
    const std::span<const fluid::Vector3> gridPhaseFractions,
    const ScenePressureCellOperatorPhaseAuditSettings& settings,
    const ScenePressureCellOperatorPhaseAuditLimits& limits) {
    if (gridPhaseFractions.empty()
        || gridPhaseFractions.size() > limits.maximumSamples) {
        throw std::length_error(
            "scene pressure-cell phase audit sample limit exceeded");
    }
    if (gridPhaseFractions.size()
        > limits.maximumOwnedBytes
            / sizeof(ScenePressureCellOperatorPhaseSample)) {
        throw std::length_error(
            "scene pressure-cell phase audit byte limit exceeded");
    }
    validatePhaseAuditSettings(settings);
    if (gridCellCount(settings.cellCounts)
        > limits.maximumGridCellsPerSample) {
        throw std::length_error(
            "scene pressure-cell phase audit grid limit exceeded");
    }
    for (std::size_t index = 0;
         index < gridPhaseFractions.size(); ++index) {
        const auto phase = gridPhaseFractions[index];
        if (!validPhaseComponent(phase.x)
            || !validPhaseComponent(phase.y)
            || !validPhaseComponent(phase.z)) {
            throw std::invalid_argument(
                "scene pressure-cell phase audit phase is invalid");
        }
        for (std::size_t earlier = 0; earlier < index; ++earlier) {
            if (gridPhaseFractions[earlier] == phase) {
                throw std::invalid_argument(
                    "scene pressure-cell phase audit phases must be unique");
            }
        }
    }

    ScenePressureCellOperatorPhaseAudit result;
    result.settings = settings;
    const Scene scene = makeScenePressureCellRefinementGeometry();
    const auto assembly = makeScenePressureCellAssembly(scene);
    const Structure structure(assembly.definition);
    result.structureDefinitionFingerprint = structure.definitionFingerprint();
    result.samples.reserve(gridPhaseFractions.size());

    for (std::size_t index = 0;
         index < gridPhaseFractions.size(); ++index) {
        const auto phase = gridPhaseFractions[index];
        ScenePressureCellOperatorPhaseSample sample;
        sample.sampleIndex = index;
        sample.gridPhaseFraction = phase;
        ScenePressureCellOperatorRefinementAuditSettings nestedSettings;
        nestedSettings.response = settings.response;
        nestedSettings.gridPhaseFraction = phase;
        ScenePressureCellOperatorRefinementAuditLimits nestedLimits;
        nestedLimits.maximumSamples = 1;
        nestedLimits.maximumGridCellsPerSample =
            limits.maximumGridCellsPerSample;
        nestedLimits.maximumOwnedBytes = limits.maximumOwnedBytes;
        nestedLimits.response = limits.response;
        const auto counts = settings.cellCounts;
        try {
            sample.acceptedAudit =
                auditScenePressureCellOperatorRefinement(
                    std::span(&counts, 1), nestedSettings, nestedLimits);
        } catch (
            const SceneFluidPressureIncompleteFaceOwnershipError& error) {
            sample.status =
                ScenePressureCellOperatorPhaseSampleStatus::
                    RejectedIncompleteFaceOwnership;
            sample.faceOwnershipRejection = error.diagnostics();
        }
        if (sample.acceptedAudit.has_value()
            && sample.acceptedAudit->structureDefinitionFingerprint
                != result.structureDefinitionFingerprint) {
            throw std::logic_error(
                "scene pressure-cell phase audit structure changed");
        }
        result.samples.push_back(std::move(sample));
        result.ownedStorageBytes = phaseAuditStorageBytes(result);
        if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
            throw std::length_error(
                "scene pressure-cell phase audit byte limit exceeded");
        }
    }
    result.statistics = phaseStatistics(result);
    result.fingerprint = phaseAuditFingerprint(result);
    validateScenePressureCellOperatorPhaseAuditIntegrity(result);
    return result;
}

void validateScenePressureCellOperatorPhaseAuditIntegrity(
    const ScenePressureCellOperatorPhaseAudit& audit) {
    if (audit.version != scenePressureCellOperatorPhaseAuditVersion
        || audit.fingerprint == 0
        || audit.structureDefinitionFingerprint == 0
        || audit.samples.empty()
        || audit.ownedStorageBytes != phaseAuditStorageBytes(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell phase audit integrity is invalid");
    }
    validatePhaseAuditSettings(audit.settings);
    for (std::size_t index = 0; index < audit.samples.size(); ++index) {
        const auto& sample = audit.samples[index];
        if (sample.sampleIndex != index
            || !validPhaseComponent(sample.gridPhaseFraction.x)
            || !validPhaseComponent(sample.gridPhaseFraction.y)
            || !validPhaseComponent(sample.gridPhaseFraction.z)) {
            throw std::invalid_argument(
                "scene pressure-cell phase audit sample is invalid");
        }
        for (std::size_t earlier = 0; earlier < index; ++earlier) {
            if (audit.samples[earlier].gridPhaseFraction
                    == sample.gridPhaseFraction) {
                throw std::invalid_argument(
                    "scene pressure-cell phase audit phases are not unique");
            }
        }
        if (sample.status
                == ScenePressureCellOperatorPhaseSampleStatus::Accepted) {
            if (sample.faceOwnershipRejection.has_value()
                || !sample.acceptedAudit.has_value()) {
                throw std::invalid_argument(
                    "scene pressure-cell accepted phase is invalid");
            }
            validateScenePressureCellOperatorRefinementAuditIntegrity(
                *sample.acceptedAudit);
            if (sample.acceptedAudit->structureDefinitionFingerprint
                    != audit.structureDefinitionFingerprint
                || sample.acceptedAudit->settings.response
                    != audit.settings.response
                || sample.acceptedAudit->settings.gridPhaseFraction
                    != sample.gridPhaseFraction
                || sample.acceptedAudit->samples.size() != 1
                || sample.acceptedAudit->samples.front().cellCounts
                    != audit.settings.cellCounts) {
                throw std::invalid_argument(
                    "scene pressure-cell accepted phase audit is mismatched");
            }
        } else if (sample.status
                       == ScenePressureCellOperatorPhaseSampleStatus::
                           RejectedIncompleteFaceOwnership) {
            if (sample.acceptedAudit.has_value()
                || !sample.faceOwnershipRejection.has_value()
                || !validIncompleteFaceOwnership(
                    *sample.faceOwnershipRejection)) {
                throw std::invalid_argument(
                    "scene pressure-cell rejected phase is invalid");
            }
        } else {
            throw std::invalid_argument(
                "scene pressure-cell phase audit status is invalid");
        }
    }
    const auto expectedStatistics = phaseStatistics(audit);
    if (audit.statistics != expectedStatistics
        || (audit.statistics.acceptedSampleCount != 0
            && (!std::isfinite(
                    audit.statistics.meanNormalizedGraphConductance)
                || !std::isfinite(
                    audit.statistics.graphCoefficientOfVariation)
                || !std::isfinite(
                    audit.statistics.meanNormalizedShadowConductance)
                || !std::isfinite(
                    audit.statistics.shadowCoefficientOfVariation)
                || !(audit.statistics.minimumNormalizedGraphConductance > 0.0)
                || !(audit.statistics.minimumNormalizedShadowConductance > 0.0)))) {
        throw std::invalid_argument(
            "scene pressure-cell phase audit statistics are invalid");
    }
    if (audit.fingerprint != phaseAuditFingerprint(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell phase audit fingerprint is invalid");
    }
}

} // namespace simwing::fsi
