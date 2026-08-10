#include "scene_pressure_cell_mimetic_conductance_phase_refinement_audit.h"

#include "scene_fluid_capped_face_partition.h"
#include "scene_fluid_cell_volume.h"
#include "scene_fluid_mimetic_condensed_trace_system.h"
#include "scene_fluid_mimetic_control_cell.h"
#include "scene_fluid_mimetic_trace_system.h"
#include "scene_fluid_opening_cap.h"
#include "scene_fluid_opening_face_crossing.h"
#include "scene_fluid_opening_patch.h"
#include "scene_fluid_opening_quadrature.h"
#include "scene_fluid_pressure_control_volume.h"
#include "scene_fluid_pressure_face_link.h"
#include "scene_fluid_region_connectivity.h"
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

class CompensatedSum final {
public:
    void add(const double value) noexcept {
        const double next = sum_ + value;
        if (std::abs(sum_) >= std::abs(value)) {
            correction_ += (sum_ - next) + value;
        } else {
            correction_ += (value - next) + sum_;
        }
        sum_ = next;
    }

    [[nodiscard]] double value() const noexcept {
        return sum_ + correction_;
    }

private:
    double sum_ = 0.0;
    double correction_ = 0.0;
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
            "scene pressure-cell mimetic conductance grid count overflows");
    }
    return xyz;
}

bool validPhaseComponent(const double value) {
    return std::isfinite(value) && value >= -0.5 && value < 0.5
        && (value != 0.0 || !std::signbit(value));
}

bool finite(const fluid::Vector3 value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

Scene translatedRefinementScene(const fluid::Vector3 translation) {
    Scene scene = makeScenePressureCellRefinementGeometry();
    if (translation == fluid::Vector3{}) {
        return scene;
    }
    for (auto& vertex : scene.vertices) {
        vertex.positionMeters.x += translation.x;
        vertex.positionMeters.y += translation.y;
        vertex.positionMeters.z += translation.z;
    }
    return scene;
}

fluid::Vector3 translatedGridLower(const fluid::Vector3 translation,
                                   const fluid::Vector3 phase,
                                   const double spacing) {
    return {
        translation.x + phase.x * spacing,
        translation.y + phase.y * spacing,
        translation.z + phase.z * spacing,
    };
}

bool validLinearConsistencyFailure(
    const fluid::MimeticLocalCellLinearConsistencyFailure& failure) {
    return failure.halfFaceCount >= 4
        && std::isfinite(failure.volumeCubicMeters)
        && failure.volumeCubicMeters > 0.0
        && std::isfinite(failure.summedAreaSquareMeters)
        && failure.summedAreaSquareMeters > 0.0
        && std::isfinite(failure.minimumFaceAreaSquareMeters)
        && failure.minimumFaceAreaSquareMeters > 0.0
        && std::isfinite(failure.maximumFaceAreaSquareMeters)
        && failure.maximumFaceAreaSquareMeters
            >= failure.minimumFaceAreaSquareMeters
        && std::isfinite(failure.maximumAreaClosureErrorSquareMeters)
        && failure.maximumAreaClosureErrorSquareMeters >= 0.0
        && std::isfinite(
            failure.maximumDivergenceTheoremErrorCubicMeters)
        && failure.maximumDivergenceTheoremErrorCubicMeters >= 0.0
        && std::isfinite(failure.consistencyGeometryConditionEstimate)
        && failure.consistencyGeometryConditionEstimate > 0.0
        && std::isfinite(failure.consistencyGramConditionEstimate)
        && failure.consistencyGramConditionEstimate > 0.0
        && std::isfinite(
            failure.stabilizationScaleInverseCubicMeters)
        && failure.stabilizationScaleInverseCubicMeters > 0.0
        && std::isfinite(failure.maximumAlgebraicConsistencyError)
        && std::isfinite(failure.algebraicConsistencyTolerance)
        && failure.algebraicConsistencyTolerance >= 0.0
        && failure.maximumAlgebraicConsistencyError
            > failure.algebraicConsistencyTolerance;
}

void validatePhases(const std::span<const fluid::Vector3> phases) {
    for (std::size_t index = 0; index < phases.size(); ++index) {
        const auto phase = phases[index];
        if (!validPhaseComponent(phase.x)
            || !validPhaseComponent(phase.y)
            || !validPhaseComponent(phase.z)) {
            throw std::invalid_argument(
                "scene pressure-cell mimetic conductance phase is invalid");
        }
        for (std::size_t earlier = 0; earlier < index; ++earlier) {
            if (phases[earlier] == phase) {
                throw std::invalid_argument(
                    "scene pressure-cell mimetic conductance phases must be unique");
            }
        }
    }
}

void validateSettings(
    const ScenePressureCellMimeticConductancePhaseRefinementAuditSettings&
        settings) {
    const auto& conductance = settings.conductance;
    const auto& solve = conductance.solve;
    if (!finite(settings.geometryTranslationMeters)
        || !std::isfinite(
            conductance.terminalIntegratedTransferPascalsMeters)
        || !(conductance.terminalIntegratedTransferPascalsMeters > 0.0)
        || !std::isfinite(
            conductance.absoluteOpeningPairAreaToleranceSquareMeters)
        || conductance.absoluteOpeningPairAreaToleranceSquareMeters < 0.0
        || !std::isfinite(
            conductance.relativeOpeningPairAreaTolerance)
        || conductance.relativeOpeningPairAreaTolerance < 0.0
        || !std::isfinite(
            solve.absoluteResidualTolerancePascalsMeters)
        || solve.absoluteResidualTolerancePascalsMeters < 0.0
        || !std::isfinite(solve.relativeResidualTolerance)
        || solve.relativeResidualTolerance < 0.0
        || !std::isfinite(
            solve.absoluteComponentCompatibilityTolerancePascalsMeters)
        || solve.absoluteComponentCompatibilityTolerancePascalsMeters < 0.0
        || solve.maximumIterations == 0) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic conductance settings are invalid");
    }
}

std::size_t storageBytes(
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& audit) {
    std::size_t phaseBytes = 0;
    std::size_t levelBytes = 0;
    std::size_t result = 0;
    if (!checkedMultiply(audit.gridPhaseFractions.size(),
                         sizeof(fluid::Vector3), phaseBytes)
        || !checkedMultiply(
            audit.levels.size(),
            sizeof(ScenePressureCellMimeticConductanceRefinementLevel),
            levelBytes)
        || !checkedAdd(phaseBytes, levelBytes, result)) {
        throw std::length_error(
            "scene pressure-cell mimetic conductance storage overflows");
    }
    for (const auto& level : audit.levels) {
        std::size_t sampleBytes = 0;
        if (!checkedMultiply(
                level.samples.size(),
                sizeof(ScenePressureCellMimeticConductancePhaseSample),
                sampleBytes)
            || !checkedAdd(result, sampleBytes, result)) {
            throw std::length_error(
                "scene pressure-cell mimetic conductance storage overflows");
        }
        for (const auto& sample : level.samples) {
            if (sample.conductanceAudit.has_value()
                && !checkedAdd(
                    result,
                    sample.conductanceAudit->ownedStorageBytes,
                    result)) {
                throw std::length_error(
                    "scene pressure-cell mimetic conductance storage overflows");
            }
        }
    }
    return result;
}

void fingerprintSettings(
    Fingerprint& fingerprint,
    const ScenePressureCellMimeticConductancePhaseRefinementAuditSettings&
        settings) {
    const auto& conductance = settings.conductance;
    fingerprint.real(settings.geometryTranslationMeters.x);
    fingerprint.real(settings.geometryTranslationMeters.y);
    fingerprint.real(settings.geometryTranslationMeters.z);
    fingerprint.real(
        conductance.terminalIntegratedTransferPascalsMeters);
    fingerprint.real(
        conductance.absoluteOpeningPairAreaToleranceSquareMeters);
    fingerprint.real(conductance.relativeOpeningPairAreaTolerance);
    fingerprint.real(
        conductance.solve.absoluteResidualTolerancePascalsMeters);
    fingerprint.real(conductance.solve.relativeResidualTolerance);
    fingerprint.real(conductance.solve
        .absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        conductance.solve.maximumIterations));
}

std::uint64_t productFingerprint(
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.structureDefinitionFingerprint);
    fingerprintSettings(fingerprint, audit.settings);
    fingerprint.integer(static_cast<std::uint64_t>(audit.ownedStorageBytes));
    fingerprint.integer(static_cast<std::uint64_t>(
        audit.gridPhaseFractions.size()));
    for (const auto phase : audit.gridPhaseFractions) {
        fingerprint.real(phase.x);
        fingerprint.real(phase.y);
        fingerprint.real(phase.z);
    }
    fingerprint.integer(static_cast<std::uint64_t>(audit.levels.size()));
    for (const auto& level : audit.levels) {
        fingerprint.integer(static_cast<std::uint64_t>(level.levelIndex));
        fingerprint.integer(static_cast<std::uint64_t>(level.cellCounts.x));
        fingerprint.integer(static_cast<std::uint64_t>(level.cellCounts.y));
        fingerprint.integer(static_cast<std::uint64_t>(level.cellCounts.z));
        fingerprint.real(level.cellSpacingMeters.x);
        fingerprint.real(level.cellSpacingMeters.y);
        fingerprint.real(level.cellSpacingMeters.z);
        fingerprint.integer(static_cast<std::uint64_t>(
            level.acceptedSampleCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            level.rejectedLocalCellLinearConsistencySampleCount));
        fingerprint.real(level.minimumNormalizedConductance);
        fingerprint.real(level.maximumNormalizedConductance);
        fingerprint.real(level.meanNormalizedConductance);
        fingerprint.real(
            level.normalizedConductanceCoefficientOfVariation);
        fingerprint.integer(static_cast<std::uint64_t>(
            level.samples.size()));
        for (const auto& sample : level.samples) {
            fingerprint.integer(static_cast<std::uint64_t>(
                sample.sampleIndex));
            fingerprint.integer(static_cast<std::uint64_t>(
                sample.phaseIndex));
            fingerprint.real(sample.gridPhaseFraction.x);
            fingerprint.real(sample.gridPhaseFraction.y);
            fingerprint.real(sample.gridPhaseFraction.z);
            fingerprint.integer(static_cast<std::uint8_t>(sample.status));
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
            fingerprint.integer(static_cast<std::uint64_t>(
                sample.openingTraceCount));
            fingerprint.real(sample.conductanceMeters);
            fingerprint.real(sample.normalizedConductance);
            fingerprint.integer(
                sample.localCellLinearConsistencyRejection.has_value()
                ? 1U : 0U);
            if (sample.localCellLinearConsistencyRejection.has_value()) {
                const auto& rejection =
                    *sample.localCellLinearConsistencyRejection;
                fingerprint.integer(static_cast<std::uint64_t>(
                    rejection.controlCellIndex));
                fingerprint.integer(rejection.controlCellStableId);
                fingerprint.integer(static_cast<std::uint64_t>(
                    rejection.gridCellIndex));
                fingerprint.integer(rejection.regionId);
                fingerprint.integer(static_cast<std::uint64_t>(
                    rejection.localCell.halfFaceCount));
                fingerprint.real(rejection.localCell.volumeCubicMeters);
                fingerprint.real(
                    rejection.localCell.summedAreaSquareMeters);
                fingerprint.real(
                    rejection.localCell.minimumFaceAreaSquareMeters);
                fingerprint.real(
                    rejection.localCell.maximumFaceAreaSquareMeters);
                fingerprint.real(rejection.localCell
                    .maximumAreaClosureErrorSquareMeters);
                fingerprint.real(rejection.localCell
                    .maximumDivergenceTheoremErrorCubicMeters);
                fingerprint.real(rejection.localCell
                    .consistencyGeometryConditionEstimate);
                fingerprint.real(rejection.localCell
                    .consistencyGramConditionEstimate);
                fingerprint.real(rejection.localCell
                    .stabilizationScaleInverseCubicMeters);
                fingerprint.real(
                    rejection.localCell
                        .maximumAlgebraicConsistencyError);
                fingerprint.real(
                    rejection.localCell.algebraicConsistencyTolerance);
            }
            fingerprint.integer(sample.conductanceAudit.has_value()
                ? sample.conductanceAudit->fingerprint : 0ULL);
        }
    }
    return fingerprint.value();
}

void summarize(
    ScenePressureCellMimeticConductanceRefinementLevel& level) {
    if (level.samples.empty()) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic conductance level is empty");
    }
    level.minimumNormalizedConductance =
        std::numeric_limits<double>::infinity();
    level.maximumNormalizedConductance = 0.0;
    level.acceptedSampleCount = 0;
    level.rejectedLocalCellLinearConsistencySampleCount = 0;
    CompensatedSum sum;
    for (const auto& sample : level.samples) {
        if (sample.status
            == ScenePressureCellMimeticConductancePhaseSampleStatus::
                RejectedLocalCellLinearConsistency) {
            ++level.rejectedLocalCellLinearConsistencySampleCount;
            continue;
        }
        if (sample.status
            != ScenePressureCellMimeticConductancePhaseSampleStatus::
                Accepted) {
            throw std::invalid_argument(
                "scene pressure-cell mimetic conductance status is invalid");
        }
        ++level.acceptedSampleCount;
        level.minimumNormalizedConductance = std::min(
            level.minimumNormalizedConductance,
            sample.normalizedConductance);
        level.maximumNormalizedConductance = std::max(
            level.maximumNormalizedConductance,
            sample.normalizedConductance);
        sum.add(sample.normalizedConductance);
    }
    if (level.acceptedSampleCount == 0) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic conductance level has no accepted samples");
    }
    level.meanNormalizedConductance = sum.value()
        / static_cast<double>(level.acceptedSampleCount);
    CompensatedSum variance;
    for (const auto& sample : level.samples) {
        if (sample.status
            != ScenePressureCellMimeticConductancePhaseSampleStatus::
                Accepted) {
            continue;
        }
        const double difference = sample.normalizedConductance
            - level.meanNormalizedConductance;
        variance.add(difference * difference);
    }
    level.normalizedConductanceCoefficientOfVariation = std::sqrt(
        std::max(0.0, variance.value()
            / static_cast<double>(level.acceptedSampleCount)))
        / level.meanNormalizedConductance;
}

ScenePressureCellMimeticConductancePhaseSample buildSample(
    const std::size_t sampleIndex,
    const std::size_t phaseIndex,
    const fluid::GridCellCounts counts,
    const fluid::Vector3 phase,
    const ScenePressureCellMimeticConductancePhaseRefinementAuditSettings&
        settings,
    const ScenePressureCellMimeticConductancePhaseRefinementAuditLimits&
        limits) {
    const Scene scene = translatedRefinementScene(
        settings.geometryTranslationMeters);
    const auto surface = assembleSceneFluidSurface(scene);
    const auto assembly = makeScenePressureCellAssembly(scene);
    Structure structure(assembly.definition);
    const SceneFluidSurfaceTransfer transfer(
        surface.definition, assembly.mappings, structure);
    const auto state = captureSceneFluidSurfaceState(
        surface.definition, assembly.mappings, structure);
    const double spacing = 4.0 / static_cast<double>(counts.x);
    const fluid::Vector3 lower = translatedGridLower(
        settings.geometryTranslationMeters, phase, spacing);
    const auto grid = makeScenePressureCellGrid(counts, lower);
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
    const auto controlCells = buildSceneFluidMimeticControlCells(
        surface.definition, state, grid, epoch, caps, openingQuadrature,
        openingPatches, pressureVolumes, faceLinks);

    ScenePressureCellMimeticConductancePhaseSample result;
    result.sampleIndex = sampleIndex;
    result.phaseIndex = phaseIndex;
    result.gridPhaseFraction = phase;
    result.gridCellCount = gridCellCount(counts);
    result.gridLowerMeters = grid.lowerMeters();
    result.cellSpacingMeters = grid.cellSpacingMeters();
    result.intakeAreaSquareMeters =
        openingQuadrature.totalAreaSquareMeters;
    result.controlVolumeCount = controlCells.controlCells.size();
    std::size_t openingHalfFaceCount = 0;
    for (const auto& halfFace : controlCells.halfFaces) {
        if (halfFace.kind
                == SceneFluidMimeticHalfFaceKind::MaterialWall
            || halfFace.otherControlVolumeIndex
                >= controlCells.controlCells.size()) {
            continue;
        }
        if (controlCells.controlCells[halfFace.controlVolumeIndex].regionId
            != controlCells.controlCells[
                halfFace.otherControlVolumeIndex].regionId) {
            ++openingHalfFaceCount;
        }
    }
    if (openingHalfFaceCount == 0 || openingHalfFaceCount % 2 != 0) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic conductance opening traces are incomplete");
    }
    result.openingTraceCount = openingHalfFaceCount / 2;

    SceneFluidMimeticTraceSystem fullTraceSystem;
    try {
        fullTraceSystem = buildSceneFluidMimeticTraceSystem(controlCells);
    } catch (
        const SceneFluidMimeticTraceLocalCellLinearConsistencyError& error) {
        result.status =
            ScenePressureCellMimeticConductancePhaseSampleStatus::
                RejectedLocalCellLinearConsistency;
        result.localCellLinearConsistencyRejection = error.diagnostics();
        return result;
    }
    const auto condensedTraceSystem =
        buildSceneFluidMimeticCondensedTraceSystem(fullTraceSystem);
    auto conductance = auditSceneFluidMimeticRegionConductance(
        controlCells, fullTraceSystem, condensedTraceSystem,
        settings.conductance, limits.conductance);
    result.fullTraceCount = fullTraceSystem.traces.size();
    result.reducedTraceCount = condensedTraceSystem.traces.size();
    if (result.openingTraceCount != conductance.openings.size()) {
        throw std::logic_error(
            "scene pressure-cell mimetic conductance opening trace count changed");
    }
    result.conductanceMeters = conductance.conductanceMeters;
    result.normalizedConductance = result.conductanceMeters
        * result.cellSpacingMeters.x / result.intakeAreaSquareMeters;
    result.conductanceAudit = std::move(conductance);
    return result;
}

} // namespace

ScenePressureCellMimeticConductancePhaseRefinementAudit
auditScenePressureCellMimeticConductancePhaseRefinement(
    const std::span<const fluid::GridCellCounts> resolutions,
    const std::span<const fluid::Vector3> gridPhaseFractions,
    const ScenePressureCellMimeticConductancePhaseRefinementAuditSettings&
        settings,
    const ScenePressureCellMimeticConductancePhaseRefinementAuditLimits&
        limits) {
    if (resolutions.empty()
        || resolutions.size() > limits.maximumResolutionCount
        || gridPhaseFractions.empty()
        || gridPhaseFractions.size() > limits.maximumPhaseCount) {
        throw std::length_error(
            "scene pressure-cell mimetic conductance sample limit exceeded");
    }
    std::size_t aggregateSamples = 0;
    if (!checkedMultiply(resolutions.size(), gridPhaseFractions.size(),
                         aggregateSamples)
        || aggregateSamples > limits.maximumAggregateSamples) {
        throw std::length_error(
            "scene pressure-cell mimetic conductance aggregate limit exceeded");
    }
    validateSettings(settings);
    validatePhases(gridPhaseFractions);

    std::size_t phaseBytes = 0;
    std::size_t levelBytes = 0;
    std::size_t sampleBytes = 0;
    std::size_t preflightBytes = 0;
    if (!checkedMultiply(gridPhaseFractions.size(),
                         sizeof(fluid::Vector3), phaseBytes)
        || !checkedMultiply(
            resolutions.size(),
            sizeof(ScenePressureCellMimeticConductanceRefinementLevel),
            levelBytes)
        || !checkedMultiply(
            aggregateSamples,
            sizeof(ScenePressureCellMimeticConductancePhaseSample),
            sampleBytes)
        || !checkedAdd(phaseBytes, levelBytes, preflightBytes)
        || !checkedAdd(preflightBytes, sampleBytes, preflightBytes)
        || preflightBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene pressure-cell mimetic conductance byte limit exceeded");
    }

    ScenePressureCellMimeticConductancePhaseRefinementAudit result;
    result.settings = settings;
    const Scene scene = translatedRefinementScene(
        settings.geometryTranslationMeters);
    const auto assembly = makeScenePressureCellAssembly(scene);
    const Structure structure(assembly.definition);
    result.structureDefinitionFingerprint =
        structure.definitionFingerprint();
    result.gridPhaseFractions.assign(
        gridPhaseFractions.begin(), gridPhaseFractions.end());
    result.levels.reserve(resolutions.size());
    std::size_t previousCellCount = 0;
    std::size_t nextSampleIndex = 0;
    for (std::size_t levelIndex = 0;
         levelIndex < resolutions.size(); ++levelIndex) {
        const auto counts = resolutions[levelIndex];
        const std::size_t cellCount = gridCellCount(counts);
        if (counts.x < 2 || counts.x != counts.y || counts.x != counts.z
            || cellCount <= previousCellCount) {
            throw std::invalid_argument(
                "scene pressure-cell mimetic conductance resolutions must be strictly increasing isotropic grids");
        }
        if (cellCount > limits.maximumGridCellsPerSample) {
            throw std::length_error(
                "scene pressure-cell mimetic conductance grid limit exceeded");
        }

        ScenePressureCellMimeticConductanceRefinementLevel level;
        level.levelIndex = levelIndex;
        level.cellCounts = counts;
        const double spacing = 4.0 / static_cast<double>(counts.x);
        level.cellSpacingMeters = {spacing, spacing, spacing};
        level.samples.reserve(gridPhaseFractions.size());
        for (std::size_t phaseIndex = 0;
             phaseIndex < gridPhaseFractions.size(); ++phaseIndex) {
            auto sample = buildSample(
                nextSampleIndex, phaseIndex, counts,
                gridPhaseFractions[phaseIndex], settings, limits);
            if (sample.conductanceAudit.has_value()
                && sample.conductanceAudit
                       ->structureDefinitionFingerprint
                    != result.structureDefinitionFingerprint) {
                throw std::logic_error(
                    "scene pressure-cell mimetic conductance structure changed");
            }
            level.samples.push_back(std::move(sample));
            ++nextSampleIndex;
        }
        summarize(level);
        result.levels.push_back(std::move(level));
        result.ownedStorageBytes = storageBytes(result);
        if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
            throw std::length_error(
                "scene pressure-cell mimetic conductance byte limit exceeded");
        }
        previousCellCount = cellCount;
    }
    result.fingerprint = productFingerprint(result);
    validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
        result);
    return result;
}

void validateScenePressureCellMimeticConductancePhaseRefinementAuditIntegrity(
    const ScenePressureCellMimeticConductancePhaseRefinementAudit& audit) {
    validateSettings(audit.settings);
    validatePhases(audit.gridPhaseFractions);
    if (audit.version
            != scenePressureCellMimeticConductancePhaseRefinementAuditVersion
        || audit.fingerprint == 0
        || audit.structureDefinitionFingerprint == 0
        || audit.gridPhaseFractions.empty() || audit.levels.empty()
        || audit.ownedStorageBytes != storageBytes(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic conductance integrity is invalid");
    }

    std::size_t previousCellCount = 0;
    std::size_t nextSampleIndex = 0;
    double intakeAreaSquareMeters = 0.0;
    for (std::size_t levelIndex = 0;
         levelIndex < audit.levels.size(); ++levelIndex) {
        const auto& level = audit.levels[levelIndex];
        const std::size_t cellCount = gridCellCount(level.cellCounts);
        const double spacing = 4.0
            / static_cast<double>(level.cellCounts.x);
        if (level.levelIndex != levelIndex
            || level.cellCounts.x < 2
            || level.cellCounts.x != level.cellCounts.y
            || level.cellCounts.x != level.cellCounts.z
            || cellCount <= previousCellCount
            || level.cellSpacingMeters
                != fluid::Vector3{spacing, spacing, spacing}
            || level.samples.size() != audit.gridPhaseFractions.size()) {
            throw std::invalid_argument(
                "scene pressure-cell mimetic conductance level is invalid");
        }
        for (std::size_t phaseIndex = 0;
             phaseIndex < level.samples.size(); ++phaseIndex) {
            const auto& sample = level.samples[phaseIndex];
            const auto phase = audit.gridPhaseFractions[phaseIndex];
            const fluid::Vector3 lower = translatedGridLower(
                audit.settings.geometryTranslationMeters, phase, spacing);
            if (sample.sampleIndex != nextSampleIndex
                || sample.phaseIndex != phaseIndex
                || sample.gridPhaseFraction != phase
                || sample.gridCellCount != cellCount
                || sample.gridLowerMeters != lower
                || sample.cellSpacingMeters != level.cellSpacingMeters
                || !std::isfinite(sample.intakeAreaSquareMeters)
                || !(sample.intakeAreaSquareMeters > 0.0)
                || (nextSampleIndex != 0
                    && sample.intakeAreaSquareMeters
                        != intakeAreaSquareMeters)
                || sample.controlVolumeCount == 0
                || sample.openingTraceCount == 0) {
                throw std::invalid_argument(
                    "scene pressure-cell mimetic conductance sample is invalid");
            }
            if (sample.status
                == ScenePressureCellMimeticConductancePhaseSampleStatus::
                    Accepted) {
                if (sample.localCellLinearConsistencyRejection.has_value()
                    || !sample.conductanceAudit.has_value()) {
                    throw std::invalid_argument(
                        "scene pressure-cell accepted mimetic conductance sample is incomplete");
                }
                const auto& conductance = *sample.conductanceAudit;
                validateSceneFluidMimeticRegionConductanceAuditIntegrity(
                    conductance);
                const double areaDifference = std::abs(
                    sample.intakeAreaSquareMeters
                    - conductance.openingAreaSquareMeters);
                const double areaTolerance = audit.settings.conductance
                        .absoluteOpeningPairAreaToleranceSquareMeters
                    + audit.settings.conductance
                        .relativeOpeningPairAreaTolerance
                        * sample.intakeAreaSquareMeters;
                if (areaDifference > areaTolerance
                    || sample.controlVolumeCount
                        != conductance.responses.size()
                    || sample.fullTraceCount == 0
                    || sample.reducedTraceCount == 0
                    || sample.reducedTraceCount > sample.fullTraceCount
                    || sample.reducedTraceCount
                        != conductance.solveDiagnostics
                            .reducedTraceSolve.traceCount
                    || sample.openingTraceCount
                        != conductance.openings.size()
                    || sample.conductanceMeters
                        != conductance.conductanceMeters
                    || sample.normalizedConductance
                        != sample.conductanceMeters
                            * spacing / sample.intakeAreaSquareMeters
                    || !std::isfinite(sample.normalizedConductance)
                    || !(sample.normalizedConductance > 0.0)
                    || conductance.structureDefinitionFingerprint
                        != audit.structureDefinitionFingerprint
                    || conductance.settings
                        != audit.settings.conductance) {
                    throw std::invalid_argument(
                        "scene pressure-cell accepted mimetic conductance sample is invalid");
                }
            } else if (sample.status
                       == ScenePressureCellMimeticConductancePhaseSampleStatus::
                           RejectedLocalCellLinearConsistency) {
                if (!sample.localCellLinearConsistencyRejection.has_value()
                    || sample.conductanceAudit.has_value()
                    || sample.fullTraceCount != 0
                    || sample.reducedTraceCount != 0
                    || sample.conductanceMeters != 0.0
                    || sample.normalizedConductance != 0.0) {
                    throw std::invalid_argument(
                        "scene pressure-cell rejected mimetic conductance sample is incomplete");
                }
                const auto& rejection =
                    *sample.localCellLinearConsistencyRejection;
                if (rejection.controlCellIndex
                        >= sample.controlVolumeCount
                    || rejection.controlCellStableId == 0
                    || rejection.regionId == invalidStableId
                    || !validLinearConsistencyFailure(
                        rejection.localCell)) {
                    throw std::invalid_argument(
                        "scene pressure-cell mimetic conductance rejection is invalid");
                }
            } else {
                throw std::invalid_argument(
                    "scene pressure-cell mimetic conductance sample status is invalid");
            }
            if (nextSampleIndex == 0) {
                intakeAreaSquareMeters = sample.intakeAreaSquareMeters;
            }
            ++nextSampleIndex;
        }
        auto expected = level;
        expected.acceptedSampleCount = 0;
        expected.rejectedLocalCellLinearConsistencySampleCount = 0;
        expected.minimumNormalizedConductance = 0.0;
        expected.maximumNormalizedConductance = 0.0;
        expected.meanNormalizedConductance = 0.0;
        expected.normalizedConductanceCoefficientOfVariation = 0.0;
        summarize(expected);
        if (level.acceptedSampleCount != expected.acceptedSampleCount
            || level.rejectedLocalCellLinearConsistencySampleCount
                != expected.rejectedLocalCellLinearConsistencySampleCount
            || level.minimumNormalizedConductance
                != expected.minimumNormalizedConductance
            || level.maximumNormalizedConductance
                != expected.maximumNormalizedConductance
            || level.meanNormalizedConductance
                != expected.meanNormalizedConductance
            || level.normalizedConductanceCoefficientOfVariation
                != expected.normalizedConductanceCoefficientOfVariation) {
            throw std::invalid_argument(
                "scene pressure-cell mimetic conductance statistics are invalid");
        }
        previousCellCount = cellCount;
    }
    if (audit.fingerprint != productFingerprint(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell mimetic conductance fingerprint is invalid");
    }
}

} // namespace simwing::fsi
