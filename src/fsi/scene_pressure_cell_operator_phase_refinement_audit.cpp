#include "scene_pressure_cell_operator_phase_refinement_audit.h"

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

std::size_t gridCellCount(const fluid::GridCellCounts counts) {
    std::size_t xy = 0;
    std::size_t xyz = 0;
    if (!checkedMultiply(counts.x, counts.y, xy)
        || !checkedMultiply(xy, counts.z, xyz)) {
        throw std::length_error(
            "scene pressure-cell phase-refinement grid count overflows");
    }
    return xyz;
}

bool validPhaseComponent(const double phase) {
    return std::isfinite(phase)
        && phase >= -0.5 && phase < 0.5
        && (phase != 0.0 || !std::signbit(phase));
}

void validatePhases(const std::span<const fluid::Vector3> phases) {
    for (std::size_t index = 0; index < phases.size(); ++index) {
        const auto phase = phases[index];
        if (!validPhaseComponent(phase.x)
            || !validPhaseComponent(phase.y)
            || !validPhaseComponent(phase.z)) {
            throw std::invalid_argument(
                "scene pressure-cell phase-refinement phase is invalid");
        }
        for (std::size_t earlier = 0; earlier < index; ++earlier) {
            if (phases[earlier] == phase) {
                throw std::invalid_argument(
                    "scene pressure-cell phase-refinement phases must be unique");
            }
        }
    }
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

std::size_t storageBytes(
    const ScenePressureCellOperatorPhaseRefinementAudit& audit) {
    std::size_t phaseBytes = 0;
    std::size_t levelBytes = 0;
    if (!checkedMultiply(audit.gridPhaseFractions.size(),
                         sizeof(fluid::Vector3), phaseBytes)
        || !checkedMultiply(audit.levels.size(),
                            sizeof(ScenePressureCellOperatorPhaseRefinementLevel),
                            levelBytes)
        || levelBytes > std::numeric_limits<std::size_t>::max() - phaseBytes) {
        throw std::length_error(
            "scene pressure-cell phase-refinement storage overflows");
    }
    std::size_t result = phaseBytes + levelBytes;
    for (const auto& level : audit.levels) {
        if (level.phaseAudit.ownedStorageBytes
            > std::numeric_limits<std::size_t>::max() - result) {
            throw std::length_error(
                "scene pressure-cell phase-refinement storage overflows");
        }
        result += level.phaseAudit.ownedStorageBytes;
    }
    return result;
}

std::uint64_t productFingerprint(
    const ScenePressureCellOperatorPhaseRefinementAudit& audit) {
    Fingerprint fingerprint;
    fingerprint.integer(audit.version);
    fingerprint.integer(audit.structureDefinitionFingerprint);
    fingerprintResponseSettings(fingerprint, audit.settings.response);
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
        fingerprint.real(level.acceptedTopologyFraction);
        fingerprint.integer(level.phaseAudit.fingerprint);
    }
    return fingerprint.value();
}

} // namespace

ScenePressureCellOperatorPhaseRefinementAudit
auditScenePressureCellOperatorPhaseRefinement(
    const std::span<const fluid::GridCellCounts> resolutions,
    const std::span<const fluid::Vector3> gridPhaseFractions,
    const ScenePressureCellOperatorPhaseRefinementAuditSettings& settings,
    const ScenePressureCellOperatorPhaseRefinementAuditLimits& limits) {
    if (resolutions.empty()
        || resolutions.size() > limits.maximumResolutionCount
        || gridPhaseFractions.empty()
        || gridPhaseFractions.size() > limits.maximumPhaseCount) {
        throw std::length_error(
            "scene pressure-cell phase-refinement sample limit exceeded");
    }
    std::size_t aggregateSampleCount = 0;
    if (!checkedMultiply(resolutions.size(), gridPhaseFractions.size(),
                         aggregateSampleCount)
        || aggregateSampleCount > limits.maximumAggregatePhaseSamples) {
        throw std::length_error(
            "scene pressure-cell phase-refinement aggregate limit exceeded");
    }
    std::size_t phaseBytes = 0;
    std::size_t levelBytes = 0;
    if (!checkedMultiply(gridPhaseFractions.size(),
                         sizeof(fluid::Vector3), phaseBytes)
        || !checkedMultiply(
            resolutions.size(),
            sizeof(ScenePressureCellOperatorPhaseRefinementLevel),
            levelBytes)
        || levelBytes > std::numeric_limits<std::size_t>::max() - phaseBytes
        || phaseBytes + levelBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene pressure-cell phase-refinement byte limit exceeded");
    }
    validatePhases(gridPhaseFractions);
    std::size_t previousGridCellCount = 0;
    for (const auto counts : resolutions) {
        const std::size_t count = gridCellCount(counts);
        if (counts.x < 2 || counts.x != counts.y || counts.x != counts.z
            || count <= previousGridCellCount) {
            throw std::invalid_argument(
                "scene pressure-cell phase-refinement resolutions must be "
                "strictly increasing isotropic grids");
        }
        if (count > limits.maximumGridCellsPerSample) {
            throw std::length_error(
                "scene pressure-cell phase-refinement grid limit exceeded");
        }
        previousGridCellCount = count;
    }

    ScenePressureCellOperatorPhaseRefinementAudit result;
    result.settings = settings;
    result.gridPhaseFractions.assign(
        gridPhaseFractions.begin(), gridPhaseFractions.end());
    result.levels.reserve(resolutions.size());
    for (std::size_t index = 0; index < resolutions.size(); ++index) {
        ScenePressureCellOperatorPhaseAuditSettings phaseSettings;
        phaseSettings.cellCounts = resolutions[index];
        phaseSettings.response = settings.response;
        ScenePressureCellOperatorPhaseAuditLimits phaseLimits;
        phaseLimits.maximumSamples = limits.maximumPhaseCount;
        phaseLimits.maximumGridCellsPerSample =
            limits.maximumGridCellsPerSample;
        phaseLimits.maximumOwnedBytes = limits.maximumOwnedBytes;
        phaseLimits.response = limits.response;
        auto phaseAudit = auditScenePressureCellOperatorGridPhases(
            gridPhaseFractions, phaseSettings, phaseLimits);

        ScenePressureCellOperatorPhaseRefinementLevel level;
        level.levelIndex = index;
        level.cellCounts = resolutions[index];
        level.cellSpacingMeters = {
            4.0 / static_cast<double>(resolutions[index].x),
            4.0 / static_cast<double>(resolutions[index].y),
            4.0 / static_cast<double>(resolutions[index].z),
        };
        level.acceptedTopologyFraction = static_cast<double>(
            phaseAudit.statistics.acceptedSampleCount)
            / static_cast<double>(gridPhaseFractions.size());
        level.phaseAudit = std::move(phaseAudit);
        if (index == 0) {
            result.structureDefinitionFingerprint =
                level.phaseAudit.structureDefinitionFingerprint;
        } else if (level.phaseAudit.structureDefinitionFingerprint
                   != result.structureDefinitionFingerprint) {
            throw std::logic_error(
                "scene pressure-cell phase-refinement structure changed");
        }
        result.levels.push_back(std::move(level));
        result.ownedStorageBytes = storageBytes(result);
        if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
            throw std::length_error(
                "scene pressure-cell phase-refinement byte limit exceeded");
        }
    }
    result.fingerprint = productFingerprint(result);
    validateScenePressureCellOperatorPhaseRefinementAuditIntegrity(result);
    return result;
}

void validateScenePressureCellOperatorPhaseRefinementAuditIntegrity(
    const ScenePressureCellOperatorPhaseRefinementAudit& audit) {
    if (audit.version
            != scenePressureCellOperatorPhaseRefinementAuditVersion
        || audit.fingerprint == 0
        || audit.structureDefinitionFingerprint == 0
        || audit.gridPhaseFractions.empty()
        || audit.levels.empty()
        || audit.ownedStorageBytes != storageBytes(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell phase-refinement integrity is invalid");
    }
    validatePhases(audit.gridPhaseFractions);
    std::size_t previousGridCellCount = 0;
    for (std::size_t index = 0; index < audit.levels.size(); ++index) {
        const auto& level = audit.levels[index];
        validateScenePressureCellOperatorPhaseAuditIntegrity(
            level.phaseAudit);
        const std::size_t count = gridCellCount(level.cellCounts);
        const fluid::Vector3 expectedSpacing{
            4.0 / static_cast<double>(level.cellCounts.x),
            4.0 / static_cast<double>(level.cellCounts.y),
            4.0 / static_cast<double>(level.cellCounts.z),
        };
        const double expectedTopologyFraction = static_cast<double>(
            level.phaseAudit.statistics.acceptedSampleCount)
            / static_cast<double>(audit.gridPhaseFractions.size());
        if (level.levelIndex != index
            || level.cellCounts.x < 2
            || level.cellCounts.x != level.cellCounts.y
            || level.cellCounts.x != level.cellCounts.z
            || count <= previousGridCellCount
            || level.cellSpacingMeters != expectedSpacing
            || level.acceptedTopologyFraction != expectedTopologyFraction
            || !std::isfinite(level.acceptedTopologyFraction)
            || level.acceptedTopologyFraction < 0.0
            || level.acceptedTopologyFraction > 1.0
            || level.phaseAudit.structureDefinitionFingerprint
                != audit.structureDefinitionFingerprint
            || level.phaseAudit.settings.cellCounts != level.cellCounts
            || level.phaseAudit.settings.response != audit.settings.response
            || level.phaseAudit.samples.size()
                != audit.gridPhaseFractions.size()) {
            throw std::invalid_argument(
                "scene pressure-cell phase-refinement level is invalid");
        }
        for (std::size_t phase = 0;
             phase < audit.gridPhaseFractions.size(); ++phase) {
            if (level.phaseAudit.samples[phase].gridPhaseFraction
                != audit.gridPhaseFractions[phase]) {
                throw std::invalid_argument(
                    "scene pressure-cell phase-refinement phases changed");
            }
        }
        previousGridCellCount = count;
    }
    if (audit.fingerprint != productFingerprint(audit)) {
        throw std::invalid_argument(
            "scene pressure-cell phase-refinement fingerprint is invalid");
    }
}

} // namespace simwing::fsi
