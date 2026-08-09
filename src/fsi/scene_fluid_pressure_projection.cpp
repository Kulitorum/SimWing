#include "scene_fluid_pressure_projection.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <ranges>
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
            byte(static_cast<std::uint8_t>(value & 0xffU));
            value >>= 8U;
        }
    }

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
    }

    void boolean(const bool value) {
        integer(static_cast<std::uint8_t>(value ? 1 : 0));
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_ == 0 ? 1 : value_;
    }

private:
    void byte(const std::uint8_t value) noexcept {
        value_ ^= value;
        value_ *= fnvPrime;
    }

    std::uint64_t value_ = fnvOffsetBasis;
};

bool checkedAdd(const std::size_t first,
                const std::size_t second,
                std::size_t& result) {
    if (second > std::numeric_limits<std::size_t>::max() - first) {
        return false;
    }
    result = first + second;
    return true;
}

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

std::size_t storageBytesForCounts(const std::size_t pressureCount,
                                  const std::size_t controlCount,
                                  const std::size_t linkCount,
                                  const std::size_t solveComponentCount) {
    std::size_t pressureBytes = 0;
    std::size_t controlBytes = 0;
    std::size_t linkBytes = 0;
    std::size_t solveComponentBytes = 0;
    std::size_t first = 0;
    std::size_t second = 0;
    std::size_t total = 0;
    if (!checkedMultiply(
            pressureCount, sizeof(double), pressureBytes)
        || !checkedMultiply(
            controlCount,
            sizeof(SceneFluidPressureProjectedControlVolume), controlBytes)
        || !checkedMultiply(
            linkCount,
            sizeof(SceneFluidPressureProjectedLink), linkBytes)
        || !checkedMultiply(
            solveComponentCount,
            sizeof(SceneFluidPressureSolveComponentDiagnostics),
            solveComponentBytes)
        || !checkedAdd(pressureBytes, controlBytes, first)
        || !checkedAdd(linkBytes, solveComponentBytes, second)
        || !checkedAdd(first, second, total)) {
        throw std::length_error(
            "scene fluid pressure-projection storage size overflows");
    }
    return total;
}

std::size_t storageBytes(const SceneFluidPressureProjection& projection) {
    return storageBytesForCounts(
        projection.pressurePascals.size(),
        projection.controlVolumes.size(),
        projection.links.size(),
        projection.diagnostics.pressureSolve.components.size());
}

void fingerprintSolveDiagnostics(
    Fingerprint& fingerprint,
    const SceneFluidPressureSolveDiagnostics& diagnostics) {
    fingerprint.boolean(diagnostics.compatible);
    fingerprint.boolean(diagnostics.converged);
    fingerprint.boolean(diagnostics.finite);
    fingerprint.integer(diagnostics.pressureOperatorFingerprint);
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.rowCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.componentCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.iterationCount));
    fingerprint.real(
        diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters);
    fingerprint.real(diagnostics.initialResidualL2PascalsMeters);
    fingerprint.real(diagnostics.finalResidualL2PascalsMeters);
    fingerprint.real(diagnostics.finalResidualMaximumPascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.components.size()));
    for (const auto& component : diagnostics.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.controlVolumeCount));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.gaugeControlVolumeIndex));
        fingerprint.real(component.rightHandSideSumPascalsMeters);
        fingerprint.real(component.compatibilityCorrectionPascalsMeters);
        fingerprint.real(component.pressureGaugeBeforePascals);
        fingerprint.real(component.pressureGaugeAfterPascals);
    }
}

std::uint64_t projectionFingerprint(
    const SceneFluidPressureProjection& projection) {
    Fingerprint fingerprint;
    fingerprint.integer(projection.version);
    fingerprint.integer(projection.pressureOperatorFingerprint);
    fingerprint.integer(projection.pressureFaceLinkFingerprint);
    fingerprint.integer(projection.pressureControlVolumeFingerprint);
    fingerprint.integer(projection.openingFluxFingerprint);
    fingerprint.integer(projection.velocityFingerprint);
    fingerprint.integer(projection.acceptedStepCount);
    fingerprint.real(projection.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(projection.cellCounts.x));
    fingerprint.integer(static_cast<std::uint64_t>(projection.cellCounts.y));
    fingerprint.integer(static_cast<std::uint64_t>(projection.cellCounts.z));
    for (const double value : {
             projection.lowerMeters.x,
             projection.lowerMeters.y,
             projection.lowerMeters.z,
             projection.upperMeters.x,
             projection.upperMeters.y,
             projection.upperMeters.z}) {
        fingerprint.real(value);
    }
    fingerprint.real(projection.settings.densityKgPerCubicMeter);
    fingerprint.real(projection.settings.timeStepSeconds);
    fingerprint.real(projection.settings
        .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond);
    fingerprint.real(
        projection.settings.relativeCorrectedVolumeRateTolerance);
    fingerprint.real(projection.settings.pressureSolve
        .absoluteResidualTolerancePascalsMeters);
    fingerprint.real(
        projection.settings.pressureSolve.relativeResidualTolerance);
    fingerprint.real(projection.settings.pressureSolve
        .absoluteComponentCompatibilityTolerancePascalsMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        projection.settings.pressureSolve.maximumIterations));
    fingerprint.integer(static_cast<std::uint64_t>(
        projection.ownedStorageBytes));
    const auto& diagnostics = projection.diagnostics;
    fingerprint.boolean(diagnostics.accepted);
    fingerprint.boolean(diagnostics.finite);
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.controlVolumeCount));
    fingerprint.integer(static_cast<std::uint64_t>(diagnostics.linkCount));
    fingerprint.integer(static_cast<std::uint64_t>(
        diagnostics.authoredOpeningLinkCount));
    fingerprint.real(
        diagnostics.predictedNetOutwardVolumeRateL2CubicMetersPerSecond);
    fingerprint.real(
        diagnostics.predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond);
    fingerprint.real(diagnostics
        .maximumPredictedComponentBalanceResidualCubicMetersPerSecond);
    fingerprint.real(
        diagnostics.correctedNetOutwardVolumeRateL2CubicMetersPerSecond);
    fingerprint.real(
        diagnostics.correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond);
    fingerprint.real(diagnostics
        .maximumCorrectedComponentBalanceResidualCubicMetersPerSecond);
    fingerprintSolveDiagnostics(fingerprint, diagnostics.pressureSolve);
    fingerprint.integer(static_cast<std::uint64_t>(
        projection.pressurePascals.size()));
    for (const double pressure : projection.pressurePascals) {
        fingerprint.real(pressure);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        projection.controlVolumes.size()));
    for (const auto& control : projection.controlVolumes) {
        fingerprint.integer(static_cast<std::uint64_t>(
            control.controlVolumeIndex));
        fingerprint.integer(control.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            control.componentIndex));
        fingerprint.real(
            control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(control.integratedRightHandSidePascalsMeters);
        fingerprint.real(
            control.correctedNetOutwardVolumeFlowRateCubicMetersPerSecond);
    }
    fingerprint.integer(static_cast<std::uint64_t>(projection.links.size()));
    for (const auto& link : projection.links) {
        fingerprint.integer(static_cast<std::uint64_t>(link.linkIndex));
        fingerprint.integer(link.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(link.faceIndex));
        fingerprint.enumeration(link.kind);
        fingerprint.integer(static_cast<std::uint64_t>(
            link.minusControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            link.plusControlVolumeIndex));
        fingerprint.integer(link.openingPatchStableId);
        fingerprint.real(
            link.predictedRelativeVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            link.pressureCorrectionVolumeFlowRateCubicMetersPerSecond);
        fingerprint.real(
            link.correctedRelativeVolumeFlowRateCubicMetersPerSecond);
    }
    return fingerprint.value();
}

void validateSettings(
    const SceneFluidPressureProjectionSettings& settings) {
    if (!std::isfinite(settings.densityKgPerCubicMeter)
        || !(settings.densityKgPerCubicMeter > 0.0)
        || !std::isfinite(settings.timeStepSeconds)
        || !(settings.timeStepSeconds > 0.0)
        || !std::isfinite(settings
            .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond)
        || settings
            .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond < 0.0
        || !std::isfinite(
            settings.relativeCorrectedVolumeRateTolerance)
        || settings.relativeCorrectedVolumeRateTolerance < 0.0
        || (settings
                .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond
                == 0.0
            && settings.relativeCorrectedVolumeRateTolerance == 0.0)
        || !std::isfinite(settings.pressureSolve
            .absoluteResidualTolerancePascalsMeters)
        || settings.pressureSolve
            .absoluteResidualTolerancePascalsMeters < 0.0
        || !std::isfinite(
            settings.pressureSolve.relativeResidualTolerance)
        || settings.pressureSolve.relativeResidualTolerance < 0.0
        || (settings.pressureSolve
                .absoluteResidualTolerancePascalsMeters == 0.0
            && settings.pressureSolve.relativeResidualTolerance == 0.0)
        || !std::isfinite(settings.pressureSolve
            .absoluteComponentCompatibilityTolerancePascalsMeters)
        || settings.pressureSolve
            .absoluteComponentCompatibilityTolerancePascalsMeters < 0.0) {
        throw std::invalid_argument(
            "scene fluid pressure-projection settings are invalid");
    }
}

double faceVelocity(const SceneFluidPressureFace& face,
                    const fluid::PeriodicCartesianGrid& grid,
                    const fluid::MacVelocityField& velocity) {
    const std::size_t index = grid.cellIndex(face.i, face.j, face.k);
    switch (face.axis) {
    case fluid::GridFaceAxis::X:
        return velocity.xFaces()[index];
    case fluid::GridFaceAxis::Y:
        return velocity.yFaces()[index];
    case fluid::GridFaceAxis::Z:
        return velocity.zFaces()[index];
    }
    throw std::invalid_argument(
        "scene fluid pressure projection has an invalid face axis");
}

double vectorL2(const std::vector<SceneFluidPressureProjectedControlVolume>&
                    controls,
                const bool corrected) {
    double squared = 0.0;
    for (const auto& control : controls) {
        const double value = corrected
            ? control.correctedNetOutwardVolumeFlowRateCubicMetersPerSecond
            : control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        squared += value * value;
    }
    return std::sqrt(squared / static_cast<double>(controls.size()));
}

double vectorMaximum(
    const std::vector<SceneFluidPressureProjectedControlVolume>& controls,
    const bool corrected) {
    double maximum = 0.0;
    for (const auto& control : controls) {
        const double value = corrected
            ? control.correctedNetOutwardVolumeFlowRateCubicMetersPerSecond
            : control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        maximum = std::max(maximum, std::abs(value));
    }
    return maximum;
}

double maximumComponentBalance(
    const SceneFluidPressureOperator& pressureOperator,
    const std::vector<SceneFluidPressureProjectedControlVolume>& controls,
    const bool corrected) {
    double maximum = 0.0;
    for (const auto& component : pressureOperator.components) {
        double balance = 0.0;
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            const std::size_t controlIndex =
                pressureOperator.componentControlVolumeIndices[
                    component.firstControlVolumeMember + offset];
            const auto& control = controls[controlIndex];
            balance += corrected
                ? control
                    .correctedNetOutwardVolumeFlowRateCubicMetersPerSecond
                : control
                    .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        }
        maximum = std::max(maximum, std::abs(balance));
    }
    return maximum;
}

void addOrientedFlow(
    std::vector<SceneFluidPressureProjectedControlVolume>& controls,
    const SceneFluidPressureProjectedLink& link,
    const double flow,
    const bool corrected) {
    auto& minus = controls[link.minusControlVolumeIndex];
    auto& plus = controls[link.plusControlVolumeIndex];
    if (corrected) {
        minus.correctedNetOutwardVolumeFlowRateCubicMetersPerSecond += flow;
        plus.correctedNetOutwardVolumeFlowRateCubicMetersPerSecond -= flow;
    } else {
        minus.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond += flow;
        plus.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond -= flow;
    }
}

} // namespace

SceneFluidPressureProjection projectSceneFluidPressureLinkFlows(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidOpeningFluxSet& openingFlux,
    const fluid::MacVelocityField& predictedVelocityMetersPerSecond,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperator& pressureOperator,
    const std::span<const double> warmPressurePascals,
    const SceneFluidPressureProjectionSettings& settings,
    const SceneFluidPressureProjectionLimits& limits) {
    validateSettings(settings);
    validateSceneFluidPressureOperator(
        pressureOperator, surface, state, grid, transfer, epoch, caps,
        openingQuadrature, openingPatches, volumes, connectivity,
        pressureVolumes, faceLinks);
    validateSceneFluidOpeningFlux(
        openingFlux, surface, state, caps, openingQuadrature,
        openingPatches, grid, predictedVelocityMetersPerSecond);
    if (pressureVolumes.controlVolumes.size() > limits.maximumControlVolumes
        || faceLinks.links.size() > limits.maximumLinks) {
        throw std::length_error(
            "scene fluid pressure projection exceeds its count limit");
    }
    const std::size_t maximumOwnedStorageBytes = storageBytesForCounts(
        pressureOperator.rows.size(), pressureVolumes.controlVolumes.size(),
        faceLinks.links.size(), pressureOperator.components.size());
    if (maximumOwnedStorageBytes > limits.maximumProjectionBytes) {
        throw std::length_error(
            "scene fluid pressure projection exceeds its byte limit");
    }
    if (warmPressurePascals.size() != pressureOperator.rows.size()
        || !std::ranges::all_of(
            warmPressurePascals,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "scene fluid pressure projection has an invalid warm start");
    }

    SceneFluidPressureProjection result;
    result.pressureOperatorFingerprint = pressureOperator.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.openingFluxFingerprint = openingFlux.fingerprint;
    result.velocityFingerprint = openingFlux.velocityFingerprint;
    result.acceptedStepCount = pressureOperator.acceptedStepCount;
    result.simulationTimeSeconds = pressureOperator.simulationTimeSeconds;
    result.cellCounts = grid.cellCounts();
    result.lowerMeters = grid.lowerMeters();
    result.upperMeters = grid.upperMeters();
    result.settings = settings;
    result.diagnostics.controlVolumeCount =
        pressureVolumes.controlVolumes.size();
    result.diagnostics.linkCount = faceLinks.links.size();
    result.controlVolumes.reserve(pressureVolumes.controlVolumes.size());
    for (const auto& source : pressureVolumes.controlVolumes) {
        result.controlVolumes.push_back({
            source.controlVolumeIndex,
            source.stableId,
            source.componentIndex,
        });
    }
    result.links.reserve(faceLinks.links.size());

    std::map<std::uint64_t, const SceneFluidOpeningFluxSample*>
        openingSamples;
    for (const auto& sample : openingFlux.samples) {
        if (sample.patchStableId == 0
            || !openingSamples.emplace(sample.patchStableId, &sample).second) {
            throw std::invalid_argument(
                "scene fluid pressure projection has duplicate opening-flux identity");
        }
    }
    std::size_t consumedOpeningSamples = 0;
    for (const auto& source : faceLinks.links) {
        const auto& face = faceLinks.faces[source.faceIndex];
        double predictedFlow = 0.0;
        if (source.kind == SceneFluidPressureFaceLinkKind::SameRegion) {
            predictedFlow = source.areaSquareMeters
                * faceVelocity(face, grid, predictedVelocityMetersPerSecond);
        } else {
            const auto found = openingSamples.find(
                source.openingPatchStableId);
            if (found == openingSamples.end()) {
                throw std::invalid_argument(
                    "scene fluid pressure projection is missing an opening-flux sample");
            }
            const auto& sample = *found->second;
            if (sample.openingId != source.openingId
                || sample.areaSquareMeters != source.areaSquareMeters) {
                throw std::invalid_argument(
                    "scene fluid pressure projection opening-flux identity is inconsistent");
            }
            if (sample.negativeSideRegionId == source.minusRegionId
                && sample.positiveSideRegionId == source.plusRegionId) {
                predictedFlow =
                    sample.relativeVolumeFlowRateCubicMetersPerSecond;
            } else if (sample.negativeSideRegionId == source.plusRegionId
                       && sample.positiveSideRegionId
                           == source.minusRegionId) {
                predictedFlow =
                    -sample.relativeVolumeFlowRateCubicMetersPerSecond;
            } else {
                throw std::invalid_argument(
                    "scene fluid pressure projection opening orientation is inconsistent");
            }
            openingSamples.erase(found);
            ++consumedOpeningSamples;
            ++result.diagnostics.authoredOpeningLinkCount;
        }
        if (!std::isfinite(predictedFlow)) {
            throw std::overflow_error(
                "scene fluid pressure projection predicted flow is not finite");
        }
        SceneFluidPressureProjectedLink link;
        link.linkIndex = source.linkIndex;
        link.stableId = source.stableId;
        link.faceIndex = source.faceIndex;
        link.kind = source.kind;
        link.minusControlVolumeIndex = source.minusControlVolumeIndex;
        link.plusControlVolumeIndex = source.plusControlVolumeIndex;
        link.openingPatchStableId = source.openingPatchStableId;
        link.predictedRelativeVolumeFlowRateCubicMetersPerSecond =
            predictedFlow;
        result.links.push_back(link);
        addOrientedFlow(
            result.controlVolumes, result.links.back(), predictedFlow,
            false);
    }
    if (consumedOpeningSamples != openingFlux.samples.size()
        || !openingSamples.empty()) {
        throw std::invalid_argument(
            "scene fluid pressure projection cannot resolve every opening-flux sample");
    }

    const double rightHandSideScale =
        -settings.densityKgPerCubicMeter / settings.timeStepSeconds;
    for (auto& control : result.controlVolumes) {
        control.integratedRightHandSidePascalsMeters =
            rightHandSideScale
            * control.predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        if (!std::isfinite(control.integratedRightHandSidePascalsMeters)) {
            throw std::overflow_error(
                "scene fluid pressure projection right-hand side overflowed");
        }
    }
    result.diagnostics.predictedNetOutwardVolumeRateL2CubicMetersPerSecond =
        vectorL2(result.controlVolumes, false);
    result.diagnostics
        .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond =
        vectorMaximum(result.controlVolumes, false);
    result.diagnostics
        .maximumPredictedComponentBalanceResidualCubicMetersPerSecond =
        maximumComponentBalance(
            pressureOperator, result.controlVolumes, false);

    std::vector<double> rightHandSide;
    rightHandSide.reserve(result.controlVolumes.size());
    for (const auto& control : result.controlVolumes) {
        rightHandSide.push_back(
            control.integratedRightHandSidePascalsMeters);
    }
    std::vector<double> candidatePressure(
        warmPressurePascals.begin(), warmPressurePascals.end());
    result.diagnostics.pressureSolve = solveSceneFluidPressureSystem(
        pressureOperator, rightHandSide, candidatePressure,
        settings.pressureSolve);
    result.diagnostics.finite = result.diagnostics.pressureSolve.finite;
    if (result.diagnostics.pressureSolve.converged
        && result.diagnostics.pressureSolve.finite) {
        const double correctionScale =
            settings.timeStepSeconds / settings.densityKgPerCubicMeter;
        bool correctedFinite = true;
        for (auto& link : result.links) {
            const auto& source = faceLinks.links[link.linkIndex];
            link.pressureCorrectionVolumeFlowRateCubicMetersPerSecond =
                correctionScale * source.geometryWeightMeters
                * (candidatePressure[link.minusControlVolumeIndex]
                   - candidatePressure[link.plusControlVolumeIndex]);
            link.correctedRelativeVolumeFlowRateCubicMetersPerSecond =
                link.predictedRelativeVolumeFlowRateCubicMetersPerSecond
                + link.pressureCorrectionVolumeFlowRateCubicMetersPerSecond;
            correctedFinite = correctedFinite
                && std::isfinite(link
                    .pressureCorrectionVolumeFlowRateCubicMetersPerSecond)
                && std::isfinite(link
                    .correctedRelativeVolumeFlowRateCubicMetersPerSecond);
            if (correctedFinite) {
                addOrientedFlow(
                    result.controlVolumes, link,
                    link.correctedRelativeVolumeFlowRateCubicMetersPerSecond,
                    true);
            }
        }
        result.diagnostics.finite = correctedFinite;
        if (correctedFinite) {
            result.diagnostics
                .correctedNetOutwardVolumeRateL2CubicMetersPerSecond =
                vectorL2(result.controlVolumes, true);
            result.diagnostics
                .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond =
                vectorMaximum(result.controlVolumes, true);
            result.diagnostics
                .maximumCorrectedComponentBalanceResidualCubicMetersPerSecond =
                maximumComponentBalance(
                    pressureOperator, result.controlVolumes, true);
            const double tolerance = std::max(
                settings
                    .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond,
                settings.relativeCorrectedVolumeRateTolerance
                    * result.diagnostics
                        .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond);
            result.diagnostics.accepted =
                result.diagnostics
                    .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                <= tolerance;
        }
        if (result.diagnostics.accepted) {
            result.pressurePascals = std::move(candidatePressure);
        } else {
            for (auto& control : result.controlVolumes) {
                control.correctedNetOutwardVolumeFlowRateCubicMetersPerSecond =
                    0.0;
            }
            for (auto& link : result.links) {
                link.pressureCorrectionVolumeFlowRateCubicMetersPerSecond =
                    0.0;
                link.correctedRelativeVolumeFlowRateCubicMetersPerSecond =
                    0.0;
            }
        }
    }

    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumProjectionBytes) {
        throw std::length_error(
            "scene fluid pressure projection exceeds its byte limit");
    }
    result.fingerprint = projectionFingerprint(result);
    validateSceneFluidPressureProjectionIntegrity(result);
    return result;
}

void validateSceneFluidPressureProjectionIntegrity(
    const SceneFluidPressureProjection& projection) {
    validateSettings(projection.settings);
    const auto finite = [](const double value) {
        return std::isfinite(value);
    };
    const auto& diagnostics = projection.diagnostics;
    if (projection.version != sceneFluidPressureProjectionVersion
        || projection.fingerprint == 0
        || projection.pressureOperatorFingerprint == 0
        || projection.pressureFaceLinkFingerprint == 0
        || projection.pressureControlVolumeFingerprint == 0
        || projection.openingFluxFingerprint == 0
        || projection.velocityFingerprint == 0
        || !finite(projection.simulationTimeSeconds)
        || !finite(projection.lowerMeters.x)
        || !finite(projection.lowerMeters.y)
        || !finite(projection.lowerMeters.z)
        || !finite(projection.upperMeters.x)
        || !finite(projection.upperMeters.y)
        || !finite(projection.upperMeters.z)
        || diagnostics.controlVolumeCount
            != projection.controlVolumes.size()
        || diagnostics.linkCount != projection.links.size()
        || projection.controlVolumes.empty()
        || !finite(diagnostics
            .predictedNetOutwardVolumeRateL2CubicMetersPerSecond)
        || !finite(diagnostics
            .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond)
        || !finite(diagnostics
            .maximumPredictedComponentBalanceResidualCubicMetersPerSecond)
        || !finite(diagnostics
            .correctedNetOutwardVolumeRateL2CubicMetersPerSecond)
        || !finite(diagnostics
            .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond)
        || !finite(diagnostics
            .maximumCorrectedComponentBalanceResidualCubicMetersPerSecond)
        || diagnostics.pressureSolve.pressureOperatorFingerprint
            != projection.pressureOperatorFingerprint
        || diagnostics.pressureSolve.rowCount
            != projection.controlVolumes.size()
        || diagnostics.pressureSolve.componentCount
            != diagnostics.pressureSolve.components.size()
        || projection.ownedStorageBytes != storageBytes(projection)
        || projection.fingerprint != projectionFingerprint(projection)) {
        throw std::invalid_argument(
            "scene fluid pressure-projection integrity is invalid");
    }
    for (std::size_t index = 0;
         index < projection.controlVolumes.size(); ++index) {
        const auto& control = projection.controlVolumes[index];
        const double expectedRightHandSide =
            -projection.settings.densityKgPerCubicMeter
            / projection.settings.timeStepSeconds
            * control
                .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond;
        if (control.controlVolumeIndex != index || control.stableId == 0
            || !finite(control
                .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond)
            || !finite(control.integratedRightHandSidePascalsMeters)
            || control.integratedRightHandSidePascalsMeters
                != expectedRightHandSide
            || !finite(control
                .correctedNetOutwardVolumeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid pressure-projection control ledger is invalid");
        }
    }
    std::size_t openingCount = 0;
    std::size_t cellCount = 0;
    std::size_t faceCount = 0;
    std::size_t xyCount = 0;
    if (!checkedMultiply(
            projection.cellCounts.x, projection.cellCounts.y, xyCount)
        || !checkedMultiply(xyCount, projection.cellCounts.z, cellCount)
        || !checkedMultiply(cellCount, std::size_t{3}, faceCount)) {
        throw std::invalid_argument(
            "scene fluid pressure-projection grid size is invalid");
    }
    for (std::size_t index = 0; index < projection.links.size(); ++index) {
        const auto& link = projection.links[index];
        if (link.linkIndex != index || link.stableId == 0
            || link.faceIndex >= faceCount
            || link.minusControlVolumeIndex
                >= projection.controlVolumes.size()
            || link.plusControlVolumeIndex
                >= projection.controlVolumes.size()
            || !finite(link
                .predictedRelativeVolumeFlowRateCubicMetersPerSecond)
            || !finite(link
                .pressureCorrectionVolumeFlowRateCubicMetersPerSecond)
            || !finite(link
                .correctedRelativeVolumeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid pressure-projection link ledger is invalid");
        }
        if (link.kind
            == SceneFluidPressureFaceLinkKind::AuthoredOpening) {
            ++openingCount;
            if (link.openingPatchStableId == 0) {
                throw std::invalid_argument(
                    "scene fluid pressure-projection opening link is invalid");
            }
        } else if (link.openingPatchStableId != 0) {
            throw std::invalid_argument(
                "scene fluid pressure-projection same-region link is invalid");
        }
    }
    if (openingCount != projection.diagnostics.authoredOpeningLinkCount
        || projection.diagnostics.accepted
            != (projection.diagnostics.pressureSolve.converged
                && projection.diagnostics.pressureSolve.finite
                && !projection.pressurePascals.empty())) {
        throw std::invalid_argument(
            "scene fluid pressure-projection acceptance is invalid");
    }
    std::vector<double> predictedLedger(
        projection.controlVolumes.size(), 0.0);
    std::vector<double> correctedLedger(
        projection.controlVolumes.size(), 0.0);
    for (const auto& link : projection.links) {
        predictedLedger[link.minusControlVolumeIndex] +=
            link.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        predictedLedger[link.plusControlVolumeIndex] -=
            link.predictedRelativeVolumeFlowRateCubicMetersPerSecond;
        if (projection.diagnostics.accepted) {
            correctedLedger[link.minusControlVolumeIndex] +=
                link.correctedRelativeVolumeFlowRateCubicMetersPerSecond;
            correctedLedger[link.plusControlVolumeIndex] -=
                link.correctedRelativeVolumeFlowRateCubicMetersPerSecond;
        }
    }
    for (std::size_t index = 0;
         index < projection.controlVolumes.size(); ++index) {
        const auto& control = projection.controlVolumes[index];
        if (predictedLedger[index]
                != control
                    .predictedNetOutwardVolumeFlowRateCubicMetersPerSecond
            || (projection.diagnostics.accepted
                && correctedLedger[index]
                    != control
                        .correctedNetOutwardVolumeFlowRateCubicMetersPerSecond)) {
            throw std::invalid_argument(
                "scene fluid pressure-projection incidence ledger is invalid");
        }
    }
    if (projection.diagnostics.accepted) {
        const double tolerance = std::max(
            projection.settings
                .absoluteCorrectedVolumeRateToleranceCubicMetersPerSecond,
            projection.settings.relativeCorrectedVolumeRateTolerance
                * projection.diagnostics
                    .predictedNetOutwardVolumeRateMaximumCubicMetersPerSecond);
        if (projection.pressurePascals.size()
                != projection.controlVolumes.size()
            || !std::ranges::all_of(projection.pressurePascals, finite)
            || projection.diagnostics
                .correctedNetOutwardVolumeRateMaximumCubicMetersPerSecond
                > tolerance) {
            throw std::invalid_argument(
                "scene fluid pressure-projection pressure is invalid");
        }
    } else {
        if (!projection.pressurePascals.empty()
            || std::ranges::any_of(
                projection.controlVolumes,
                [](const auto& control) {
                    return control
                        .correctedNetOutwardVolumeFlowRateCubicMetersPerSecond
                        != 0.0;
                })
            || std::ranges::any_of(
                projection.links,
                [](const auto& link) {
                    return link
                            .pressureCorrectionVolumeFlowRateCubicMetersPerSecond
                            != 0.0
                        || link
                            .correctedRelativeVolumeFlowRateCubicMetersPerSecond
                            != 0.0;
                })) {
            throw std::invalid_argument(
                "scene fluid rejected pressure projection exposes corrected state");
        }
    }
}

} // namespace simwing::fsi
