#include "scene_fluid_pressure_operator.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

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

std::size_t storageBytesForCounts(const std::size_t rowCount,
                                  const std::size_t entryCount,
                                  const std::size_t componentCount,
                                  const std::size_t memberCount) {
    std::size_t rowBytes = 0;
    std::size_t entryBytes = 0;
    std::size_t componentBytes = 0;
    std::size_t memberBytes = 0;
    std::size_t firstTotal = 0;
    std::size_t secondTotal = 0;
    std::size_t total = 0;
    if (!checkedMultiply(rowCount,
                         sizeof(SceneFluidPressureOperatorRow), rowBytes)
        || !checkedMultiply(entryCount,
                            sizeof(SceneFluidPressureOperatorEntry),
                            entryBytes)
        || !checkedMultiply(componentCount,
                            sizeof(SceneFluidPressureOperatorComponent),
                            componentBytes)
        || !checkedMultiply(memberCount, sizeof(std::size_t), memberBytes)
        || !checkedAdd(rowBytes, entryBytes, firstTotal)
        || !checkedAdd(componentBytes, memberBytes, secondTotal)
        || !checkedAdd(firstTotal, secondTotal, total)) {
        throw std::length_error(
            "scene fluid pressure-operator storage size overflows");
    }
    return total;
}

std::size_t storageBytes(
    const SceneFluidPressureOperator& pressureOperator) {
    return storageBytesForCounts(
        pressureOperator.rows.size(), pressureOperator.entries.size(),
        pressureOperator.components.size(),
        pressureOperator.componentControlVolumeIndices.size());
}

std::uint64_t entryStableId(
    const SceneFluidPressureFaceLink& link,
    const SceneFluidPressureControlVolume& source,
    const SceneFluidPressureControlVolume& destination) {
    Fingerprint fingerprint;
    fingerprint.integer(std::uint64_t{0x70726573736f7065ULL});
    fingerprint.integer(link.stableId);
    fingerprint.integer(source.stableId);
    fingerprint.integer(destination.stableId);
    return fingerprint.value();
}

std::uint64_t operatorFingerprint(
    const SceneFluidPressureOperator& pressureOperator) {
    Fingerprint fingerprint;
    fingerprint.integer(pressureOperator.version);
    for (const std::uint64_t value : {
             pressureOperator.pressureControlVolumeFingerprint,
             pressureOperator.pressureFaceLinkFingerprint,
             pressureOperator.surfaceDefinitionFingerprint,
             pressureOperator.surfaceStateFingerprint,
             pressureOperator.structureDefinitionFingerprint,
             pressureOperator.acceptedStepCount}) {
        fingerprint.integer(value);
    }
    fingerprint.real(pressureOperator.simulationTimeSeconds);
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.ownedStorageBytes));
    fingerprint.real(pressureOperator.totalGeometryWeightMeters);
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.rows.size()));
    for (const auto& row : pressureOperator.rows) {
        fingerprint.integer(static_cast<std::uint64_t>(row.rowIndex));
        fingerprint.integer(row.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            row.controlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(row.componentIndex));
        fingerprint.integer(static_cast<std::uint8_t>(row.isGauge));
        fingerprint.integer(static_cast<std::uint64_t>(row.firstEntry));
        fingerprint.integer(static_cast<std::uint64_t>(row.entryCount));
        fingerprint.real(row.diagonalGeometryWeightMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.entries.size()));
    for (const auto& entry : pressureOperator.entries) {
        fingerprint.integer(static_cast<std::uint64_t>(entry.entryIndex));
        fingerprint.integer(entry.stableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            entry.sourceFaceLinkIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            entry.columnControlVolumeIndex));
        fingerprint.real(entry.geometryWeightMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.components.size()));
    for (const auto& component : pressureOperator.components) {
        fingerprint.integer(static_cast<std::uint64_t>(
            component.componentIndex));
        fingerprint.integer(component.gaugeRegionId);
        fingerprint.integer(static_cast<std::uint64_t>(
            component.gaugeControlVolumeIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.firstControlVolumeMember));
        fingerprint.integer(static_cast<std::uint64_t>(
            component.controlVolumeCount));
        fingerprint.real(component.totalVolumeCubicMeters);
        fingerprint.real(component.totalGeometryWeightMeters);
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        pressureOperator.componentControlVolumeIndices.size()));
    for (const auto index :
         pressureOperator.componentControlVolumeIndices) {
        fingerprint.integer(static_cast<std::uint64_t>(index));
    }
    return fingerprint.value();
}

void validatePayload(
    const SceneFluidPressureOperator& pressureOperator) {
    if (pressureOperator.version != sceneFluidPressureOperatorVersion
        || pressureOperator.fingerprint == 0
        || pressureOperator.pressureControlVolumeFingerprint == 0
        || pressureOperator.pressureFaceLinkFingerprint == 0
        || pressureOperator.surfaceDefinitionFingerprint == 0
        || pressureOperator.surfaceStateFingerprint == 0
        || pressureOperator.structureDefinitionFingerprint == 0
        || !std::isfinite(pressureOperator.simulationTimeSeconds)
        || !std::isfinite(pressureOperator.totalGeometryWeightMeters)
        || pressureOperator.totalGeometryWeightMeters < 0.0
        || pressureOperator.rows.empty()
        || pressureOperator.components.empty()
        || pressureOperator.ownedStorageBytes
            != storageBytes(pressureOperator)
        || pressureOperator.fingerprint
            != operatorFingerprint(pressureOperator)) {
        throw std::invalid_argument(
            "scene fluid pressure-operator identity is invalid");
    }
    std::size_t expectedFirstEntry = 0;
    for (std::size_t rowIndex = 0;
         rowIndex < pressureOperator.rows.size(); ++rowIndex) {
        const auto& row = pressureOperator.rows[rowIndex];
        std::size_t rowEnd = 0;
        if (row.rowIndex != rowIndex
            || row.controlVolumeIndex != rowIndex
            || row.stableId == 0
            || row.componentIndex >= pressureOperator.components.size()
            || row.firstEntry != expectedFirstEntry
            || !checkedAdd(row.firstEntry, row.entryCount, rowEnd)
            || rowEnd > pressureOperator.entries.size()
            || !std::isfinite(row.diagonalGeometryWeightMeters)
            || row.diagonalGeometryWeightMeters < 0.0) {
            throw std::invalid_argument(
                "scene fluid pressure-operator row is invalid");
        }
        double diagonal = 0.0;
        std::size_t previousColumn = 0;
        std::size_t previousLink = 0;
        bool hasPrevious = false;
        for (std::size_t entryIndex = row.firstEntry;
             entryIndex < rowEnd; ++entryIndex) {
            const auto& entry = pressureOperator.entries[entryIndex];
            if (entry.entryIndex != entryIndex
                || entry.stableId == 0
                || entry.columnControlVolumeIndex
                    >= pressureOperator.rows.size()
                || entry.columnControlVolumeIndex == rowIndex
                || !std::isfinite(entry.geometryWeightMeters)
                || !(entry.geometryWeightMeters > 0.0)
                || (hasPrevious
                    && (entry.columnControlVolumeIndex < previousColumn
                        || (entry.columnControlVolumeIndex == previousColumn
                            && entry.sourceFaceLinkIndex < previousLink)))) {
                throw std::invalid_argument(
                    "scene fluid pressure-operator entry is invalid");
            }
            diagonal += entry.geometryWeightMeters;
            previousColumn = entry.columnControlVolumeIndex;
            previousLink = entry.sourceFaceLinkIndex;
            hasPrevious = true;
        }
        if (diagonal != row.diagonalGeometryWeightMeters) {
            throw std::invalid_argument(
                "scene fluid pressure-operator row diagonal is invalid");
        }
        expectedFirstEntry = rowEnd;
    }
    if (expectedFirstEntry != pressureOperator.entries.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-operator entry ranges are invalid");
    }
    std::size_t expectedFirstMember = 0;
    for (std::size_t componentIndex = 0;
         componentIndex < pressureOperator.components.size();
         ++componentIndex) {
        const auto& component = pressureOperator.components[componentIndex];
        std::size_t componentEnd = 0;
        if (component.componentIndex != componentIndex
            || component.gaugeRegionId == invalidStableId
            || component.gaugeControlVolumeIndex
                >= pressureOperator.rows.size()
            || component.firstControlVolumeMember != expectedFirstMember
            || component.controlVolumeCount == 0
            || !checkedAdd(component.firstControlVolumeMember,
                           component.controlVolumeCount, componentEnd)
            || componentEnd
                > pressureOperator.componentControlVolumeIndices.size()
            || !std::isfinite(component.totalVolumeCubicMeters)
            || !(component.totalVolumeCubicMeters > 0.0)
            || !std::isfinite(component.totalGeometryWeightMeters)
            || component.totalGeometryWeightMeters < 0.0
            || !pressureOperator.rows[
                    component.gaugeControlVolumeIndex].isGauge
            || pressureOperator.rows[
                    component.gaugeControlVolumeIndex].componentIndex
                != componentIndex) {
            throw std::invalid_argument(
                "scene fluid pressure-operator component is invalid");
        }
        for (std::size_t member = component.firstControlVolumeMember;
             member < componentEnd; ++member) {
            const std::size_t rowIndex =
                pressureOperator.componentControlVolumeIndices[member];
            if (rowIndex >= pressureOperator.rows.size()
                || pressureOperator.rows[rowIndex].componentIndex
                    != componentIndex) {
                throw std::invalid_argument(
                    "scene fluid pressure-operator component member is invalid");
            }
        }
        expectedFirstMember = componentEnd;
    }
    if (expectedFirstMember
        != pressureOperator.componentControlVolumeIndices.size()
        || expectedFirstMember != pressureOperator.rows.size()) {
        throw std::invalid_argument(
            "scene fluid pressure-operator component ranges are incomplete");
    }
}

SceneFluidPressureOperator buildOperator(
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperatorLimits& limits) {
    const std::size_t unresolvedTopologyCount =
        faceLinks.unresolvedActiveFaceCount
        + faceLinks.unresolvedCappedFaceCount
        + faceLinks.unresolvedAmbiguousFaceCount
        + faceLinks.unresolvedOpeningFaceCount
        + faceLinks.unresolvedEmbeddedOpeningPatchCount;
    if (unresolvedTopologyCount != 0
        || faceLinks.resolvedFullFaceCount
            + faceLinks.resolvedPartitionFaceCount
            + faceLinks.resolvedOpeningFaceCount
            != faceLinks.faces.size()) {
        throw SceneFluidPressureIncompleteFaceOwnershipError({
            faceLinks.faces.size(),
            faceLinks.resolvedFullFaceCount,
            faceLinks.resolvedPartitionFaceCount,
            faceLinks.resolvedOpeningFaceCount,
            faceLinks.unresolvedActiveFaceCount,
            faceLinks.unresolvedCappedFaceCount,
            faceLinks.unresolvedAmbiguousFaceCount,
            faceLinks.unresolvedOpeningFaceCount,
            faceLinks.unresolvedEmbeddedOpeningPatchCount,
        });
    }
    if (pressureVolumes.fingerprint
            != faceLinks.pressureControlVolumeFingerprint
        || pressureVolumes.surfaceDefinitionFingerprint
            != faceLinks.surfaceDefinitionFingerprint
        || pressureVolumes.surfaceStateFingerprint
            != faceLinks.surfaceStateFingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure-operator source identity is invalid");
    }

    const std::size_t rowCount = pressureVolumes.controlVolumes.size();
    std::size_t entryCount = 0;
    if (!checkedMultiply(faceLinks.links.size(), std::size_t{2}, entryCount)) {
        throw std::length_error(
            "scene fluid pressure-operator entry count overflows");
    }
    const std::size_t componentCount = pressureVolumes.components.size();
    if (rowCount > limits.maximumRows) {
        throw std::length_error(
            "scene fluid pressure operator exceeds its row limit");
    }
    if (entryCount > limits.maximumEntries) {
        throw std::length_error(
            "scene fluid pressure operator exceeds its entry limit");
    }
    if (componentCount > limits.maximumComponents) {
        throw std::length_error(
            "scene fluid pressure operator exceeds its component limit");
    }
    const std::size_t ownedStorageBytes = storageBytesForCounts(
        rowCount, entryCount, componentCount, rowCount);
    if (ownedStorageBytes > limits.maximumOperatorBytes) {
        throw std::length_error(
            "scene fluid pressure operator exceeds its byte limit");
    }

    SceneFluidPressureOperator result;
    result.pressureControlVolumeFingerprint = pressureVolumes.fingerprint;
    result.pressureFaceLinkFingerprint = faceLinks.fingerprint;
    result.surfaceDefinitionFingerprint =
        pressureVolumes.surfaceDefinitionFingerprint;
    result.surfaceStateFingerprint = pressureVolumes.surfaceStateFingerprint;
    result.structureDefinitionFingerprint =
        pressureVolumes.structureDefinitionFingerprint;
    result.acceptedStepCount = pressureVolumes.acceptedStepCount;
    result.simulationTimeSeconds = pressureVolumes.simulationTimeSeconds;
    result.ownedStorageBytes = ownedStorageBytes;
    result.rows.resize(rowCount);
    result.entries.resize(entryCount);
    result.components.reserve(componentCount);
    result.componentControlVolumeIndices.reserve(rowCount);

    std::vector<std::size_t> degrees(rowCount, 0);
    for (const auto& link : faceLinks.links) {
        if (link.minusControlVolumeIndex >= rowCount
            || link.plusControlVolumeIndex >= rowCount
            || link.minusControlVolumeIndex
                == link.plusControlVolumeIndex
            || pressureVolumes.controlVolumes[
                    link.minusControlVolumeIndex].componentIndex
                != pressureVolumes.controlVolumes[
                    link.plusControlVolumeIndex].componentIndex) {
            throw std::invalid_argument(
                "scene fluid pressure face link cannot enter the operator");
        }
        ++degrees[link.minusControlVolumeIndex];
        ++degrees[link.plusControlVolumeIndex];
        result.totalGeometryWeightMeters += link.geometryWeightMeters;
    }

    std::size_t nextEntry = 0;
    for (std::size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        const auto& control = pressureVolumes.controlVolumes[rowIndex];
        auto& row = result.rows[rowIndex];
        row.rowIndex = rowIndex;
        row.stableId = control.stableId;
        row.controlVolumeIndex = rowIndex;
        row.componentIndex = control.componentIndex;
        row.isGauge = pressureVolumes.components[
            control.componentIndex].gaugeControlVolumeIndex == rowIndex;
        row.firstEntry = nextEntry;
        row.entryCount = degrees[rowIndex];
        if (!checkedAdd(nextEntry, row.entryCount, nextEntry)) {
            throw std::length_error(
                "scene fluid pressure-operator row range overflows");
        }
    }
    if (nextEntry != entryCount) {
        throw std::logic_error(
            "scene fluid pressure-operator entry count is inconsistent");
    }

    std::vector<std::size_t> cursors(rowCount, 0);
    for (std::size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        cursors[rowIndex] = result.rows[rowIndex].firstEntry;
    }
    std::unordered_set<std::uint64_t> entryStableIds;
    entryStableIds.reserve(entryCount);
    const auto appendEntry = [&](const SceneFluidPressureFaceLink& link,
                                 const std::size_t sourceIndex,
                                 const std::size_t destinationIndex) {
        auto& entry = result.entries[cursors[sourceIndex]++];
        entry.stableId = entryStableId(
            link, pressureVolumes.controlVolumes[sourceIndex],
            pressureVolumes.controlVolumes[destinationIndex]);
        entry.sourceFaceLinkIndex = link.linkIndex;
        entry.columnControlVolumeIndex = destinationIndex;
        entry.geometryWeightMeters = link.geometryWeightMeters;
        if (!entryStableIds.insert(entry.stableId).second) {
            throw std::invalid_argument(
                "scene fluid pressure-operator entry stable ID collides");
        }
    };
    for (const auto& link : faceLinks.links) {
        appendEntry(link, link.minusControlVolumeIndex,
                    link.plusControlVolumeIndex);
        appendEntry(link, link.plusControlVolumeIndex,
                    link.minusControlVolumeIndex);
    }
    for (auto& row : result.rows) {
        auto first = result.entries.begin()
            + static_cast<std::ptrdiff_t>(row.firstEntry);
        auto last = first + static_cast<std::ptrdiff_t>(row.entryCount);
        std::sort(first, last,
                  [](const SceneFluidPressureOperatorEntry& firstEntry,
                     const SceneFluidPressureOperatorEntry& secondEntry) {
                      if (firstEntry.columnControlVolumeIndex
                          != secondEntry.columnControlVolumeIndex) {
                          return firstEntry.columnControlVolumeIndex
                              < secondEntry.columnControlVolumeIndex;
                      }
                      return firstEntry.sourceFaceLinkIndex
                          < secondEntry.sourceFaceLinkIndex;
                  });
        for (auto entry = first; entry != last; ++entry) {
            entry->entryIndex = static_cast<std::size_t>(
                entry - result.entries.begin());
            row.diagonalGeometryWeightMeters +=
                entry->geometryWeightMeters;
        }
    }

    constexpr std::size_t unassigned =
        std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> graphComponents(rowCount, unassigned);
    std::vector<std::size_t> queue;
    std::size_t graphComponentCount = 0;
    for (std::size_t seed = 0; seed < rowCount; ++seed) {
        if (graphComponents[seed] != unassigned) continue;
        queue.clear();
        queue.push_back(seed);
        graphComponents[seed] = graphComponentCount;
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) {
            const std::size_t rowIndex = queue[cursor];
            const auto& row = result.rows[rowIndex];
            for (std::size_t offset = 0; offset < row.entryCount; ++offset) {
                const std::size_t neighbour = result.entries[
                    row.firstEntry + offset].columnControlVolumeIndex;
                if (result.rows[neighbour].componentIndex
                    != row.componentIndex) {
                    throw std::invalid_argument(
                        "pressure-operator link crosses authored components");
                }
                if (graphComponents[neighbour] == unassigned) {
                    graphComponents[neighbour] = graphComponentCount;
                    queue.push_back(neighbour);
                }
            }
        }
        ++graphComponentCount;
    }
    if (graphComponentCount != componentCount) {
        throw std::invalid_argument(
            "scene fluid pressure component is disconnected in its face graph");
    }
    std::vector<std::size_t> graphForPressureComponent(
        componentCount, unassigned);
    std::vector<std::size_t> pressureComponentForGraph(
        graphComponentCount, unassigned);
    for (std::size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        const std::size_t pressureComponent =
            result.rows[rowIndex].componentIndex;
        const std::size_t graphComponent = graphComponents[rowIndex];
        auto& expectedGraph =
            graphForPressureComponent[pressureComponent];
        auto& expectedPressure =
            pressureComponentForGraph[graphComponent];
        if ((expectedGraph != unassigned
             && expectedGraph != graphComponent)
            || (expectedPressure != unassigned
                && expectedPressure != pressureComponent)) {
            throw std::invalid_argument(
                "scene fluid pressure component does not match its face graph");
        }
        expectedGraph = graphComponent;
        expectedPressure = pressureComponent;
    }

    for (const auto& source : pressureVolumes.components) {
        SceneFluidPressureOperatorComponent component;
        component.componentIndex = source.componentIndex;
        component.gaugeRegionId = source.gaugeRegionId;
        component.gaugeControlVolumeIndex =
            source.gaugeControlVolumeIndex;
        component.firstControlVolumeMember =
            result.componentControlVolumeIndices.size();
        component.controlVolumeCount = source.controlVolumeCount;
        component.totalVolumeCubicMeters = source.totalVolumeCubicMeters;
        for (std::size_t offset = 0;
             offset < source.controlVolumeCount; ++offset) {
            const std::size_t controlIndex =
                pressureVolumes.componentControlVolumeIndices[
                    source.firstControlVolumeMember + offset];
            result.componentControlVolumeIndices.push_back(controlIndex);
            component.totalGeometryWeightMeters +=
                0.5 * result.rows[controlIndex]
                    .diagonalGeometryWeightMeters;
        }
        result.components.push_back(component);
    }
    result.fingerprint = operatorFingerprint(result);
    validatePayload(result);
    return result;
}

void validateSolveSettings(
    const SceneFluidPressureSolveSettings& settings) {
    if (!std::isfinite(
            settings.absoluteResidualTolerancePascalsMeters)
        || settings.absoluteResidualTolerancePascalsMeters < 0.0
        || !std::isfinite(settings.relativeResidualTolerance)
        || settings.relativeResidualTolerance < 0.0
        || (settings.absoluteResidualTolerancePascalsMeters == 0.0
            && settings.relativeResidualTolerance == 0.0)
        || !std::isfinite(
            settings.absoluteComponentCompatibilityTolerancePascalsMeters)
        || settings.absoluteComponentCompatibilityTolerancePascalsMeters
            < 0.0) {
        throw std::invalid_argument(
            "scene fluid pressure-solve settings are invalid");
    }
}

bool applyOperatorUnchecked(
    const SceneFluidPressureOperator& pressureOperator,
    const std::span<const double> pressureValues,
    std::vector<double>& result) {
    result.assign(pressureOperator.rows.size(), 0.0);
    for (const auto& row : pressureOperator.rows) {
        for (std::size_t offset = 0; offset < row.entryCount; ++offset) {
            const auto& entry = pressureOperator.entries[
                row.firstEntry + offset];
            result[row.rowIndex] += entry.geometryWeightMeters
                * (pressureValues[row.rowIndex]
                   - pressureValues[entry.columnControlVolumeIndex]);
        }
        if (!std::isfinite(result[row.rowIndex])) return false;
    }
    return true;
}

double vectorDot(const std::span<const double> first,
                 const std::span<const double> second) {
    double result = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        result += first[index] * second[index];
    }
    return result;
}

double vectorL2(const std::span<const double> values) {
    return std::sqrt(
        vectorDot(values, values) / static_cast<double>(values.size()));
}

double vectorMaximumAbsolute(const std::span<const double> values) {
    double result = 0.0;
    for (const double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

void subtractComponentMeans(
    const SceneFluidPressureOperator& pressureOperator,
    std::vector<double>& values) {
    for (const auto& component : pressureOperator.components) {
        double sum = 0.0;
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            sum += values[pressureOperator.componentControlVolumeIndices[
                component.firstControlVolumeMember + offset]];
        }
        const double mean = sum
            / static_cast<double>(component.controlVolumeCount);
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            values[pressureOperator.componentControlVolumeIndices[
                component.firstControlVolumeMember + offset]] -= mean;
        }
    }
}

void shiftComponentGauges(
    const SceneFluidPressureOperator& pressureOperator,
    std::vector<double>& pressureValues) {
    for (const auto& component : pressureOperator.components) {
        const double gauge = pressureValues[
            component.gaugeControlVolumeIndex];
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            pressureValues[pressureOperator.componentControlVolumeIndices[
                component.firstControlVolumeMember + offset]] -= gauge;
        }
        pressureValues[component.gaugeControlVolumeIndex] = 0.0;
    }
}

void validateOperatorInputs(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks) {
    validateSceneFluidGridEpoch(epoch, surface, state, grid, transfer);
    validateSceneFluidOpeningGridPatches(
        openingPatches, surface, state, caps, openingQuadrature, grid);
    validateSceneFluidCellVolumes(
        volumes, surface, state, grid, transfer, epoch);
    validateSceneFluidPressureControlVolumes(
        pressureVolumes, surface, volumes, connectivity);
    validateSceneFluidPressureFaceLinkIntegrity(faceLinks);
    if (faceLinks.surfaceDefinitionFingerprint != surface.fingerprint
        || faceLinks.surfaceStateFingerprint != state.fingerprint
        || faceLinks.gridEpochFingerprint != epoch.fingerprint
        || faceLinks.openingPatchFingerprint != openingPatches.fingerprint
        || faceLinks.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure-operator face-link identity is invalid");
    }
}

} // namespace

SceneFluidPressureOperator buildSceneFluidPressureOperator(
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks,
    const SceneFluidPressureOperatorLimits& limits) {
    validateOperatorInputs(
        surface, state, grid, transfer, epoch, caps, openingQuadrature,
        openingPatches, volumes, connectivity, pressureVolumes, faceLinks);
    auto result = buildOperator(pressureVolumes, faceLinks, limits);
    validateSceneFluidPressureOperator(
        result, surface, state, grid, transfer, epoch, caps,
        openingQuadrature, openingPatches, volumes, connectivity,
        pressureVolumes, faceLinks);
    return result;
}

void validateSceneFluidPressureOperator(
    const SceneFluidPressureOperator& pressureOperator,
    const SceneFluidSurfaceDefinition& surface,
    const SceneFluidSurfaceState& state,
    const fluid::PeriodicCartesianGrid& grid,
    const SceneFluidSurfaceTransfer& transfer,
    const SceneFluidGridEpoch& epoch,
    const SceneFluidOpeningCapSet& caps,
    const SceneFluidOpeningQuadratureSet& openingQuadrature,
    const SceneFluidOpeningGridPatchSet& openingPatches,
    const SceneFluidCellVolumeSet& volumes,
    const SceneFluidRegionConnectivity& connectivity,
    const SceneFluidPressureControlVolumeSet& pressureVolumes,
    const SceneFluidPressureFaceLinkSet& faceLinks) {
    validateOperatorInputs(
        surface, state, grid, transfer, epoch, caps, openingQuadrature,
        openingPatches, volumes, connectivity, pressureVolumes, faceLinks);
    validatePayload(pressureOperator);
    if (pressureOperator.pressureControlVolumeFingerprint
            != pressureVolumes.fingerprint
        || pressureOperator.pressureFaceLinkFingerprint
            != faceLinks.fingerprint) {
        throw std::invalid_argument(
            "scene fluid pressure-operator source identity is invalid");
    }
    const SceneFluidPressureOperatorLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildOperator(
        pressureVolumes, faceLinks, unlimited);
    if (pressureOperator != expected) {
        throw std::invalid_argument(
            "scene fluid pressure-operator payload is invalid");
    }
}

std::vector<double> applySceneFluidPressureOperator(
    const SceneFluidPressureOperator& pressureOperator,
    const std::span<const double> pressureValues) {
    validatePayload(pressureOperator);
    if (pressureValues.size() != pressureOperator.rows.size()
        || !std::ranges::all_of(pressureValues, [](const double value) {
               return std::isfinite(value);
           })) {
        throw std::invalid_argument(
            "scene fluid pressure-operator input is invalid");
    }
    std::vector<double> result;
    if (!applyOperatorUnchecked(
            pressureOperator, pressureValues, result)) {
        throw std::overflow_error(
            "scene fluid pressure-operator application overflowed");
    }
    return result;
}

SceneFluidPressureSolveDiagnostics solveSceneFluidPressureSystem(
    const SceneFluidPressureOperator& pressureOperator,
    const std::span<const double> integratedRightHandSidePascalsMeters,
    std::vector<double>& pressurePascals,
    const SceneFluidPressureSolveSettings& settings) {
    validatePayload(pressureOperator);
    validateSolveSettings(settings);
    if (integratedRightHandSidePascalsMeters.size()
            != pressureOperator.rows.size()
        || pressurePascals.size() != pressureOperator.rows.size()
        || !std::ranges::all_of(
            integratedRightHandSidePascalsMeters,
            [](const double value) { return std::isfinite(value); })
        || !std::ranges::all_of(
            pressurePascals,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "scene fluid pressure-solve fields are invalid");
    }

    SceneFluidPressureSolveDiagnostics diagnostics;
    diagnostics.finite = true;
    diagnostics.pressureOperatorFingerprint = pressureOperator.fingerprint;
    diagnostics.rowCount = pressureOperator.rows.size();
    diagnostics.componentCount = pressureOperator.components.size();
    diagnostics.components.reserve(pressureOperator.components.size());
    std::vector<double> rightHandSide(
        integratedRightHandSidePascalsMeters.begin(),
        integratedRightHandSidePascalsMeters.end());
    diagnostics.compatible = true;
    for (const auto& component : pressureOperator.components) {
        SceneFluidPressureSolveComponentDiagnostics componentDiagnostics;
        componentDiagnostics.componentIndex = component.componentIndex;
        componentDiagnostics.controlVolumeCount =
            component.controlVolumeCount;
        componentDiagnostics.gaugeControlVolumeIndex =
            component.gaugeControlVolumeIndex;
        componentDiagnostics.pressureGaugeBeforePascals =
            pressurePascals[component.gaugeControlVolumeIndex];
        componentDiagnostics.pressureGaugeAfterPascals =
            componentDiagnostics.pressureGaugeBeforePascals;
        for (std::size_t offset = 0;
             offset < component.controlVolumeCount; ++offset) {
            componentDiagnostics.rightHandSideSumPascalsMeters +=
                rightHandSide[
                    pressureOperator.componentControlVolumeIndices[
                        component.firstControlVolumeMember + offset]];
        }
        diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters =
            std::max(
                diagnostics.maximumAbsoluteComponentCompatibilityPascalsMeters,
                std::abs(componentDiagnostics
                    .rightHandSideSumPascalsMeters));
        if (!std::isfinite(
                componentDiagnostics.rightHandSideSumPascalsMeters)) {
            diagnostics.finite = false;
            diagnostics.compatible = false;
        } else if (std::abs(
                       componentDiagnostics.rightHandSideSumPascalsMeters)
                   > settings
                       .absoluteComponentCompatibilityTolerancePascalsMeters) {
            diagnostics.compatible = false;
        }
        componentDiagnostics.compatibilityCorrectionPascalsMeters =
            componentDiagnostics.rightHandSideSumPascalsMeters
            / static_cast<double>(component.controlVolumeCount);
        diagnostics.components.push_back(componentDiagnostics);
    }
    if (!diagnostics.compatible) return diagnostics;

    for (const auto& component : diagnostics.components) {
        const auto& owner = pressureOperator.components[
            component.componentIndex];
        for (std::size_t offset = 0;
             offset < owner.controlVolumeCount; ++offset) {
            rightHandSide[
                pressureOperator.componentControlVolumeIndices[
                    owner.firstControlVolumeMember + offset]]
                -= component.compatibilityCorrectionPascalsMeters;
        }
    }

    std::vector<double> candidatePressure = pressurePascals;
    shiftComponentGauges(pressureOperator, candidatePressure);
    std::vector<double> operatorPressure;
    if (!applyOperatorUnchecked(
            pressureOperator, candidatePressure, operatorPressure)) {
        diagnostics.finite = false;
        return diagnostics;
    }
    std::vector<double> residual(pressureOperator.rows.size(), 0.0);
    for (std::size_t row = 0; row < residual.size(); ++row) {
        residual[row] = rightHandSide[row] - operatorPressure[row];
    }
    subtractComponentMeans(pressureOperator, residual);
    if (!std::ranges::all_of(
            residual,
            [](const double value) { return std::isfinite(value); })) {
        diagnostics.finite = false;
        return diagnostics;
    }
    diagnostics.initialResidualL2PascalsMeters = vectorL2(residual);
    diagnostics.finalResidualL2PascalsMeters =
        diagnostics.initialResidualL2PascalsMeters;
    diagnostics.finalResidualMaximumPascalsMeters =
        vectorMaximumAbsolute(residual);
    if (!std::isfinite(diagnostics.initialResidualL2PascalsMeters)
        || !std::isfinite(
            diagnostics.finalResidualMaximumPascalsMeters)) {
        diagnostics.finite = false;
        return diagnostics;
    }
    const double convergenceThreshold = std::max(
        settings.absoluteResidualTolerancePascalsMeters,
        settings.relativeResidualTolerance
            * diagnostics.initialResidualL2PascalsMeters);

    std::vector<double> direction = residual;
    std::vector<double> operatorDirection;
    double residualSquared = vectorDot(residual, residual);
    diagnostics.converged =
        diagnostics.finalResidualL2PascalsMeters <= convergenceThreshold;
    while (!diagnostics.converged
           && diagnostics.iterationCount < settings.maximumIterations) {
        if (!applyOperatorUnchecked(
                pressureOperator, direction, operatorDirection)) {
            diagnostics.finite = false;
            break;
        }
        const double denominator = vectorDot(direction, operatorDirection);
        if (!std::isfinite(denominator) || !(denominator > 0.0)
            || !std::isfinite(residualSquared)
            || !(residualSquared > 0.0)) {
            diagnostics.finite = false;
            break;
        }
        const double alpha = residualSquared / denominator;
        if (!std::isfinite(alpha)) {
            diagnostics.finite = false;
            break;
        }
        for (std::size_t row = 0; row < residual.size(); ++row) {
            candidatePressure[row] += alpha * direction[row];
            residual[row] -= alpha * operatorDirection[row];
        }
        ++diagnostics.iterationCount;
        shiftComponentGauges(pressureOperator, candidatePressure);
        subtractComponentMeans(pressureOperator, residual);
        const double nextResidualSquared = vectorDot(residual, residual);
        diagnostics.finalResidualL2PascalsMeters = std::sqrt(
            nextResidualSquared / static_cast<double>(residual.size()));
        diagnostics.finalResidualMaximumPascalsMeters =
            vectorMaximumAbsolute(residual);
        diagnostics.converged = std::isfinite(nextResidualSquared)
            && std::isfinite(diagnostics.finalResidualL2PascalsMeters)
            && diagnostics.finalResidualL2PascalsMeters
                <= convergenceThreshold;
        if (diagnostics.converged) {
            residualSquared = nextResidualSquared;
            break;
        }
        if (!std::isfinite(nextResidualSquared)
            || !(residualSquared > 0.0)) {
            diagnostics.finite = false;
            break;
        }
        const double beta = nextResidualSquared / residualSquared;
        if (!std::isfinite(beta)) {
            diagnostics.finite = false;
            break;
        }
        for (std::size_t row = 0; row < residual.size(); ++row) {
            direction[row] = residual[row] + beta * direction[row];
        }
        subtractComponentMeans(pressureOperator, direction);
        residualSquared = nextResidualSquared;
    }

    if (diagnostics.converged) {
        if (!applyOperatorUnchecked(
                pressureOperator, candidatePressure, operatorPressure)) {
            diagnostics.finite = false;
            diagnostics.converged = false;
        } else {
            for (std::size_t row = 0; row < residual.size(); ++row) {
                residual[row] = rightHandSide[row] - operatorPressure[row];
            }
            subtractComponentMeans(pressureOperator, residual);
            diagnostics.finalResidualL2PascalsMeters = vectorL2(residual);
            diagnostics.finalResidualMaximumPascalsMeters =
                vectorMaximumAbsolute(residual);
            diagnostics.converged = std::isfinite(
                diagnostics.finalResidualL2PascalsMeters)
                && diagnostics.finalResidualL2PascalsMeters
                    <= convergenceThreshold;
        }
    }
    if (!diagnostics.converged) return diagnostics;

    shiftComponentGauges(pressureOperator, candidatePressure);
    for (auto& component : diagnostics.components) {
        component.pressureGaugeAfterPascals = candidatePressure[
            component.gaugeControlVolumeIndex];
    }
    pressurePascals = std::move(candidatePressure);
    return diagnostics;
}

} // namespace simwing::fsi
