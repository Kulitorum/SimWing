#include "scene_fluid_mimetic_trace_system.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {

SceneFluidMimeticTraceLocalCellLinearConsistencyError::
    SceneFluidMimeticTraceLocalCellLinearConsistencyError(
        SceneFluidMimeticTraceLocalCellLinearConsistencyFailure diagnostics)
    : std::invalid_argument(
        "scene fluid mimetic trace local cell failed linear consistency"),
      diagnostics_(std::move(diagnostics)) {}

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

    template<typename Enumeration>
    void enumeration(const Enumeration value) {
        using Underlying = std::underlying_type_t<Enumeration>;
        using Unsigned = std::make_unsigned_t<Underlying>;
        integer(static_cast<Unsigned>(value));
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

std::uint64_t materialWallTraceStableId(
    const std::uint64_t halfFaceStableId) {
    Fingerprint fingerprint;
    fingerprint.integer(std::uint64_t{0x6d696d6577616c6cULL});
    fingerprint.integer(halfFaceStableId);
    return fingerprint.value();
}

void fingerprintSettings(
    Fingerprint& fingerprint,
    const SceneFluidMimeticTraceSystemSettings& settings) {
    const auto& local = settings.localCell;
    for (const double value : {
             local.absoluteAreaClosureToleranceSquareMeters,
             local.absoluteDivergenceTheoremToleranceCubicMeters,
             local.relativeGeometryTolerance,
             local.unitNormalTolerance,
             local.algebraicConsistencyTolerance}) {
        fingerprint.real(value);
    }
    fingerprint.integer(
        static_cast<std::uint64_t>(local.maximumHalfFaces));
    fingerprint.integer(
        static_cast<std::uint64_t>(local.maximumOperatorBytes));
}

std::uint64_t productFingerprint(
    const SceneFluidMimeticTraceSystem& system) {
    Fingerprint fingerprint;
    fingerprint.integer(system.version);
    fingerprint.integer(system.mimeticControlCellFingerprint);
    fingerprint.integer(system.structureDefinitionFingerprint);
    fingerprint.integer(system.acceptedStepCount);
    fingerprint.real(system.simulationTimeSeconds);
    fingerprintSettings(fingerprint, system.settings);
    for (const std::size_t value : {
             system.ownedStorageBytes,
             system.localOperatorStorageBytes,
             system.componentCount,
             system.sharedTraceCount,
             system.materialWallTraceCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.real(system.minimumPositiveOperatorDiagonal);
    fingerprint.real(system.maximumOperatorDiagonal);
    fingerprint.integer(static_cast<std::uint64_t>(system.traces.size()));
    for (const auto& trace : system.traces) {
        fingerprint.integer(static_cast<std::uint64_t>(trace.traceIndex));
        fingerprint.integer(trace.stableId);
        fingerprint.enumeration(trace.kind);
        fingerprint.integer(trace.sourceStableId);
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.componentIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.firstIncidence));
        fingerprint.integer(static_cast<std::uint64_t>(
            trace.incidenceCount));
        fingerprint.real(trace.operatorDiagonal);
        fingerprint.integer(static_cast<std::uint8_t>(trace.isGauge));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        system.incidences.size()));
    for (const auto& incidence : system.incidences) {
        fingerprint.integer(static_cast<std::uint64_t>(
            incidence.incidenceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            incidence.traceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            incidence.halfFaceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            incidence.controlCellIndex));
        fingerprint.integer(static_cast<std::uint64_t>(
            incidence.localHalfFaceIndex));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        system.halfFaceTraceIndices.size()));
    for (const std::size_t value : system.halfFaceTraceIndices) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        system.componentGaugeTraceIndices.size()));
    for (const std::size_t value : system.componentGaugeTraceIndices) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        system.localOperators.size()));
    for (const auto& local : system.localOperators) {
        fingerprint.integer(local.fingerprint);
    }
    return fingerprint.value();
}

std::size_t storageBytes(
    const SceneFluidMimeticTraceSystem& system) {
    std::size_t total = 0;
    const auto addVector = [&](const std::size_t count,
                               const std::size_t elementSize) {
        std::size_t bytes = 0;
        if (!checkedMultiply(count, elementSize, bytes)
            || !checkedAdd(total, bytes, total)) {
            throw std::length_error(
                "scene fluid mimetic trace-system storage overflows");
        }
    };
    addVector(system.traces.size(), sizeof(SceneFluidMimeticTrace));
    addVector(system.incidences.size(),
              sizeof(SceneFluidMimeticTraceIncidence));
    addVector(system.halfFaceTraceIndices.size(), sizeof(std::size_t));
    addVector(system.componentGaugeTraceIndices.size(), sizeof(std::size_t));
    addVector(system.localOperators.size(),
              sizeof(fluid::MimeticLocalCellOperator));
    if (!checkedAdd(total, system.localOperatorStorageBytes, total)) {
        throw std::length_error(
            "scene fluid mimetic trace-system storage overflows");
    }
    return total;
}

struct TraceAssembly {
    SceneFluidMimeticHalfFaceKind kind =
        SceneFluidMimeticHalfFaceKind::CartesianTrace;
    std::uint64_t sourceStableId = 0;
    std::size_t componentIndex = 0;
    std::vector<std::size_t> halfFaceIndices;
};

class DisjointControlCells final {
public:
    explicit DisjointControlCells(const std::size_t count)
        : parents_(count) {
        for (std::size_t index = 0; index < count; ++index) {
            parents_[index] = index;
        }
    }

    std::size_t root(std::size_t index) {
        while (parents_[index] != index) {
            parents_[index] = parents_[parents_[index]];
            index = parents_[index];
        }
        return index;
    }

    void join(const std::size_t first, const std::size_t second) {
        const std::size_t firstRoot = root(first);
        const std::size_t secondRoot = root(second);
        if (firstRoot != secondRoot) parents_[secondRoot] = firstRoot;
    }

private:
    std::vector<std::size_t> parents_;
};

double localOperatorDiagonal(
    const fluid::MimeticLocalCellOperator& localOperator,
    const std::size_t face,
    const std::vector<double>& inverseFluxDiagonal,
    const std::vector<double>& inverseFluxArea,
    const double conservationDenominator) {
    const double area = localOperator.faceAreasSquareMeters[face];
    const double coupling = area * inverseFluxArea[face];
    double result = area * area * inverseFluxDiagonal[face]
        - coupling * coupling / conservationDenominator;
    const double scale = std::max(
        area * area * inverseFluxDiagonal[face],
        coupling * coupling / conservationDenominator);
    const double tolerance = 512.0
        * std::numeric_limits<double>::epsilon() * scale;
    if (result < 0.0 && result >= -tolerance) result = 0.0;
    if (!std::isfinite(result) || result < 0.0) {
        throw std::invalid_argument(
            "scene fluid mimetic trace diagonal is invalid");
    }
    return result;
}

SceneFluidMimeticTraceSystem buildSystem(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystemSettings& settings,
    const SceneFluidMimeticTraceSystemLimits& limits) {
    validateSceneFluidMimeticControlCellIntegrity(controlCells);
    const std::size_t controlCount = controlCells.controlCells.size();
    const std::size_t halfFaceCount = controlCells.halfFaces.size();
    if (controlCount == 0
        || controlCells.readyControlCellCount != controlCount
        || controlCells.incompleteTopologyControlCellCount != 0
        || controlCells.nonclosingControlCellCount != 0) {
        throw std::invalid_argument(
            "scene fluid mimetic trace system requires complete closed cells");
    }
    if (controlCount > limits.maximumLocalOperators) {
        throw std::length_error(
            "scene fluid mimetic trace-system local-operator limit exceeded");
    }
    if (halfFaceCount > limits.maximumIncidences) {
        throw std::length_error(
            "scene fluid mimetic trace-system incidence limit exceeded");
    }

    std::size_t componentCount = 0;
    for (const auto& cell : controlCells.controlCells) {
        if (cell.controlCellIndex >= controlCount
            || cell.controlVolumeIndex != cell.controlCellIndex
            || cell.componentIndex >= controlCount
            || !cell.readyForMimeticOperator) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-system control cell is invalid");
        }
        componentCount = std::max(componentCount, cell.componentIndex + 1);
    }
    std::vector<bool> componentSeen(componentCount, false);
    for (const auto& cell : controlCells.controlCells) {
        componentSeen[cell.componentIndex] = true;
    }
    if (componentCount == 0
        || std::ranges::find(componentSeen, false) != componentSeen.end()) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-system components are not contiguous");
    }

    std::map<std::uint64_t, TraceAssembly> traceAssemblies;
    for (const auto& halfFace : controlCells.halfFaces) {
        if (halfFace.controlVolumeIndex >= controlCount) {
            throw std::invalid_argument(
                "scene fluid mimetic trace-system half-face owner is invalid");
        }
        const auto& cell =
            controlCells.controlCells[halfFace.controlVolumeIndex];
        const bool materialWall = halfFace.kind
            == SceneFluidMimeticHalfFaceKind::MaterialWall;
        const std::uint64_t traceId = materialWall
            ? materialWallTraceStableId(halfFace.stableId)
            : halfFace.traceStableId;
        auto [iterator, inserted] = traceAssemblies.try_emplace(traceId);
        if (inserted && traceAssemblies.size() > limits.maximumTraces) {
            throw std::length_error(
                "scene fluid mimetic trace-system trace limit exceeded");
        }
        auto& trace = iterator->second;
        if (inserted) {
            trace.kind = halfFace.kind;
            trace.sourceStableId = halfFace.sourceStableId;
            trace.componentIndex = cell.componentIndex;
        } else if (materialWall
                   || trace.kind != halfFace.kind
                   || trace.sourceStableId != halfFace.sourceStableId
                   || trace.componentIndex != cell.componentIndex) {
            throw std::invalid_argument(
                "scene fluid mimetic trace stable ID collides");
        }
        trace.halfFaceIndices.push_back(halfFace.halfFaceIndex);
    }
    SceneFluidMimeticTraceSystem result;
    result.mimeticControlCellFingerprint = controlCells.fingerprint;
    result.structureDefinitionFingerprint =
        controlCells.structureDefinitionFingerprint;
    result.acceptedStepCount = controlCells.acceptedStepCount;
    result.simulationTimeSeconds = controlCells.simulationTimeSeconds;
    result.settings = settings;
    result.componentCount = componentCount;
    result.traces.reserve(traceAssemblies.size());
    result.incidences.reserve(halfFaceCount);
    result.halfFaceTraceIndices.assign(
        halfFaceCount, invalidSceneFluidMimeticTraceIndex);
    result.componentGaugeTraceIndices.assign(
        componentCount, invalidSceneFluidMimeticTraceIndex);

    for (const auto& [stableId, source] : traceAssemblies) {
        const bool wall = source.kind
            == SceneFluidMimeticHalfFaceKind::MaterialWall;
        if ((wall && source.halfFaceIndices.size() != 1)
            || (!wall && source.halfFaceIndices.size() != 2)) {
            throw std::invalid_argument(
                "scene fluid mimetic trace incidence count is invalid");
        }
        SceneFluidMimeticTrace trace;
        trace.traceIndex = result.traces.size();
        trace.stableId = stableId;
        trace.kind = source.kind;
        trace.sourceStableId = source.sourceStableId;
        trace.componentIndex = source.componentIndex;
        trace.firstIncidence = result.incidences.size();
        trace.incidenceCount = source.halfFaceIndices.size();
        if (wall) {
            ++result.materialWallTraceCount;
        } else {
            ++result.sharedTraceCount;
        }
        for (const std::size_t halfFaceIndex : source.halfFaceIndices) {
            if (halfFaceIndex >= halfFaceCount
                || result.halfFaceTraceIndices[halfFaceIndex]
                    != invalidSceneFluidMimeticTraceIndex) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace incidence is duplicated");
            }
            const auto& halfFace = controlCells.halfFaces[halfFaceIndex];
            const auto& cell = controlCells.controlCells[
                halfFace.controlVolumeIndex];
            if (halfFaceIndex < cell.firstHalfFace
                || halfFaceIndex >= cell.firstHalfFace + cell.halfFaceCount) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace local incidence is invalid");
            }
            SceneFluidMimeticTraceIncidence incidence;
            incidence.incidenceIndex = result.incidences.size();
            incidence.traceIndex = trace.traceIndex;
            incidence.halfFaceIndex = halfFaceIndex;
            incidence.controlCellIndex = cell.controlCellIndex;
            incidence.localHalfFaceIndex =
                halfFaceIndex - cell.firstHalfFace;
            result.halfFaceTraceIndices[halfFaceIndex] = trace.traceIndex;
            result.incidences.push_back(incidence);
        }
        if (!wall) {
            const auto& first = controlCells.halfFaces[
                source.halfFaceIndices[0]];
            const auto& second = controlCells.halfFaces[
                source.halfFaceIndices[1]];
            if (first.otherControlVolumeIndex
                    != second.controlVolumeIndex
                || second.otherControlVolumeIndex
                    != first.controlVolumeIndex
                || first.areaSquareMeters != second.areaSquareMeters
                || first.outwardUnitNormal.x
                    != -second.outwardUnitNormal.x
                || first.outwardUnitNormal.y
                    != -second.outwardUnitNormal.y
                || first.outwardUnitNormal.z
                    != -second.outwardUnitNormal.z) {
                throw std::invalid_argument(
                    "scene fluid mimetic shared trace does not pair exactly");
            }
        }
        result.traces.push_back(trace);
    }
    if (std::ranges::find(
            result.halfFaceTraceIndices,
            invalidSceneFluidMimeticTraceIndex)
        != result.halfFaceTraceIndices.end()) {
        throw std::logic_error(
            "scene fluid mimetic trace mapping is incomplete");
    }

    DisjointControlCells connectedControls(controlCount);
    for (const auto& trace : result.traces) {
        if (trace.incidenceCount == 2) {
            connectedControls.join(
                result.incidences[trace.firstIncidence].controlCellIndex,
                result.incidences[trace.firstIncidence + 1]
                    .controlCellIndex);
        }
    }
    std::vector<std::size_t> componentRoots(
        componentCount, invalidSceneFluidMimeticTraceIndex);
    for (const auto& cell : controlCells.controlCells) {
        const std::size_t root = connectedControls.root(
            cell.controlCellIndex);
        auto& componentRoot = componentRoots[cell.componentIndex];
        if (componentRoot == invalidSceneFluidMimeticTraceIndex) {
            componentRoot = root;
        } else if (componentRoot != root) {
            throw std::invalid_argument(
                "scene fluid mimetic trace component is disconnected");
        }
    }

    result.localOperators.reserve(controlCount);
    for (const auto& cell : controlCells.controlCells) {
        fluid::MimeticLocalCellGeometry geometry;
        geometry.volumeCubicMeters = cell.volumeCubicMeters;
        geometry.centroidMeters = cell.centroidMeters;
        geometry.halfFaces.reserve(cell.halfFaceCount);
        for (std::size_t offset = 0; offset < cell.halfFaceCount; ++offset) {
            const auto& halfFace = controlCells.halfFaces[
                cell.firstHalfFace + offset];
            geometry.halfFaces.push_back({
                halfFace.areaSquareMeters,
                halfFace.centroidMeters,
                halfFace.outwardUnitNormal,
            });
        }
        fluid::MimeticLocalCellOperator local;
        try {
            local = fluid::buildMimeticLocalCellOperator(
                geometry, settings.localCell);
        } catch (
            const fluid::MimeticLocalCellLinearConsistencyError& error) {
            throw SceneFluidMimeticTraceLocalCellLinearConsistencyError({
                cell.controlCellIndex,
                cell.stableId,
                cell.cellIndex,
                cell.regionId,
                error.diagnostics().maximumAlgebraicConsistencyError,
                error.diagnostics().algebraicConsistencyTolerance,
            });
        }
        if (!checkedAdd(result.localOperatorStorageBytes,
                        local.ownedStorageBytes,
                        result.localOperatorStorageBytes)
            || result.localOperatorStorageBytes
                > limits.maximumLocalOperatorBytes) {
            throw std::length_error(
                "scene fluid mimetic trace-system local storage limit exceeded");
        }
        const auto inverseFluxArea =
            fluid::applyMimeticInverseFluxInnerProduct(
                local, local.faceAreasSquareMeters);
        const auto inverseFluxDiagonal =
            fluid::mimeticInverseFluxInnerProductDiagonal(local);
        double conservationDenominator = 0.0;
        for (std::size_t face = 0; face < local.halfFaceCount; ++face) {
            conservationDenominator +=
                local.faceAreasSquareMeters[face] * inverseFluxArea[face];
        }
        if (!std::isfinite(conservationDenominator)
            || conservationDenominator <= 0.0) {
            throw std::invalid_argument(
                "scene fluid mimetic trace conservation denominator is invalid");
        }
        for (std::size_t face = 0; face < local.halfFaceCount; ++face) {
            const std::size_t traceIndex = result.halfFaceTraceIndices[
                cell.firstHalfFace + face];
            result.traces[traceIndex].operatorDiagonal +=
                localOperatorDiagonal(
                    local, face, inverseFluxDiagonal, inverseFluxArea,
                    conservationDenominator);
        }
        result.localOperators.push_back(std::move(local));
    }

    result.minimumPositiveOperatorDiagonal =
        std::numeric_limits<double>::infinity();
    for (auto& trace : result.traces) {
        if (!std::isfinite(trace.operatorDiagonal)
            || trace.operatorDiagonal <= 0.0) {
            throw std::invalid_argument(
                "scene fluid mimetic global trace has no positive diagonal");
        }
        result.minimumPositiveOperatorDiagonal = std::min(
            result.minimumPositiveOperatorDiagonal,
            trace.operatorDiagonal);
        result.maximumOperatorDiagonal = std::max(
            result.maximumOperatorDiagonal, trace.operatorDiagonal);
        auto& gauge = result.componentGaugeTraceIndices[
            trace.componentIndex];
        if (gauge == invalidSceneFluidMimeticTraceIndex) {
            gauge = trace.traceIndex;
            trace.isGauge = true;
        }
    }
    if (!std::isfinite(result.minimumPositiveOperatorDiagonal)
        || result.maximumOperatorDiagonal <= 0.0
        || std::ranges::find(
            result.componentGaugeTraceIndices,
            invalidSceneFluidMimeticTraceIndex)
            != result.componentGaugeTraceIndices.end()) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-system gauge ownership is invalid");
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid mimetic trace-system byte limit exceeded");
    }
    result.fingerprint = productFingerprint(result);
    return result;
}

} // namespace

SceneFluidMimeticTraceSystem buildSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticControlCellSet& controlCells,
    const SceneFluidMimeticTraceSystemSettings& settings,
    const SceneFluidMimeticTraceSystemLimits& limits) {
    auto result = buildSystem(controlCells, settings, limits);
    validateSceneFluidMimeticTraceSystem(result, controlCells);
    return result;
}

void validateSceneFluidMimeticTraceSystemIntegrity(
    const SceneFluidMimeticTraceSystem& system) {
    if (system.version != sceneFluidMimeticTraceSystemVersion
        || system.fingerprint == 0
        || system.mimeticControlCellFingerprint == 0
        || system.structureDefinitionFingerprint == 0
        || system.componentCount == 0
        || system.traces.empty()
        || system.incidences.empty()
        || system.localOperators.empty()
        || system.incidences.size()
            != system.halfFaceTraceIndices.size()
        || system.localOperators.size()
            > system.incidences.size()
        || system.componentGaugeTraceIndices.size()
            != system.componentCount) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-system integrity is invalid");
    }
    std::size_t localBytes = 0;
    for (const auto& local : system.localOperators) {
        fluid::validateMimeticLocalCellOperator(local);
        if (!checkedAdd(localBytes, local.ownedStorageBytes, localBytes)) {
            throw std::length_error(
                "scene fluid mimetic trace-system local storage overflows");
        }
    }
    std::size_t sharedCount = 0;
    std::size_t wallCount = 0;
    double minimumDiagonal = std::numeric_limits<double>::infinity();
    double maximumDiagonal = 0.0;
    std::vector<bool> gaugeSeen(system.componentCount, false);
    std::uint64_t previousStableId = 0;
    for (const auto& trace : system.traces) {
        const bool wall = trace.kind
            == SceneFluidMimeticHalfFaceKind::MaterialWall;
        if (trace.traceIndex >= system.traces.size()
            || trace.traceIndex != &trace - system.traces.data()
            || trace.stableId == 0
            || (trace.traceIndex != 0
                && trace.stableId <= previousStableId)
            || trace.sourceStableId == 0
            || trace.componentIndex >= system.componentCount
            || trace.firstIncidence > system.incidences.size()
            || trace.incidenceCount
                > system.incidences.size() - trace.firstIncidence
            || (wall ? trace.incidenceCount != 1
                     : trace.incidenceCount != 2)
            || !std::isfinite(trace.operatorDiagonal)
            || trace.operatorDiagonal <= 0.0
            || (trace.isGauge && gaugeSeen[trace.componentIndex])) {
            throw std::invalid_argument(
                "scene fluid mimetic trace record is invalid");
        }
        previousStableId = trace.stableId;
        if (wall) ++wallCount;
        else ++sharedCount;
        if (trace.isGauge) gaugeSeen[trace.componentIndex] = true;
        minimumDiagonal = std::min(
            minimumDiagonal, trace.operatorDiagonal);
        maximumDiagonal = std::max(
            maximumDiagonal, trace.operatorDiagonal);
        for (std::size_t offset = 0;
             offset < trace.incidenceCount; ++offset) {
            const auto& incidence = system.incidences[
                trace.firstIncidence + offset];
            if (incidence.incidenceIndex
                    != trace.firstIncidence + offset
                || incidence.traceIndex != trace.traceIndex
                || incidence.halfFaceIndex
                    >= system.halfFaceTraceIndices.size()
                || system.halfFaceTraceIndices[incidence.halfFaceIndex]
                    != trace.traceIndex
                || incidence.controlCellIndex
                    >= system.localOperators.size()
                || incidence.localHalfFaceIndex
                    >= system.localOperators[incidence.controlCellIndex]
                        .halfFaceCount) {
                throw std::invalid_argument(
                    "scene fluid mimetic trace incidence is invalid");
            }
        }
    }
    for (std::size_t component = 0;
         component < system.componentCount; ++component) {
        const std::size_t gauge =
            system.componentGaugeTraceIndices[component];
        if (!gaugeSeen[component]
            || gauge >= system.traces.size()
            || !system.traces[gauge].isGauge
            || system.traces[gauge].componentIndex != component) {
            throw std::invalid_argument(
                "scene fluid mimetic trace gauge is invalid");
        }
    }
    if (localBytes != system.localOperatorStorageBytes
        || sharedCount != system.sharedTraceCount
        || wallCount != system.materialWallTraceCount
        || minimumDiagonal != system.minimumPositiveOperatorDiagonal
        || maximumDiagonal != system.maximumOperatorDiagonal
        || system.ownedStorageBytes != storageBytes(system)
        || productFingerprint(system) != system.fingerprint) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-system summary is invalid");
    }
}

void validateSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticTraceSystem& system,
    const SceneFluidMimeticControlCellSet& controlCells) {
    validateSceneFluidMimeticControlCellIntegrity(controlCells);
    validateSceneFluidMimeticTraceSystemIntegrity(system);
    if (system.mimeticControlCellFingerprint != controlCells.fingerprint
        || system.structureDefinitionFingerprint
            != controlCells.structureDefinitionFingerprint
        || system.acceptedStepCount != controlCells.acceptedStepCount
        || system.simulationTimeSeconds
            != controlCells.simulationTimeSeconds) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-system identity is invalid");
    }
    const SceneFluidMimeticTraceSystemLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
    };
    const auto expected = buildSystem(
        controlCells, system.settings, unlimited);
    if (system != expected) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-system payload is invalid");
    }
}

SceneFluidMimeticTraceEvaluation evaluateSceneFluidMimeticTraceSystem(
    const SceneFluidMimeticTraceSystem& system,
    const std::span<const double> traceScalars,
    const std::span<const double> integratedCellSources) {
    validateSceneFluidMimeticTraceSystemIntegrity(system);
    if (traceScalars.size() != system.traces.size()
        || integratedCellSources.size() != system.localOperators.size()
        || !std::ranges::all_of(
            traceScalars,
            [](const double value) { return std::isfinite(value); })
        || !std::ranges::all_of(
            integratedCellSources,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "scene fluid mimetic trace-system field is invalid");
    }
    SceneFluidMimeticTraceEvaluation result;
    result.cellScalars.resize(system.localOperators.size(), 0.0);
    result.halfFaceIntegratedOutwardFluxes.resize(
        system.halfFaceTraceIndices.size(), 0.0);
    result.traceIntegratedOutwardFluxSums.resize(
        system.traces.size(), 0.0);
    std::size_t firstHalfFace = 0;
    for (std::size_t cell = 0;
         cell < system.localOperators.size(); ++cell) {
        const auto& local = system.localOperators[cell];
        std::vector<double> localTraces(local.halfFaceCount, 0.0);
        for (std::size_t face = 0; face < local.halfFaceCount; ++face) {
            localTraces[face] = traceScalars[
                system.halfFaceTraceIndices[firstHalfFace + face]];
        }
        const auto balance = fluid::balanceMimeticLocalCell(
            local, localTraces, integratedCellSources[cell]);
        result.cellScalars[cell] = balance.cellScalar;
        result.maximumCellConservationResidual = std::max(
            result.maximumCellConservationResidual,
            std::abs(balance.conservationResidual));
        for (std::size_t face = 0; face < local.halfFaceCount; ++face) {
            const std::size_t halfFaceIndex = firstHalfFace + face;
            const std::size_t traceIndex =
                system.halfFaceTraceIndices[halfFaceIndex];
            const double flux = balance.integratedOutwardFluxes[face];
            result.halfFaceIntegratedOutwardFluxes[halfFaceIndex] = flux;
            result.traceIntegratedOutwardFluxSums[traceIndex] += flux;
        }
        firstHalfFace += local.halfFaceCount;
    }
    if (firstHalfFace != system.halfFaceTraceIndices.size()) {
        throw std::logic_error(
            "scene fluid mimetic trace-system local ranges do not close");
    }
    for (const double flux : result.traceIntegratedOutwardFluxSums) {
        result.maximumTraceFluxImbalance = std::max(
            result.maximumTraceFluxImbalance, std::abs(flux));
    }
    return result;
}

std::vector<double> applySceneFluidMimeticTraceOperator(
    const SceneFluidMimeticTraceSystem& system,
    const std::span<const double> traceScalars) {
    const std::vector<double> sources(system.localOperators.size(), 0.0);
    auto result = evaluateSceneFluidMimeticTraceSystem(
        system, traceScalars, sources).traceIntegratedOutwardFluxSums;
    for (double& value : result) value = -value;
    return result;
}

std::vector<double> buildSceneFluidMimeticTraceRightHandSide(
    const SceneFluidMimeticTraceSystem& system,
    const std::span<const double> integratedCellSources) {
    const std::vector<double> traces(system.traces.size(), 0.0);
    return evaluateSceneFluidMimeticTraceSystem(
        system, traces, integratedCellSources)
        .traceIntegratedOutwardFluxSums;
}

} // namespace simwing::fsi
