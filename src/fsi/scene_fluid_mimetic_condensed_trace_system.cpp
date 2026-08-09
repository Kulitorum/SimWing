#include "scene_fluid_mimetic_condensed_trace_system.h"

#include <algorithm>
#include <bit>
#include <cmath>
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

std::uint64_t productFingerprint(
    const SceneFluidMimeticCondensedTraceSystem& system) {
    Fingerprint fingerprint;
    fingerprint.integer(system.version);
    fingerprint.integer(system.fullTraceSystemFingerprint);
    fingerprint.integer(system.structureDefinitionFingerprint);
    fingerprint.integer(system.acceptedStepCount);
    fingerprint.real(system.simulationTimeSeconds);
    for (const std::size_t value : {
             system.ownedStorageBytes,
             system.localCondensationStorageBytes,
             system.componentCount,
             system.eliminatedMaterialWallTraceCount}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.real(system.minimumPositiveOperatorDiagonal);
    fingerprint.real(system.maximumOperatorDiagonal);
    fingerprint.integer(static_cast<std::uint64_t>(system.traces.size()));
    for (const auto& trace : system.traces) {
        fingerprint.integer(static_cast<std::uint64_t>(trace.traceIndex));
        fingerprint.integer(static_cast<std::uint64_t>(trace.fullTraceIndex));
        fingerprint.integer(trace.stableId);
        fingerprint.enumeration(trace.kind);
        fingerprint.integer(static_cast<std::uint64_t>(trace.componentIndex));
        fingerprint.real(trace.operatorDiagonal);
        fingerprint.integer(static_cast<std::uint8_t>(trace.isGauge));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        system.fullTraceReducedIndices.size()));
    for (const std::size_t value : system.fullTraceReducedIndices) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        system.componentGaugeTraceIndices.size()));
    for (const std::size_t value : system.componentGaugeTraceIndices) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.integer(static_cast<std::uint64_t>(
        system.localCondensations.size()));
    for (const auto& local : system.localCondensations) {
        fingerprint.integer(local.fingerprint);
    }
    return fingerprint.value();
}

std::size_t storageBytes(
    const SceneFluidMimeticCondensedTraceSystem& system) {
    std::size_t total = 0;
    const auto add = [&](const std::size_t count,
                         const std::size_t elementSize) {
        std::size_t bytes = 0;
        if (!checkedMultiply(count, elementSize, bytes)
            || !checkedAdd(total, bytes, total)) {
            throw std::length_error(
                "scene fluid condensed trace-system storage overflows");
        }
    };
    add(system.traces.size(), sizeof(SceneFluidMimeticCondensedTrace));
    add(system.fullTraceReducedIndices.size(), sizeof(std::size_t));
    add(system.componentGaugeTraceIndices.size(), sizeof(std::size_t));
    add(system.localCondensations.size(),
        sizeof(fluid::MimeticWallCondensation));
    if (!checkedAdd(total, system.localCondensationStorageBytes, total)) {
        throw std::length_error(
            "scene fluid condensed trace-system storage overflows");
    }
    return total;
}

bool finiteField(const std::span<const double> values) {
    return std::ranges::all_of(
        values, [](const double value) { return std::isfinite(value); });
}

SceneFluidMimeticCondensedTraceSystem buildSystem(
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystemLimits& limits) {
    validateSceneFluidMimeticTraceSystemIntegrity(fullSystem);
    if (fullSystem.sharedTraceCount > limits.maximumReducedTraces
        || fullSystem.localOperators.size()
            > limits.maximumLocalCondensations) {
        throw std::length_error(
            "scene fluid condensed trace-system count limit exceeded");
    }

    SceneFluidMimeticCondensedTraceSystem result;
    result.fullTraceSystemFingerprint = fullSystem.fingerprint;
    result.structureDefinitionFingerprint =
        fullSystem.structureDefinitionFingerprint;
    result.acceptedStepCount = fullSystem.acceptedStepCount;
    result.simulationTimeSeconds = fullSystem.simulationTimeSeconds;
    result.componentCount = fullSystem.componentCount;
    result.eliminatedMaterialWallTraceCount =
        fullSystem.materialWallTraceCount;
    result.traces.reserve(fullSystem.sharedTraceCount);
    result.fullTraceReducedIndices.assign(
        fullSystem.traces.size(),
        invalidSceneFluidMimeticCondensedTraceIndex);
    result.componentGaugeTraceIndices.assign(
        result.componentCount,
        invalidSceneFluidMimeticCondensedTraceIndex);
    for (const auto& full : fullSystem.traces) {
        if (full.kind == SceneFluidMimeticHalfFaceKind::MaterialWall) {
            continue;
        }
        SceneFluidMimeticCondensedTrace trace;
        trace.traceIndex = result.traces.size();
        trace.fullTraceIndex = full.traceIndex;
        trace.stableId = full.stableId;
        trace.kind = full.kind;
        trace.componentIndex = full.componentIndex;
        result.fullTraceReducedIndices[full.traceIndex] = trace.traceIndex;
        auto& gauge = result.componentGaugeTraceIndices[
            trace.componentIndex];
        if (gauge == invalidSceneFluidMimeticCondensedTraceIndex) {
            gauge = trace.traceIndex;
            trace.isGauge = true;
        }
        result.traces.push_back(trace);
    }
    if (result.traces.size() != fullSystem.sharedTraceCount
        || std::ranges::find(
            result.componentGaugeTraceIndices,
            invalidSceneFluidMimeticCondensedTraceIndex)
            != result.componentGaugeTraceIndices.end()) {
        throw std::invalid_argument(
            "scene fluid condensed trace-system component has no shared trace");
    }

    result.localCondensations.reserve(fullSystem.localOperators.size());
    std::size_t firstHalfFace = 0;
    for (const auto& localOperator : fullSystem.localOperators) {
        std::vector<std::uint8_t> wallMask(
            localOperator.halfFaceCount, 0);
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            wallMask[face] = fullSystem.traces[fullTrace].kind
                == SceneFluidMimeticHalfFaceKind::MaterialWall;
        }
        auto local = fluid::buildMimeticWallCondensation(
            localOperator, wallMask, limits.local);
        if (!checkedAdd(result.localCondensationStorageBytes,
                        local.ownedStorageBytes,
                        result.localCondensationStorageBytes)
            || result.localCondensationStorageBytes
                > limits.maximumLocalCondensationBytes) {
            throw std::length_error(
                "scene fluid condensed trace-system local byte limit exceeded");
        }
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            if (wallMask[face] != 0) continue;
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            const std::size_t reducedTrace =
                result.fullTraceReducedIndices[fullTrace];
            if (reducedTrace >= result.traces.size()) {
                throw std::logic_error(
                    "scene fluid condensed trace mapping is incomplete");
            }
            result.traces[reducedTrace].operatorDiagonal +=
                local.condensedOperatorDiagonal[face];
        }
        firstHalfFace += localOperator.halfFaceCount;
        result.localCondensations.push_back(std::move(local));
    }
    if (firstHalfFace != fullSystem.halfFaceTraceIndices.size()) {
        throw std::logic_error(
            "scene fluid condensed trace local ranges do not close");
    }

    result.minimumPositiveOperatorDiagonal =
        std::numeric_limits<double>::infinity();
    for (const auto& trace : result.traces) {
        if (!std::isfinite(trace.operatorDiagonal)
            || trace.operatorDiagonal <= 0.0) {
            throw std::invalid_argument(
                "scene fluid condensed trace has no positive diagonal");
        }
        result.minimumPositiveOperatorDiagonal = std::min(
            result.minimumPositiveOperatorDiagonal,
            trace.operatorDiagonal);
        result.maximumOperatorDiagonal = std::max(
            result.maximumOperatorDiagonal, trace.operatorDiagonal);
    }
    result.ownedStorageBytes = storageBytes(result);
    if (result.ownedStorageBytes > limits.maximumOwnedBytes) {
        throw std::length_error(
            "scene fluid condensed trace-system byte limit exceeded");
    }
    result.fingerprint = productFingerprint(result);
    return result;
}

void validatePayload(
    const SceneFluidMimeticCondensedTraceSystem& system,
    const SceneFluidMimeticTraceSystem& fullSystem) {
    validateSceneFluidMimeticTraceSystemIntegrity(fullSystem);
    if (system.version != sceneFluidMimeticCondensedTraceSystemVersion
        || system.fingerprint == 0
        || system.fullTraceSystemFingerprint != fullSystem.fingerprint
        || system.structureDefinitionFingerprint
            != fullSystem.structureDefinitionFingerprint
        || system.acceptedStepCount != fullSystem.acceptedStepCount
        || system.simulationTimeSeconds != fullSystem.simulationTimeSeconds
        || system.componentCount != fullSystem.componentCount
        || system.eliminatedMaterialWallTraceCount
            != fullSystem.materialWallTraceCount
        || system.traces.size() != fullSystem.sharedTraceCount
        || system.fullTraceReducedIndices.size()
            != fullSystem.traces.size()
        || system.componentGaugeTraceIndices.size()
            != system.componentCount
        || system.localCondensations.size()
            != fullSystem.localOperators.size()) {
        throw std::invalid_argument(
            "scene fluid condensed trace-system integrity is invalid");
    }
    std::size_t localBytes = 0;
    for (std::size_t local = 0;
         local < system.localCondensations.size(); ++local) {
        fluid::validateMimeticWallCondensation(
            system.localCondensations[local],
            fullSystem.localOperators[local]);
        if (!checkedAdd(localBytes,
                        system.localCondensations[local].ownedStorageBytes,
                        localBytes)) {
            throw std::length_error(
                "scene fluid condensed trace local storage overflows");
        }
    }
    std::vector<bool> gauges(system.componentCount, false);
    double minimumDiagonal = std::numeric_limits<double>::infinity();
    double maximumDiagonal = 0.0;
    std::uint64_t previousStableId = 0;
    for (const auto& trace : system.traces) {
        if (trace.traceIndex != &trace - system.traces.data()
            || trace.fullTraceIndex >= fullSystem.traces.size()
            || trace.stableId == 0
            || (trace.traceIndex != 0
                && trace.stableId <= previousStableId)
            || trace.kind == SceneFluidMimeticHalfFaceKind::MaterialWall
            || trace.componentIndex >= system.componentCount
            || !std::isfinite(trace.operatorDiagonal)
            || trace.operatorDiagonal <= 0.0
            || system.fullTraceReducedIndices[trace.fullTraceIndex]
                != trace.traceIndex) {
            throw std::invalid_argument(
                "scene fluid condensed trace record is invalid");
        }
        const auto& full = fullSystem.traces[trace.fullTraceIndex];
        if (trace.stableId != full.stableId
            || trace.kind != full.kind
            || trace.componentIndex != full.componentIndex
            || (trace.isGauge && gauges[trace.componentIndex])) {
            throw std::invalid_argument(
                "scene fluid condensed trace source is invalid");
        }
        previousStableId = trace.stableId;
        gauges[trace.componentIndex] =
            gauges[trace.componentIndex] || trace.isGauge;
        minimumDiagonal = std::min(
            minimumDiagonal, trace.operatorDiagonal);
        maximumDiagonal = std::max(
            maximumDiagonal, trace.operatorDiagonal);
    }
    for (const auto& full : fullSystem.traces) {
        const std::size_t reduced =
            system.fullTraceReducedIndices[full.traceIndex];
        if ((full.kind == SceneFluidMimeticHalfFaceKind::MaterialWall
             && reduced != invalidSceneFluidMimeticCondensedTraceIndex)
            || (full.kind != SceneFluidMimeticHalfFaceKind::MaterialWall
                && (reduced >= system.traces.size()
                    || system.traces[reduced].fullTraceIndex
                        != full.traceIndex))) {
            throw std::invalid_argument(
                "scene fluid condensed trace reverse mapping is invalid");
        }
    }
    for (std::size_t component = 0;
         component < system.componentCount; ++component) {
        const std::size_t gauge =
            system.componentGaugeTraceIndices[component];
        if (!gauges[component] || gauge >= system.traces.size()
            || !system.traces[gauge].isGauge
            || system.traces[gauge].componentIndex != component) {
            throw std::invalid_argument(
                "scene fluid condensed trace gauge is invalid");
        }
    }
    if (localBytes != system.localCondensationStorageBytes
        || minimumDiagonal != system.minimumPositiveOperatorDiagonal
        || maximumDiagonal != system.maximumOperatorDiagonal
        || system.ownedStorageBytes != storageBytes(system)
        || productFingerprint(system) != system.fingerprint) {
        throw std::invalid_argument(
            "scene fluid condensed trace-system summary is invalid");
    }
}

} // namespace

SceneFluidMimeticCondensedTraceSystem
buildSceneFluidMimeticCondensedTraceSystem(
    const SceneFluidMimeticTraceSystem& fullSystem,
    const SceneFluidMimeticCondensedTraceSystemLimits& limits) {
    auto result = buildSystem(fullSystem, limits);
    validateSceneFluidMimeticCondensedTraceSystem(result, fullSystem);
    return result;
}

void validateSceneFluidMimeticCondensedTraceSystem(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem) {
    validatePayload(condensedSystem, fullSystem);
    const SceneFluidMimeticCondensedTraceSystemLimits unlimited{
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max(),
        {std::numeric_limits<std::size_t>::max(),
         std::numeric_limits<std::size_t>::max()},
    };
    const auto expected = buildSystem(fullSystem, unlimited);
    if (condensedSystem != expected) {
        throw std::invalid_argument(
            "scene fluid condensed trace-system payload is invalid");
    }
}

std::vector<double> applySceneFluidMimeticCondensedTraceOperator(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const std::span<const double> reducedTraceScalars) {
    validatePayload(condensedSystem, fullSystem);
    if (reducedTraceScalars.size() != condensedSystem.traces.size()
        || !finiteField(reducedTraceScalars)) {
        throw std::invalid_argument(
            "scene fluid condensed trace field is invalid");
    }
    std::vector<double> result(condensedSystem.traces.size(), 0.0);
    std::size_t firstHalfFace = 0;
    for (std::size_t cell = 0;
         cell < fullSystem.localOperators.size(); ++cell) {
        const auto& localOperator = fullSystem.localOperators[cell];
        std::vector<double> localTraces(localOperator.halfFaceCount, 0.0);
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            const std::size_t reducedTrace =
                condensedSystem.fullTraceReducedIndices[fullTrace];
            if (reducedTrace != invalidSceneFluidMimeticCondensedTraceIndex) {
                localTraces[face] = reducedTraceScalars[reducedTrace];
            }
        }
        const auto localAction =
            fluid::applyMimeticWallCondensedTraceOperator(
                condensedSystem.localCondensations[cell],
                localOperator, localTraces);
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            const std::size_t reducedTrace =
                condensedSystem.fullTraceReducedIndices[fullTrace];
            if (reducedTrace != invalidSceneFluidMimeticCondensedTraceIndex) {
                result[reducedTrace] += localAction[face];
            }
        }
        firstHalfFace += localOperator.halfFaceCount;
    }
    if (!finiteField(result)) {
        throw std::overflow_error(
            "scene fluid condensed trace action overflowed");
    }
    return result;
}

std::vector<double> condenseSceneFluidMimeticTraceRightHandSide(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const std::span<const double> fullTraceRightHandSide) {
    validatePayload(condensedSystem, fullSystem);
    if (fullTraceRightHandSide.size() != fullSystem.traces.size()
        || !finiteField(fullTraceRightHandSide)) {
        throw std::invalid_argument(
            "scene fluid condensed trace RHS is invalid");
    }
    std::vector<double> result(condensedSystem.traces.size(), 0.0);
    for (const auto& trace : condensedSystem.traces) {
        result[trace.traceIndex] =
            fullTraceRightHandSide[trace.fullTraceIndex];
    }
    std::size_t firstHalfFace = 0;
    for (std::size_t cell = 0;
         cell < fullSystem.localOperators.size(); ++cell) {
        const auto& localOperator = fullSystem.localOperators[cell];
        std::vector<double> localRightHandSide(
            localOperator.halfFaceCount, 0.0);
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            if (condensedSystem.fullTraceReducedIndices[fullTrace]
                == invalidSceneFluidMimeticCondensedTraceIndex) {
                localRightHandSide[face] =
                    fullTraceRightHandSide[fullTrace];
            }
        }
        const auto correction =
            fluid::condenseMimeticWallTraceRightHandSide(
                condensedSystem.localCondensations[cell],
                localOperator, localRightHandSide);
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            const std::size_t reducedTrace =
                condensedSystem.fullTraceReducedIndices[fullTrace];
            if (reducedTrace != invalidSceneFluidMimeticCondensedTraceIndex) {
                result[reducedTrace] += correction[face];
            }
        }
        firstHalfFace += localOperator.halfFaceCount;
    }
    if (!finiteField(result)) {
        throw std::overflow_error(
            "scene fluid condensed trace RHS overflowed");
    }
    return result;
}

std::vector<double> reconstructSceneFluidMimeticFullTraces(
    const SceneFluidMimeticCondensedTraceSystem& condensedSystem,
    const SceneFluidMimeticTraceSystem& fullSystem,
    const std::span<const double> fullTraceRightHandSide,
    const std::span<const double> reducedTraceScalars) {
    validatePayload(condensedSystem, fullSystem);
    if (fullTraceRightHandSide.size() != fullSystem.traces.size()
        || reducedTraceScalars.size() != condensedSystem.traces.size()
        || !finiteField(fullTraceRightHandSide)
        || !finiteField(reducedTraceScalars)) {
        throw std::invalid_argument(
            "scene fluid full-trace reconstruction field is invalid");
    }
    std::vector<double> result(fullSystem.traces.size(), 0.0);
    std::vector<bool> wallWritten(fullSystem.traces.size(), false);
    for (const auto& trace : condensedSystem.traces) {
        result[trace.fullTraceIndex] = reducedTraceScalars[trace.traceIndex];
    }
    std::size_t firstHalfFace = 0;
    for (std::size_t cell = 0;
         cell < fullSystem.localOperators.size(); ++cell) {
        const auto& localOperator = fullSystem.localOperators[cell];
        std::vector<double> localRightHandSide(
            localOperator.halfFaceCount, 0.0);
        std::vector<double> localActiveTraces(
            localOperator.halfFaceCount, 0.0);
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            const std::size_t reducedTrace =
                condensedSystem.fullTraceReducedIndices[fullTrace];
            if (reducedTrace == invalidSceneFluidMimeticCondensedTraceIndex) {
                localRightHandSide[face] =
                    fullTraceRightHandSide[fullTrace];
            } else {
                localActiveTraces[face] =
                    reducedTraceScalars[reducedTrace];
            }
        }
        const auto localTraces = fluid::reconstructMimeticWallTraces(
            condensedSystem.localCondensations[cell],
            localOperator, localRightHandSide, localActiveTraces);
        for (std::size_t face = 0;
             face < localOperator.halfFaceCount; ++face) {
            const std::size_t fullTrace =
                fullSystem.halfFaceTraceIndices[firstHalfFace + face];
            if (condensedSystem.fullTraceReducedIndices[fullTrace]
                == invalidSceneFluidMimeticCondensedTraceIndex) {
                if (wallWritten[fullTrace]) {
                    throw std::logic_error(
                        "scene fluid wall trace reconstruction is duplicated");
                }
                result[fullTrace] = localTraces[face];
                wallWritten[fullTrace] = true;
            }
        }
        firstHalfFace += localOperator.halfFaceCount;
    }
    for (const auto& full : fullSystem.traces) {
        if (full.kind == SceneFluidMimeticHalfFaceKind::MaterialWall
            && !wallWritten[full.traceIndex]) {
            throw std::logic_error(
                "scene fluid wall trace reconstruction is incomplete");
        }
    }
    if (!finiteField(result)) {
        throw std::overflow_error(
            "scene fluid full-trace reconstruction overflowed");
    }
    return result;
}

} // namespace simwing::fsi
