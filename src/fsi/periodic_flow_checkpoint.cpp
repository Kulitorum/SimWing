#include "periodic_flow_case.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::array<std::uint8_t, 8> checkpointMagic{
    'S', 'W', 'P', 'C', 'K', 'P', 'T', '1'};
constexpr std::size_t checkpointEnvelopeBytes = 28;
constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

bool fail(PeriodicFlowCaseCheckpointError* error,
          const PeriodicFlowCaseCheckpointErrorCode code,
          std::string message) {
    if (error != nullptr) {
        error->code = code;
        error->message = std::move(message);
    }
    return false;
}

void clearError(PeriodicFlowCaseCheckpointError* error) {
    if (error != nullptr) {
        *error = {};
    }
}

std::uint64_t checksum(const std::span<const std::uint8_t> bytes) noexcept {
    std::uint64_t result = fnvOffsetBasis;
    for (const std::uint8_t value : bytes) {
        result ^= value;
        result *= fnvPrime;
    }
    return result;
}

class BufferWriter final {
public:
    BufferWriter(std::vector<std::uint8_t>& bytes, const std::size_t limit)
        : bytes_(bytes), limit_(limit) {}

    bool u8(const std::uint8_t value) {
        return raw(&value, sizeof(value));
    }

    bool u16(const std::uint16_t value) {
        return unsignedInteger(value);
    }

    bool u32(const std::uint32_t value) {
        return unsignedInteger(value);
    }

    bool u64(const std::uint64_t value) {
        return unsignedInteger(value);
    }

    bool size(const std::size_t value) {
        if (value > std::numeric_limits<std::uint64_t>::max()) {
            return false;
        }
        return u64(static_cast<std::uint64_t>(value));
    }

    bool finiteDouble(const double value) {
        if (!std::isfinite(value)) {
            return false;
        }
        return u64(std::bit_cast<std::uint64_t>(value));
    }

    bool boolean(const bool value) {
        return u8(value ? 1U : 0U);
    }

    bool vector(const fluid::Vector3 value) {
        return finiteDouble(value.x)
            && finiteDouble(value.y)
            && finiteDouble(value.z);
    }

    [[nodiscard]] bool exceeded() const noexcept { return exceeded_; }

private:
    template<typename Unsigned>
    bool unsignedInteger(Unsigned value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        std::array<std::uint8_t, sizeof(Unsigned)> bytes{};
        for (std::size_t index = 0; index < bytes.size(); ++index) {
            bytes[index] = static_cast<std::uint8_t>(value & 0xffU);
            value >>= 8U;
        }
        return raw(bytes.data(), bytes.size());
    }

    bool raw(const void* data, const std::size_t size) {
        if (size > limit_ || bytes_.size() > limit_ - size) {
            exceeded_ = true;
            return false;
        }
        const auto* begin = static_cast<const std::uint8_t*>(data);
        bytes_.insert(bytes_.end(), begin, begin + size);
        return true;
    }

    std::vector<std::uint8_t>& bytes_;
    std::size_t limit_ = 0;
    bool exceeded_ = false;
};

class BufferReader final {
public:
    explicit BufferReader(const std::span<const std::uint8_t> bytes)
        : bytes_(bytes) {}

    bool u8(std::uint8_t& value) {
        return unsignedInteger(value);
    }

    bool u16(std::uint16_t& value) {
        return unsignedInteger(value);
    }

    bool u32(std::uint32_t& value) {
        return unsignedInteger(value);
    }

    bool u64(std::uint64_t& value) {
        return unsignedInteger(value);
    }

    bool size(std::size_t& value) {
        std::uint64_t encoded = 0;
        if (!u64(encoded)
            || encoded > std::numeric_limits<std::size_t>::max()) {
            return false;
        }
        value = static_cast<std::size_t>(encoded);
        return true;
    }

    bool finiteDouble(double& value) {
        std::uint64_t bits = 0;
        if (!u64(bits)) {
            return false;
        }
        value = std::bit_cast<double>(bits);
        if (!std::isfinite(value)) {
            return false;
        }
        return true;
    }

    bool boolean(bool& value) {
        std::uint8_t encoded = 0;
        if (!u8(encoded)) {
            return false;
        }
        if (encoded > 1U) {
            return false;
        }
        value = encoded != 0U;
        return true;
    }

    bool vector(fluid::Vector3& value) {
        return finiteDouble(value.x)
            && finiteDouble(value.y)
            && finiteDouble(value.z);
    }

    [[nodiscard]] bool atEnd() const noexcept {
        return position_ == bytes_.size();
    }

    [[nodiscard]] bool truncated() const noexcept { return truncated_; }
    [[nodiscard]] bool limitExceeded() const noexcept {
        return limitExceeded_;
    }
    void markLimitExceeded() noexcept { limitExceeded_ = true; }

private:
    template<typename Unsigned>
    bool unsignedInteger(Unsigned& value) {
        static_assert(std::is_unsigned_v<Unsigned>);
        if (sizeof(Unsigned) > bytes_.size() - position_) {
            truncated_ = true;
            return false;
        }
        value = 0;
        for (std::size_t index = 0; index < sizeof(Unsigned); ++index) {
            value |= static_cast<Unsigned>(bytes_[position_ + index])
                << (8U * index);
        }
        position_ += sizeof(Unsigned);
        return true;
    }

    std::span<const std::uint8_t> bytes_;
    std::size_t position_ = 0;
    bool truncated_ = false;
    bool limitExceeded_ = false;
};

template<typename Enum>
bool writeEnum(BufferWriter& writer, const Enum value) {
    return writer.u8(static_cast<std::uint8_t>(value));
}

template<typename Enum>
bool readEnum(BufferReader& reader, Enum& value, const Enum maximum) {
    std::uint8_t encoded = 0;
    if (!reader.u8(encoded)
        || encoded > static_cast<std::uint8_t>(maximum)) {
        return false;
    }
    value = static_cast<Enum>(encoded);
    return true;
}

bool writeDiffusion(
    BufferWriter& writer,
    const fluid::PeriodicMacDiffusionDiagnostics& value) {
    return writer.u32(value.version)
        && writer.finiteDouble(value.densityKgPerCubicMeter)
        && writer.finiteDouble(
            value.kinematicViscositySquareMetersPerSecond)
        && writer.finiteDouble(value.timeStepSeconds)
        && writer.vector(value.directionalDiffusionNumbers)
        && writer.finiteDouble(value.totalDiffusionNumber)
        && writer.finiteDouble(value.maximumAcceptedDiffusionNumber)
        && writer.vector(value.momentumBeforeNewtonSeconds)
        && writer.vector(value.momentumAfterNewtonSeconds)
        && writer.vector(value.momentumResidualNewtonSeconds)
        && writer.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.dissipatedKineticEnergyJoules)
        && writer.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && writer.boolean(value.stable)
        && writer.boolean(value.finite)
        && writer.boolean(value.accepted);
}

bool readDiffusion(
    BufferReader& reader,
    fluid::PeriodicMacDiffusionDiagnostics& value) {
    return reader.u32(value.version)
        && value.version == fluid::periodicMacDiffusionVersion
        && reader.finiteDouble(value.densityKgPerCubicMeter)
        && reader.finiteDouble(
            value.kinematicViscositySquareMetersPerSecond)
        && reader.finiteDouble(value.timeStepSeconds)
        && reader.vector(value.directionalDiffusionNumbers)
        && reader.finiteDouble(value.totalDiffusionNumber)
        && reader.finiteDouble(value.maximumAcceptedDiffusionNumber)
        && reader.vector(value.momentumBeforeNewtonSeconds)
        && reader.vector(value.momentumAfterNewtonSeconds)
        && reader.vector(value.momentumResidualNewtonSeconds)
        && reader.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.dissipatedKineticEnergyJoules)
        && reader.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && reader.boolean(value.stable)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeDiffusionSspRk2(
    BufferWriter& writer,
    const fluid::PeriodicMacDiffusionSspRk2Diagnostics& value) {
    return writer.u32(value.version)
        && writeDiffusion(writer, value.firstEulerStage)
        && writeDiffusion(writer, value.secondEulerStage)
        && writer.vector(value.momentumBeforeNewtonSeconds)
        && writer.vector(value.momentumAfterNewtonSeconds)
        && writer.vector(value.momentumResidualNewtonSeconds)
        && writer.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.dissipatedKineticEnergyJoules)
        && writer.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && writer.boolean(value.finite)
        && writer.boolean(value.accepted);
}

bool readDiffusionSspRk2(
    BufferReader& reader,
    fluid::PeriodicMacDiffusionSspRk2Diagnostics& value) {
    return reader.u32(value.version)
        && value.version == fluid::periodicMacDiffusionSspRk2Version
        && readDiffusion(reader, value.firstEulerStage)
        && readDiffusion(reader, value.secondEulerStage)
        && reader.vector(value.momentumBeforeNewtonSeconds)
        && reader.vector(value.momentumAfterNewtonSeconds)
        && reader.vector(value.momentumResidualNewtonSeconds)
        && reader.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.dissipatedKineticEnergyJoules)
        && reader.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeVariableAdvection(
    BufferWriter& writer,
    const fluid::VariableMacAdvectionDiagnostics& value) {
    return writer.u32(value.version)
        && writeEnum(writer, value.reconstruction)
        && writer.finiteDouble(value.densityKgPerCubicMeter)
        && writer.finiteDouble(value.timeStepSeconds)
        && writer.finiteDouble(value.maximumLocalOutgoingCourantNumber)
        && writer.finiteDouble(
            value.maximumAcceptedLocalOutgoingCourantNumber)
        && writer.finiteDouble(value.maximumAdvectingDivergencePerSecond)
        && writer.finiteDouble(value.maximumControlVolumeDivergencePerSecond)
        && writer.finiteDouble(value.acceptedDivergenceTolerancePerSecond)
        && writer.vector(value.momentumBeforeNewtonSeconds)
        && writer.vector(value.momentumAfterNewtonSeconds)
        && writer.vector(value.momentumResidualNewtonSeconds)
        && writer.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.numericalKineticEnergyLossJoules)
        && writer.vector(value.componentMinimumBeforeMetersPerSecond)
        && writer.vector(value.componentMaximumBeforeMetersPerSecond)
        && writer.vector(value.componentMinimumAfterMetersPerSecond)
        && writer.vector(value.componentMaximumAfterMetersPerSecond)
        && writer.finiteDouble(value.maximumBoundViolationMetersPerSecond)
        && writer.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && writer.boolean(value.energyCriterionEnabled)
        && writer.boolean(value.energyNonIncreasing)
        && writer.boolean(value.uniformAdvector)
        && writer.boolean(value.divergenceCompatible)
        && writer.boolean(value.stable)
        && writer.boolean(value.bounded)
        && writer.boolean(value.finite)
        && writer.boolean(value.accepted);
}

bool readVariableAdvection(
    BufferReader& reader,
    fluid::VariableMacAdvectionDiagnostics& value) {
    return reader.u32(value.version)
        && value.version == fluid::variableMacAdvectionVersion
        && readEnum(
            reader, value.reconstruction,
            fluid::VariableMacReconstruction::MonotonizedCentral)
        && reader.finiteDouble(value.densityKgPerCubicMeter)
        && reader.finiteDouble(value.timeStepSeconds)
        && reader.finiteDouble(value.maximumLocalOutgoingCourantNumber)
        && reader.finiteDouble(
            value.maximumAcceptedLocalOutgoingCourantNumber)
        && reader.finiteDouble(value.maximumAdvectingDivergencePerSecond)
        && reader.finiteDouble(value.maximumControlVolumeDivergencePerSecond)
        && reader.finiteDouble(value.acceptedDivergenceTolerancePerSecond)
        && reader.vector(value.momentumBeforeNewtonSeconds)
        && reader.vector(value.momentumAfterNewtonSeconds)
        && reader.vector(value.momentumResidualNewtonSeconds)
        && reader.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.numericalKineticEnergyLossJoules)
        && reader.vector(value.componentMinimumBeforeMetersPerSecond)
        && reader.vector(value.componentMaximumBeforeMetersPerSecond)
        && reader.vector(value.componentMinimumAfterMetersPerSecond)
        && reader.vector(value.componentMaximumAfterMetersPerSecond)
        && reader.finiteDouble(value.maximumBoundViolationMetersPerSecond)
        && reader.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && reader.boolean(value.energyCriterionEnabled)
        && reader.boolean(value.energyNonIncreasing)
        && reader.boolean(value.uniformAdvector)
        && reader.boolean(value.divergenceCompatible)
        && reader.boolean(value.stable)
        && reader.boolean(value.bounded)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeProjection(BufferWriter& writer,
                     const fluid::ProjectionDiagnostics& value) {
    return writer.boolean(value.converged)
        && writer.size(value.iterationCount)
        && writer.finiteDouble(value.compatibilityDivergencePerSecond)
        && writer.finiteDouble(value.initialResidualPascalsPerSquareMeter)
        && writer.finiteDouble(value.finalResidualPascalsPerSquareMeter)
        && writer.finiteDouble(value.divergenceL2BeforePerSecond)
        && writer.finiteDouble(value.divergenceL2AfterPerSecond)
        && writer.finiteDouble(value.divergenceMaximumBeforePerSecond)
        && writer.finiteDouble(value.divergenceMaximumAfterPerSecond)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.pressureMeanPascals)
        && writer.size(value.pressureJumpFaceCount)
        && writer.finiteDouble(
            value.pressureJumpSourceCompatibilityPascalsPerSquareMeter);
}

bool readProjection(BufferReader& reader,
                    fluid::ProjectionDiagnostics& value) {
    return reader.boolean(value.converged)
        && reader.size(value.iterationCount)
        && reader.finiteDouble(value.compatibilityDivergencePerSecond)
        && reader.finiteDouble(value.initialResidualPascalsPerSquareMeter)
        && reader.finiteDouble(value.finalResidualPascalsPerSquareMeter)
        && reader.finiteDouble(value.divergenceL2BeforePerSecond)
        && reader.finiteDouble(value.divergenceL2AfterPerSecond)
        && reader.finiteDouble(value.divergenceMaximumBeforePerSecond)
        && reader.finiteDouble(value.divergenceMaximumAfterPerSecond)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.pressureMeanPascals)
        && reader.size(value.pressureJumpFaceCount)
        && reader.finiteDouble(
            value.pressureJumpSourceCompatibilityPascalsPerSquareMeter);
}

bool writeProjectedAdvection(
    BufferWriter& writer,
    const fluid::ProjectedMacAdvectionSspRk2Diagnostics& value) {
    return writer.u32(value.version)
        && writeEnum(writer, value.reconstruction)
        && writeVariableAdvection(writer, value.firstAdvection)
        && writeProjection(writer, value.firstProjection)
        && writeVariableAdvection(writer, value.secondAdvection)
        && writeProjection(writer, value.secondProjection)
        && writer.vector(value.momentumBeforeNewtonSeconds)
        && writer.vector(value.momentumAfterNewtonSeconds)
        && writer.vector(value.momentumResidualNewtonSeconds)
        && writer.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.totalKineticEnergyLossJoules)
        && writer.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && writer.finiteDouble(value.initialDivergenceL2PerSecond)
        && writer.finiteDouble(value.finalDivergenceL2PerSecond)
        && writeEnum(writer, value.failureStage)
        && writer.boolean(value.finite)
        && writer.boolean(value.accepted);
}

bool readProjectedAdvection(
    BufferReader& reader,
    fluid::ProjectedMacAdvectionSspRk2Diagnostics& value) {
    return reader.u32(value.version)
        && value.version == fluid::projectedMacAdvectionSspRk2Version
        && readEnum(
            reader, value.reconstruction,
            fluid::VariableMacReconstruction::MonotonizedCentral)
        && readVariableAdvection(reader, value.firstAdvection)
        && readProjection(reader, value.firstProjection)
        && readVariableAdvection(reader, value.secondAdvection)
        && readProjection(reader, value.secondProjection)
        && reader.vector(value.momentumBeforeNewtonSeconds)
        && reader.vector(value.momentumAfterNewtonSeconds)
        && reader.vector(value.momentumResidualNewtonSeconds)
        && reader.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.totalKineticEnergyLossJoules)
        && reader.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && reader.finiteDouble(value.initialDivergenceL2PerSecond)
        && reader.finiteDouble(value.finalDivergenceL2PerSecond)
        && readEnum(
            reader, value.failureStage,
            fluid::ProjectedMacAdvectionFailureStage::Conservation)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeStrang(
    BufferWriter& writer,
    const fluid::PeriodicFlowStrangSspRk2Diagnostics& value) {
    return writer.u32(value.version)
        && writeDiffusionSspRk2(writer, value.firstHalfDiffusion)
        && writeProjectedAdvection(writer, value.projectedAdvection)
        && writeDiffusionSspRk2(writer, value.secondHalfDiffusion)
        && writer.vector(value.momentumBeforeNewtonSeconds)
        && writer.vector(value.momentumAfterNewtonSeconds)
        && writer.vector(value.momentumResidualNewtonSeconds)
        && writer.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.firstHalfViscousEnergyLossJoules)
        && writer.finiteDouble(value.transportProjectionEnergyLossJoules)
        && writer.finiteDouble(value.secondHalfViscousEnergyLossJoules)
        && writer.finiteDouble(value.totalEnergyLossJoules)
        && writer.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && writer.finiteDouble(value.initialDivergenceL2PerSecond)
        && writer.finiteDouble(value.finalDivergenceL2PerSecond)
        && writeEnum(writer, value.failureStage)
        && writer.boolean(value.finite)
        && writer.boolean(value.accepted);
}

bool readStrang(
    BufferReader& reader,
    fluid::PeriodicFlowStrangSspRk2Diagnostics& value) {
    return reader.u32(value.version)
        && value.version == fluid::periodicFlowStrangSspRk2Version
        && readDiffusionSspRk2(reader, value.firstHalfDiffusion)
        && readProjectedAdvection(reader, value.projectedAdvection)
        && readDiffusionSspRk2(reader, value.secondHalfDiffusion)
        && reader.vector(value.momentumBeforeNewtonSeconds)
        && reader.vector(value.momentumAfterNewtonSeconds)
        && reader.vector(value.momentumResidualNewtonSeconds)
        && reader.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.firstHalfViscousEnergyLossJoules)
        && reader.finiteDouble(value.transportProjectionEnergyLossJoules)
        && reader.finiteDouble(value.secondHalfViscousEnergyLossJoules)
        && reader.finiteDouble(value.totalEnergyLossJoules)
        && reader.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && reader.finiteDouble(value.initialDivergenceL2PerSecond)
        && reader.finiteDouble(value.finalDivergenceL2PerSecond)
        && readEnum(
            reader, value.failureStage,
            fluid::PeriodicFlowStrangFailureStage::Conservation)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeSubcycling(
    BufferWriter& writer,
    const fluid::PeriodicFlowStrangSubcyclingDiagnostics& value) {
    if (!writer.u32(value.version)
        || !writer.finiteDouble(value.requestedIntervalSeconds)
        || !writer.finiteDouble(value.substepSeconds)
        || !writer.size(value.plannedSubstepCount)
        || !writer.size(value.completedSubstepCount)
        || !writer.size(value.stabilityRetryCount)
        || !writer.size(value.failedSubstepIndex)
        || !writer.size(value.substeps.size())) {
        return false;
    }
    for (const auto& substep : value.substeps) {
        if (!writeStrang(writer, substep)) {
            return false;
        }
    }
    return writeStrang(writer, value.failedSubstep)
        && writer.finiteDouble(
            value.maximumObservedOutgoingCourantNumber)
        && writer.finiteDouble(value.maximumObservedDiffusionNumber)
        && writer.vector(value.momentumBeforeNewtonSeconds)
        && writer.vector(value.momentumAfterNewtonSeconds)
        && writer.vector(value.momentumResidualNewtonSeconds)
        && writer.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && writer.finiteDouble(value.kineticEnergyBeforeJoules)
        && writer.finiteDouble(value.kineticEnergyAfterJoules)
        && writer.finiteDouble(value.totalEnergyLossJoules)
        && writer.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && writer.finiteDouble(value.initialDivergenceL2PerSecond)
        && writer.finiteDouble(value.finalDivergenceL2PerSecond)
        && writeEnum(writer, value.failureStage)
        && writer.boolean(value.finite)
        && writer.boolean(value.accepted);
}

bool readSubcycling(
    BufferReader& reader,
    fluid::PeriodicFlowStrangSubcyclingDiagnostics& value,
    const std::uint64_t maximumSubsteps) {
    std::size_t substepCount = 0;
    if (!reader.u32(value.version)
        || value.version != fluid::periodicFlowStrangSubcyclingVersion
        || !reader.finiteDouble(value.requestedIntervalSeconds)
        || !reader.finiteDouble(value.substepSeconds)
        || !reader.size(value.plannedSubstepCount)
        || !reader.size(value.completedSubstepCount)
        || !reader.size(value.stabilityRetryCount)
        || !reader.size(value.failedSubstepIndex)
        || !reader.size(substepCount)) {
        return false;
    }
    if (substepCount > maximumSubsteps) {
        reader.markLimitExceeded();
        return false;
    }
    value.substeps.resize(substepCount);
    for (auto& substep : value.substeps) {
        if (!readStrang(reader, substep)) {
            return false;
        }
    }
    return readStrang(reader, value.failedSubstep)
        && reader.finiteDouble(
            value.maximumObservedOutgoingCourantNumber)
        && reader.finiteDouble(value.maximumObservedDiffusionNumber)
        && reader.vector(value.momentumBeforeNewtonSeconds)
        && reader.vector(value.momentumAfterNewtonSeconds)
        && reader.vector(value.momentumResidualNewtonSeconds)
        && reader.finiteDouble(value.momentumResidualNormNewtonSeconds)
        && reader.finiteDouble(value.kineticEnergyBeforeJoules)
        && reader.finiteDouble(value.kineticEnergyAfterJoules)
        && reader.finiteDouble(value.totalEnergyLossJoules)
        && reader.finiteDouble(value.maximumVelocityChangeMetersPerSecond)
        && reader.finiteDouble(value.initialDivergenceL2PerSecond)
        && reader.finiteDouble(value.finalDivergenceL2PerSecond)
        && readEnum(
            reader, value.failureStage,
            fluid::PeriodicFlowStrangSubcyclingFailureStage::Conservation)
        && reader.boolean(value.finite)
        && reader.boolean(value.accepted);
}

bool writeField(BufferWriter& writer, const std::span<const double> values) {
    if (!writer.size(values.size())) {
        return false;
    }
    for (const double value : values) {
        if (!writer.finiteDouble(value)) {
            return false;
        }
    }
    return true;
}

bool readField(BufferReader& reader,
               const std::span<double> values,
               const std::uint64_t maximumScalarSamples) {
    std::size_t count = 0;
    if (!reader.size(count) || count != values.size()) {
        return false;
    }
    if (count > maximumScalarSamples) {
        reader.markLimitExceeded();
        return false;
    }
    for (double& value : values) {
        if (!reader.finiteDouble(value)) {
            return false;
        }
    }
    return true;
}

PeriodicFlowCaseCheckpointErrorCode readerErrorCode(
    const BufferReader& reader) noexcept {
    if (reader.limitExceeded()) {
        return PeriodicFlowCaseCheckpointErrorCode::LimitExceeded;
    }
    if (reader.truncated()) {
        return PeriodicFlowCaseCheckpointErrorCode::Truncated;
    }
    return PeriodicFlowCaseCheckpointErrorCode::InvalidData;
}

} // namespace

bool serializePeriodicFlowCaseCheckpoint(
    const PeriodicFlowCaseCheckpoint& checkpoint,
    std::vector<std::uint8_t>& bytes,
    PeriodicFlowCaseCheckpointError* error,
    const PeriodicFlowCaseCheckpointLimits& limits) {
    clearError(error);
    bytes.clear();
    if (limits.maximumBytes < checkpointEnvelopeBytes
        || limits.maximumBytes > std::numeric_limits<std::size_t>::max()) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                    "configured checkpoint byte limit is invalid");
    }
    try {
        PeriodicFlowCase validator;
        validator.restore(checkpoint);
        if (checkpoint.scalarSampleCount > limits.maximumScalarSamples
            || checkpoint.detail->diagnostics.substeps.size()
                > limits.maximumSubsteps) {
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                        "periodic-flow checkpoint exceeds configured limits");
        }

        std::vector<std::uint8_t> payload;
        BufferWriter writer(
            payload,
            static_cast<std::size_t>(limits.maximumBytes)
                - checkpointEnvelopeBytes);
        const auto& detail = *checkpoint.detail;
        const bool wrote = writer.u32(checkpoint.version)
            && writer.u64(checkpoint.caseDefinitionFingerprint)
            && writer.size(checkpoint.cellCounts.x)
            && writer.size(checkpoint.cellCounts.y)
            && writer.size(checkpoint.cellCounts.z)
            && writer.vector(checkpoint.lowerMeters)
            && writer.vector(checkpoint.upperMeters)
            && writer.size(checkpoint.scalarSampleCount)
            && writer.u64(checkpoint.acceptedStepCount)
            && writer.finiteDouble(checkpoint.simulationTimeSeconds)
            && writeField(writer, detail.velocityMetersPerSecond.xFaces())
            && writeField(writer, detail.velocityMetersPerSecond.yFaces())
            && writeField(writer, detail.velocityMetersPerSecond.zFaces())
            && writeField(writer, detail.pressurePascals.values())
            && writeSubcycling(writer, detail.diagnostics);
        if (!wrote) {
            if (writer.exceeded()) {
                return fail(
                    error,
                    PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                    "serialized periodic-flow checkpoint exceeds byte limit");
            }
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::InvalidData,
                        "periodic-flow checkpoint contains invalid data");
        }

        bytes.reserve(checkpointEnvelopeBytes + payload.size());
        bytes.insert(
            bytes.end(), checkpointMagic.begin(), checkpointMagic.end());
        BufferWriter envelope(
            bytes, static_cast<std::size_t>(limits.maximumBytes));
        if (!envelope.u16(periodicFlowCaseCheckpointProtocolVersion)
            || !envelope.u16(0)
            || !envelope.size(payload.size())
            || !envelope.u64(checksum(payload))) {
            bytes.clear();
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                        "periodic-flow checkpoint envelope exceeds byte limit");
        }
        bytes.insert(bytes.end(), payload.begin(), payload.end());
        return true;
    } catch (const std::invalid_argument&) {
        bytes.clear();
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::InvalidData,
                    "periodic-flow checkpoint is not restorable");
    } catch (const std::bad_alloc&) {
        bytes.clear();
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                    "unable to allocate periodic-flow checkpoint buffer");
    } catch (const std::length_error&) {
        bytes.clear();
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                    "periodic-flow checkpoint exceeds platform limits");
    }
}

bool deserializePeriodicFlowCaseCheckpoint(
    const std::span<const std::uint8_t> bytes,
    PeriodicFlowCaseCheckpoint& checkpoint,
    PeriodicFlowCaseCheckpointError* error,
    const PeriodicFlowCaseCheckpointLimits& limits) {
    clearError(error);
    if (limits.maximumBytes < checkpointEnvelopeBytes
        || limits.maximumBytes > std::numeric_limits<std::size_t>::max()
        || bytes.size() > limits.maximumBytes) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                    "periodic-flow checkpoint exceeds configured byte limit");
    }
    if (bytes.size() < checkpointEnvelopeBytes) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::Truncated,
                    "periodic-flow checkpoint envelope is truncated");
    }
    if (!std::ranges::equal(
            checkpointMagic,
            bytes.first(checkpointMagic.size()))) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::InvalidMagic,
                    "periodic-flow checkpoint magic is invalid");
    }

    BufferReader envelope(bytes.subspan(checkpointMagic.size()));
    std::uint16_t protocolVersion = 0;
    std::uint16_t reserved = 0;
    std::size_t payloadSize = 0;
    std::uint64_t expectedChecksum = 0;
    if (!envelope.u16(protocolVersion)
        || !envelope.u16(reserved)
        || !envelope.size(payloadSize)
        || !envelope.u64(expectedChecksum)) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::Truncated,
                    "periodic-flow checkpoint envelope is truncated");
    }
    if (protocolVersion != periodicFlowCaseCheckpointProtocolVersion) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::UnsupportedVersion,
                    "periodic-flow checkpoint protocol version is unsupported");
    }
    if (reserved != 0) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::InvalidData,
                    "periodic-flow checkpoint reserved bits are nonzero");
    }
    const std::size_t availablePayload = bytes.size() - checkpointEnvelopeBytes;
    if (payloadSize > availablePayload) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::Truncated,
                    "periodic-flow checkpoint payload is truncated");
    }
    if (payloadSize < availablePayload) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::TrailingData,
                    "periodic-flow checkpoint has trailing data");
    }
    const auto payload = bytes.subspan(checkpointEnvelopeBytes, payloadSize);
    if (checksum(payload) != expectedChecksum) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::ChecksumMismatch,
                    "periodic-flow checkpoint checksum does not match");
    }

    try {
        BufferReader reader(payload);
        PeriodicFlowCaseCheckpoint candidate;
        std::size_t countX = 0;
        std::size_t countY = 0;
        std::size_t countZ = 0;
        if (!reader.u32(candidate.version)
            || !reader.u64(candidate.caseDefinitionFingerprint)
            || !reader.size(countX)
            || !reader.size(countY)
            || !reader.size(countZ)
            || !reader.vector(candidate.lowerMeters)
            || !reader.vector(candidate.upperMeters)
            || !reader.size(candidate.scalarSampleCount)
            || !reader.u64(candidate.acceptedStepCount)
            || !reader.finiteDouble(candidate.simulationTimeSeconds)) {
            return fail(error, readerErrorCode(reader),
                        "periodic-flow checkpoint metadata is invalid");
        }
        if (candidate.version != periodicFlowCaseCheckpointVersion) {
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::UnsupportedVersion,
                        "periodic-flow checkpoint version is unsupported");
        }
        if (candidate.caseDefinitionFingerprint
            != periodicFlowCaseDefinitionFingerprint) {
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::InvalidData,
                        "periodic-flow checkpoint identity is invalid");
        }
        if (candidate.scalarSampleCount > limits.maximumScalarSamples) {
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                        "periodic-flow checkpoint sample count exceeds limit");
        }
        candidate.cellCounts = {countX, countY, countZ};
        const fluid::PeriodicCartesianGrid grid(
            candidate.cellCounts,
            candidate.lowerMeters,
            candidate.upperMeters);
        if (grid.cellCount() != candidate.scalarSampleCount) {
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::InvalidData,
                        "periodic-flow checkpoint grid size is inconsistent");
        }

        fluid::MacVelocityField velocity(grid);
        fluid::CellScalarField pressure(grid);
        fluid::PeriodicFlowStrangSubcyclingDiagnostics diagnostics;
        if (!readField(
                reader, velocity.xFaces(), limits.maximumScalarSamples)
            || !readField(
                reader, velocity.yFaces(), limits.maximumScalarSamples)
            || !readField(
                reader, velocity.zFaces(), limits.maximumScalarSamples)
            || !readField(
                reader, pressure.values(), limits.maximumScalarSamples)
            || !readSubcycling(
                reader, diagnostics, limits.maximumSubsteps)) {
            return fail(error, readerErrorCode(reader),
                        "periodic-flow checkpoint payload is invalid");
        }
        if (!reader.atEnd()) {
            return fail(error,
                        PeriodicFlowCaseCheckpointErrorCode::TrailingData,
                        "periodic-flow checkpoint payload has trailing data");
        }
        candidate.detail =
            std::make_shared<PeriodicFlowCaseCheckpoint::Detail>(
                PeriodicFlowCaseCheckpoint::Detail{
                    candidate.cellCounts,
                    candidate.lowerMeters,
                    candidate.upperMeters,
                    std::move(velocity),
                    std::move(pressure),
                    std::move(diagnostics),
                    candidate.acceptedStepCount,
                    candidate.simulationTimeSeconds,
                });
        PeriodicFlowCase validator;
        validator.restore(candidate);
        checkpoint = std::move(candidate);
        return true;
    } catch (const std::invalid_argument&) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::InvalidData,
                    "periodic-flow checkpoint cannot restore into this case");
    } catch (const std::bad_alloc&) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                    "unable to allocate periodic-flow checkpoint state");
    } catch (const std::length_error&) {
        return fail(error,
                    PeriodicFlowCaseCheckpointErrorCode::LimitExceeded,
                    "periodic-flow checkpoint exceeds platform limits");
    }
}

} // namespace simwing::fsi
