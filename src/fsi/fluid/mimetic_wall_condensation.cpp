#include "mimetic_wall_condensation.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace simwing::fsi::fluid {
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

std::size_t ownedBytes(const std::size_t count) {
    constexpr std::size_t bytesPerFace = sizeof(std::uint8_t)
        + 2 * sizeof(double);
    if (count > std::numeric_limits<std::size_t>::max() / bytesPerFace) {
        throw std::length_error(
            "mimetic wall-condensation storage overflows");
    }
    return count * bytesPerFace;
}

std::uint64_t productFingerprint(
    const MimeticWallCondensation& condensation) {
    Fingerprint fingerprint;
    fingerprint.integer(condensation.version);
    fingerprint.integer(condensation.localOperatorFingerprint);
    for (const std::size_t value : {
             condensation.halfFaceCount,
             condensation.wallHalfFaceCount,
             condensation.activeHalfFaceCount,
             condensation.ownedStorageBytes}) {
        fingerprint.integer(static_cast<std::uint64_t>(value));
    }
    fingerprint.real(condensation.conservationDenominator);
    fingerprint.integer(static_cast<std::uint64_t>(
        condensation.wallMask.size()));
    for (const std::uint8_t value : condensation.wallMask) {
        fingerprint.integer(value);
    }
    for (const auto& values : {
             &condensation.conservationCoupling,
             &condensation.condensedOperatorDiagonal}) {
        fingerprint.integer(static_cast<std::uint64_t>(values->size()));
        for (const double value : *values) fingerprint.real(value);
    }
    for (const double value : condensation.inverseWoodburyCore) {
        fingerprint.real(value);
    }
    for (const double value : condensation.wallSchurMetric) {
        fingerprint.real(value);
    }
    return fingerprint.value();
}

template<std::size_t Size>
std::array<double, Size * Size> invertEquilibrated(
    const std::array<double, Size * Size>& matrix,
    const char* message) {
    std::array<double, Size> scales{};
    for (std::size_t row = 0; row < Size; ++row) {
        double maximum = 0.0;
        for (std::size_t column = 0; column < Size; ++column) {
            maximum = std::max(
                maximum, std::abs(matrix[row * Size + column]));
        }
        if (!std::isfinite(maximum) || maximum <= 0.0) {
            throw std::invalid_argument(message);
        }
        scales[row] = std::sqrt(maximum);
    }
    std::array<double, Size * 2 * Size> augmented{};
    double maximumScaled = 0.0;
    for (std::size_t row = 0; row < Size; ++row) {
        for (std::size_t column = 0; column < Size; ++column) {
            const double value = matrix[row * Size + column]
                / (scales[row] * scales[column]);
            augmented[row * 2 * Size + column] = value;
            maximumScaled = std::max(maximumScaled, std::abs(value));
        }
        augmented[row * 2 * Size + Size + row] = 1.0;
    }
    const double tolerance = 4096.0
        * std::numeric_limits<double>::epsilon() * maximumScaled;
    for (std::size_t column = 0; column < Size; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < Size; ++row) {
            if (std::abs(augmented[row * 2 * Size + column])
                > std::abs(augmented[pivot * 2 * Size + column])) {
                pivot = row;
            }
        }
        const double pivotValue =
            augmented[pivot * 2 * Size + column];
        if (!std::isfinite(pivotValue)
            || std::abs(pivotValue) <= tolerance) {
            throw std::invalid_argument(message);
        }
        if (pivot != column) {
            for (std::size_t entry = 0; entry < 2 * Size; ++entry) {
                std::swap(augmented[column * 2 * Size + entry],
                          augmented[pivot * 2 * Size + entry]);
            }
        }
        const double diagonal =
            augmented[column * 2 * Size + column];
        for (std::size_t entry = 0; entry < 2 * Size; ++entry) {
            augmented[column * 2 * Size + entry] /= diagonal;
        }
        for (std::size_t row = 0; row < Size; ++row) {
            if (row == column) continue;
            const double factor = augmented[row * 2 * Size + column];
            for (std::size_t entry = 0; entry < 2 * Size; ++entry) {
                augmented[row * 2 * Size + entry] -= factor
                    * augmented[column * 2 * Size + entry];
            }
        }
    }
    std::array<double, Size * Size> result{};
    for (std::size_t row = 0; row < Size; ++row) {
        for (std::size_t column = 0; column < Size; ++column) {
            result[row * Size + column] = augmented[
                row * 2 * Size + Size + column]
                / (scales[row] * scales[column]);
        }
    }
    return result;
}

std::array<double, 7> lowRankRow(
    const MimeticLocalCellOperator& localOperator,
    const MimeticWallCondensation& condensation,
    const std::size_t face) {
    const double area = localOperator.faceAreasSquareMeters[face];
    std::array<double, 7> result{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        result[axis] = area
            * localOperator.normalRows[face * 3 + axis];
        result[3 + axis] = area
            * localOperator.consistencyRows[face * 3 + axis];
    }
    result[6] = condensation.conservationCoupling[face];
    return result;
}

std::array<double, 7> multiplySignedCore(
    const MimeticLocalCellOperator& localOperator,
    const MimeticWallCondensation& condensation,
    const std::array<double, 7>& value) {
    std::array<double, 7> result{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            result[row] += localOperator.inverseConsistencyGeometry[
                row * 3 + column] * value[column];
            result[3 + row] -=
                localOperator.stabilizationScaleInverseCubicMeters
                * localOperator.inverseConsistencyGram[
                    row * 3 + column] * value[3 + column];
        }
    }
    result[6] = -value[6] / condensation.conservationDenominator;
    return result;
}

double quadratic7(const std::array<double, 7>& value,
                  const std::array<double, 49>& matrix) {
    CompensatedSum result;
    for (std::size_t row = 0; row < 7; ++row) {
        for (std::size_t column = 0; column < 7; ++column) {
            result.add(value[row] * matrix[row * 7 + column]
                       * value[column]);
        }
    }
    return result.value();
}

std::vector<double> applyFullLocalTraceOperator(
    const MimeticLocalCellOperator& localOperator,
    const std::span<const double> traces) {
    auto result = balanceMimeticLocalCell(
        localOperator, traces, 0.0).integratedOutwardFluxes;
    for (double& value : result) value = -value;
    return result;
}

std::vector<double> solveWallBlock(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    const std::span<const double> rightHandSide) {
    std::array<CompensatedSum, 7> transposeScaledRightHandSide;
    for (std::size_t face = 0;
         face < condensation.halfFaceCount; ++face) {
        if (condensation.wallMask[face] == 0) continue;
        const double area = localOperator.faceAreasSquareMeters[face];
        const double inverseDiagonal = 1.0
            / (localOperator.stabilizationScaleInverseCubicMeters
               * area * area);
        const auto row = lowRankRow(localOperator, condensation, face);
        for (std::size_t mode = 0; mode < 7; ++mode) {
            transposeScaledRightHandSide[mode].add(
                row[mode] * inverseDiagonal * rightHandSide[face]);
        }
    }
    std::array<double, 7> coefficient{};
    for (std::size_t row = 0; row < 7; ++row) {
        for (std::size_t column = 0; column < 7; ++column) {
            coefficient[row] += condensation.inverseWoodburyCore[
                row * 7 + column]
                * transposeScaledRightHandSide[column].value();
        }
    }
    std::vector<double> result(condensation.halfFaceCount, 0.0);
    for (std::size_t face = 0;
         face < condensation.halfFaceCount; ++face) {
        if (condensation.wallMask[face] == 0) continue;
        const double area = localOperator.faceAreasSquareMeters[face];
        const double inverseDiagonal = 1.0
            / (localOperator.stabilizationScaleInverseCubicMeters
               * area * area);
        const auto row = lowRankRow(localOperator, condensation, face);
        CompensatedSum correction;
        for (std::size_t mode = 0; mode < 7; ++mode) {
            correction.add(row[mode] * coefficient[mode]);
        }
        result[face] = inverseDiagonal
            * (rightHandSide[face] - correction.value());
        if (!std::isfinite(result[face])) {
            throw std::overflow_error(
                "mimetic wall-block solve overflowed");
        }
    }
    return result;
}

void validateField(const MimeticWallCondensation& condensation,
                   const std::span<const double> values,
                   const bool requireZeroWall) {
    if (values.size() != condensation.halfFaceCount
        || !std::ranges::all_of(values, [](const double value) {
               return std::isfinite(value);
           })) {
        throw std::invalid_argument(
            "mimetic wall-condensation field is invalid");
    }
    if (requireZeroWall) {
        for (std::size_t face = 0; face < values.size(); ++face) {
            if (condensation.wallMask[face] != 0 && values[face] != 0.0) {
                throw std::invalid_argument(
                    "mimetic wall-condensation active field has wall data");
            }
        }
    }
}

} // namespace

MimeticWallCondensation buildMimeticWallCondensation(
    const MimeticLocalCellOperator& localOperator,
    const std::span<const std::uint8_t> wallMask,
    const MimeticWallCondensationSettings& settings) {
    validateMimeticLocalCellOperator(localOperator);
    if (wallMask.size() != localOperator.halfFaceCount) {
        throw std::invalid_argument(
            "mimetic wall-condensation mask size is invalid");
    }
    if (localOperator.halfFaceCount > settings.maximumHalfFaces) {
        throw std::length_error(
            "mimetic wall-condensation face limit exceeded");
    }
    MimeticWallCondensation result;
    result.localOperatorFingerprint = localOperator.fingerprint;
    result.halfFaceCount = localOperator.halfFaceCount;
    result.wallMask.assign(wallMask.begin(), wallMask.end());
    for (const std::uint8_t value : result.wallMask) {
        if (value > 1) {
            throw std::invalid_argument(
                "mimetic wall-condensation mask is invalid");
        }
        result.wallHalfFaceCount += value != 0;
    }
    result.activeHalfFaceCount = result.halfFaceCount
        - result.wallHalfFaceCount;
    if (result.activeHalfFaceCount == 0) {
        throw std::invalid_argument(
            "mimetic wall condensation cannot eliminate every half-face");
    }
    result.ownedStorageBytes = ownedBytes(result.halfFaceCount);
    if (result.ownedStorageBytes > settings.maximumOwnedBytes) {
        throw std::length_error(
            "mimetic wall-condensation byte limit exceeded");
    }

    const auto inverseFluxArea = applyMimeticInverseFluxInnerProduct(
        localOperator, localOperator.faceAreasSquareMeters);
    result.conservationCoupling.resize(result.halfFaceCount, 0.0);
    CompensatedSum denominator;
    for (std::size_t face = 0; face < result.halfFaceCount; ++face) {
        result.conservationCoupling[face] =
            localOperator.faceAreasSquareMeters[face]
            * inverseFluxArea[face];
        denominator.add(result.conservationCoupling[face]);
    }
    result.conservationDenominator = denominator.value();
    if (!std::isfinite(result.conservationDenominator)
        || result.conservationDenominator <= 0.0) {
        throw std::invalid_argument(
            "mimetic wall-condensation denominator is invalid");
    }

    std::array<double, 49> gram{};
    if (result.wallHalfFaceCount != 0) {
        for (std::size_t face = 0; face < result.halfFaceCount; ++face) {
            if (result.wallMask[face] == 0) continue;
            const double area = localOperator.faceAreasSquareMeters[face];
            const double inverseDiagonal = 1.0
                / (localOperator.stabilizationScaleInverseCubicMeters
                   * area * area);
            const auto row = lowRankRow(localOperator, result, face);
            for (std::size_t first = 0; first < 7; ++first) {
                for (std::size_t second = 0; second < 7; ++second) {
                    gram[first * 7 + second] +=
                        row[first] * inverseDiagonal * row[second];
                }
            }
        }
        const auto geometry = invertEquilibrated<3>(
            localOperator.inverseConsistencyGeometry,
            "mimetic wall-condensation geometry inverse is singular");
        const auto consistencyGram = invertEquilibrated<3>(
            localOperator.inverseConsistencyGram,
            "mimetic wall-condensation consistency inverse is singular");
        std::array<double, 49> woodburyCore = gram;
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t column = 0; column < 3; ++column) {
                woodburyCore[row * 7 + column] +=
                    geometry[row * 3 + column];
                woodburyCore[(3 + row) * 7 + 3 + column] -=
                    consistencyGram[row * 3 + column]
                    / localOperator.stabilizationScaleInverseCubicMeters;
            }
        }
        woodburyCore[6 * 7 + 6] -= result.conservationDenominator;
        result.inverseWoodburyCore = invertEquilibrated<7>(
            woodburyCore,
            "mimetic wall-condensation wall block is singular");

        std::array<double, 49> gramInverseCore{};
        for (std::size_t row = 0; row < 7; ++row) {
            for (std::size_t column = 0; column < 7; ++column) {
                for (std::size_t middle = 0; middle < 7; ++middle) {
                    gramInverseCore[row * 7 + column] +=
                        gram[row * 7 + middle]
                        * result.inverseWoodburyCore[
                            middle * 7 + column];
                }
            }
        }
        for (std::size_t row = 0; row < 7; ++row) {
            for (std::size_t column = 0; column < 7; ++column) {
                double value = gram[row * 7 + column];
                for (std::size_t middle = 0; middle < 7; ++middle) {
                    value -= gramInverseCore[row * 7 + middle]
                        * gram[middle * 7 + column];
                }
                result.wallSchurMetric[row * 7 + column] = value;
            }
        }
    }

    const auto inverseFluxDiagonal =
        mimeticInverseFluxInnerProductDiagonal(localOperator);
    result.condensedOperatorDiagonal.resize(result.halfFaceCount, 0.0);
    for (std::size_t face = 0; face < result.halfFaceCount; ++face) {
        if (result.wallMask[face] != 0) continue;
        const double area = localOperator.faceAreasSquareMeters[face];
        const double fullDiagonal = area * area
            * inverseFluxDiagonal[face]
            - result.conservationCoupling[face]
                * result.conservationCoupling[face]
                / result.conservationDenominator;
        double correction = 0.0;
        if (result.wallHalfFaceCount != 0) {
            correction = quadratic7(
                multiplySignedCore(
                    localOperator, result,
                    lowRankRow(localOperator, result, face)),
                result.wallSchurMetric);
        }
        double diagonal = fullDiagonal - correction;
        const double tolerance = 8192.0
            * std::numeric_limits<double>::epsilon()
            * std::max({1.0, std::abs(fullDiagonal), std::abs(correction)});
        if (diagonal < 0.0 && diagonal >= -tolerance) diagonal = 0.0;
        if (!std::isfinite(diagonal) || diagonal < 0.0) {
            throw std::invalid_argument(
                "mimetic wall-condensation diagonal is invalid");
        }
        result.condensedOperatorDiagonal[face] = diagonal;
    }
    result.fingerprint = productFingerprint(result);
    validateMimeticWallCondensation(result, localOperator);
    return result;
}

void validateMimeticWallCondensation(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator) {
    validateMimeticLocalCellOperator(localOperator);
    if (condensation.version != mimeticWallCondensationVersion
        || condensation.fingerprint == 0
        || condensation.localOperatorFingerprint != localOperator.fingerprint
        || condensation.halfFaceCount != localOperator.halfFaceCount
        || condensation.wallMask.size() != condensation.halfFaceCount
        || condensation.conservationCoupling.size()
            != condensation.halfFaceCount
        || condensation.condensedOperatorDiagonal.size()
            != condensation.halfFaceCount
        || condensation.wallHalfFaceCount
                + condensation.activeHalfFaceCount
            != condensation.halfFaceCount
        || condensation.activeHalfFaceCount == 0
        || condensation.ownedStorageBytes
            != ownedBytes(condensation.halfFaceCount)
        || !std::isfinite(condensation.conservationDenominator)
        || condensation.conservationDenominator <= 0.0) {
        throw std::invalid_argument(
            "mimetic wall-condensation integrity is invalid");
    }
    std::size_t wallCount = 0;
    for (std::size_t face = 0;
         face < condensation.halfFaceCount; ++face) {
        const std::uint8_t wall = condensation.wallMask[face];
        const double coupling = condensation.conservationCoupling[face];
        const double diagonal =
            condensation.condensedOperatorDiagonal[face];
        if (wall > 1 || !std::isfinite(coupling)
            || !std::isfinite(diagonal) || diagonal < 0.0
            || (wall != 0 && diagonal != 0.0)) {
            throw std::invalid_argument(
                "mimetic wall-condensation face data is invalid");
        }
        wallCount += wall != 0;
    }
    if (wallCount != condensation.wallHalfFaceCount
        || !std::ranges::all_of(
            condensation.inverseWoodburyCore,
            [](const double value) { return std::isfinite(value); })
        || !std::ranges::all_of(
            condensation.wallSchurMetric,
            [](const double value) { return std::isfinite(value); })
        || productFingerprint(condensation) != condensation.fingerprint) {
        throw std::invalid_argument(
            "mimetic wall-condensation fingerprint is invalid");
    }
}

std::vector<double> applyMimeticWallCondensedTraceOperator(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    const std::span<const double> activeTraceValues) {
    validateMimeticWallCondensation(condensation, localOperator);
    validateField(condensation, activeTraceValues, true);
    auto first = applyFullLocalTraceOperator(
        localOperator, activeTraceValues);
    if (condensation.wallHalfFaceCount == 0) return first;
    const auto wallSolution = solveWallBlock(
        condensation, localOperator, first);
    std::vector<double> correction(condensation.halfFaceCount, 0.0);
    for (std::size_t face = 0; face < correction.size(); ++face) {
        if (condensation.wallMask[face] != 0) {
            correction[face] = -wallSolution[face];
        }
    }
    const auto second = applyFullLocalTraceOperator(
        localOperator, correction);
    for (std::size_t face = 0; face < first.size(); ++face) {
        first[face] = condensation.wallMask[face] == 0
            ? first[face] + second[face] : 0.0;
        if (!std::isfinite(first[face])) {
            throw std::overflow_error(
                "mimetic wall-condensed action overflowed");
        }
    }
    return first;
}

std::vector<double> condenseMimeticWallTraceRightHandSide(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    const std::span<const double> fullRightHandSide) {
    validateMimeticWallCondensation(condensation, localOperator);
    validateField(condensation, fullRightHandSide, false);
    if (condensation.wallHalfFaceCount == 0) {
        return std::vector<double>(
            fullRightHandSide.begin(), fullRightHandSide.end());
    }
    const auto wallSolution = solveWallBlock(
        condensation, localOperator, fullRightHandSide);
    std::vector<double> correction(condensation.halfFaceCount, 0.0);
    for (std::size_t face = 0; face < correction.size(); ++face) {
        if (condensation.wallMask[face] != 0) {
            correction[face] = -wallSolution[face];
        }
    }
    const auto action = applyFullLocalTraceOperator(
        localOperator, correction);
    std::vector<double> result(condensation.halfFaceCount, 0.0);
    for (std::size_t face = 0; face < result.size(); ++face) {
        if (condensation.wallMask[face] == 0) {
            result[face] = fullRightHandSide[face] + action[face];
            if (!std::isfinite(result[face])) {
                throw std::overflow_error(
                    "mimetic wall-condensed right-hand side overflowed");
            }
        }
    }
    return result;
}

std::vector<double> reconstructMimeticWallTraces(
    const MimeticWallCondensation& condensation,
    const MimeticLocalCellOperator& localOperator,
    const std::span<const double> fullRightHandSide,
    const std::span<const double> activeTraceValues) {
    validateMimeticWallCondensation(condensation, localOperator);
    validateField(condensation, fullRightHandSide, false);
    validateField(condensation, activeTraceValues, true);
    if (condensation.wallHalfFaceCount == 0) {
        return std::vector<double>(
            activeTraceValues.begin(), activeTraceValues.end());
    }
    const auto action = applyFullLocalTraceOperator(
        localOperator, activeTraceValues);
    std::vector<double> wallRightHandSide(
        condensation.halfFaceCount, 0.0);
    for (std::size_t face = 0; face < wallRightHandSide.size(); ++face) {
        if (condensation.wallMask[face] != 0) {
            wallRightHandSide[face] = fullRightHandSide[face]
                - action[face];
        }
    }
    const auto wallSolution = solveWallBlock(
        condensation, localOperator, wallRightHandSide);
    std::vector<double> result(
        activeTraceValues.begin(), activeTraceValues.end());
    for (std::size_t face = 0; face < result.size(); ++face) {
        if (condensation.wallMask[face] != 0) {
            result[face] = wallSolution[face];
        }
    }
    return result;
}

} // namespace simwing::fsi::fluid
