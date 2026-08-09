#include "mimetic_local_cell.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace simwing::fsi::fluid {
namespace {

using Matrix3 = std::array<double, 9>;

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

class Fingerprint final {
public:
    void integer(const std::uint64_t value) {
        for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
            value_ ^= static_cast<std::uint8_t>(value >> (8 * byte));
            value_ *= fnvPrime;
        }
    }

    void real(const double value) {
        integer(std::bit_cast<std::uint64_t>(value));
    }

    [[nodiscard]] std::uint64_t value() const noexcept {
        return value_;
    }

private:
    std::uint64_t value_ = fnvOffsetBasis;
};

bool finite(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double component(const Vector3& value, const std::size_t axis) {
    return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

double dot(const Vector3& first, const Vector3& second) {
    return first.x * second.x
        + first.y * second.y
        + first.z * second.z;
}

Vector3 subtract(const Vector3& first, const Vector3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z,
    };
}

Matrix3 symmetricInverse(const Matrix3& matrix,
                         const char* singularMessage) {
    Matrix3 lower{};
    const double maximumDiagonal = std::max({
        matrix[0], matrix[4], matrix[8],
    });
    if (!std::isfinite(maximumDiagonal) || maximumDiagonal <= 0.0) {
        throw std::invalid_argument(singularMessage);
    }
    const double pivotTolerance = maximumDiagonal
        * 128.0 * std::numeric_limits<double>::epsilon();
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
            double value = matrix[row * 3 + column];
            for (std::size_t inner = 0; inner < column; ++inner) {
                value -= lower[row * 3 + inner]
                    * lower[column * 3 + inner];
            }
            if (row == column) {
                if (!std::isfinite(value) || value <= pivotTolerance) {
                    throw std::invalid_argument(singularMessage);
                }
                lower[row * 3 + column] = std::sqrt(value);
            } else {
                lower[row * 3 + column] = value
                    / lower[column * 3 + column];
            }
        }
    }

    Matrix3 inverse{};
    for (std::size_t rightHandSide = 0;
         rightHandSide < 3; ++rightHandSide) {
        std::array<double, 3> intermediate{};
        for (std::size_t row = 0; row < 3; ++row) {
            double value = row == rightHandSide ? 1.0 : 0.0;
            for (std::size_t inner = 0; inner < row; ++inner) {
                value -= lower[row * 3 + inner] * intermediate[inner];
            }
            intermediate[row] = value / lower[row * 3 + row];
        }
        std::array<double, 3> solution{};
        for (std::size_t reverse = 3; reverse > 0; --reverse) {
            const std::size_t row = reverse - 1;
            double value = intermediate[row];
            for (std::size_t inner = row + 1; inner < 3; ++inner) {
                value -= lower[inner * 3 + row] * solution[inner];
            }
            solution[row] = value / lower[row * 3 + row];
        }
        for (std::size_t row = 0; row < 3; ++row) {
            inverse[row * 3 + rightHandSide] = solution[row];
        }
    }
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = row + 1; column < 3; ++column) {
            const double value = 0.5 * (
                inverse[row * 3 + column]
                + inverse[column * 3 + row]);
            inverse[row * 3 + column] = value;
            inverse[column * 3 + row] = value;
        }
    }
    return inverse;
}

double bilinear3(const std::array<double, 3>& first,
                 const Matrix3& matrix,
                 const std::array<double, 3>& second) {
    double result = 0.0;
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            result += first[row] * matrix[row * 3 + column]
                * second[column];
        }
    }
    return result;
}

void validateSettings(const MimeticLocalCellSettings& settings) {
    if (!std::isfinite(
            settings.absoluteAreaClosureToleranceSquareMeters)
        || settings.absoluteAreaClosureToleranceSquareMeters < 0.0
        || !std::isfinite(
            settings.absoluteDivergenceTheoremToleranceCubicMeters)
        || settings.absoluteDivergenceTheoremToleranceCubicMeters < 0.0
        || !std::isfinite(settings.relativeGeometryTolerance)
        || settings.relativeGeometryTolerance < 0.0
        || !std::isfinite(settings.unitNormalTolerance)
        || settings.unitNormalTolerance < 0.0
        || !std::isfinite(settings.algebraicConsistencyTolerance)
        || settings.algebraicConsistencyTolerance < 0.0
        || settings.maximumHalfFaces < 4
        || settings.maximumOperatorBytes == 0) {
        throw std::invalid_argument(
            "invalid mimetic local-cell settings");
    }
}

std::size_t checkedOperatorBytes(const std::size_t halfFaceCount) {
    if (halfFaceCount > std::numeric_limits<std::size_t>::max()
            / halfFaceCount) {
        throw std::length_error(
            "mimetic local-cell matrix size overflows");
    }
    const std::size_t matrixEntries = halfFaceCount * halfFaceCount;
    const std::size_t maximumEntries =
        std::numeric_limits<std::size_t>::max() / sizeof(double);
    if (halfFaceCount > maximumEntries
        || matrixEntries > maximumEntries - halfFaceCount) {
        throw std::length_error(
            "mimetic local-cell matrix size overflows");
    }
    return (matrixEntries + halfFaceCount) * sizeof(double);
}

void validateInputScalars(const MimeticLocalCellOperator& localOperator,
                          const double cellScalar,
                          const std::span<const double> faceTraceScalars) {
    validateMimeticLocalCellOperator(localOperator);
    if (faceTraceScalars.size() != localOperator.halfFaceCount
        || !std::isfinite(cellScalar)
        || !std::ranges::all_of(
            faceTraceScalars,
            [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "invalid mimetic local-cell scalar field");
    }
}

std::uint64_t operatorFingerprint(
    const MimeticLocalCellOperator& localOperator) {
    Fingerprint fingerprint;
    fingerprint.integer(localOperator.version);
    fingerprint.integer(localOperator.halfFaceCount);
    fingerprint.integer(localOperator.ownedStorageBytes);
    fingerprint.real(localOperator.volumeCubicMeters);
    fingerprint.real(
        localOperator.stabilizationScaleInverseCubicMeters);
    fingerprint.real(
        localOperator.maximumAreaClosureErrorSquareMeters);
    fingerprint.real(
        localOperator.maximumDivergenceTheoremErrorCubicMeters);
    fingerprint.real(
        localOperator.maximumAlgebraicConsistencyError);
    fingerprint.integer(localOperator.faceAreasSquareMeters.size());
    for (const double area : localOperator.faceAreasSquareMeters) {
        fingerprint.real(area);
    }
    fingerprint.integer(localOperator.inverseFluxInnerProduct.size());
    for (const double value : localOperator.inverseFluxInnerProduct) {
        fingerprint.real(value);
    }
    return fingerprint.value();
}

} // namespace

MimeticLocalCellOperator buildMimeticLocalCellOperator(
    const MimeticLocalCellGeometry& geometry,
    const MimeticLocalCellSettings& settings) {
    validateSettings(settings);
    const std::size_t count = geometry.halfFaces.size();
    if (count < 4) {
        throw std::invalid_argument(
            "mimetic local cell requires at least four half-faces");
    }
    if (count > settings.maximumHalfFaces) {
        throw std::length_error(
            "mimetic local-cell half-face limit exceeded");
    }
    const std::size_t ownedStorageBytes = checkedOperatorBytes(count);
    if (ownedStorageBytes > settings.maximumOperatorBytes) {
        throw std::length_error(
            "mimetic local-cell matrix byte limit exceeded");
    }
    if (!std::isfinite(geometry.volumeCubicMeters)
        || geometry.volumeCubicMeters <= 0.0
        || !finite(geometry.centroidMeters)) {
        throw std::invalid_argument(
            "invalid mimetic local-cell volume or centroid");
    }

    std::vector<std::array<double, 3>> consistencyRows(count);
    std::vector<std::array<double, 3>> normalRows(count);
    std::vector<double> areas(count);
    Vector3 areaClosure;
    double summedArea = 0.0;
    for (std::size_t face = 0; face < count; ++face) {
        const auto& source = geometry.halfFaces[face];
        if (!std::isfinite(source.areaSquareMeters)
            || source.areaSquareMeters <= 0.0
            || !finite(source.centroidMeters)
            || !finite(source.outwardUnitNormal)) {
            throw std::invalid_argument(
                "invalid mimetic half-face geometry");
        }
        const double normalLength = std::sqrt(
            dot(source.outwardUnitNormal, source.outwardUnitNormal));
        if (!std::isfinite(normalLength)
            || std::abs(normalLength - 1.0)
                > settings.unitNormalTolerance) {
            throw std::invalid_argument(
                "mimetic half-face normal is not unit length");
        }
        areas[face] = source.areaSquareMeters;
        summedArea += source.areaSquareMeters;
        areaClosure.x += source.areaSquareMeters
            * source.outwardUnitNormal.x;
        areaClosure.y += source.areaSquareMeters
            * source.outwardUnitNormal.y;
        areaClosure.z += source.areaSquareMeters
            * source.outwardUnitNormal.z;
        const Vector3 offset = subtract(
            source.centroidMeters, geometry.centroidMeters);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            consistencyRows[face][axis] = source.areaSquareMeters
                * component(offset, axis);
            normalRows[face][axis] = component(
                source.outwardUnitNormal, axis);
        }
    }
    if (!std::isfinite(summedArea)) {
        throw std::invalid_argument(
            "non-finite mimetic local-cell surface area");
    }
    const double maximumAreaClosureError = std::max({
        std::abs(areaClosure.x),
        std::abs(areaClosure.y),
        std::abs(areaClosure.z),
    });
    const double areaClosureTolerance =
        settings.absoluteAreaClosureToleranceSquareMeters
        + settings.relativeGeometryTolerance * summedArea;
    if (maximumAreaClosureError > areaClosureTolerance) {
        throw std::invalid_argument(
            "mimetic local-cell half-faces do not close");
    }

    Matrix3 normalTransposeConsistency{};
    Matrix3 consistencyTransposeConsistency{};
    for (std::size_t face = 0; face < count; ++face) {
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t column = 0; column < 3; ++column) {
                normalTransposeConsistency[row * 3 + column] +=
                    normalRows[face][row]
                    * consistencyRows[face][column];
                consistencyTransposeConsistency[row * 3 + column] +=
                    consistencyRows[face][row]
                    * consistencyRows[face][column];
            }
        }
    }
    double maximumDivergenceTheoremError = 0.0;
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            const double expected = row == column
                ? geometry.volumeCubicMeters : 0.0;
            maximumDivergenceTheoremError = std::max(
                maximumDivergenceTheoremError,
                std::abs(normalTransposeConsistency[row * 3 + column]
                         - expected));
        }
    }
    const double volumeTolerance =
        settings.absoluteDivergenceTheoremToleranceCubicMeters
        + settings.relativeGeometryTolerance
            * geometry.volumeCubicMeters;
    if (maximumDivergenceTheoremError > volumeTolerance) {
        throw std::invalid_argument(
            "mimetic local-cell geometry violates the divergence theorem");
    }

    Matrix3 symmetricGeometry{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            symmetricGeometry[row * 3 + column] = 0.5 * (
                normalTransposeConsistency[row * 3 + column]
                + normalTransposeConsistency[column * 3 + row]);
        }
    }
    const Matrix3 inverseGeometry = symmetricInverse(
        symmetricGeometry,
        "singular mimetic local-cell consistency geometry");
    const Matrix3 inverseConsistencyGram = symmetricInverse(
        consistencyTransposeConsistency,
        "rank-deficient mimetic local-cell consistency rows");

    std::vector<double> consistentMatrix(count * count, 0.0);
    std::vector<double> nullspaceProjector(count * count, 0.0);
    double consistentTrace = 0.0;
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t column = 0; column < count; ++column) {
            const double consistent = bilinear3(
                normalRows[row], inverseGeometry, normalRows[column]);
            const double projected = (row == column ? 1.0 : 0.0)
                - bilinear3(
                    consistencyRows[row], inverseConsistencyGram,
                    consistencyRows[column]);
            consistentMatrix[row * count + column] = consistent;
            nullspaceProjector[row * count + column] = projected;
        }
        consistentTrace += consistentMatrix[row * count + row];
    }
    const double stabilizationScale = consistentTrace
        / static_cast<double>(count - 3);
    if (!std::isfinite(stabilizationScale)
        || stabilizationScale <= 0.0) {
        throw std::invalid_argument(
            "invalid mimetic local-cell stabilization scale");
    }

    std::vector<double> matrix(count * count, 0.0);
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t column = row; column < count; ++column) {
            const double forward = consistentMatrix[row * count + column]
                + stabilizationScale
                    * nullspaceProjector[row * count + column];
            const double reverse = consistentMatrix[column * count + row]
                + stabilizationScale
                    * nullspaceProjector[column * count + row];
            const double value = 0.5 * (forward + reverse);
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "non-finite mimetic local-cell matrix");
            }
            matrix[row * count + column] = value;
            matrix[column * count + row] = value;
        }
    }

    double maximumAlgebraicConsistencyError = 0.0;
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            double value = 0.0;
            for (std::size_t column = 0; column < count; ++column) {
                value += matrix[row * count + column]
                    * consistencyRows[column][axis];
            }
            maximumAlgebraicConsistencyError = std::max(
                maximumAlgebraicConsistencyError,
                std::abs(value - normalRows[row][axis]));
        }
    }
    if (maximumAlgebraicConsistencyError
        > settings.algebraicConsistencyTolerance) {
        throw std::invalid_argument(
            "mimetic local-cell matrix failed linear consistency");
    }

    MimeticLocalCellOperator result;
    result.halfFaceCount = count;
    result.ownedStorageBytes = ownedStorageBytes;
    result.volumeCubicMeters = geometry.volumeCubicMeters;
    result.stabilizationScaleInverseCubicMeters = stabilizationScale;
    result.maximumAreaClosureErrorSquareMeters =
        maximumAreaClosureError;
    result.maximumDivergenceTheoremErrorCubicMeters =
        maximumDivergenceTheoremError;
    result.maximumAlgebraicConsistencyError =
        maximumAlgebraicConsistencyError;
    result.faceAreasSquareMeters = std::move(areas);
    result.inverseFluxInnerProduct = std::move(matrix);
    result.fingerprint = operatorFingerprint(result);
    validateMimeticLocalCellOperator(result);
    return result;
}

void validateMimeticLocalCellOperator(
    const MimeticLocalCellOperator& localOperator) {
    const std::size_t count = localOperator.halfFaceCount;
    if (localOperator.version != mimeticLocalCellVersion
        || localOperator.fingerprint == 0
        || count < 4
        || localOperator.faceAreasSquareMeters.size() != count
        || count > std::numeric_limits<std::size_t>::max() / count
        || localOperator.inverseFluxInnerProduct.size() != count * count
        || localOperator.ownedStorageBytes != checkedOperatorBytes(count)
        || !std::isfinite(localOperator.volumeCubicMeters)
        || localOperator.volumeCubicMeters <= 0.0
        || !std::isfinite(
            localOperator.stabilizationScaleInverseCubicMeters)
        || localOperator.stabilizationScaleInverseCubicMeters <= 0.0
        || !std::isfinite(
            localOperator.maximumAreaClosureErrorSquareMeters)
        || localOperator.maximumAreaClosureErrorSquareMeters < 0.0
        || !std::isfinite(
            localOperator.maximumDivergenceTheoremErrorCubicMeters)
        || localOperator.maximumDivergenceTheoremErrorCubicMeters < 0.0
        || !std::isfinite(
            localOperator.maximumAlgebraicConsistencyError)
        || localOperator.maximumAlgebraicConsistencyError < 0.0) {
        throw std::invalid_argument(
            "invalid mimetic local-cell operator");
    }
    for (const double area : localOperator.faceAreasSquareMeters) {
        if (!std::isfinite(area) || area <= 0.0) {
            throw std::invalid_argument(
                "invalid mimetic local-cell operator area");
        }
    }
    for (std::size_t row = 0; row < count; ++row) {
        if (localOperator.inverseFluxInnerProduct[row * count + row]
            <= 0.0) {
            throw std::invalid_argument(
                "invalid mimetic local-cell matrix diagonal");
        }
        for (std::size_t column = 0; column < count; ++column) {
            const double value =
                localOperator.inverseFluxInnerProduct[
                    row * count + column];
            if (!std::isfinite(value)
                || value != localOperator.inverseFluxInnerProduct[
                    column * count + row]) {
                throw std::invalid_argument(
                    "invalid mimetic local-cell symmetric matrix");
            }
        }
    }
    if (operatorFingerprint(localOperator) != localOperator.fingerprint) {
        throw std::invalid_argument(
            "mimetic local-cell operator fingerprint mismatch");
    }
}

std::vector<double> applyMimeticLocalNormalFlux(
    const MimeticLocalCellOperator& localOperator,
    const double cellScalar,
    const std::span<const double> faceTraceScalars) {
    validateInputScalars(
        localOperator, cellScalar, faceTraceScalars);
    const std::size_t count = localOperator.halfFaceCount;
    std::vector<double> weightedDifferences(count, 0.0);
    for (std::size_t face = 0; face < count; ++face) {
        weightedDifferences[face] =
            localOperator.faceAreasSquareMeters[face]
            * (faceTraceScalars[face] - cellScalar);
    }
    std::vector<double> result(count, 0.0);
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t column = 0; column < count; ++column) {
            result[row] -= localOperator.inverseFluxInnerProduct[
                row * count + column] * weightedDifferences[column];
        }
        if (!std::isfinite(result[row])) {
            throw std::invalid_argument(
                "non-finite mimetic local-cell normal flux");
        }
    }
    return result;
}

MimeticLocalCellBalance balanceMimeticLocalCell(
    const MimeticLocalCellOperator& localOperator,
    const std::span<const double> faceTraceScalars,
    const double integratedSource) {
    validateInputScalars(localOperator, 0.0, faceTraceScalars);
    if (!std::isfinite(integratedSource)) {
        throw std::invalid_argument(
            "non-finite mimetic local-cell integrated source");
    }
    const std::size_t count = localOperator.halfFaceCount;
    std::vector<double> matrixArea(count, 0.0);
    std::vector<double> matrixWeightedTrace(count, 0.0);
    for (std::size_t row = 0; row < count; ++row) {
        for (std::size_t column = 0; column < count; ++column) {
            const double matrix =
                localOperator.inverseFluxInnerProduct[
                    row * count + column];
            matrixArea[row] += matrix
                * localOperator.faceAreasSquareMeters[column];
            matrixWeightedTrace[row] += matrix
                * localOperator.faceAreasSquareMeters[column]
                * faceTraceScalars[column];
        }
    }
    double denominator = 0.0;
    double traceCoupling = 0.0;
    for (std::size_t face = 0; face < count; ++face) {
        const double area = localOperator.faceAreasSquareMeters[face];
        denominator += area * matrixArea[face];
        traceCoupling += area * matrixWeightedTrace[face];
    }
    if (!std::isfinite(denominator) || denominator <= 0.0
        || !std::isfinite(traceCoupling)) {
        throw std::invalid_argument(
            "invalid mimetic local-cell conservation denominator");
    }
    MimeticLocalCellBalance result;
    result.cellScalar = (integratedSource + traceCoupling)
        / denominator;
    result.requestedIntegratedSource = integratedSource;
    result.outwardNormalFluxes = applyMimeticLocalNormalFlux(
        localOperator, result.cellScalar, faceTraceScalars);
    result.integratedOutwardFluxes.resize(count);
    for (std::size_t face = 0; face < count; ++face) {
        result.integratedOutwardFluxes[face] =
            localOperator.faceAreasSquareMeters[face]
            * result.outwardNormalFluxes[face];
        result.integratedOutwardFluxSum +=
            result.integratedOutwardFluxes[face];
    }
    result.conservationResidual =
        result.integratedOutwardFluxSum - integratedSource;
    if (!std::isfinite(result.cellScalar)
        || !std::isfinite(result.integratedOutwardFluxSum)
        || !std::isfinite(result.conservationResidual)) {
        throw std::invalid_argument(
            "non-finite mimetic local-cell balance");
    }
    return result;
}

} // namespace simwing::fsi::fluid
