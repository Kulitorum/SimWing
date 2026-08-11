#include "softwing/membrane.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace softwing {
namespace {

bool finite(const Matrix2& matrix) {
    return std::isfinite(matrix.m00) && std::isfinite(matrix.m01) &&
           std::isfinite(matrix.m10) && std::isfinite(matrix.m11);
}

bool finite(const SymmetricMatrix3& matrix) {
    return std::isfinite(matrix.xx) && std::isfinite(matrix.yy) &&
           std::isfinite(matrix.zz) && std::isfinite(matrix.xy) &&
           std::isfinite(matrix.xz) && std::isfinite(matrix.yz);
}

double maximumMagnitude(const SymmetricMatrix3& matrix) {
    return std::max({std::abs(matrix.xx),
                     std::abs(matrix.yy),
                     std::abs(matrix.zz),
                     std::abs(matrix.xy),
                     std::abs(matrix.xz),
                     std::abs(matrix.yz)});
}

Vec3 solveScaledPositiveDefinite(
    const SymmetricMatrix3& matrix,
    const Vec3& rightHandSide) {
    const double scaleX = std::sqrt(matrix.xx);
    const double scaleY = std::sqrt(matrix.yy);
    const double scaleZ = std::sqrt(matrix.zz);
    if (!(scaleX > 0.0) || !(scaleY > 0.0) || !(scaleZ > 0.0)
        || !std::isfinite(scaleX) || !std::isfinite(scaleY)
        || !std::isfinite(scaleZ)) {
        throw std::invalid_argument(
            "SymmetricMatrix3 is singular or ill-conditioned");
    }

    const double l10 = matrix.xy / scaleX / scaleY;
    const double l20 = matrix.xz / scaleX / scaleZ;
    const double firstPivotSquared = 1.0 - l10 * l10;
    constexpr double minimumScaledPivotSquared =
        4096.0 * std::numeric_limits<double>::epsilon();
    if (!std::isfinite(firstPivotSquared)
        || !(firstPivotSquared > minimumScaledPivotSquared)) {
        throw std::invalid_argument(
            "SymmetricMatrix3 is singular or ill-conditioned");
    }
    const double l11 = std::sqrt(firstPivotSquared);
    const double l21 =
        (matrix.yz / scaleY / scaleZ - l20 * l10) / l11;
    const double secondPivotSquared =
        1.0 - l20 * l20 - l21 * l21;
    if (!std::isfinite(l21) || !std::isfinite(secondPivotSquared)
        || !(secondPivotSquared > minimumScaledPivotSquared)) {
        throw std::invalid_argument(
            "SymmetricMatrix3 is singular or ill-conditioned");
    }
    const double l22 = std::sqrt(secondPivotSquared);

    const Vec3 scaledRightHandSide{
        rightHandSide.x / scaleX,
        rightHandSide.y / scaleY,
        rightHandSide.z / scaleZ,
    };
    const double forward0 = scaledRightHandSide.x;
    const double forward1 =
        (scaledRightHandSide.y - l10 * forward0) / l11;
    const double forward2 =
        (scaledRightHandSide.z - l20 * forward0 - l21 * forward1)
        / l22;
    const double scaledResultZ = forward2 / l22;
    const double scaledResultY =
        (forward1 - l21 * scaledResultZ) / l11;
    const double scaledResultX =
        forward0 - l10 * scaledResultY - l20 * scaledResultZ;
    const Vec3 result{
        scaledResultX / scaleX,
        scaledResultY / scaleY,
        scaledResultZ / scaleZ,
    };
    if (!std::isfinite(result.x) || !std::isfinite(result.y)
        || !std::isfinite(result.z)) {
        throw std::invalid_argument(
            "SymmetricMatrix3 solve result is non-finite");
    }
    const Vec3 residual = matrix * result - rightHandSide;
    const double rowScaleX =
        std::abs(matrix.xx * result.x)
        + std::abs(matrix.xy * result.y)
        + std::abs(matrix.xz * result.z)
        + std::abs(rightHandSide.x);
    const double rowScaleY =
        std::abs(matrix.xy * result.x)
        + std::abs(matrix.yy * result.y)
        + std::abs(matrix.yz * result.z)
        + std::abs(rightHandSide.y);
    const double rowScaleZ =
        std::abs(matrix.xz * result.x)
        + std::abs(matrix.yz * result.y)
        + std::abs(matrix.zz * result.z)
        + std::abs(rightHandSide.z);
    constexpr double maximumRelativeResidual =
        4096.0 * std::numeric_limits<double>::epsilon();
    if (!std::isfinite(residual.x) || !std::isfinite(residual.y)
        || !std::isfinite(residual.z)
        || !std::isfinite(rowScaleX) || !std::isfinite(rowScaleY)
        || !std::isfinite(rowScaleZ)
        || (rowScaleX == 0.0 ? residual.x != 0.0
                             : std::abs(residual.x)
                                 > maximumRelativeResidual * rowScaleX)
        || (rowScaleY == 0.0 ? residual.y != 0.0
                             : std::abs(residual.y)
                                 > maximumRelativeResidual * rowScaleY)
        || (rowScaleZ == 0.0 ? residual.z != 0.0
                             : std::abs(residual.z)
                                 > maximumRelativeResidual * rowScaleZ)) {
        throw std::invalid_argument(
            "SymmetricMatrix3 is singular or ill-conditioned");
    }
    return result;
}

} // namespace

Matrix2 checkedInverse(const Matrix2& matrix) {
    if (!finite(matrix)) {
        throw std::invalid_argument("Matrix2 must be finite");
    }
    const double scale = std::max({std::abs(matrix.m00),
                                   std::abs(matrix.m01),
                                   std::abs(matrix.m10),
                                   std::abs(matrix.m11)});
    const double determinant = matrix.determinant();
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() * scale * scale;
    if (!(scale > 0.0) || !std::isfinite(determinant) ||
        std::abs(determinant) <= tolerance) {
        throw std::invalid_argument("Matrix2 is singular or ill-conditioned");
    }
    const double inverseDeterminant = 1.0 / determinant;
    return {matrix.m11 * inverseDeterminant,
            -matrix.m01 * inverseDeterminant,
            -matrix.m10 * inverseDeterminant,
            matrix.m00 * inverseDeterminant};
}

double SymmetricMatrix3::determinant() const {
    return xx * (yy * zz - yz * yz) -
           xy * (xy * zz - yz * xz) +
           xz * (xy * yz - yy * xz);
}

Vec3 SymmetricMatrix3::operator*(const Vec3& value) const {
    return {xx * value.x + xy * value.y + xz * value.z,
            xy * value.x + yy * value.y + yz * value.z,
            xz * value.x + yz * value.y + zz * value.z};
}

SymmetricMatrix3 operator+(const SymmetricMatrix3& lhs,
                           const SymmetricMatrix3& rhs) {
    return {lhs.xx + rhs.xx,
            lhs.yy + rhs.yy,
            lhs.zz + rhs.zz,
            lhs.xy + rhs.xy,
            lhs.xz + rhs.xz,
            lhs.yz + rhs.yz};
}

SymmetricMatrix3 operator*(double scalar, const SymmetricMatrix3& matrix) {
    return {scalar * matrix.xx,
            scalar * matrix.yy,
            scalar * matrix.zz,
            scalar * matrix.xy,
            scalar * matrix.xz,
            scalar * matrix.yz};
}

bool isPositiveDefinite(const SymmetricMatrix3& matrix) {
    if (!finite(matrix) || !(matrix.xx > 0.0)
        || !(matrix.yy > 0.0) || !(matrix.zz > 0.0)) {
        return false;
    }
    const double scaleX = std::sqrt(matrix.xx);
    const double scaleY = std::sqrt(matrix.yy);
    const double scaleZ = std::sqrt(matrix.zz);
    const double l10 = matrix.xy / scaleX / scaleY;
    const double l20 = matrix.xz / scaleX / scaleZ;
    const double firstPivotSquared = 1.0 - l10 * l10;
    if (!std::isfinite(firstPivotSquared)
        || !(firstPivotSquared > 0.0)) {
        return false;
    }
    const double l11 = std::sqrt(firstPivotSquared);
    const double l21 =
        (matrix.yz / scaleY / scaleZ - l20 * l10) / l11;
    const double secondPivotSquared =
        1.0 - l20 * l20 - l21 * l21;
    return std::isfinite(secondPivotSquared)
        && secondPivotSquared > 0.0;
}

SymmetricMatrix3 checkedInverse(const SymmetricMatrix3& matrix) {
    if (!finite(matrix)) {
        throw std::invalid_argument("SymmetricMatrix3 must be finite");
    }
    const double scale = maximumMagnitude(matrix);
    const double determinant = matrix.determinant();
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() *
                             scale * scale * scale;
    if (!(scale > 0.0) || !std::isfinite(determinant) ||
        std::abs(determinant) <= tolerance) {
        throw std::invalid_argument(
            "SymmetricMatrix3 is singular or ill-conditioned");
    }

    const double inverseDeterminant = 1.0 / determinant;
    return {(matrix.yy * matrix.zz - matrix.yz * matrix.yz) *
                inverseDeterminant,
            (matrix.xx * matrix.zz - matrix.xz * matrix.xz) *
                inverseDeterminant,
            (matrix.xx * matrix.yy - matrix.xy * matrix.xy) *
                inverseDeterminant,
            (matrix.xz * matrix.yz - matrix.xy * matrix.zz) *
                inverseDeterminant,
            (matrix.xy * matrix.yz - matrix.xz * matrix.yy) *
                inverseDeterminant,
            (matrix.xy * matrix.xz - matrix.xx * matrix.yz) *
                inverseDeterminant};
}

Vec3 checkedSolve(const SymmetricMatrix3& matrix,
                  const Vec3& rightHandSide) {
    if (!std::isfinite(rightHandSide.x) || !std::isfinite(rightHandSide.y) ||
        !std::isfinite(rightHandSide.z)) {
        throw std::invalid_argument("Symmetric solve right-hand side must be finite");
    }
    try {
        // Preserve the historical arithmetic exactly for every system that
        // the explicit inverse already accepts.
        return checkedInverse(matrix) * rightHandSide;
    } catch (const std::invalid_argument&) {
        if (!finite(matrix) || !isPositiveDefinite(matrix)) {
            throw;
        }
        // Tiny authored membrane charts can create a large unscaled dynamic
        // range while the normalized SPD system remains well conditioned.
        // Diagonal equilibration followed by Cholesky solves that rejected
        // subset without weakening the inverse API or perturbing old solves.
        return solveScaledPositiveDefinite(matrix, rightHandSide);
    }
}

SymmetricMatrix3 OrthotropicMembraneMaterial::stiffnessMatrix() const {
    return {warpStiffness,
            weftStiffness,
            shearStiffness,
            couplingStiffness,
            0.0,
            0.0};
}

SymmetricMatrix3 OrthotropicMembraneMaterial::complianceMatrix() const {
    validateOrthotropicMembraneMaterial(*this);
    return checkedInverse(stiffnessMatrix());
}

SymmetricMatrix3 effectiveMembraneStiffness(
    const OrthotropicMembraneMaterial& material,
    const Vec3& greenStrain) {
    const SymmetricMatrix3 stiffness = material.stiffnessMatrix();
    // This branch is deliberately exact. Existing generic membrane users
    // retain both the same matrix and the same downstream arithmetic.
    if (material.compressionStiffnessRatio == 1.0) {
        return stiffness;
    }

    const double retained = std::sqrt(material.compressionStiffnessRatio);
    const double warpScale = greenStrain.x < 0.0 ? retained : 1.0;
    const double weftScale = greenStrain.y < 0.0 ? retained : 1.0;
    // Engineering shear belongs to both material directions. Its D entry is
    // the geometric mean of their normal entries: one compressed direction
    // softens shear partially, both soften it by the full retained factor.
    const double shearScale = std::sqrt(warpScale * weftScale);
    return {stiffness.xx * warpScale * warpScale,
            stiffness.yy * weftScale * weftScale,
            stiffness.zz * shearScale * shearScale,
            stiffness.xy * warpScale * weftScale,
            stiffness.xz * warpScale * shearScale,
            stiffness.yz * weftScale * shearScale};
}

void validateOrthotropicMembraneMaterial(
    const OrthotropicMembraneMaterial& material) {
    const double values[]{material.warpStiffness,
                          material.weftStiffness,
                          material.couplingStiffness,
                          material.shearStiffness,
                          material.warpPreTension,
                          material.weftPreTension,
                          material.dampingTime,
                          material.compressionStiffnessRatio};
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Membrane material values must be finite");
        }
    }
    if (!(material.warpStiffness > 0.0) ||
        !(material.weftStiffness > 0.0) ||
        !(material.shearStiffness > 0.0) ||
        material.couplingStiffness < 0.0 ||
        !(material.couplingStiffness * material.couplingStiffness <
          material.warpStiffness * material.weftStiffness) ||
        material.warpPreTension < 0.0 || material.warpPreTension >= 0.25 ||
        material.weftPreTension < 0.0 || material.weftPreTension >= 0.25 ||
        material.dampingTime < 0.0 ||
        !(material.compressionStiffnessRatio > 0.0) ||
        material.compressionStiffnessRatio > 1.0 ||
        !isPositiveDefinite(material.stiffnessMatrix())) {
        throw std::invalid_argument("Invalid orthotropic membrane material");
    }
}

} // namespace softwing
