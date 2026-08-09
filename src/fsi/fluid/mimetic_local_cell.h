#pragma once

#include "grid.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace simwing::fsi::fluid {

inline constexpr std::uint32_t mimeticLocalCellVersion = 1;

struct MimeticHalfFaceGeometry {
    double areaSquareMeters = 0.0;
    Vector3 centroidMeters;
    Vector3 outwardUnitNormal;

    bool operator==(const MimeticHalfFaceGeometry&) const = default;
};

struct MimeticLocalCellGeometry {
    double volumeCubicMeters = 0.0;
    Vector3 centroidMeters;
    std::vector<MimeticHalfFaceGeometry> halfFaces;

    bool operator==(const MimeticLocalCellGeometry&) const = default;
};

struct MimeticLocalCellSettings {
    double absoluteAreaClosureToleranceSquareMeters = 1.0e-12;
    double absoluteDivergenceTheoremToleranceCubicMeters = 1.0e-12;
    double relativeGeometryTolerance = 1.0e-10;
    double unitNormalTolerance = 1.0e-10;
    double algebraicConsistencyTolerance = 1.0e-10;
    std::size_t maximumHalfFaces = 4096;
    std::size_t maximumOperatorBytes = 256ULL * 1024ULL * 1024ULL;

    bool operator==(const MimeticLocalCellSettings&) const = default;
};

// One isotropic mixed-hybrid mimetic cell. The row-major matrix W is symmetric
// positive definite and maps area-weighted face-trace differences to outward
// normal flux:
//
//     u = -W diag(area) (lambda - p_cell 1).
//
// Its consistency identity is W R = N, where row f of R is
// area_f * (centroid_f - centroid_cell) and row f of N is the outward unit
// normal. Half-face centroids must use one unwrapped local periodic image.
// This kernel owns no global trace unknowns or production pressure topology.
struct MimeticLocalCellOperator {
    std::uint32_t version = mimeticLocalCellVersion;
    std::uint64_t fingerprint = 0;
    std::size_t halfFaceCount = 0;
    std::size_t ownedStorageBytes = 0;
    double volumeCubicMeters = 0.0;
    double stabilizationScaleInverseCubicMeters = 0.0;
    double maximumAreaClosureErrorSquareMeters = 0.0;
    double maximumDivergenceTheoremErrorCubicMeters = 0.0;
    double maximumAlgebraicConsistencyError = 0.0;
    std::vector<double> faceAreasSquareMeters;
    std::vector<double> inverseFluxInnerProduct;

    bool operator==(const MimeticLocalCellOperator&) const = default;
};

struct MimeticLocalCellBalance {
    double cellScalar = 0.0;
    double requestedIntegratedSource = 0.0;
    double integratedOutwardFluxSum = 0.0;
    double conservationResidual = 0.0;
    std::vector<double> outwardNormalFluxes;
    std::vector<double> integratedOutwardFluxes;

    bool operator==(const MimeticLocalCellBalance&) const = default;
};

[[nodiscard]] MimeticLocalCellOperator buildMimeticLocalCellOperator(
    const MimeticLocalCellGeometry& geometry,
    const MimeticLocalCellSettings& settings = {});

void validateMimeticLocalCellOperator(
    const MimeticLocalCellOperator& localOperator);

[[nodiscard]] std::vector<double> applyMimeticLocalNormalFlux(
    const MimeticLocalCellOperator& localOperator,
    double cellScalar,
    std::span<const double> faceTraceScalars);

// Eliminates the one cell scalar against exact integrated conservation while
// retaining the face traces. The returned outward integrated fluxes sum to the
// requested source up to reported floating-point residual only.
[[nodiscard]] MimeticLocalCellBalance balanceMimeticLocalCell(
    const MimeticLocalCellOperator& localOperator,
    std::span<const double> faceTraceScalars,
    double integratedSource);

} // namespace simwing::fsi::fluid
